#!/usr/bin/env python3
"""Does shared larval host identity predict real adult site-year co-detection?

Independent POSTHOC descriptive NICHE/FILTERING question. Original adult counts
and fixed eligible pair correlations are held constant across 199 butterfly-host
bipartite swaps. No inferred causal competition or butterfly recruitment.

Protocol: provenance/reviewer_defenses/protocols/SBUS_ILLINOIS_COOCCURRENCE_HOST_DEGREE_NULL_20261026.md
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from analyze_sbus_illinois_niche_counts_20261025 import (
    PROGRAMS, MIN_EVENTS, MIN_SEASON_BINS, frozen_focal, load_source, read_events,
)
from analyze_sbus_illinois_host_degree_null_20261025 import (
    accepted_edges, butterfly_families, draw_graph
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DOI = "10.6084/m9.figshare.27934602.v1"
N_SIM = 199
SEED = 20261026


def field_panel(events):
    siteyears = defaultdict(list)
    for (program,event),row in events.items():
        if program!=PROGRAMS[0] or row["bad_duration"]:
            continue
        siteyears[(row["site"],row["year"])].append(row)
    good={}
    for (site,year),rows in siteyears.items():
        if (len(rows)<MIN_EVENTS or
            len({(r["doy"]-1)//31 for r in rows})<MIN_SEASON_BINS):
            continue
        effort=sum(r["duration"] for r in rows)/60.0
        if effort<=0 or not math.isfinite(effort):
            continue
        counts=defaultdict(float)
        for row in rows:
            for butterfly,number in row["counts"].items():
                counts[butterfly]+=number
        good[(site,year)]=(effort,counts)
    yearset=defaultdict(set)
    for site,year in good:
        yearset[site].add(year)
    keep={site for site,years in yearset.items()
          if any(y+1 in years and y+2 in years for y in years)}
    good={key:row for key,row in good.items() if key[0] in keep}
    sitekeys=sorted(good)
    if len(keep)!=121 or len(sitekeys)!=851:
        raise ValueError(f"Original 121 sites / 851 site-years not reproduced: {len(keep)}, {len(sitekeys)}")
    species=sorted({sp for effort,counts in good.values()
                    for sp,n in counts.items() if n>0})
    if len(species)!=47:
        raise ValueError(f"Original 47 adult butterflies not reproduced: {len(species)}")
    vector=np.zeros((len(sitekeys),len(species)),dtype=np.int32)
    loc={sp:i for i,sp in enumerate(species)}
    for i,key in enumerate(sitekeys):
        for sp,n in good[key][1].items():
            if n>0:
                vector[i,loc[sp]]=1
    events_audit={
        "qualified_sites":len(keep),
        "qualified_site_years":len(sitekeys),
        "matched_239_butterflies":len(species),
        "source_presence_cells":int(np.sum(vector)),
        "eligible_survey_method":"Illinois Butterfly Monitoring Network, 5+ events, 3+ seasonal 31-day bins, 3 consecutive eligible years",
    }
    return good,sitekeys,species,vector,events_audit


def pair_correlations(good,sitekeys,species,B):
    presence=B.sum(axis=0)
    n=B.shape[0]
    keep=(presence>=10)&(presence<=n-10)
    ids=np.flatnonzero(keep).tolist()
    if len(ids)<15:
        return {"status":"NOT_EVALUABLE_TOO_FEW_ORIGINAL_BUTTERFLIES","n":len(ids)}
    candidate_names=[species[i] for i in ids]
    original=B[:,ids].astype(np.float64)
    ncols=len(ids)
    P=original.sum(axis=0)
    z=(original-P[None,:]/n)
    cov=z.T@z
    diag=np.sqrt(np.maximum(0,np.diag(cov)))
    with np.errstate(divide="ignore",invalid="ignore"):
        phi=cov/np.outer(diag,diag)
    observed_both=original.T@original
    union=P[:,None]+P[None,:]-observed_both
    with np.errstate(divide="ignore",invalid="ignore"):
        jac=observed_both/union

    # Secondary: same-site adult count/hour temporal anomaly correlation.
    # Every site-specific species series is demeaned before pooling across sites.
    loc={sp:i for i,sp in enumerate(candidate_names)}
    temporal_xy=np.zeros((ncols,ncols),dtype=float)
    temporal_x2=np.zeros((ncols,ncols),dtype=float)
    temporal_y2=np.zeros((ncols,ncols),dtype=float)
    temporal_yrs=np.zeros((ncols,ncols),dtype=int)
    temporal_sites=np.zeros((ncols,ncols),dtype=int)
    site_entries=defaultdict(list)
    for st,y in sitekeys:
        site_entries[st].append(y)
    for st,years in site_entries.items():
        years=sorted(years)
        mat=np.zeros((len(years),ncols),dtype=float)
        had=np.zeros(ncols,dtype=bool)
        for i,yy in enumerate(years):
            effort,counts=good[(st,yy)]
            for sp,value in counts.items():
                idx=loc.get(sp)
                if idx is not None and value>0:
                    had[idx]=True
                    mat[i,idx]=math.log1p(value/effort)
        mat=mat-mat.mean(axis=0,keepdims=True)
        candidates=np.flatnonzero(had)
        for ii,i in enumerate(candidates):
            for j in candidates[ii+1:]:
                temporal_xy[i,j]+=float(mat[:,i]@mat[:,j])
                temporal_x2[i,j]+=float(mat[:,i]@mat[:,i])
                temporal_y2[i,j]+=float(mat[:,j]@mat[:,j])
                temporal_yrs[i,j]+=len(years)
                temporal_sites[i,j]+=1
    pairs=[]
    temporals=[]
    for i in range(ncols):
        for j in range(i+1,ncols):
            r=float(phi[i,j])
            v=float(jac[i,j])
            if not math.isfinite(r) or not math.isfinite(v):
                raise ValueError("Binary phi/Jaccard undefined after 10-site-year prevalence eligibility")
            pairs.append((candidate_names[i],candidate_names[j],r,v))
            if (temporal_sites[i,j]>=5 and temporal_yrs[i,j]>=30 and
                temporal_x2[i,j]>1e-10 and temporal_y2[i,j]>1e-10):
                val=temporal_xy[i,j]/math.sqrt(temporal_x2[i,j]*temporal_y2[i,j])
                if abs(val)>1+1e-9:
                    raise ValueError("Impossible within-site residual correlation")
                temporals.append((candidate_names[i],candidate_names[j],float(np.clip(val,-1,1))))
    return {
        "status":"FIXED_SPECIES_PAIRS_READY",
        "eligible_butterflies":len(ids),
        "positive_presence_range":[int(min(presence[ids])),int(max(presence[ids]))],
        "eligible_pair_count":len(pairs),
        "pairs":pairs,
        "within_site_temporal_candidate_pairs":temporals,
    }


def meancontrast(records,hostmap,index=2):
    share=[];non=[]
    for pair in records:
        a,b=pair[:2]
        (share if hostmap[a]&hostmap[b] else non).append(pair[index])
    if not share or not non:
        return {"status":"NO_VALID_BOTH_HOST_CLASSIFICATION","shared":len(share),"nonshared":len(non)}
    return {
        "status":"OK",
        "shared":len(share),
        "nonshared":len(non),
        "mean_shared":float(np.mean(share)),
        "mean_nonshared":float(np.mean(non)),
        "difference":float(np.mean(share)-np.mean(non)),
    }


def null_test(records, observed, base_graph, edge_pairs, all_groups, viable,
              plantfam, families,metric_index,seed_shift=0):
    if observed["status"]!="OK" or observed["shared"]<20 or observed["nonshared"]<20:
        return {"status":"NOT_EVALUABLE_OBSERVED_CLASS_COVERAGE","observed":observed}
    samples=[]
    failure_reasons=Counter()
    diagnostics=[]
    for i in range(N_SIM):
        hostgraph,diag=draw_graph(
            base_graph,edge_pairs,all_groups,viable,plantfam,families,i,
            seed_override=SEED+seed_shift+7919*i,
        )
        r=meancontrast(records,hostgraph,metric_index)
        reasons=[]
        if diag["accepted"]<100:reasons.append("ACCEPTED_SWAP_LT100")
        if diag["moved_viable_edges_fraction"]<.1:reasons.append("MOVED_LT10PCT")
        if r["status"]!="OK" or min(r["shared"],r["nonshared"])<20:
            reasons.append("CLASS_LT20")
        if reasons:
            for reason in reasons:failure_reasons[reason]+=1
        else:
            samples.append(r["difference"])
        diagnostics.append({
            "replicate":i,"shared":r["shared"],"nonshared":r["nonshared"],
            "swaps":diag["accepted"],"moved_viable":diag["moved_viable_edges_fraction"],
            "invalid":reasons,
        })
    eligible=len(samples)>=190 and sum(failure_reasons.values())<=9
    quantiles=np.quantile(samples,[.025,.5,.975]).tolist() if eligible else None
    return {
        "status":("EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED"
                  if eligible else "NOT_EVALUABLE_NETWORK_COVERAGE_OR_SWAP_GATE"),
        "observed":observed,
        "null_draws_requested":N_SIM,
        "null_draws_valid":len(samples),
        "null_invalid_reasons":dict(failure_reasons),
        "null_difference_q025_median_q975":quantiles,
        "observed_minus_null_median":float(observed["difference"]-quantiles[1])
                                    if eligible else None,
        "right_tail_fraction_exploratory":(
            (1+sum(x>=observed["difference"] for x in samples))/(len(samples)+1)
            if eligible else None),
        "left_tail_fraction_exploratory":(
            (1+sum(x<=observed["difference"] for x in samples))/(len(samples)+1)
            if eligible else None),
        "null_difference_history":samples if eligible else None,
        "n_shared_pairs_range":[min(x["shared"] for x in diagnostics),
                                max(x["shared"] for x in diagnostics)],
        "range_accepted_swaps":[min(x["swaps"] for x in diagnostics),
                                max(x["swaps"] for x in diagnostics)],
        "range_changed_viable_edge_fraction":[min(x["moved_viable"] for x in diagnostics),
                                               max(x["moved_viable"] for x in diagnostics)],
        "pair_measurements_and_sample_identical_across_all_null_states":True,
        "not_a_causal_field_filter_or_competition_effect":True,
    }


def run(hostcsv):
    focal=frozen_focal()
    families=butterfly_families(focal)
    hostmap,edges,groups,viable,plantfam=accepted_edges(hostcsv,focal,families)
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/"source.zip"
        load_source(path)
        events,_=read_events(path,focal)
    good,keys,species,vector,field=field_panel(events)
    pairdata=pair_correlations(good,keys,species,vector)
    out={
        "schema":"chocho_illinois_niche_coded_field_codetection_v01_20261026",
        "source_doi":SOURCE_DOI,
        "source_file_md5":"3161b1c8802a0d67771b2044e18da928",
        "host_source_git_commit":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
        "wcvp_git_commit":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
        "original_focal":len(focal),
        "adult_field_panel":field,
        "original_network_edges":len(edges),
        "swappable_host_network_edges":sum(len(v) for v in viable.values()),
        "repeated_site_year_imputed_non_detection_not_validated_absence":True,
        "source_site_year_or_species_pair_raw_values_not_saved":True,
        "source_legal":"Henry et al. (2025), CC BY 4.0; preserve full author attribution in any paper",
        "limitation":[
            "Adult positive counts under standardized survey effort are species detection; no local plant checks",
            "Lack of species in reported survey events is treated as a detection-zero proxy, not observed absence",
            "Host overlap from HOSTS literature is not verified plant feeding at source sites",
            "Host-network null holds degrees within broad plant/butterfly families, but not climatic niche or phylogeny",
            "Co-detection associations and family-conditioned source null do not prove competition, facilitation or causal global spread",
        ],
    }
    if pairdata["status"]!="FIXED_SPECIES_PAIRS_READY":
        out["status"]=pairdata["status"]
        out["pair_source_gate"]=pairdata
        return out
    pairs=pairdata["pairs"]
    ph=meancontrast(pairs,hostmap,2)
    jac=meancontrast(pairs,hostmap,3)
    out["sample"]={
        "eligible_butterflies":pairdata["eligible_butterflies"],
        "positive_detection_site_year_range":pairdata["positive_presence_range"],
        "fixed_species_pairs":len(pairs),
        "contemporaneous_phi":ph,
        "secondary_jaccard":jac,
        "within_site_temporal_pair_count":len(pairdata["within_site_temporal_candidate_pairs"]),
    }
    out["null_phi"]=null_test(pairs,ph,hostmap,edges,groups,viable,plantfam,families,2)
    # Jaccard on the SAME graphs as primary in an independent deterministic
    # replay using the same graph seeds; secondary only and not a rescue test.
    if jac["status"]=="OK":
        out["null_jaccard"]=null_test(pairs,jac,hostmap,edges,groups,viable,plantfam,families,3)
    timed=pairdata["within_site_temporal_candidate_pairs"]
    out["within_site_temporal"]={
        "status":"NOT_EVALUABLE_BELOW_PREDECLARED_PAIR_COVERAGE",
        "fixed_pairs":len(timed),
    }
    if timed:
        te=meancontrast(timed,hostmap,2)
        out["within_site_temporal"]["observed"]=te
        if te["status"]=="OK" and te["shared"]>=20 and te["nonshared"]>=20:
            out["within_site_temporal"]=null_test(timed,te,hostmap,edges,groups,
                                                 viable,plantfam,families,2)
    out["status"]="POSTHOC_FIELD_CODETECTION_NICHE_SOURCE_AUDIT_COMPLETE"
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--accepted-host-csv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=run(a.accepted_host_csv)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                     encoding="utf-8")
    brief={
        "status":result["status"],
        "panel":result["adult_field_panel"],
        "sample":result.get("sample"),
        "primary_null":{k:result["null_phi"].get(k) for k in (
            "status","null_draws_valid","null_difference_q025_median_q975",
            "right_tail_fraction_exploratory")} if "null_phi" in result else None,
        "within_site_temporal":{k:result["within_site_temporal"].get(k) for k in (
            "status","fixed_pairs","null_difference_q025_median_q975")} if
            "within_site_temporal" in result else None,
    }
    print(json.dumps(brief,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()
