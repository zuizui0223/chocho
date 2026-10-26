#!/usr/bin/env python3
"""Posthoc independent-seed and site/pair fragility of annual adult co-fluctuation.

Reuses exact archived host and original US transect source. Raw site rows and
pairs are local only; only source aggregate diagnostics and 199 null scalar
histories are saved. Does not infer competition.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from analyze_sbus_illinois_niche_counts_20261025 import (
    frozen_focal,load_source,read_events
)
from analyze_sbus_illinois_host_degree_null_20261025 import (
    accepted_edges,butterfly_families,draw_graph
)
from analyze_sbus_illinois_adult_codetection_hostnull_20261026 import (
    field_panel,pair_correlations,meancontrast
)

ROOT=Path(__file__).resolve().parents[1]
SEED=20261027
DRAWS=199
ORIG_DELTA=0.14323064389607087
ORIG_SHARED=62
ORIG_NONSHARED=578
ORIG_TEMPORAL_PAIRS=640
ORIG_SITES=121
ORIGINAL_RECEIPT=ROOT/"provenance/reviewer_defenses/results/sbus_illinois_adult_codetection_hostnull_20261026.json"


def per_site_temporal_contributions(good,sitekeys,pair_records):
    names=sorted(set(sp for p in pair_records for sp in p[:2]))
    by_site=defaultdict(list)
    for site,year in sitekeys:
        by_site[site].append(year)
    sites=sorted(by_site)
    pairs=pair_records
    m=len(pairs)
    contributions=np.zeros((len(sites),m,4),dtype=float)
    site_counts=np.zeros(m,dtype=int)
    loc={sp:i for i,sp in enumerate(names)}
    pairs_ij=[(loc[p[0]],loc[p[1]]) for p in pairs]
    for sitei,site in enumerate(sites):
        years=sorted(by_site[site])
        mat=np.zeros((len(years),len(names)),dtype=float)
        had=np.zeros(len(names),dtype=bool)
        for t,year in enumerate(years):
            effort,counts=good[(site,year)]
            for sp,cnt in counts.items():
                i=loc.get(sp)
                if i is not None and cnt>0:
                    had[i]=True
                    mat[t,i]=math.log1p(cnt/effort)
        mat-=mat.mean(axis=0,keepdims=True)
        for pidx,(i,j) in enumerate(pairs_ij):
            if not (had[i] and had[j]):
                continue
            x=mat[:,i]
            y=mat[:,j]
            contributions[sitei,pidx]=[
                float(x@y),float(x@x),float(y@y),float(len(years))
            ]
            site_counts[pidx]+=1
    total=np.sum(contributions,axis=0)
    if np.any(total[:,1]<=1e-10) or np.any(total[:,2]<=1e-10):
        raise ValueError("An original eligible temporal pair lost all variance")
    if np.any(site_counts<5) or np.any(total[:,3]<30):
        raise ValueError("Original temporal pair cohort not supported")
    correlations=total[:,0]/np.sqrt(total[:,1]*total[:,2])
    raw=np.asarray([p[2] for p in pair_records],dtype=float)
    if not np.allclose(correlations,raw,rtol=0,atol=1e-9):
        raise ValueError("Independent temporal site components don't reproduce original cofluctuations")
    return sites,contributions,total,site_counts


def group_average(values,weights,share):
    mask=np.asarray(share,dtype=bool)
    if mask.sum()<20 or (~mask).sum()<20:
        return None
    if weights is None:
        return float(np.mean(values[mask])-np.mean(values[~mask]))
    return float(np.average(values[mask],weights=weights[mask])-
                 np.average(values[~mask],weights=weights[~mask]))


def null_summary(draws,obs):
    values=np.asarray(draws,dtype=float)
    if len(values)<190:
        return {"status":"NOT_EVALUABLE_FEWER_THAN_190_GRAPH_SAMPLES","n":len(values)}
    q=np.quantile(values,[.025,.5,.975],method="linear").tolist()
    return {
        "status":"SOURCE_CONDITIONED_GRAPH_NULL_QUALIFIED",
        "n_graphs":len(values),
        "null_q025_median_q975":q,
        "observed":obs,
        "observed_minus_null_median":obs-q[1],
        "left_tail_fraction":(1+sum(x<=obs for x in draws))/(len(draws)+1),
        "right_tail_fraction":(1+sum(x>=obs for x in draws))/(len(draws)+1),
        "history":list(map(float,draws)),
    }


def run(accepted_hosts):
    first=json.loads(ORIGINAL_RECEIPT.read_text(encoding="utf-8"))
    if(first["null_phi"]["status"]!="EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED"
       or first["sample"]["within_site_temporal_pair_count"]!=ORIG_TEMPORAL_PAIRS
       or first["within_site_temporal"]["observed"]["shared"]!=ORIG_SHARED
       or abs(first["within_site_temporal"]["observed"]["difference"]-ORIG_DELTA)>1e-12):
        raise ValueError("First-seed original source receipts drifted")
    focal=frozen_focal()
    families=butterfly_families(focal)
    hostmap,edges,groups,viable,plantfam=accepted_edges(accepted_hosts,focal,families)
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/"source.zip"
        load_source(path)
        events,_=read_events(path,focal)
    good,keys,species,vector,field=field_panel(events)
    pairs=pair_correlations(good,keys,species,vector)
    if(pairs["status"]!="FIXED_SPECIES_PAIRS_READY" or
       pairs["eligible_butterflies"]!=39 or
       len(pairs["within_site_temporal_candidate_pairs"])!=ORIG_TEMPORAL_PAIRS):
        raise ValueError("Original repeated-site species and 640 temporal-pair cohort drifted")
    records=pairs["within_site_temporal_candidate_pairs"]
    observed=meancontrast(records,hostmap,index=2)
    if(observed["shared"]!=ORIG_SHARED or observed["nonshared"]!=ORIG_NONSHARED or
       abs(observed["difference"]-ORIG_DELTA)>1e-10):
        raise ValueError("Original time-synchronous cofluctuation contrast drifted")
    sites,parts,total,site_counts=per_site_temporal_contributions(good,keys,records)
    weights=total[:,3]
    data_corr=np.asarray([p[2] for p in records],dtype=float)
    H=np.asarray([bool(hostmap[a]&hostmap[b]) for a,b,*rest in records],dtype=bool)
    if len(sites)!=ORIG_SITES:
        raise ValueError("Independent Illinois site sample drifted")
    orig_weighted=group_average(data_corr,weights,H)
    if orig_weighted is None:
        raise ValueError("Original weighted pair categories too small")

    # Fixed-pair leave-one-site-out resampling. No new pair eligibility testing.
    site_jacks=[]
    jack_fail=0
    for i in range(len(sites)):
        remaining=total-parts[i]
        if np.any(remaining[:,1]<=1e-10) or np.any(remaining[:,2]<=1e-10):
            jack_fail+=1
            continue
        r=remaining[:,0]/np.sqrt(remaining[:,1]*remaining[:,2])
        val=group_average(r,None,H)
        if val is None or not math.isfinite(val):
            jack_fail+=1
        else:
            site_jacks.append(val)
    excluded_species=sorted(set(sp for p in records for sp in p[:2]))
    species_jacks=[]
    for sp in excluded_species:
        sel=np.asarray([sp not in p[:2] for p in records])
        if (H&sel).sum()<20 or ((~H)&sel).sum()<20:
            continue
        species_jacks.append(group_average(data_corr[sel],None,H[sel]))
    if not site_jacks or not species_jacks:
        raise ValueError("Jackknife source support entirely absent")

    # IMPORTANT: Fixed 640 sampled pairs and original pair correlations.
    # Only the reported HOSTS network assignments are randomized independently.
    equal_null=[]
    weighted_null=[]
    fails=Counter()
    metrics=[]
    for rep in range(DRAWS):
        graph,diag=draw_graph(hostmap,edges,groups,viable,plantfam,families,rep,
                             seed_override=SEED+7919*rep)
        labels=np.asarray([bool(graph[a]&graph[b]) for a,b,*rest in records])
        problems=[]
        if diag["accepted"]<100:problems.append("LESS_THAN_100_ACCEPTED_SWAPS")
        if diag["moved_viable_edges_fraction"]<.1:problems.append("LESS_THAN_10_PERCENT_VIABLE_MOVED")
        e=group_average(data_corr,None,labels)
        w=group_average(data_corr,weights,labels)
        if e is None or w is None:problems.append("PAIR_CLASS_BELOW_20")
        if problems:
            for p in problems:fails[p]+=1
        else:
            equal_null.append(e)
            weighted_null.append(w)
        metrics.append({
            "changed_viable":diag["moved_viable_edges_fraction"],
            "accepted":diag["accepted"],
            "n_shared":int(labels.sum()),
        })
    sources_valid=len(equal_null)>=190 and len(weighted_null)>=190 and sum(fails.values())<=9
    equal=null_summary(equal_null,observed["difference"]) if sources_valid else {
        "status":"NOT_EVALUABLE_SOURCE_NULL_GATE"
    }
    weighted=null_summary(weighted_null,orig_weighted) if sources_valid else {
        "status":"NOT_EVALUABLE_SOURCE_NULL_GATE"
    }
    first_null=first["within_site_temporal"]
    # 3 outcomes were observed: φ, Jaccard and within-site yearly co-fluctuation.
    # The exact observed first right-tail .03 is not a corrected inference.
    original_p=first_null.get("right_tail_fraction_exploratory")
    if original_p is None:
        original_p=first_null.get("right_tail_fraction_exploratory_only")
    meta={
        "schema":"chocho_illinois_host_sharing_temporal_source_fragility_20261026_v01",
        "status":"POSTHOC_SECOND_SEED_AND_FIXED_PAIR_SITE_JACKKNIFE_COMPLETED",
        "original_data_doi":"10.6084/m9.figshare.27934602.v1",
        "focal_original_239_names":239,
        "accepted_hosts_commit":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
        "wcvp_commit":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
        "fixed_field_panel":field,
        "fixed_adult_pair_cohort":len(records),
        "host_sharing_pairs":ORIG_SHARED,
        "nonsharing_pairs":ORIG_NONSHARED,
        "observed_equal_weighted_temporal_contrast":ORIG_DELTA,
        "observed_siteyear_weighted_temporal_contrast":orig_weighted,
        "first_seed_null":{
            "q025_median_q975":first_null.get("null_difference_q025_median_q975"),
            "right_tail_exploratory":original_p,
            "seed":20261026,
        },
        "independent_seed_null_equal_pairs":equal,
        "independent_seed_null_weighted_pairs":weighted,
        "site_jackknife":{
            "sites_total":len(sites),
            "sites_evaluable":len(site_jacks),
            "sites_not_evaluable":jack_fail,
            "fixed_species_pairs_every_replicate":len(records),
            "min_median_max_contrast":np.quantile(site_jacks,[0,.5,1]).tolist(),
            "fraction_positive":float(np.mean(np.asarray(site_jacks)>0)),
        },
        "leave_one_butterfly_out":{
            "source_butterflies":len(excluded_species),
            "models_evaluable":len(species_jacks),
            "min_median_max_contrast":np.quantile(species_jacks,[0,.5,1]).tolist(),
            "fraction_positive":float(np.mean(np.asarray(species_jacks)>0)),
        },
        "network_quality":{
            "requested":DRAWS,
            "accepted_model_draws":len(equal_null),
            "invalid_reason_counts":dict(fails),
            "accept_swap_min_max":[min(d["accepted"] for d in metrics),
                                   max(d["accepted"] for d in metrics)],
            "viable_graph_movement_min_max":[min(d["changed_viable"] for d in metrics),
                                             max(d["changed_viable"] for d in metrics)],
            "randomized_shared_pair_count_range":[min(d["n_shared"] for d in metrics),
                                                   max(d["n_shared"] for d in metrics)],
            "all_null_graphs_same_observed_adult_site_years_and_fixed_pair_scores":True,
        },
        "multiple_endpoint_note":"Original exploratory φ/Jaccard/temporal comparisons all examined; raw tail .03 does not pass Bonferroni alpha .0167 across three metrics, and source-graph randomization is not a causal field p value.",
        "not_proven":[
            "Direct local host-plant consumption, local larval biomass, actual competition",
            "Host-plant introduction caused butterfly growth",
            "Phylogenetic or environmental independence of these related adult butterfly pairs",
            "True zero adult abundance on non-detection in standard transect count data",
        ],
    }
    return meta


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--accepted-host-csv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    r=run(a.accepted_host_csv)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":r["status"],
        "observed_equal":r["observed_equal_weighted_temporal_contrast"],
        "observed_weighted":r["observed_siteyear_weighted_temporal_contrast"],
        "first":r["first_seed_null"],
        "second_equal":{k:v for k,v in r["independent_seed_null_equal_pairs"].items() if k!="history"},
        "second_weighted":{k:v for k,v in r["independent_seed_null_weighted_pairs"].items() if k!="history"},
        "sites":r["site_jackknife"],
        "species":r["leave_one_butterfly_out"],
        "graph":r["network_quality"],
    },ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()
