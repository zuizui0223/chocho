#!/usr/bin/env python3
"""Audit wing size versus ORIGINAL 224-species source-derived first-index spatial continuity.

Inputs: previously verified *derived* 205-species × distance-bin source summary,
not private raw 83,573-row Parquet. 2011–2022 first INDEXED island records,
strict preceding-year original geographic source, 2001–2024 background.

This does NOT measure actual flight distance, movement mode, port transport,
colonization or arrival; it cannot replace the withheld island×species×year
Parquet. Test species-equal outcomes, not 12,028 independent record events.
"""
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

COHORTS=("all","prior_opportunity","both_period_opportunity")
SOURCE_TOTALS={"all":(205,12028),"prior_opportunity":(194,3813),
               "both_period_opportunity":(193,3136)}
BANDS=("10 km未満","10–100 km","100–500 km","500 km以上")
SEED=20261103
PERM=4999

def read_csv(path):
    with path.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h))
def write_csv(path,rr):
    if not rr:raise RuntimeError("Empty result")
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(rr[0]))
        w.writeheader();w.writerows(rr)
def calculate(source,wing,host,trait):
    taxa={r["taxon_key"]:r["species"] for r in source}
    if len(taxa)!=225:raise RuntimeError("Original speciesKey bridge not unique 225")
    original={r["species"]:r for r in host}
    if len(original)!=239:raise RuntimeError("GEB 239 host analysis panel drift")
    families=defaultdict(set)
    for r in trait:
        sp=r["Species"];fam=r["Family"]
        if fam:families[sp].add(fam)
    x={}
    for r in wing:
        c,k=r["cohort"],r["speciesKey"]
        if c not in COHORTS:raise RuntimeError("Unexpected original island cohort")
        x.setdefault((c,k),{"species":taxa.get(k),"cohort":c,"speciesKey":k,
                         "total_events":int(r["total_events"]),"by":{}})
        if x[(c,k)]["total_events"]!=int(r["total_events"]):raise RuntimeError("Coherent event denominators lost")
        if r["distance_band"] in x[(c,k)]["by"]:
            raise RuntimeError("Duplicate distance band")
        x[(c,k)]["by"][r["distance_band"]]=int(r["events"])
    result=[]
    for cohort in COHORTS:
        now=[v for (c,k),v in x.items() if c==cohort]
        N,total=SOURCE_TOTALS[cohort]
        if len(now)!=N or sum(v["total_events"] for v in now)!=total:
            raise RuntimeError(f"Original source sample changed {cohort}")
        for v in now:
            if set(v["by"])!=set(BANDS):raise RuntimeError("Four distance bins not complete")
            if sum(v["by"].values())!=v["total_events"]:raise RuntimeError("Distance counts do not partition original events")
            sp=v["species"]
            if sp is None or sp not in original or not sp in families or len(families[sp])!=1:
                raise RuntimeError(f"Missing or contradictory speciesKey family/host source {sp}")
            ww=float(original[sp]["wing_size_proxy"])
            if not math.isfinite(ww) or ww<=0:raise RuntimeError("Invalid continuous wing proxy "+sp)
            h=original[sp]
            n=v["total_events"]
            far100=v["by"]["100–500 km"]+v["by"]["500 km以上"]
            far500=v["by"]["500 km以上"]
            result.append({
             "cohort":cohort,
             "speciesKey":v["speciesKey"],"species":sp,
             "family":next(iter(families[sp])),
             "source_island_first_record_events":n,
             "prior_nearest_record_lt10km":v["by"]["10 km未満"],
             "prior_nearest_record_10to100km":v["by"]["10–100 km"],
             "prior_nearest_record_100to500km":v["by"]["100–500 km"],
             "prior_nearest_record_ge500km":far500,
             "prior_nearest_record_ge100km":far100,
             "share_prior_nearest_ge500km":far500/n,
             "share_prior_nearest_ge100km":far100/n,
             "continuous_wing_proxy":ww,
             "native_resource_units":int(h["native_resource_units"]),
             "host_family_count":float(h["host_family_count"]),
             "GEB_introduced_added_resource_units":int(h["introduced_added_units"]),
             "GEB_native_only_resource_baseline":int(h["native_resource_units"]),
             "original_index_events_not_colonizations":1
            })
    return result

def partial_test(rows,outcome,mode):
    n=len(rows)
    if n<28:return {"status":"TOO_FEW_ELIGIBLE_BUTTERFLY_SPECIES","n":n}
    wings=np.array([r["continuous_wing_proxy"] for r in rows],float)
    counts=np.array([r["source_island_first_record_events"] for r in rows],float)
    y_num=np.array([r[outcome] for r in rows],float)
    yy=y_num/counts
    # y = square-root/arcsine transformed species-level fraction, Jeffreys
    # half-event stabilization of exact zero distance events.
    y=np.arcsin(np.sqrt((y_num+.5)/(counts+1)))
    x=np.log(wings)
    family=sorted(set(r["family"] for r in rows))
    zz=[
        np.ones(n),
        np.log1p(counts),
        np.log1p(np.array([r["native_resource_units"] for r in rows],float)),
        np.log1p(np.array([r["host_family_count"] for r in rows],float))
    ]
    for fam in family[1:]:
        zz.append(np.array([float(r["family"]==fam) for r in rows],float))
    Z=np.column_stack(zz)
    if np.linalg.matrix_rank(Z)!=Z.shape[1]:
        return {"status":"RANK_DEFICIENT_COVARIATES","n":n}
    # Partial regression on species-equal y and log wingspan.
    Q,_=np.linalg.qr(Z,mode="reduced")
    yr=y-Q@(Q.T@y)
    def slope(xx):
        z=xx-Q@(Q.T@xx)
        d=float(z@z)
        return float(z@yr/d) if d>1e-12 else None
    observed=slope(x)
    raw=spearmanr(wings,yy)
    # Residualized rank correlation sensitivity: continuous wing proxy
    # and event-ratio rank share independently partialed against Z.
    def ranks(xx):
        from scipy.stats import rankdata
        return rankdata(xx).astype(float)
    zr=ranks(wings);wr=ranks(yy)
    zr=zr-Q@(Q.T@zr);wr=wr-Q@(Q.T@wr)
    residual_rank_r=float(zr@wr/np.sqrt((zr@zr)*(wr@wr)))
    group=defaultdict(list)
    if mode=="global":
        group["all"]=list(range(n))
    elif mode=="family":
        for i,r in enumerate(rows):group[r["family"]].append(i)
    elif mode=="family_event_intensity_halves":
        famvalues=defaultdict(list)
        for i,r in enumerate(rows):famvalues[r["family"]].append(counts[i])
        cutoff={k:float(np.median(v)) for k,v in famvalues.items()}
        for i,r in enumerate(rows):
            group[r["family"]+"|"+("small" if counts[i]<=cutoff[r["family"]] else "large")].append(i)
    else:raise ValueError("unknown permutation design")
    rng=random.Random(SEED)
    movers=[arr for arr in group.values() if len(arr)>=2 and len(set(x[arr]))>1]
    null=[]
    for it in range(PERM):
        xp=x.copy()
        for inds in movers:
            order=inds.copy();rng.shuffle(order)
            xp[inds]=x[order]
        beta=slope(xp)
        if beta is not None:null.append(beta)
    if not null:
        return {"status":"NO_EXCHANGEABLE_PERMUTATIONS","n":n}
    mu=float(np.mean(null))
    obs=float(observed)
    p2=(1+sum(abs(b-mu)>=abs(obs-mu)-1e-12 for b in null))/(1+len(null))
    ph=(1+sum(b>=obs-1e-12 for b in null))/(1+len(null))
    # Wing quartiles provide a non-regression description across species.
    edges=np.quantile(wings,[.25,.5,.75])
    quartiles=[]
    for q in range(4):
        ids=[i for i,v in enumerate(wings)
             if (v>edges[q-1] if q>0 else True) and (v<=edges[q] if q<3 else True)]
        quartiles.append({
           "wing_quartile":q+1,"species":len(ids),
           "mean_species_equal_far_event_fraction":float(np.mean(yy[ids])) if ids else None,
           "event_weighted_far_fraction":
               float(np.sum(y_num[ids])/np.sum(counts[ids])) if ids else None,
           "wing_proxy_median":float(np.median(wings[ids])) if ids else None,
        })
    return {
     "status":"SPECIES_EQUAL_FIRST_INDEX_DISTANCE_PROXY_ANALYSIS",
     "n_butterfly_species":n,
     "species_event_minimum":min(counts),
     "total_original_first_index_events":int(counts.sum()),
     "outcome":"first indexed island records whose earlier indexed nearest coordinate was 100+ or 500+ km; NOT butterfly physical dispersal",
     "source_distance_threshold_km":100 if outcome.endswith("100km") else 500,
     "raw_spearman_rho_wing_vs_species_far_fraction":float(raw.statistic),
     "raw_spearman_p_naive":float(raw.pvalue),
     "partial_spearman_rank_control":residual_rank_r,
     "family_count":len(family),
     "partial_regression_arcsin_sqrt_beta_per_ln_wing":obs,
     "nuisance":["log1p_n_original_source_first_records","log1p_GEB_native_resource_units",
                 "log1p_HOSTS_host_family_breadth","LepTraits_family_fixed_effect"],
     "permute_wing_within":mode,
     "null":{"valid_permutations":len(null),"exchangeable_species":len({i for inds in movers for i in inds}),
             "exchangeable_blocks":len(movers),"null_mean":mu,
             "q025":float(np.quantile(null,.025)),"q975":float(np.quantile(null,.975)),
             "p_two_sided_centered":p2,"p_high":ph},
     "wing_quartiles":quartiles,
     "limits":[
       "Species-equal event-ratio outcome is not island first-colonization probability, presence, absence, or dispersal distance.",
       "Wingspan is a limited proxy for butterfly flight ability; natural rafting, wind, horticultural transit, shipping and air movement are not distinguished.",
       "The 205 species have source-selected 2011–2022 indexed first-record island events and earlier-year points; no full island-year denominator or unrecorded species.",
       "Family adjustment does not remove within-family phylogenetic dependence and source reporting heterogeneity.",
       "The 3 original nested cohorts are not independent; tests and thresholds are exploratory and multiplicity exposed.",
       "No geographic port accessibility, modern climate, host introductions or vessel routes are contained in this original event distance summary.",
     ]
    }

def main():
    p=argparse.ArgumentParser()
    for v in ("species-distance-summary","species-key-bridge","geb-species-metrics","leptraits","out-dir"):
        p.add_argument("--"+v,required=True,type=Path)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    species_data=calculate(read_csv(a.species_key_bridge),
                           read_csv(a.species_distance_summary),
                           read_csv(a.geb_species_metrics),
                           read_csv(a.leptraits))
    write_csv(a.out_dir/"original205_wing_distance_species_source.csv",species_data)
    analyzed={}
    for cohort in COHORTS:
        analyzed[cohort]={}
        for threshold in ("prior_nearest_record_ge100km","prior_nearest_record_ge500km"):
            analyzed[cohort][threshold]={}
            for min_n in (1,10):
                inp=[r for r in species_data if r["cohort"]==cohort and
                     r["source_island_first_record_events"]>=min_n]
                analyzed[cohort][threshold][str(min_n)]={
                  key:partial_test(inp,threshold,mode)
                  for key,mode in (
                    ("global","global"),("family","family"),
                    ("family_and_event_intensity","family_event_intensity_halves")
                  )
                }
    summary={
     "schema":"chocho_original205_island_first_index_distance_vs_continuous_wingspan_v0.1",
     "status":"SOURCE_DERIVED_FROM_PRIVATE_PARQUET_PREVIOUSLY_VERIFIED_AGGREGATES",
     "source":{
       "private_original_224_species_83573_island_years_sha256":"8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f",
       "source_subdataset":"provenance/reviewer_defenses/results/original_island_continuity_20261030/global_original_continuity_species_counts.csv",
       "source_original_event_checksum":"e374df514117e6091cd3426af59a1ba5401199b1ebf5e4a0b994a45886a2d09f",
       "joined_species_name_key":"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv",
       "joined_plant_resource_wings":"data/frozen/figure_sources/anthropogenic_species_metrics.csv",
       "joined_family":"data/external/leptraits_consensus_v1.0.csv",
       "local_raw_224species_Parquet_NOT_read_this_run":True
     },
     "cohort_event_audits":SOURCE_TOTALS,
     "prespecified_primary":"both_period_opportunity / >=100 km / species >=10 first-index events; family x event intensity half source permutation",
     "exploratory_multi_cohort_sensitivities":analyzed,
     "no_flight_or_shipment_claim":True
    }
    (a.out_dir/"original205_wing_distance_audit.json").write_text(
       json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    p0=analyzed["both_period_opportunity"]["prior_nearest_record_ge100km"]["10"]["family_and_event_intensity"]
    print("PRIMARY_ORIGINAL205_WING_VS_DISTANCE_INDEX",json.dumps(p0,ensure_ascii=False,indent=2))
    print("ALL_COHORTS_SOURCE_LOCKS",json.dumps({k:{"n":len([r for r in species_data if r["cohort"]==k]),
                                               "events":sum(r["source_island_first_record_events"] for r in species_data if r["cohort"]==k)}
                                         for k in COHORTS},ensure_ascii=False))
if __name__=="__main__":main()
