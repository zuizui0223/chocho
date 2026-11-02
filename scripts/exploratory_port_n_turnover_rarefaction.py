#!/usr/bin/env python3
"""Source-locked butterfly host vs other plant N: observation-effort rarefaction.

Input: previously downloaded GH Actions artifacts, no live GBIF requests.
Recorded plant speciesKey occurrences are a GBIF indexing proxy, not
independent plants, known abundance, confirmed host use, or colonization.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

PERIODS = ("early_2001_2010", "later_2011_2024")
SIDES = ("port", "nonport")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_inputs(raw_path, mapping_path, aggregate_path):
    raw = [json.loads(x) for x in Path(raw_path).read_text(encoding="utf-8").splitlines() if x.strip()]
    mapping = read_json(mapping_path)
    agg = read_json(aggregate_path)
    if len(raw)!=312 or agg["site_period_facets"]!=312:
        raise ValueError("Original WPI78 source does not cover 312 site periods")
    if (agg["predeclared_source_pairs"], agg["independent_islands"])!=(78,35):
        raise ValueError("Unexpected WPI78 geometry or island support")
    if any(v !=0 for v in agg["drift_counts"].values()):
        raise ValueError("GBIF index drift already observed")
    rawdict={}
    for r in raw:
        k=(str(r["island_id"]),int(r["port_id"]),r["side"],r["period"])
        if k in rawdict:raise ValueError("Duplicated raw site-period species facets")
        rawdict[k]=r["facet"]
    obs={}
    for r in agg["observations"]:
        k=(str(r["island_id"]),int(r["port_id"]),r["side"],r["period"])
        if k in obs: raise ValueError("Duplicated summary site-period")
        obs[k]=r
    if rawdict.keys()!=obs.keys():
        raise ValueError("Raw and summary facet keys do not match")
    for k,r in rawdict.items():
        if int(r["indexed_plant_records"])!=obs[k]["plant_records"]:
            raise ValueError(f"GBIF record count changed {k}")
        if int(r["indexed_species_keys"])!=obs[k]["all_plant_species_keys"]:
            raise ValueError(f"GBIF species-key count changed {k}")
    bykey={int(k):v for k,v in mapping["accepted_unambiguous_by_speciesKey"].items()}
    focal={int(k) for k in mapping["focal_host_by_speciesKey"]}
    if not focal.issubset(bykey):
        raise ValueError("Focal host keys not present in trait mapping")
    return rawdict,obs,bykey,focal


def source_eligible_fixed_pairs(observations):
    bypair=defaultdict(dict)
    for (island,port,side,period),r in observations.items():
        if (r["plant_records"]>=30 and r["N_scored_candidate_host_species"]>=2
            and r["N_scored_other_plant_species"]>=5):
            bypair[(island,port)][(side,period)]=r
    return sorted((k for k,v in bypair.items() if len(v)==4),
                  key=lambda x:(int(x[0]),x[1]))


def exact_set_turnover(raw,trait,focal,pairs):
    result=defaultdict(list)
    for island,port in pairs:
        for side in SIDES:
            a=raw[(island,port,side,PERIODS[0])]["faceted_specieskey_counts"]
            b=raw[(island,port,side,PERIODS[1])]["faceted_specieskey_counts"]
            ka={int(x) for x,v in a.items() if int(v)>0 and int(x) in trait}
            kb={int(x) for x,v in b.items() if int(v)>0 and int(x) in trait}
            for g,subset in (("focal_host_N",focal),
                             ("other_plant_N",set(trait)-focal)):
                before=ka&subset
                after=kb&subset
                shared=before&after
                union=before|after
                result[(island,side,g)].append({
                    "island_id":island,"port_id":port,"side":side,"group":g,
                    "old_species":len(before),"late_species":len(after),
                    "shared_species":len(shared),
                    "newly_indexed_species":len(after-before),
                    "lost_species":len(before-after),
                    "jaccard":len(shared)/len(union) if union else None,
                    "fraction_late_newly_indexed":len(after-before)/len(after) if after else None,
                    "fraction_early_not_in_late":len(before-after)/len(before) if before else None
                })
    allrows=[r for q in result.values() for r in q]
    summary={}
    for g in ("focal_host_N","other_plant_N"):
        summary[g]={}
        for side in SIDES:
            subset=[r for r in allrows if r["group"]==g and r["side"]==side]
            byisland=defaultdict(list)
            for r in subset:
                if r["fraction_late_newly_indexed"] is not None:
                    byisland[r["island_id"]].append(r["fraction_late_newly_indexed"])
            vals=[statistics.mean(v) for v in byisland.values()]
            summary[g][side]={
                "pairs":len(subset),
                "independent_islands_with_late_records":len(vals),
                "island_equal_mean_late_only_species_fraction":statistics.mean(vals) if vals else None,
                "island_equal_median_late_only_species_fraction":statistics.median(vals) if vals else None,
                "no_biological_first_arrival_inference":True
            }
    return allrows,summary


def sampled_keys_from_records(counts,rng,depth):
    """Exact simple random sample without replacement from GBIF indexed records.

    No expansion of the entire speciesKey-by-count vector is required.
    """
    ordered=sorted((int(k),int(v)) for k,v in counts.items() if int(v)>0)
    keys=[];cum=[];total=0
    for key,count in ordered:
        keys.append(key)
        total+=count
        cum.append(total)
    if depth>total:
        raise ValueError("Not enough speciesKey-indexed records to rarefy")
    draws=rng.sample(range(total),depth)
    return {keys[bisect.bisect_right(cum,index)] for index in draws}


def mean_group_N(keys,traits,focal):
    host=[traits[k]["N"] for k in keys if k in focal and k in traits]
    other=[traits[k]["N"] for k in keys if k not in focal and k in traits]
    if len(host)<2 or len(other)<5:
        return None
    return statistics.mean(host)-statistics.mean(other)


def island_temporal_rarefaction(raw,traits,focal,pairs,*,draws=399,depth_cap=120,seed=312211):
    if depth_cap<5 or draws<1:raise ValueError("Invalid fixed rarefaction design")
    template=[]
    for island,port in pairs:
        rawsites=[raw[(island,port,side,period)]["faceted_specieskey_counts"]
                  for side in SIDES for period in PERIODS]
        n_per_site=[sum(int(v) for v in site.values()) for site in rawsites]
        common=min(depth_cap,*n_per_site)
        if common<30:
            template.append({"island_id":island,"port_id":port,
                             "records_per_site":n_per_site,"depth":common,
                             "eligible_min_30":False})
        else:
            template.append({"island_id":island,"port_id":port,
                             "records_per_site":n_per_site,"depth":common,
                             "eligible_min_30":True})
    rng=random.Random(seed)
    statistics_by_draw=[]
    for _ in range(draws):
        values=defaultdict(list)
        evaluated=0
        for item in template:
            if not item["eligible_min_30"]:continue
            island,port,depth=item["island_id"],item["port_id"],item["depth"]
            gaps={}
            valid=True
            for side in SIDES:
                for period in PERIODS:
                    site=raw[(island,port,side,period)]["faceted_specieskey_counts"]
                    detected=sampled_keys_from_records(site,rng,depth)
                    v=mean_group_N(detected,traits,focal)
                    if v is None:
                        valid=False
                        break
                    gaps[(side,period)]=v
                if not valid:break
            if not valid:continue
            d=((gaps[("port",PERIODS[1])]-gaps[("nonport",PERIODS[1])])-
               (gaps[("port",PERIODS[0])]-gaps[("nonport",PERIODS[0])]))
            values[island].append(d)
            evaluated+=1
        islandmeans=[statistics.mean(x) for x in values.values()]
        statistics_by_draw.append({
            "eligible_port_pairs":evaluated,
            "eligible_independent_islands":len(values),
            "island_equal_mean":statistics.mean(islandmeans) if islandmeans else None,
            "all_island_means_positive":bool(islandmeans) and all(x>0 for x in islandmeans)
        })
    supported=[x for x in statistics_by_draw if x["eligible_independent_islands"]==6]
    all_pos=[x for x in statistics_by_draw if x["island_equal_mean"] is not None]
    def quantile(data,q):
        z=sorted(data)
        if not z:return None
        return z[round((len(z)-1)*q)]
    return {
        "design":"with-four-cell-within-pair exact random sample without replacement of GBIF indexed record counts",
        "draws":draws,"random_seed":seed,
        "per_site_record_cap":depth_cap,
        "pairs_that_can_reach_30_specieskey_records":sum(x["eligible_min_30"] for x in template),
        "original_eligible_pairs":len(template),
        "common_per_site_depth_min_median_max":[
            min(x["depth"] for x in template),
            statistics.median(x["depth"] for x in template),
            max(x["depth"] for x in template)
        ],
        "simulations_with_all_six_islands_supported":len(supported),
        "simulations_with_at_least_one_estimate":len(all_pos),
        "fraction_positive_mean_given_at_least_one_island":(
            sum(x["island_equal_mean"]>0 for x in all_pos)/len(all_pos)
            if all_pos else None),
        "mean_of_island_mean_estimates_when_all_six_survive":(
            statistics.mean(x["island_equal_mean"] for x in supported)
            if supported else None),
        "fraction_all_six_islands_positive_given_full_support":(
            sum(x["all_island_means_positive"] for x in supported)/len(supported)
            if supported else None),
        "quantiles_island_count":[quantile([x["eligible_independent_islands"] for x in statistics_by_draw],q)
                                  for q in (0,.25,.5,.75,1)],
        "quantiles_paired_port_count":[quantile([x["eligible_port_pairs"] for x in statistics_by_draw],q)
                                      for q in (0,.25,.5,.75,1)],
        "no_inference_if_sampling_loses_required_trait_groups":True,
        "rarefaction_is_only_detection_sensitivity_not_physical_plant_abundance":True,
        "individual_draws":statistics_by_draw
    }


def run(raw_path,mapping_path,agg_path,*,draws,depth_cap,seed):
    raw,observations,trait,focal=load_inputs(raw_path,mapping_path,agg_path)
    pairs=source_eligible_fixed_pairs(observations)
    if len(pairs)!=27 or len({i for i,_ in pairs})!=6:
        raise ValueError("Original exact 27 port-pair 6 independent-island cohort changed")
    turn,turn_summary=exact_set_turnover(raw,trait,focal,pairs)
    rare=island_temporal_rarefaction(raw,trait,focal,pairs,draws=draws,
                                   depth_cap=depth_cap,seed=seed)
    return {
        "schema":"chocho_N_trait_27pairs_six_Europe_islands_turnover_and_rarefaction_v1",
        "status":"DETECTION_EFFORT_SENSITIVITY_NOT_BIOLOGICAL_COLONIZATION",
        "raw_specieskey_source_sha256":hashlib.sha256(Path(raw_path).read_bytes()).hexdigest(),
        "GBIF_trait_mapping_sha256":hashlib.sha256(Path(mapping_path).read_bytes()).hexdigest(),
        "specieskey_facet_records_verified":len(raw),
        "comparison_port_pairs":len(pairs),
        "comparison_independent_islands":6,
        "original_6_island_observed_temporal_N_DiD":.1493596770601339,
        "facets_used_first_second_periods":"2001-2010_vs_2011-2024",
        "species_turnover":turn_summary,
        "rarefaction":rare,
        "turnover_site_rows":turn,
        "limitations":[
          "GBIF records need not represent independent plant individuals or identical collection protocols.",
          "Newly recorded plant species in later data are not confirmed introductions or colonists.",
          "The two observation intervals differ in length (10 versus 14 years).",
          "Only six European independent islands, with UK island contributing 18 of the 27 matched port pairs.",
          "Missing nutrient-scored taxa in draws are unknown, not absence.",
          "A positive bootstrap or signflip from the original comparison is exploratory and does not establish a causal nutrient gate."
        ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-facets",required=True,type=Path)
    ap.add_argument("--mapping",required=True,type=Path)
    ap.add_argument("--aggregates",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--draws",default=399,type=int)
    ap.add_argument("--depth-cap",default=120,type=int)
    ap.add_argument("--seed",default=312211,type=int)
    p=ap.parse_args()
    out=run(p.raw_facets,p.mapping,p.aggregates,
            draws=p.draws,depth_cap=p.depth_cap,seed=p.seed)
    p.output.parent.mkdir(parents=True,exist_ok=True)
    p.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("ALL_FLORA_N_TURNOVER_AND_RAREFACTION",json.dumps({
        "cohort":out["comparison_port_pairs"],
        "turnover":out["species_turnover"],
        "rarefaction":{k:v for k,v in out["rarefaction"].items() if k!="individual_draws"}
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
