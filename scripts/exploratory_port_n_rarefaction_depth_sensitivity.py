#!/usr/bin/env python3
"""Record-vs-species rarefaction sensitivity of six-island plant N time contrast.

All inputs are frozen GH Actions provenance. Record counts are GBIF catalogue
intensity, NOT independent individual plants. No fresh GBIF calls.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from pathlib import Path
from collections import defaultdict

SRC=Path(__file__).resolve().parent
sys.path.insert(0,str(SRC))
from exploratory_port_n_turnover_rarefaction import (
    PERIODS,SIDES,load_inputs,source_eligible_fixed_pairs,
    sampled_keys_from_records,mean_group_N
)


def cell_support_for_pair(raw,pair,unit,cap):
    island,port=pair
    stats=[]
    for side in SIDES:
        for period in PERIODS:
            f=raw[(island,port,side,period)]["faceted_specieskey_counts"]
            records=sum(int(v) for v in f.values())
            species=sum(int(v)>0 for v in f.values())
            stats.append(records if unit=="record" else species)
    return min(int(cap),*stats)


def simulate_mode(raw,trait,focal,pairs,*,unit,cap,n_draws,seed):
    assert unit in ("record","species")
    rng=random.Random(seed)
    pair_depth={pair:cell_support_for_pair(raw,pair,unit,cap) for pair in pairs}
    tested=[(pair,d) for pair,d in pair_depth.items() if d>=10]
    source_spec={}
    for pair,_ in tested:
        island,port=pair
        for side in SIDES:
            for period in PERIODS:
                species=raw[(island,port,side,period)]["faceted_specieskey_counts"]
                positive={int(k):int(v) for k,v in species.items() if int(v)>0}
                source_spec[(island,port,side,period)]=positive
    stats=[]
    for _ in range(n_draws):
        island_effects=defaultdict(list)
        for (island,port),depth in tested:
            contrast={}
            usable=True
            for side in SIDES:
                for period in PERIODS:
                    site=source_spec[(island,port,side,period)]
                    if unit=="record":
                        observed=sampled_keys_from_records(site,rng,depth)
                    else:
                        observed=set(rng.sample(sorted(site),depth))
                    score=mean_group_N(observed,trait,focal)
                    if score is None:
                        usable=False;break
                    contrast[(side,period)]=score
                if not usable:break
            if usable:
                change=((contrast[("port",PERIODS[1])]-contrast[("nonport",PERIODS[1])])-
                        (contrast[("port",PERIODS[0])]-contrast[("nonport",PERIODS[0])]))
                island_effects[island].append(change)
        means=[statistics.mean(v) for v in island_effects.values()]
        stats.append({"eligible_pairs":sum(len(v) for v in island_effects.values()),
                      "eligible_islands":len(means),
                      "mean_island_delta":statistics.mean(means) if means else None,
                      "every_island_positive":len(means)==6 and all(v>0 for v in means)})
    total=len(stats)
    full=[s for s in stats if s["eligible_islands"]==6 and s["mean_island_delta"] is not None]
    any_ok=[s for s in stats if s["mean_island_delta"] is not None]
    median=lambda values:statistics.median(values) if values else None
    quant=lambda values,p:sorted(values)[round((len(values)-1)*p)] if values else None
    full_effects=[s["mean_island_delta"] for s in full]
    return {
        "unit":unit,"per_cell_cap":cap,"draws":total,"seed":seed,
        "original_pairs":len(pairs),"pairs_eligible_for_min_10_sampling":len(tested),
        "depth_median":median(list(pair_depth.values())),
        "depth_min":min(pair_depth.values()),"depth_max":max(pair_depth.values()),
        "median_islands_retained_per_draw":median([s["eligible_islands"] for s in stats]),
        "median_pairs_retained_per_draw":median([s["eligible_pairs"] for s in stats]),
        "full_six_islands_probability":len(full)/total,
        "mean_if_all_six_islands_supported":statistics.mean(full_effects) if full_effects else None,
        "fraction_positive_if_all_six_supported":sum(v>0 for v in full_effects)/len(full_effects) if full_effects else None,
        "fraction_every_island_positive_if_all_six_supported":sum(s["every_island_positive"] for s in full)/len(full) if full else None,
        "empirical_sampling_only_95_if_all_six_supported":[quant(full_effects,.025),quant(full_effects,.975)] if full else None,
        "fraction_positive_if_any_island_supported":sum(s["mean_island_delta"]>0 for s in any_ok)/len(any_ok) if any_ok else None,
        "conditional_on_trait_coverage_and_rarity":True,
        "interpretation":"Sensitivity to GBIF record counts vs uniformly indexed distinct plant species; not abundance/field-population evidence."
    }


def run(raw_path,map_path,aggregate_path,*,draws,seed):
    raw,obs,traits,focal=load_inputs(raw_path,map_path,aggregate_path)
    pairs=source_eligible_fixed_pairs(obs)
    if len(pairs)!=27 or len({i for i,_ in pairs})!=6:
        raise ValueError("Fixed eligible source cohort is not 27 ports / 6 islands")
    results=[]
    settings=[("record",30),("record",60),("record",120),("record",250),
              ("record",500),("record",1000),
              ("species",20),("species",40),("species",80),("species",160)]
    for i,(unit,cap) in enumerate(settings):
        z=simulate_mode(raw,traits,focal,pairs,unit=unit,cap=cap,n_draws=draws,seed=seed+101*i)
        results.append(z)
        print("GBIF_RECORD_VS_SPECIES_RAREFACTION",json.dumps(z,sort_keys=True),flush=True)
    return {
        "schema":"chocho_N_port_27pair_europe_record_species_depth_sensitivity_v1",
        "status":"DETECTION_EFFORT_ESTIMAND_SENSITIVITY_NOT_TRUE_PLANT_ABUNDANCE",
        "input_312_facets_verified":len(raw),
        "focal_exactly_matched_pairs":len(pairs),
        "island_count":6,
        "original_species_presence_weighted_signal":.1493596770601339,
        "settings":results,
        "limitations":[
          "GBIF duplicate and field-reporting procedures make occurrence records neither independent plants nor unbiased sampling intensity.",
          "Record-weighted and species-uniform rarefaction answer different questions and are not estimators of plant biomass.",
          "The first and later record windows have different lengths (10 vs 14 years).",
          "Only six European islands and candidate host identities defined by a global catalog rather than island breeding data.",
          "Missing eligible trait species in a simulated cell are excluded, not interpreted as biological absence.",
          "All ten combinations of sampling-unit and depth are exploratory post-hoc sensitivity diagnostics."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw",type=Path,required=True)
    p.add_argument("--map",type=Path,required=True)
    p.add_argument("--aggregate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--draws",type=int,default=299)
    p.add_argument("--seed",type=int,default=549311)
    a=p.parse_args()
    z=run(a.raw,a.map,a.aggregate,draws=a.draws,seed=a.seed)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RECORD_VS_SPECIES_DESIGN_COMPLETE",len(z["settings"]),flush=True)


if __name__=="__main__":
    main()
