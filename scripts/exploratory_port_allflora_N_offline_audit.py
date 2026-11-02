#!/usr/bin/env python3
"""Independent offline sensitivity audit for frozen matched WPI78 all-plant N data.

All estimates are exploratory GBIF-index composition contrasts, not nutrient
loading, local butterfly host use, growth, or plant colonization.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


def check_source(data: dict):
    assert data["predeclared_source_pairs"]==78
    assert data["independent_islands"]==35
    assert data["site_period_facets"]==312
    assert len(data["observations"])==312
    assert all(v==0 for v in data["drift_counts"].values())
    keys=[(str(x["island_id"]),int(x["port_id"]),x["side"],x["period"]) for x in data["observations"]]
    if len(keys)!=len(set(keys)):raise ValueError("Duplicated site period")
    for item in data["observations"]:
        if item["N_scored_plant_species"]!=(
                item["N_scored_candidate_host_species"]+
                item["N_scored_other_plant_species"]):
            raise ValueError("Incomplete partition of N-scored plant species")
        if item["plant_records"]<item["all_plant_species_keys"]:
            # GBIF sometimes may index imperfect species facets, but not
            # more distinct named species than records in an accepted query.
            raise ValueError("Too many species keys relative to GBIF records")


def paired(data, *, min_records=30,min_target=2,min_other=5,min_coverage=0.0,
           min_effort_ratio=0.0):
    sites=defaultdict(dict)
    for r in data["observations"]:
        k=(str(r["island_id"]),int(r["port_id"]),r["period"])
        if r["side"] in sites[k]:raise ValueError("Duplicate site")
        sites[k][r["side"]]=r
    out={}
    for k,site in sites.items():
        if set(site)!={"port","nonport"}:raise ValueError("Missing side")
        p,c=site["port"],site["nonport"]
        mx=max(p["plant_records"],c["plant_records"])
        effort=min(p["plant_records"],c["plant_records"])/mx if mx else 0
        good=all(x["plant_records"]>=min_records and
                 x["N_scored_candidate_host_species"]>=min_target and
                 x["N_scored_other_plant_species"]>=min_other and
                 (x["N_scored_fraction_of_all_indexed_plant_species"] or 0)>=min_coverage
                 for x in (p,c)) and effort>=min_effort_ratio
        if not good:continue
        target=p["mean_N_candidate_hosts"]-c["mean_N_candidate_hosts"]
        background=p["mean_N_other_scored_plants"]-c["mean_N_other_scored_plants"]
        out[k]={"island_id":k[0],"port_id":k[1],"period":k[2],
                "target":target,"background":background,
                "target_minus_background":target-background,
                "port_records":p["plant_records"],
                "nonport_records":c["plant_records"],
                "recording_effort_min_over_max":effort,
                "N_coverage_min":min(p["N_scored_fraction_of_all_indexed_plant_species"],
                                     c["N_scored_fraction_of_all_indexed_plant_species"])}
    return out


def island_means(rows,key):
    by=defaultdict(list)
    for r in rows:
        by[r["island_id"]].append(r[key])
    return {k:statistics.mean(v) for k,v in by.items()}


def interval(values,seed=63872,draws=4999):
    if not values:return None
    rng=random.Random(seed);n=len(values)
    samples=sorted(statistics.mean(rng.choices(values,k=n)) for _ in range(draws))
    return [samples[round((len(samples)-1)*.025)],samples[round((len(samples)-1)*.975)]]


def signflip_islands(values):
    """Exact two-sided sign-flip enumeration, descriptive only, not a confirmatory p."""
    n=len(values)
    if not n or n>18:return None
    observed=abs(statistics.mean(values))
    hits=sum(abs(sum(s*x for s,x in zip(signs,values))/n)>=observed-1e-12
             for signs in itertools.product((-1,1),repeat=n))
    return hits/(2**n)


def stats(rows,key):
    group=island_means(rows,key)
    vals=list(group.values())
    if not vals:return {"pairs":0,"islands":0,"island_mean":None,"interval95":None}
    overall=statistics.mean(vals)
    loo=[statistics.mean([v for k,v in group.items() if k!=island])
         for island in group] if len(group)>1 else []
    return {"pairs":len(rows),"islands":len(group),
            "island_mean":overall,
            "island_median":statistics.median(vals),
            "interval95":interval(vals),
            "island_leave_one_out_mean_range":[min(loo),max(loo)] if loo else None,
            "positive_islands":sum(x>0 for x in vals),
            "negative_islands":sum(x<0 for x in vals),
            "two_sided_exact_island_signflip_descriptive":signflip_islands(vals),
            "island_details":[{"island_id":k,"mean":v,
                               "n_site_pairs":sum(r["island_id"]==k for r in rows)}
                              for k,v in sorted(group.items(),key=lambda q:int(q[0]))]}


def temporal_same_pairs(records):
    early={k[:2]:v for k,v in records.items() if k[2]=="early_2001_2010"}
    late={k[:2]:v for k,v in records.items() if k[2]=="later_2011_2024"}
    shared=sorted(set(early)&set(late),key=lambda v:(int(v[0]),v[1]))
    rows=[]
    for k in shared:
        row={"island_id":k[0],"port_id":k[1]}
        for key in ("target","background","target_minus_background"):
            row[key]=late[k][key]-early[k][key]
        rows.append(row)
    return rows


def evaluate(data):
    check_source(data)
    settings=[
        {"name":"base","min_records":30,"min_target":2,"min_other":5},
        {"name":"max_5fold_site_effort_imbalance","min_records":30,"min_target":2,"min_other":5,"min_effort_ratio":.2},
        {"name":"min_10_percent_N_trait_coverage_each_side","min_records":30,"min_target":2,"min_other":5,"min_coverage":.10},
        {"name":"min_20_percent_N_trait_coverage_each_side","min_records":30,"min_target":2,"min_other":5,"min_coverage":.20},
        {"name":"min_100_all_plant_records_each_side","min_records":100,"min_target":2,"min_other":5},
        {"name":"min_3_hosts_10_other_scored_plants","min_records":30,"min_target":3,"min_other":10}
    ]
    outputs={}
    for conf in settings:
        name=conf["name"]
        pair=paired(data,**{k:v for k,v in conf.items() if k!="name"})
        years={}
        for per in ("early_2001_2010","later_2011_2024"):
            elig=[v for k,v in pair.items() if k[2]==per]
            years[per]={k:stats(elig,k) for k in ("target","background","target_minus_background")}
        same=temporal_same_pairs(pair)
        years["same_exact_port_pairs_period_change"]={
            k:stats(same,k) for k in ("target","background","target_minus_background")
        }
        outputs[name]=years
    base=outputs["base"]
    assert base["early_2001_2010"]["target_minus_background"]["pairs"]==27
    assert base["later_2011_2024"]["target_minus_background"]["pairs"]==43
    assert base["same_exact_port_pairs_period_change"]["target_minus_background"]["pairs"]==27
    assert abs(base["same_exact_port_pairs_period_change"]["target_minus_background"]["island_mean"]-.1493596770601339)<1e-9
    return {
       "schema":"chocho_WPI78_N_target_vs_other_temporal_island_exact_audit_v1",
       "status":"EXPLORATORY_NO_BIOLOGICAL_CAUSAL_INFERENCE",
       "source_GBIF_312_site_period_facet_receipt":"https://github.com/zuizui0223/chocho/actions/runs/41407370295",
       "source_original_78_geom_no_overlap":True,
       "source_aggregate_response_unchanged_under_GBIF_reindexing":True,
       "results":outputs,
       "limitations":[
          "Six eligible independent islands in the identical pair temporal contrast, not 27 independent experimental replicates.",
          "A positive descriptive island bootstrap interval and sign-flip statistic do not prove a port mechanism after multiple ecological exploration.",
          "Background 'other' N-scored flora may include many hosts of other butterfly species.",
          "German Ellenberg N is an affinity score not soil or wastewater N measured at paired sites.",
          "Later/earlier are GBIF index windows and change in recorded taxa, not direct biological colonization.",
          "Different records, collection sources, and observation effort can change apparent site composition.",
          "Focus group includes globally catalogued non-native hosts, not historically introduced hosts at those islands.",
          "Low trait coverage may be selective and vary geographically."
       ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    data=json.loads(args.source.read_text(encoding="utf-8"))
    result=evaluate(data)
    result["source_sha256"]=hashlib.sha256(args.source.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    for setting,o in result["results"].items():
        print("ALL_FLORA_N_ISLAND_SENSITIVITY",setting,json.dumps({
            k:{m:{"pairs":d["pairs"],"islands":d["islands"],
                 "mean":d["island_mean"],"CI":d["interval95"],
                 "loo":d.get("island_leave_one_out_mean_range"),
                 "signflip":d.get("two_sided_exact_island_signflip_descriptive"),
                 "positive_islands":d.get("positive_islands")}
                 for m,d in v.items()}
            for k,v in o.items()
        },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
