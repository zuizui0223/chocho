#!/usr/bin/env python3
"""Audit whether static 2021 urban land cover is separable from port identity
in the source-locked 78 matched island GBIF butterfly-host N study.

This does not estimate shipping causal effect. Class 50 built area is a
2021 land-cover proxy (500m around site), whereas plants were observed
in GBIF 5km buffers pooled in 2001-2010 and 2011-2024. Different grain
and time must not be conflated.
"""
from __future__ import annotations
import argparse
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


PERIODS=("early_2001_2010","later_2011_2024")


def ranks(values):
    sorted_indices=sorted(range(len(values)),key=lambda i:values[i])
    ret=[None]*len(values);i=0
    while i<len(values):
        j=i+1
        while j<len(values) and values[sorted_indices[j]]==values[sorted_indices[i]]:
            j+=1
        rank=(i+j+1)/2
        for k in range(i,j):ret[sorted_indices[k]]=rank
        i=j
    return ret


def rho(x,y):
    if len(x)!=len(y) or len(x)<4:return None
    a,b=ranks(x),ranks(y)
    ma=statistics.mean(a);mb=statistics.mean(b)
    numerator=sum((u-ma)*(v-mb) for u,v in zip(a,b))
    den=math.sqrt(sum((u-ma)**2 for u in a)*sum((v-mb)**2 for v in b))
    return numerator/den if den else None


def aggregate_island(rows,metric):
    groups=defaultdict(list)
    for r in rows:
        if r.get(metric) is not None:
            groups[r["island_id"]].append(r)
    return [
        {"island_id":island,"country":v[0]["country"],"n_port_pairs":len(v),
         metric:statistics.mean(x[metric] for x in v),
         "built_gap":statistics.mean(x["built_gap"] for x in v),
         "crop_gap":statistics.mean(x["crop_gap"] for x in v)}
        for island,v in sorted(groups.items(),key=lambda t:int(t[0]))
    ]


def n_contrast(plant,site_period):
    p=plant.get(site_period)
    if p is None:
        raise ValueError("Source plant site missing "+str(site_period))
    if (p["plant_records"]<30 or p["N_scored_candidate_host_species"]<2 or
        p["N_scored_other_plant_species"]<5):
        return None
    return p["mean_N_candidate_hosts"]-p["mean_N_other_scored_plants"]


def classify_land_use(p,n):
    # Avoid treating "urban" as a social/economic classification.
    a,b=p["built_fraction_among_land_pixels"],n["built_fraction_among_land_pixels"]
    if a is None or b is None:return "UNCLASSIFIED"
    if a>=.20 and b>=.20:return "BOTH_BUILT20"
    if a>=.20 and b<.05:return "PORT_BUILT20_REMOTE_UNDER5"
    if a<.05 and b<.05:return "BOTH_UNDER5"
    return "MIXED_OR_MODERATE"


def run(urban,flora):
    if (urban["cohort"],urban["source_year"],urban["landcover_buffer_radius_m"])!=("78",2021,500):
        raise ValueError("2021 ESA urban study must be exact source 78 x 500m")
    if (flora["predeclared_source_pairs"],flora["site_period_facets"])!=(78,312):
        raise ValueError("Plant index source wrong cohort")
    site={(str(x["island_id"]),int(x["port_id"]),x["side"]):x for x in urban["site_buffers"]}
    plant={(str(x["island_id"]),int(x["port_id"]),x["side"],x["period"]):x for x in flora["observations"]}
    if len(site)!=156 or len(plant)!=312:raise ValueError("Duplicate site identities")
    pair=[]
    for x in urban["paired_sites"]:
        island=str(x["island_id"]);pid=int(x["port_id"])
        p=site[(island,pid,"port")];c=site[(island,pid,"nonport")]
        valid=p["quality"]=="COMPLETE" and c["quality"]=="COMPLETE"
        entry={
            "island_id":island,"port_id":pid,"country":x["country"],
            "urban_source_complete":valid,
            "landcover_group":classify_land_use(p,c),
            "port_built":p["built_fraction_among_land_pixels"] if valid else None,
            "nonport_built":c["built_fraction_among_land_pixels"] if valid else None,
            "built_gap":p["built_fraction_among_land_pixels"]-c["built_fraction_among_land_pixels"] if valid else None,
            "crop_gap":p["cropland_fraction_among_land_pixels"]-c["cropland_fraction_among_land_pixels"] if valid else None,
            "port_terrestrial_pixel_count":p["count_terrestrial_nonwater_pixels"],
            "nonport_terrestrial_pixel_count":c["count_terrestrial_nonwater_pixels"],
        }
        for period in PERIODS:
            a=n_contrast(plant,(island,pid,"port",period))
            b=n_contrast(plant,(island,pid,"nonport",period))
            entry[period+"_N_gap"] = a-b if valid and a is not None and b is not None else None
        if entry[PERIODS[0]+"_N_gap"] is not None and entry[PERIODS[1]+"_N_gap"] is not None:
            entry["same_pair_late_minus_early_N_gap"]=entry[PERIODS[1]+"_N_gap"]-entry[PERIODS[0]+"_N_gap"]
        else:
            entry["same_pair_late_minus_early_N_gap"]=None
        pair.append(entry)
    if len(pair)!=78:raise ValueError("Expected all 78 port/coast pairs")
    complete=[x for x in pair if x["urban_source_complete"]]
    return {
        "schema":"chocho_worldcover2021_urban_vs_port_N_common_support_audit_v1",
        "status":"STATIC_2021_URBAN_PROXY_CONFOUNDING_SCREEN_ONLY_NOT_PORT_EFFECT",
        "source_urban_valid_pairs":len(complete),"source_botanical_pairs":len(pair),
        "2021_urban_raster_is_post_2001_2010_period":True,
        "site_500m_not_identical_to_plant_5km":True,
        "common_support":{
           "nonport_2021_built_ge20pct":sum(x["nonport_built"] is not None and x["nonport_built"]>=.2 for x in complete),
           "nonport_2021_built_lt5pct":sum(x["nonport_built"] is not None and x["nonport_built"]<.05 for x in complete),
           "port_2021_built_ge20pct":sum(x["port_built"] is not None and x["port_built"]>=.2 for x in complete),
           "port_2021_built_lt5pct":sum(x["port_built"] is not None and x["port_built"]<.05 for x in complete),
           "both_urban_ge20pct_pairs":sum(x["landcover_group"]=="BOTH_BUILT20" for x in complete),
           "both_urban_ge20pct_islands":len({x["island_id"] for x in complete if x["landcover_group"]=="BOTH_BUILT20"}),
           "port_urban20_remote_rural5_pairs":sum(x["landcover_group"]=="PORT_BUILT20_REMOTE_UNDER5" for x in complete),
           "both_rural_under5_pairs":sum(x["landcover_group"]=="BOTH_UNDER5" for x in complete),
           "urban_control_gate":"NEEDS_10_BOTH_URBAN_PAIRS_IN_5_ISLANDS"
        },
        "results":summaries(pair),
        "pair_records":pair,
        "limitations":[
          "ESA WorldCover class 50 represents a satellite land-cover estimate of built-up in 2021; not population density, port shipping, or waste discharge.",
          "2021 built fraction cannot explain the 2001-2010 botanical record composition temporally.",
          "500m terrestrial landcover pixel composition is not identical to 5km GBIF site botanical radius.",
          "Ports and nearby urbanization are correlated; without urban nonport common support the independent port coefficient is not identified.",
          "GBIF observation, planted/introduced status, actual butterflies and nutrient exposure remain unresolved.",
          "Within-island correlation only descriptive if a few islands have >1 botanical eligible pair."
        ]
    }


def summaries(pairs):
    keys=[(PERIODS[0]+"_N_gap","early"),(PERIODS[1]+"_N_gap","late"),
          ("same_pair_late_minus_early_N_gap","same_pair_temporal")]
    out={}
    for key,label in keys:
        rows=[x for x in pairs if x[key] is not None]
        by_island=aggregate_island(rows,key)
        valid_geo=[v for v in by_island if v[key] is not None]
        out[label]={
            "eligible_pairs":len(rows),"independent_islands":len(by_island),
            "island_weighted_mean_N_gap":statistics.mean(v[key] for v in by_island) if by_island else None,
            "pair_level_spearman_N_gap_vs_port_minus_remote_built":rho(
                [r["built_gap"] for r in rows],[r[key] for r in rows]),
            "island_level_spearman_N_gap_vs_port_minus_remote_built":rho(
                [r["built_gap"] for r in by_island],[r[key] for r in by_island]),
            "pair_level_spearman_N_gap_vs_port_minus_remote_cropland":rho(
                [r["crop_gap"] for r in rows],[r[key] for r in rows]),
            "island_mean_details":valid_geo,
            "urban_both_ge20pct_eligible_pairs":sum(r["landcover_group"]=="BOTH_BUILT20" for r in rows),
            "landcover_port_built_nonport_rural_pairs":sum(r["landcover_group"]=="PORT_BUILT20_REMOTE_UNDER5" for r in rows)
        }
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--urban",type=Path,required=True)
    ap.add_argument("--allflora",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    urban=json.loads(args.urban.read_text(encoding="utf-8"))
    flora=json.loads(args.allflora.read_text(encoding="utf-8"))
    r=run(urban,flora)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("WORLDCOVER_2021_URBAN_N_COMMON_SUPPORT",json.dumps({
        "support":r["common_support"],
        "eligible":{k:{k2:v2 for k2,v2 in value.items() if k2!="island_mean_details"}
                    for k,value in r["results"].items()}
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
