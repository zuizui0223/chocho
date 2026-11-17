#!/usr/bin/env python3
"""Quality gate for actual 2005 NOAA old OSCAR globally sampled coastal starts.

Reuses immutable 3-date NOAA source-receipt WITHOUT re-running currents, and
assesses how many islands remain when permitted distance between an original
GSHHG ocean-shore geometry terminal and the NOAA velocity sampling point is
restricted from 200km to a more defensible 25/50/80/100/150km.
No observed seed launch nor directed transport inferred.
"""
from __future__ import annotations
import argparse,json,statistics
from collections import Counter,defaultdict
from pathlib import Path

THRESHOLDS=(25,50,80,100,150,200)
DATES=3
FULL_ORIGINAL_ISLANDS=126
GSHHG_ELIGIBLE=104

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    data=json.loads(a.source.read_text())
    assert data["status"]=="GLOBAL_NOAA_OLD_OSCAR_ACTUAL_UV_SOURCE_COVERAGE_NOT_126_ISLAND_PARTICLE_TRANSPORT"
    assert data["NOAA_not_NASA_daily_v2"] and data["source_GSHHG_shore_graph_valid_islands"]==GSHHG_ELIGIBLE
    rows=data["126_source_region_release_eligibility"]
    if len(rows)!=FULL_ORIGINAL_ISLANDS or len({x["region"] for x in rows})!=FULL_ORIGINAL_ISLANDS:
        raise RuntimeError("Original WCVP island sample changed")
    if not (sum(r["fixed_season_validities"]==3 for r in rows)==93 and
            sum(r["fixed_season_validities"]>=1 for r in rows)==94):
        raise RuntimeError("Original 2005 NOAA global actual u/v validity changed")
    by_threshold={}
    per_island=[]
    for r in rows:
        dist=r.get("nearest_valid_offshore_source_km_per_season") or [None]*3
        if len(dist)!=DATES:raise RuntimeError("NOAA release-date source original 3 times missing")
        status=r["noaa_spring_summer_autumn_launch_status"]
        is_coastal=status!="ORIGINAL_GSHHG_SHORE_ANCHOR_NOT_AVAILABLE"
        per_island.append({
         "WGSRPD3":r["region"],
         "GSHHG_original_source_shore_anchor_available":is_coastal,
         "noaa_quality_status_at_200km":status,
         "source_3days_nearest_valid_NOAA_1deg_start_offset_km":dist,
         "offshore_offset_min_of_available_km":min((float(x) for x in dist if x is not None),default=None),
         "offshore_offset_max_of_available_km":max((float(x) for x in dist if x is not None),default=None),
         "qualifies_3dates_at50km":all(x is not None and float(x)<=50 for x in dist),
         "qualifies_any_date_at50km":any(x is not None and float(x)<=50 for x in dist),
        })
    byregion={r["WGSRPD3"]:r for r in per_island}
    for km in THRESHOLDS:
        allthree=[
          c for c in per_island if all(
            d is not None and float(d)<=km
            for d in c["source_3days_nearest_valid_NOAA_1deg_start_offset_km"])
        ]
        anyday=[
          c for c in per_island if any(
            d is not None and float(d)<=km
            for d in c["source_3days_nearest_valid_NOAA_1deg_start_offset_km"])
        ]
        n=len(allthree)
        by_threshold[str(km)]={
          "max_allowed_distance_from_original_GSHHG_shore_km":km,
          "source_regions_all_three_dates_valid":n,
          "source_regions_any_one_date_valid":len(anyday),
          "potential_unordered_island_pairs_both_3date_source_eligible":n*(n-1)//2,
          "potential_unordered_island_pairs_both_any_date_source_eligible":len(anyday)*(len(anyday)-1)//2,
          "eligible_WGSRPD3_region_codes":sorted(x["WGSRPD3"] for x in allthree),
          "actual_real_particle_trajectories_calculated":False
        }
    assert (by_threshold["200"]["source_regions_all_three_dates_valid"],
            by_threshold["200"]["source_regions_any_one_date_valid"])==(93,94)
    assert all(by_threshold[str(x)]["source_regions_all_three_dates_valid"]<=
               by_threshold[str(y)]["source_regions_all_three_dates_valid"]
               for x,y in zip(THRESHOLDS,THRESHOLDS[1:]))
    elig=[r for r in per_island if r["noaa_quality_status_at_200km"]=="ALL_THREE_SOURCE_DATES_VALID"]
    original_min=[r["offshore_offset_min_of_available_km"] for r in elig]
    original_max=[r["offshore_offset_max_of_available_km"] for r in elig]
    report={
      "schema":"chocho_original_NOAA2005_126island_global_one_degree_coastal_start_offset_quality_v1",
      "status":"ORIGINAL_REAL_NOAA_UV_COVERAGE_ONLY_WITH_OFFSHORE_DISTANCE_SENSITIVITY_NOT_OCEAN_DISPERSAL",
      "original_source_CI":44011612861,
      "original_source_release_times":data["noaa_actual_source_release_dates"],
      "NOAA_native_original_dataset":data["NOAA_dataset"],
      "NOAA_grid_source":"1deg subsample of original 1/3degree legacy 2005 OSCAR actual u/v",
      "GSHHG_source_SHA256":data["source_GSHHG_graph_SHA256"],
      "source_global_u_v_original_NetCDF_subset_SHA256":data["NOAA_source_NetCDF_original_subset_sha256"],
      "original_source_verified_total_126_islands":len(rows),
      "original_GSHHG_source_resolvable_104_islands":GSHHG_ELIGIBLE,
      "original_NOAA2005_200km_three_date_complete_islands":len(elig),
      "original_NOAA2005_start_offshore_offset_3date_min_median_km":statistics.median(original_min),
      "original_NOAA2005_start_offshore_offset_3date_max_median_km":statistics.median(original_max),
      "distance_threshold_source_quality":by_threshold,
      "per_original_region":per_island,
      "limits":[
        "Original NOAA source grid was subsampled to ~1 degree and 15m water currents; no 1/3-degree coastal trajectories are yet available.",
        "A valid NOAA sea-current sample located up to 200km from its frozen GSHHG coast source is not a biological seed departure location.",
        "Distances are source availability constraints; no historical propagule arrival or viability and no downstream recipient coastline detection were measured.",
        "The number of pair combinations is combinatorial capacity, not a verified 126×126 directional transition matrix or any successful sea transport event.",
        "Do not construct a claimed directional particle-access matrix from these three independent source dates without continuous historical currents and land-contact censoring."
      ]
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GLOBAL_NOAA126_SHORE_OFFSHORE_OFFSET_SENSITIVITY",json.dumps({
      "status":report["status"],"median_best_start_km":report["original_NOAA2005_start_offshore_offset_3date_min_median_km"],
      "median_worst_start_km":report["original_NOAA2005_start_offshore_offset_3date_max_median_km"],
      "thresholds":{k:{q:v[q] for q in ("source_regions_all_three_dates_valid","source_regions_any_one_date_valid","potential_unordered_island_pairs_both_3date_source_eligible")} for k,v in by_threshold.items()}
    }),flush=True)

if __name__=="__main__":
    main()
