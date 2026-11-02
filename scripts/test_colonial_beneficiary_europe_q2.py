#!/usr/bin/env python3
"""Q2-G exploratory: do historically colonized islands disproportionately add
potential resources for butterflies with Europe-centered native host geography?

Careful: 'Europe-centered native HOST RESOURCE' is not native butterfly range.
No maritime route, shipment, colonization event, host suitability or occupancy.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter
from pathlib import Path
import numpy as np
from test_colonial_duration_resource_q2 import within_permutation

def readcsv(path):
    with path.open(newline="",encoding="utf-8") as h:
        yield from csv.DictReader(h)

def main():
    p=argparse.ArgumentParser()
    for x in ("beneficiary-audit","beneficiary-islands","island-area","out"):
        p.add_argument("--"+x,required=True,type=Path)
    args=p.parse_args();args.out.parent.mkdir(parents=True,exist_ok=True)
    base=json.loads(args.beneficiary_audit.read_text(encoding="utf-8"))
    if base["population"]["island_added_resource_cells"]!=5674 or base["population"]["eligible_butterflies"]!=239:
        raise RuntimeError("Q2 source-derived beneficiary population changed")
    a={r["code"]:r for r in readcsv(args.island_area)}
    data=[]
    for r in readcsv(args.beneficiary_islands):
        code=r["island_code"]
        if code not in a:raise RuntimeError("Missing official island area "+code)
        if float(r["added_resource_butterflies"])<1:raise RuntimeError("Zero denominator")
        data.append({
           "code":code,"level1":r["level1"],
           "area_km2":float(a[code]["area_km2"]),
           "coverage":float(r["lenzner_polygon_coverage"]),
           "maxyears":float(r["max_four_european_empire_years"]),
           "log_legacy":math.log1p(max(0,float(r["max_four_european_empire_years"]))),
           "dominant_europe_share":float(r["europe_dominant_resource_fraction"]),
           "continuous_europe_share":float(r["mean_share_native_resources_in_europe"]),
           "n_cells":int(r["added_resource_butterflies"])
        })
    if len(data)!=122 or sum(r["n_cells"] for r in data)!=5674:
        raise RuntimeError("Region sum does not recover fixed 5674 resource cells")
    variants={
      "primary_cover80":lambda r:r["coverage"]>=.8,
      "non_europe_cover80":lambda r:r["coverage"]>=.8 and r["level1"]!="1",
      "cover50":lambda r:r["coverage"]>=.5,
      "positive_50year_cover80":lambda r:r["coverage"]>=.8 and r["maxyears"]>=50,
    }
    output={}
    for label,condition in variants.items():
        rows=[x for x in data if condition(x)]
        if len(rows)<30:
            output[label]={"n":len(rows),"status":"SMALL_COHORT"}
            continue
        median_area=np.median([r["area_km2"] for r in rows])
        by_area=[r["level1"]+"|"+("small" if r["area_km2"]<=median_area else "large")
                 for r in rows]
        tests={}
        for outcome in ("dominant_europe_share","continuous_europe_share"):
            y=[{**r,"share":r[outcome]} for r in rows]
            tests[outcome]={
              "observed_mean":float(np.mean([x["share"] for x in y])),
              "by_level1":within_permutation(y,[r["level1"] for r in y],draws=4999),
              "by_level1_and_area_half":within_permutation(y,by_area,draws=4999),
            }
        output[label]={
          "n":len(rows),
          "sum_new_resource_cells":sum(x["n_cells"] for x in rows),
          "known_four_empire_50yr":sum(r["maxyears"]>=50 for r in rows),
          "tests":tests
        }
    result={
      "schema":"chocho_colonial_europe_native_resource_beneficiary_v0.1",
      "status":"EXPLORATORY_NONCAUSAL_NARROW_EUROPE_SOURCE_CONTRAST",
      "historical_exposure":"log1p(max 4 Lenzner former empire area-weighted years)",
      "responses":[
         "proportion of butterflies newly acquiring island resource opportunity whose WCVP/HOSTS native host-resource footprint is Europe-dominant",
         "mean Europe fraction of butterflies' native host-resource footprints, across new butterfly-island opportunity cells"
      ],
      "adjustment":"OLS level1 fixed effects + ln WGSRPD3 region area; permutation of historical exposure within level1 and level1×area-half",
      "weights":"Each WGSRPD3 island equally weighted; opportunity counts recorded separately",
      "results":output,
      "caution":[
        "European native host-resource footprint is not the butterfly's ancestral or original geographic range.",
        "Post-hoc exploratory extension after earlier colonial association screens; multiple-testing limitations apply.",
        "Historic 4-empire mapping is incomplete, particularly for islands with smaller polygons and other colonial powers.",
        "Positive source correlation, if any, is not transport pathway, plant arrival, or butterflies successfully colonizing."
      ]
    }
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
