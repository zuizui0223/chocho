#!/usr/bin/env python3
"""Run only the five not-yet-recovered original GBIF HTTP400 historical
plant × island pairs. Reuse the validated split-polygon query method,
with ≤20 respectful API requests and no automatic imputation.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from recover_q4_gbif_multipolygon_years import (
    records as readcsv,savecsv,geo_parts,run_case,MAX_CALLS
)
from audit_q4_gbif_plant_first_year_pilot import build_geometries

def main():
    p=argparse.ArgumentParser()
    for name in ("retry-source","previous-audit","wgsrpd3-geojson","out-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    old=json.loads(a.previous_audit.read_text(encoding="utf-8"))
    if not old["positive_control_match"] or old["fully_recovered_original_400_pairs"]!=2:
        raise RuntimeError("Positive control and first two complete-region recoveries not source-verified")
    source=readcsv(a.retry_source)
    failed=[r for r in source if r["status"]=="HTTP_400"]
    if len(source)!=12 or len(failed)!=7:
        raise RuntimeError("Original source failure status drift")
    previous={(r["plant_id"],r["region"]) for r in old["per_case"]
              if r["kind"]=="RESCUE" and r["complete_full_WGSRPD_region"]}
    if previous!={("2559676","MOL"),("2394644","TUA")}:
        raise RuntimeError("Unexpected previously recovered source pairs")
    todo=[r for r in failed if (r["plant_id"],r["island_code"]) not in previous]
    if len(todo)!=5:raise RuntimeError("Expected five pending source HTTP400 cases")
    geo=build_geometries(a.wgsrpd3_geojson,{r["island_code"] for r in todo})
    todo.sort(key=lambda r:(len(geo_parts(r["island_code"],geo)),r["island_code"],r["plant_id"]))
    shared={"calls":0,"consecutive429":0,"stopped":False}
    outcomes=[];chunks=[]
    for row in todo:
        if shared["stopped"]:break
        result,details=run_case(row,geo,shared,"RESCUE_REMAINDER")
        outcomes.append(result)
        chunks+=details
    savecsv(a.out_dir/"remaining_five_split_wgsrpd_years.csv",outcomes)
    savecsv(a.out_dir/"remaining_polygon_query_statuses.csv",chunks)
    complete=[r for r in outcomes if r["complete_full_WGSRPD_region"]]
    audit={
      "schema":"chocho_GBIF_Q4_split_polygon_remaining_five_v0.1",
      "status":"COMPLETE_COHORT_PARTITION" if len(complete)==5 else "INCOMPLETE_REQUEST_OR_RATE_LIMIT",
      "prior_source_audit":str(a.previous_audit),
      "previous_full_recoveries":2,
      "remaining_original_http400_pairs":len(todo),
      "attempted":len(outcomes),"fully_recovered":len(complete),
      "new_dated_pairs":sum(r["status"]=="DATED_OCCURRENCE_FOUND" for r in complete),
      "new_no_indexed_record_pairs":sum(r["status"]=="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT" for r in complete),
      "GBIF_API_requests":shared["calls"],"request_limit":MAX_CALLS,
      "stopped_after_two_429":shared["stopped"],
      "per_case":outcomes,
      "limits":[
        "The GBIF first indexed plant year is a historical documentation date, NOT actual introduction time.",
        "Complete source WKT partition is required before setting a first year or confirming no GBIF record.",
        "Record counts might duplicate observations exactly on touching polygon boundaries.",
        "These five are from the original source error cohort and must not be treated as an independent new species sample.",
        "No butterfly settlement, flight, shipping arrival year or climatic threshold tested."
      ],
    }
    (a.out_dir/"remaining_five_polygon_recovery_audit.json").write_text(
        json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(audit,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
