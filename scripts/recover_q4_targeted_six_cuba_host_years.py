#!/usr/bin/env python3
"""Recover the six source-locked historical botanical GBIF HTTP400 failures
in a 47-literature-case targeted introduced-host first-record audit. All six
were Cuba WGSRPD3 polygons, not six different islands. Reuse separately
validated repeatable POLYGON year-facet queries and full-region source WKT.
This is case-conditioned documentary chronology, never proof of plant arrival.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from recover_q4_gbif_multipolygon_years import run_case,geo_parts,savecsv,MAX_CALLS
from audit_q4_gbif_plant_first_year_pilot import build_geometries

def read(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser()
    for n in ("targeted-source","split-control-audit","geojson","out-dir"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    audit=json.loads(a.split_control_audit.read_text(encoding="utf8"))
    if audit.get("positive_control_match") is not True:
        raise RuntimeError("Source reference count/year for GBIF split POLYGON was not previously validated")
    source=read(a.targeted_source)
    if len(source)!=84:raise RuntimeError("Original 84 case-selected plant-island first-year sample changed")
    missing=[r for r in source if r["status"]=="GBIF_QUERY_FAILED"]
    if len(missing)!=6 or {r["wgsrpd3_code"] for r in missing}!={"CUB"} or any(
        "HTTP Error 400" not in r.get("query_failure","") for r in missing):
        raise RuntimeError("Unexpected source missingness; only six rejected Cuban queries authorized")
    transformed=[]
    for r in missing:
        if not r["gbif_species_key"]:raise RuntimeError("GBIF accepted taxon missing")
        transformed.append({
            "plant_id":r["plant_id"],"plant_name":r["wcvp_plant_name"],
            "island_code":r["wgsrpd3_code"],"gbif_species_key":r["gbif_species_key"],
            "status":r["status"],"original_status":r["status"]})
    geo=build_geometries(a.geojson,{"CUB"})
    if len(geo_parts("CUB",geo))!=26:raise RuntimeError("Original Cuban region geometry differs")
    shared={"calls":0,"consecutive429":0,"stopped":False}
    rows=[];chunks=[]
    for case in sorted(transformed,key=lambda r:r["plant_id"]):
        if shared["stopped"]:break
        result,querylog=run_case(case,geo,shared,"LITERATURE_SELECTED_CUBAN_SOURCE")
        rows.append(result);chunks+=querylog
    savecsv(a.out_dir/"six_cuba_targeted_plant_year_results.csv",rows)
    savecsv(a.out_dir/"six_cuba_polygon_query_chunks.csv",chunks)
    completed=[r for r in rows if r["complete_full_WGSRPD_region"]]
    summary={
        "schema":"chocho_q4_literature_selected_six_cuba_gbif_year_recovery_v0.1",
        "status":"CASE_CONDITIONED_SOURCE_DATE_RETRIEVAL_NOT_CAUSAL_COLONIZATION",
        "original_literature_selected_plant_island_pairs":84,
        "original_cuban_http400_queries":6,
        "positive_control_11_part_MRQ_6records_2004_previously_verified":True,
        "original_cuba_polygon_components":26,
        "attempted":len(rows),"fully_covered":len(completed),
        "new_dated_pairs":sum(r["status"]=="DATED_OCCURRENCE_FOUND" for r in completed),
        "new_no_indexed_records":sum(r["status"]=="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT" for r in completed),
        "new_gbif_queries":shared["calls"],"max_queries":MAX_CALLS,
        "GBIF_429_halt":shared["stopped"],
        "per_plant":rows,
        "limits":[
           "All six are targeted on published island butterfly observations, so they cannot be pooled with the response-blind 500-pair pilot as confirmatory sampling.",
           "An earliest GBIF plant year is the earliest indexed botanical record, not introduction or establishment.",
           "The Cuban botanical regional polygon is a macroregion. Historical larval feeding and butterfly first arrival are not independently dated.",
           "No butterfly risk-set, 3-gate causal ordering or wing-versus-ship mechanism is inferred.",
           "Any partial or failed polygon query remains unknown; neither absence nor earliest source year is imputed.",
           "GBIF 429 rate limiting halts further queries after two consecutive responses."
        ]
    }
    (a.out_dir/"six_cuba_targeted_host_date_recovery_audit.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
