#!/usr/bin/env python3
"""Query first GBIF-indexed years of 3 response-blind Canary plants after CNY island correction.

This is an additive pilot, NOT a rerun or re-selection of the original 500.
It reuses the original GBIF year facet and polygon validation logic. A plant
record year is not an introduction year; a missing record is not true absence.
"""
import argparse,csv,json,datetime
from pathlib import Path
from audit_q4_gbif_plant_first_year_pilot import (build_geometries,key_from_gbif,facet_year,
                                                   STATUS_NO_EVIDENCE)
def rows(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def main():
    p=argparse.ArgumentParser()
    for k in ("CNY-three-pairs","WCVP-three-names","wgsrpd-geojson","q4-v02-source-audit","out-dir"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    data=rows(a.CNY_three_pairs)
    names=rows(a.WCVP_three_names)
    meta=json.loads(a.q4_v02_source_audit.read_text(encoding="utf8"))
    if len(data)!=3 or len(names)!=3 or any(r["island_code"]!="CNY" for r in data):
        raise RuntimeError("Not exactly three Canary plant pairs")
    ids={r["plant_id"] for r in data}
    source={r["plant_id"]:r for r in names}
    if len(ids)!=3 or set(source)!=ids or meta["fixed_pilot"]["new_response_blind_CNY_host_pairs"]!=3:
        raise RuntimeError("Fixed 503 sample geographic source changed")
    geom=build_geometries(a.wgsrpd_geojson,{"CNY"})["CNY"]
    report=[]
    statuses={}
    for i,pid in enumerate(sorted(ids),1):
        r=source[pid]
        if r["gbif_species_query_eligible"]!="True":
            status={"status":"NOT_WCVP_ACCEPTED_SPECIES_RANK","gbif_record_count":"",
                    "earliest_indexed_year":"","detail":"Taxonomic rank excluded"}
        elif geom["status"]!="GEOMETRY_READY":
            status={"status":"GEOMETRY_NOT_QUERYABLE","gbif_record_count":"",
                    "earliest_indexed_year":"","detail":geom["status"]}
        else:
            match=key_from_gbif(r["wcvp_accepted_name"])
            if match["status"]!="EXACT_SPECIES_MATCH":
                status={"status":match["status"],"gbif_record_count":"",
                        "earliest_indexed_year":"","detail":match.get("detail","")}
            else:
                try:
                    response=facet_year(match["key"],geom["wkt"])
                except Exception as e:
                    response={"status":"GBIF_QUERY_FAILED","gbif_record_count":"",
                              "earliest_indexed_year":"","detail":str(e)}
                status={**response,"detail":response.get("detail",""),"gbif_species_key":match["key"]}
        statuses[pid]=status["status"]
        report.append({
            "plant_id":pid,"wcvp_accepted_name":r["wcvp_accepted_name"],
            "wcvp_taxon_rank":r["wcvp_taxon_rank"],
            "island_code":"CNY","GBIF_query_status":status["status"],
            "gbif_species_key":status.get("gbif_species_key",""),
            "earliest_gbif_plant_indexed_year":status.get("earliest_indexed_year",""),
            "gbif_plant_total_records":status.get("gbif_record_count",""),
            "pre1900_record_count":status.get("pre1900_records",""),
            "pre2001_record_count":status.get("pre2001_records",""),
            "distinct_years":status.get("distinct_recorded_years",""),
            "error":status.get("detail","")
        })
        print("CNY_HISTORY",i,pid,status["status"],status.get("earliest_indexed_year",""),flush=True)
    effort={"status":"GEOMETRY_NOT_QUERYABLE"}
    if geom["status"]=="GEOMETRY_READY":
        try:effort=facet_year(6,geom["wkt"])
        except Exception as e:effort={"status":"GBIF_QUERY_FAILED","detail":str(e)}
    with (a.out_dir/"CNY_three_plant_years.csv").open("w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=list(report[0]));wr.writeheader();wr.writerows(report)
    audit={
      "schema":"chocho_geo_corrected_Q4_canary_three_plant_firstyear_index_v0.2",
      "status":"SOURCE_DATED_BOTANICAL_FIRST_RECORDS_NOT_INTRODUCTION",
      "source":{"q4_v02_population":meta["corrected_population"],
                "CNY_response_blind_new_pilot_n":3,"original_500_pilot_unchanged":True,
                "WCVP_git_SHA":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
                "WGSRPD_git_SHA":"52da7828aba9d461dd133c27b3bd7a4407161f54"},
      "geometry":{"status":geom["status"],"fidelity":geom.get("geometry_fidelity")},
      "plant_statuses":statuses,
      "plants_dated":sum(r["GBIF_query_status"]=="DATED_OCCURRENCE_FOUND" for r in report),
      "plant_first_years":[int(r["earliest_gbif_plant_indexed_year"]) for r in report if r["GBIF_query_status"]=="DATED_OCCURRENCE_FOUND"],
      "all_botanical_indexed_effort":{"status":effort["status"],
                                      "first_year":effort.get("earliest_indexed_year"),
                                      "total_indexed_records":effort.get("gbif_record_count"),
                                      "pre1900":effort.get("pre1900_records")},
      "last_checked_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "limits":["GBIF year facets depend on digitized specimens and sampling intensity, not actual introductions.",
                "CNY is an archipelago; this single WCVP region does not establish that a host and a butterfly appeared on the same physical island.",
                "Current WCVP alien-status does not identify arrival year, host quality or local larval reproduction.",
                "No butterfly first-record or wing×port effect is tested in this plant inventory."]
    }
    (a.out_dir/"CNY_three_plant_firstyear_audit.json").write_text(
        json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("CNY_FIRST_YEAR_RESULT",json.dumps(audit,ensure_ascii=False))
if __name__=="__main__":main()
