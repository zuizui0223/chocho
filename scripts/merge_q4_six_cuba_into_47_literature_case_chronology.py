#!/usr/bin/env python3
"""Rebuild fixed 84 case-targeted plant dates and 47 literature chronology rows
after six full-geometry Cuba GBIF index-date recoveries. This does NOT estimate
arrival, butterfly colonization, or historical consumer lag.
"""
import argparse,csv,json
from pathlib import Path
from collections import Counter

def csvread(p):
    with p.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def csvsave(p,data):
    if not data:raise RuntimeError("Empty source")
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]))
        w.writeheader();w.writerows(data)

def main():
    p=argparse.ArgumentParser()
    for name in ("source-targeted84","source-events47","source-audit","cuba-rescue","cuba-audit","out-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    host=csvread(a.source_targeted84)
    literature=csvread(a.source_events47)
    base=json.loads(a.source_audit.read_text(encoding="utf-8"))
    cuba=csvread(a.cuba_rescue)
    receipt=json.loads(a.cuba_audit.read_text(encoding="utf-8"))
    if len(host)!=84 or len(literature)!=47 or len(cuba)!=6 or base["plant_pairs_with_first_index_year"]!=67:
        raise RuntimeError("Original fixed literature case/host input drift")
    if receipt["original_cuban_http400_queries"]!=6 or receipt["new_dated_pairs"]!=6 or receipt["fully_covered"]!=6:
        raise RuntimeError("Cuba WKT source not completely recovered")
    patches={}
    for rr in cuba:
        if rr["region"]!="CUB" or rr["status"]!="DATED_OCCURRENCE_FOUND" or rr["complete_full_WGSRPD_region"].lower()!="true":
            raise RuntimeError("Source polygon coverage failure")
        k=(rr["plant_id"],rr["region"])
        if k in patches:raise RuntimeError("Duplicate Cuba source host")
        patches[k]=rr
    repaired=set();modified=[]
    for src in host:
        r={**src}
        key=(r["plant_id"],r["wgsrpd3_code"])
        if key in patches:
            rr=patches[key]
            if r["status"]!="GBIF_QUERY_FAILED" or "HTTP Error 400" not in r["query_failure"]:
                raise RuntimeError("Original source outcome cannot be overwritten")
            if not r["gbif_species_key"] or not rr["earliest_indexed_year_if_complete"]:
                raise RuntimeError("Missing source taxonomy or year")
            r["status"]="DATED_OCCURRENCE_FOUND"
            r["source_query_mode"]="RECOVERED_FULL_CUBA_GBIF_MULTIPART_POLYGON"
            r["first_gbif_indexed_plant_year"]=rr["earliest_indexed_year_if_complete"]
            r["gbif_indexed_plant_record_count"]=rr["total_records_if_complete"]
            r["gbif_pre1900_count"]=""
            r["gbif_pre2001_count"]=""
            r["query_failure"]=""
            repaired.add(key)
        modified.append(r)
    if repaired!=set(patches) or len({(r["plant_id"],r["wgsrpd3_code"]) for r in modified})!=84:
        raise RuntimeError("Cuba source mapping or fixed host cohort mismatch")
    newcase=[];n_cuba=0
    for src in literature:
        r={**src}
        if (r["species"],r["physical_island"],r["wgsrpd3_code"])==("Papilio demoleus","Cuba","CUB"):
            n_cuba+=1
            if int(r["first_butterfly_documented_year_min"])!=2006 or int(r["modern_introduced_host_candidates"])!=6:
                raise RuntimeError("Cuba first record/host count drift")
            plants=[x for x in modified if (x["plant_id"],"CUB") in patches]
            if len(plants)!=6:raise RuntimeError("Cuba plants not all resolved")
            years=[int(x["first_gbif_indexed_plant_year"]) for x in plants]
            previous=[y for y in years if y<2006]
            if len(previous)!=5 or min(previous)!=1902:
                raise RuntimeError("Expected five original host dates preceding 2006")
            r["host_plant_region_pairs_with_GBIF_firstyear"]=6
            r["host_plant_region_pairs_documented_before_butterfly"]=5
            r["earliest_documented_candidate_host_year"]=min(years)
            r["earliest_preliterature_candidate_host_year"]=min(previous)
        newcase.append(r)
    if n_cuba!=1:raise RuntimeError("Must affect one literature butterfly case")
    statuses=Counter(r["status"] for r in modified)
    if statuses!=Counter({"DATED_OCCURRENCE_FOUND":73,"NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT":10,
                         "WCVP_NON_SPECIES_OR_NONACCEPTED_RANK":1}):
        raise RuntimeError(f"Unreconciled 84 plant results {statuses}")
    hasdated=sum(int(r["host_plant_region_pairs_with_GBIF_firstyear"])>0 for r in newcase)
    prior=sum(int(r["host_plant_region_pairs_documented_before_butterfly"])>0 for r in newcase)
    cells={(r["species"],r["wgsrpd3_code"]) for r in newcase if r["wgsrpd3_code"] and int(r["host_plant_region_pairs_documented_before_butterfly"])>0}
    if (hasdated,prior,len(cells))!=(25,22,12):
        raise RuntimeError(f"Source case-date coverage drift: {hasdated,prior,len(cells)}")
    report={
      "schema":"chocho_Q4_case_selected_Cuba_plant_index_full_geometry_recovery_v0.2",
      "status":"SOURCE_INDEX_RECORD_ORDER_NOT_REAL_ESTABLISHMENT_OR_LAG",
      "fixed_original_plants":84,"fixed_literature_butterfly_cases":47,
      "new_source_dated_Cuba_plants":6,
      "dated_plant_pairs_before":67,"dated_plant_pairs_after":73,
      "plant_statuses":dict(statuses),
      "literature_cases_with_dated_candidate_host":hasdated,
      "literature_cases_with_documented_prior_candidate_host":prior,
      "distinct_butterfly_x_WGSRPD3_with_prior_documented_host":len(cells),
      "Papilio_demoleus_Cuba":{
         "first_butterfly_literature_record_year":2006,
         "GBIF_dated_candidate_host_pairs":6,
         "host_pairs_indexed_before_first_butterfly_record":5,
         "oldest_candidate_host_GBIF_index_year":1902,
         "difference_between_record_dates_NOT_colonization_lag_years":104,
      },
      "all_six_original_failed_gbif_queries_recovered_from_complete_official_Cuba_polygon":True,
      "limitations":[
         "Host targets were selected because butterflies were already in literature. This is not a prospective risk set.",
         "First plant GBIF indexed year and first butterfly literature year are observation dates, not introduction or island-colonization years.",
         "Modern WCVP host opportunity is not verified breeding or feeding on the same physical island.",
         "The Cuban botanical source year can be shifted by historical collecting effort and taxonomy.",
         "Do not use a 104-year record-order interval as biological butterfly lag, plant arrival time or access gate duration."
      ]
    }
    csvsave(a.out_dir/"targeted_84_plant_firstyears_v02.csv",modified)
    csvsave(a.out_dir/"literature_47_host_record_order_v02.csv",newcase)
    (a.out_dir/"Cuba_source_documentation_correction_audit.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
