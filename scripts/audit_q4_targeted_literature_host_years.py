#!/usr/bin/env python3
"""Targeted historical GBIF botanical coverage for pre-existing butterfly island literature cases.

Only the fixed 47 literature leads and the audited Canary-corrected Q2 host
edges determine which host×WGSRPD3 region pairs are queried. Results are a
CASE-SELECTED evidence enrichment, not an unbiased estimate of historical
host detection or colonization. Reuses source-locked 500-pair audit wherever
it had a valid date or an explicitly successful no-index search.
"""
from __future__ import annotations
import argparse,csv,json,hashlib,datetime,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import defaultdict,Counter
from pathlib import Path

from audit_q4_gbif_plant_first_year_pilot import (
    build_geometries, key_from_gbif, facet_year, STATUS_NO_EVIDENCE,
)

def csvrows(p):
    with p.open(newline="",encoding="utf-8") as fh:
        return list(csv.DictReader(fh))
def savecsv(p, rr):
    if not rr:raise RuntimeError("No rows")
    with p.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rr[0]))
        w.writeheader()
        w.writerows(rr)

def main():
    p=argparse.ArgumentParser()
    for x in ("literature","q2-edges","q2-audit","wcvp-names-v02",
              "previous-pilot","previous-taxon","wgsrpd-geojson","out-dir"):
        p.add_argument("--"+x,type=Path,required=True)
    p.add_argument("--workers",type=int,default=2)
    a=p.parse_args()
    if not 1<=a.workers<=3:raise RuntimeError("Overconcurrent GBIF targeted recovery")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    literature=csvrows(a.literature)
    edges=csvrows(a.q2_edges)
    audit=json.loads(a.q2_audit.read_text())
    names={r["plant_id"]:r for r in csvrows(a.wcvp_names_v02)}
    previous=csvrows(a.previous_pilot)
    prior={(r["plant_id"],r["island_code"]):r for r in previous}
    previous_tax={r["plant_id"]:r for r in csvrows(a.previous_taxon)}
    if len(literature)!=47 or len(edges)!=11573 or len(names)!=564:
        raise RuntimeError("Core fixed target population changed")
    if audit["scope"]["added_strict_island_cells"]!=5741 or len(previous)!=500:
        raise RuntimeError("Q2 v02 or Q4 500-pair source cohort changed")
    q2=defaultdict(set)
    for row in edges:
        q2[(row["species"],row["island_code"])].add(row["plant_id"])
    if len(q2)!=5741:raise RuntimeError("Q2 butterfly×island unique count drift")
    eligible=[r for r in literature if (r["species"],r["wgsrpd3_code"]) in q2]
    if len(eligible)<15:
        raise RuntimeError("Insufficient fixed literature overlaps")
    cells=sorted({(r["species"],r["wgsrpd3_code"]) for r in eligible})
    targets=sorted({(pid,region) for species,region in cells for pid in q2[(species,region)]})
    if len(targets)>450:raise RuntimeError("Unexpected >450 targeted plant×island queries")
    regions={region for _,region in targets}
    geos=build_geometries(a.wgsrpd_geojson,regions)
    targetplants={x for x,_ in targets}
    old_success={"DATED_OCCURRENCE_FOUND",STATUS_NO_EVIDENCE}
    cached={(pid,area):prior[(pid,area)] for pid,area in targets
        if (pid,area) in prior and prior[(pid,area)]["q4_pilot_status"] in old_success}
    taxonomy={}
    for pid in sorted(targetplants):
        obj=names[pid]
        if obj["gbif_species_query_eligible"]!="True":
            taxonomy[pid]={"status":"WCVP_NON_SPECIES_OR_NONACCEPTED_RANK",
                            "key":"","canonical":obj["wcvp_accepted_name"],
                            "confidence":"","match_type":"","detail":"Source species-only screen"}
        elif pid in previous_tax and previous_tax[pid]["status"]=="EXACT_SPECIES_MATCH":
            taxonomy[pid]={"status":"EXACT_SPECIES_MATCH",
                           "key":previous_tax[pid]["key"],
                           "canonical":previous_tax[pid]["canonical"],
                           "confidence":previous_tax[pid]["confidence"],
                           "match_type":previous_tax[pid]["match_type"],
                           "detail":"source_frozen_500_pilot_taxonomic_match"}
    needs=[p for p in sorted(targetplants) if p not in taxonomy]
    def match_one(pid):
        return pid,key_from_gbif(names[pid]["wcvp_accepted_name"])
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(match_one,pid) for pid in needs]):
            pid,d=f.result()
            taxonomy[pid]=d
    if set(taxonomy)!=targetplants:
        raise RuntimeError("Incomplete taxonomy lookup")
    savecsv(a.out_dir/"targeted_plant_taxon_source.csv",[
       {"plant_id":p,"accepted_WCVP_name":names[p]["wcvp_accepted_name"],
        "WCVP_rank":names[p]["wcvp_taxon_rank"],
        "WCVP_accepted_status":names[p]["wcvp_taxon_status"],
        **taxonomy[p]} for p in sorted(targetplants)])

    toquery=[]
    for pid,region in targets:
        if (pid,region) in cached:continue
        tax=taxonomy[pid]
        g=geos[region]
        if tax["status"]!="EXACT_SPECIES_MATCH" or g["status"]!="GEOMETRY_READY":continue
        toquery.append(("PLANT",pid,region,str(tax["key"]),g["wkt"]))
    # Also query independent annual botanical-index effort for the same regions.
    # 'taxonKey=6' is all vascular/plantae kingdom as in the preceding Q4 pilot,
    # a reporting-effort proxy, not an unbiased botanical completeness estimator.
    for region in sorted(regions):
        if geos[region]["status"]=="GEOMETRY_READY":
            toquery.append(("EFFORT","6",region,"6",geos[region]["wkt"]))
    result={}
    def get_year(task):
        typ,pid,reg,key,wkt=task
        try:
            dat=facet_year(key,wkt)
            dat["query_type"]=typ
        except Exception as ex:
            dat={"status":"GBIF_QUERY_FAILED","error":type(ex).__name__+": "+str(ex)[:220],
                 "gbif_record_count":"","earliest_indexed_year":"","query_type":typ}
        return (typ,pid,reg),dat
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for num, future in enumerate(as_completed([pool.submit(get_year,task) for task in toquery]),1):
            key,dat=future.result()
            result[key]=dat
            if num%25==0:
                print("Q4_CASE_TARGETED_YEAR_FACETS",num,len(toquery),dict(Counter(
                    row["status"] for row in result.values())),flush=True)
    reportrows=[]
    for pid,reg in targets:
        old=cached.get((pid,reg))
        if old:
            data={"status":old["q4_pilot_status"],
                  "earliest_indexed_year":old["earliest_gbif_indexed_year"],
                  "gbif_record_count":old["gbif_plant_record_count"],
                  "pre1900_records":old["gbif_plant_pre1900_count"],
                  "pre2001_records":old["gbif_plant_pre2001_count"]}
            mode="REUSED_FROZEN_PILOT_SUCCESS"
        else:
            data=result.get(("PLANT",pid,reg))
            if data is None:
                status=(taxonomy[pid]["status"] if taxonomy[pid]["status"]!="EXACT_SPECIES_MATCH" else geos[reg]["status"])
                data={"status":status,"earliest_indexed_year":"",
                      "gbif_record_count":"","error":""}
            mode="NEW_TARGETED_CASE_QUERY"
        reportrows.append({
          "plant_id":pid,"wcvp_plant_name":names[pid]["wcvp_accepted_name"],
          "wgsrpd3_code":reg,
          "gbif_species_key":taxonomy[pid].get("key",""),
          "gbif_taxon_match_status":taxonomy[pid]["status"],
          "source_query_mode":mode,
          "status":data.get("status"),
          "first_gbif_indexed_plant_year":data.get("earliest_indexed_year",""),
          "gbif_indexed_plant_record_count":data.get("gbif_record_count",""),
          "gbif_pre1900_count":data.get("pre1900_records",""),
          "gbif_pre2001_count":data.get("pre2001_records",""),
          "wgsrpd3_polygon_query_quality":geos[reg].get("geometry_fidelity"),
          "query_failure":data.get("error","")
        })
    savecsv(a.out_dir/"targeted_introduced_plant_island_years.csv",reportrows)
    bypair={(r["plant_id"],r["wgsrpd3_code"]):r for r in reportrows}
    byregion={}
    for reg in sorted(regions):
        data=result.get(("EFFORT","6",reg))
        if data is None:
            data={"status":geos[reg]["status"],"gbif_record_count":"","earliest_indexed_year":"",
                 "pre1900_records":"","pre2001_records":""}
        byregion[reg]=data
    effrows=[{"wgsrpd3_code":reg,"index_status":dat["status"],
              "all_plant_records":dat.get("gbif_record_count",""),
              "all_plant_earliest_index_year":dat.get("earliest_indexed_year",""),
              "all_plant_pre1900_record_count":dat.get("pre1900_records",""),
              "all_plant_pre2001_record_count":dat.get("pre2001_records",""),
              "error":dat.get("error","")} for reg,dat in sorted(byregion.items())]
    savecsv(a.out_dir/"case_regions_botanical_reporting_effort.csv",effrows)
    out_events=[]
    for c in literature:
        sp,reg=c["species"],c["wgsrpd3_code"]
        hosts=q2.get((sp,reg),())
        plants=[bypair[(p,reg)] for p in hosts] if hosts else []
        dated=[p for p in plants if p["status"]=="DATED_OCCURRENCE_FOUND"]
        firstyear=int(c["first_evidence_earliest_year"])
        prior=[p for p in dated if int(p["first_gbif_indexed_plant_year"])<firstyear]
        out_events.append({
          "species":sp,
          "physical_island":c["physical_island"],
          "wgsrpd3_code":reg,
          "first_butterfly_documented_year_min":firstyear,
          "first_butterfly_documented_year_max":int(c["first_evidence_latest_year"]),
          "documented_establishment_by_year":c["established_by_year"],
          "resource_added_by_modern_WCVP_HOSTS":int(bool(hosts)),
          "modern_introduced_host_candidates":len(hosts),
          "host_plant_region_pairs_with_GBIF_firstyear":len(dated),
          "host_plant_region_pairs_documented_before_butterfly":len(prior),
          "earliest_documented_candidate_host_year":min(
             (int(p["first_gbif_indexed_plant_year"]) for p in dated),default=""),
          "earliest_preliterature_candidate_host_year":min(
             (int(p["first_gbif_indexed_plant_year"]) for p in prior),default=""),
          "same_physical_island_plant_record_confirmed":0,
          "observation_not_causal_arrival_order":1,
          "indexed_botanical_effort_status":byregion.get(reg,{}).get("status","UNKNOWN_REGION"),
          "indexed_botanical_pre2001_records":byregion.get(reg,{}).get("pre2001_records",""),
        })
    savecsv(a.out_dir/"literature_case_host_year_order.csv",out_events)
    positives={r["species"]+"|"+r["wgsrpd3_code"] for r in out_events
               if r["host_plant_region_pairs_documented_before_butterfly"] and r["wgsrpd3_code"]}
    stat={
      "schema":"chocho_Q4_targeted_literature_host_index_dates_v0.1",
      "status":"EXPLORATORY_LITERATURE_SELECTED_YEAR_OVERLAPS_NOT_COLONIZATION",
      "source":{"candidate_literature_events":len(literature),
                 "unique_eligible_butterfly_region_cells":len(cells),
                 "Q2_corrected_island_resource_cells":len(q2),
                 "Q2_corrected_source_edges":len(edges),
                 "GBIF_earliest_species_index_range":"1700-2026",
                 "recorded_on_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()},
      "targets":{"distinct_plant_island_pairs":len(targets),
                 "distinct_plants":len(targetplants),"regions":len(regions),
                 "reused_previous_GBIF_query_success":len(cached),
                 "new_GBIF_queries_attempted":sum(x[0]=="PLANT" for x in toquery),
                 "botanical_effort_region_queries":sum(x[0]=="EFFORT" for x in toquery)},
      "plant_query_status":dict(Counter(r["status"] for r in reportrows)),
      "plant_pairs_with_first_index_year":sum(r["status"]=="DATED_OCCURRENCE_FOUND" for r in reportrows),
      "literature_events_Q2_modern_host_resource_gain":sum(r["resource_added_by_modern_WCVP_HOSTS"] for r in out_events),
      "literature_events_with_dated_compatible_plant_index":sum(r["host_plant_region_pairs_with_GBIF_firstyear"]>0 for r in out_events),
      "literature_events_with_compatible_plant_index_before_butterfly_first_report":sum(
        r["host_plant_region_pairs_documented_before_butterfly"]>0 for r in out_events),
      "distinct_butterfly_WGSRPD3_cells_with_prior_host_index":len(positives),
      "same_physical_island_host_prior_verified":0,
      "botanical_effort_status":dict(Counter(r["index_status"] for r in effrows)),
      "limits":[
       "Target pairs are chosen conditional on a published butterfly first report: no natural control cohort, not a confirmatory test of lag or access.",
       "GBIF first indexed plant year is not date of arrival, planting, establishment or naturalization; retrospective date gaps are only record-order gaps.",
       "A broad WGSRPD3 region can contain multiple physical islands, including the Canary Islands; regional index order does not prove physical-island temporal precedence.",
       "The same butterfly species can have multiple physical-island literature cases within a single region: these are non-independent.",
       "No butterfly absence or island-first-colonization risk set, no species climate suitability, no port vessel arrivals and no causal last-gate labels.",
       "Some historical botanical first records are controlled by taxonomic and detection biases; raw documentation dates are not comparable without effort by era.",
      ],
    }
    (a.out_dir/"targeted_literature_hostyear_audit.json").write_text(
        json.dumps(stat,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("Q4_LITERATURE_TARGET_GBIF",json.dumps(stat,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
