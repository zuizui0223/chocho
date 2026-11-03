#!/usr/bin/env python3
"""Resolve exact historic GBIF plant indexed years for 13 butterfly DOI-indexed
island resources and 39 response-matched nonindexed *DOI* comparison cells.

No butterfly absence, biological establishment, plant import year, flight vs
ship mode, or causation is inferred. This is a source-dated documentation audit.
"""
from __future__ import annotations
import argparse,csv,json,time
from collections import defaultdict,Counter
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import audit_q4_gbif_plant_first_year_pilot as gbif

def csv_read(path):
    with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def csv_save(path,values):
    if not values:raise RuntimeError("no values")
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(values[0]))
        w.writeheader();w.writerows(values)
def exact_year(s):
    if not str(s or "").strip():return None
    k=int(str(s))
    if not 1700<=k<=2026:raise RuntimeError("bad indexed year "+str(s))
    return k

def main():
    p=argparse.ArgumentParser()
    for name in ("matched-host-pairs","matched-species-islands","q2-edges",
                 "wcvp-plant-names","pilot-plant-years","targeted-plant-years",
                 "canary-plant-years","q3-positive","wgsrpd3-geojson","output-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    pairs=csv_read(a.matched_host_pairs)
    cases=csv_read(a.matched_species_islands)
    taxa=csv_read(a.wcvp_plant_names)
    pilot=csv_read(a.pilot_plant_years)
    targeted=csv_read(a.targeted_plant_years)
    canary=csv_read(a.canary_plant_years)
    pos=csv_read(a.q3_positive)
    q2_edges=csv_read(a.q2_edges)
    if len(q2_edges)!=11573:raise RuntimeError("Q2 corrected source host edges drift")
    pair_membership=defaultdict(set)
    for e in q2_edges:
        pair_membership[(e["species"],e["island_code"])].add(e["plant_id"])
    if len(pair_membership)!=5741:raise RuntimeError("Q2 host resource membership drift")
    if (len(pairs),len(cases),len(taxa),len(pilot),len(targeted),len(canary),len(pos))!=(100,52,44,500,84,3,13):
        raise RuntimeError("frozen source populations changed")
    pids={r["plant_id"] for r in pairs}
    if len({(r["plant_id"],r["island_code"]) for r in pairs})!=100 or pids!={r["plant_id"] for r in taxa}:
        raise RuntimeError("not 100 unique host-island and 44 distinct plant IDs")
    names={r["plant_id"]:r for r in taxa}
    selection=defaultdict(set)
    for r in pairs:selection[r["island_code"]].add(r["plant_id"])
    geos=gbif.build_geometries(a.wgsrpd3_geojson,set(selection))
    # Source cache is an indexing optimization. Do NOT treat source queries
    # selected by literature butterflies as an independent control sample.
    cached={}
    keys=defaultdict(set)
    def add(plant,code,status,year,gbifkey,source):
        if plant not in pids:return
        key=(plant,code)
        if status not in ("DATED_OCCURRENCE_FOUND",gbif.STATUS_NO_EVIDENCE):
            return
        val=exact_year(year) if status=="DATED_OCCURRENCE_FOUND" else None
        if status=="DATED_OCCURRENCE_FOUND" and val is None:
            raise RuntimeError("Date flagged but empty indexed year")
        if key in cached and (cached[key]["year"]!=val or cached[key]["status"]!=status):
            raise RuntimeError(f"Contradictory cached GBIF result for {key}, query sources need reconciliation")
        cached[key]={"status":status,"year":val,"source":source}
        if gbifkey:keys[plant].add(str(gbifkey))
    for r in pilot:add(r["plant_id"],r["island_code"],r["q4_pilot_status"],
        r["earliest_gbif_indexed_year"],r["gbif_species_key"],"Q4_BLIND_500")
    for r in targeted:add(r["plant_id"],r["wgsrpd3_code"],r["status"],
        r["first_gbif_indexed_plant_year"],r["gbif_species_key"],"Q4_LITERATURE_TARGETED_84")
    for r in canary:add(r["plant_id"],r["island_code"],r["GBIF_query_status"],
        r["earliest_gbif_plant_indexed_year"],r["gbif_species_key"],"Q4_CNY_PRESELECTED_3")

    matched={}
    for plant in sorted(pids):
        tax=names[plant]
        if tax["gbif_species_query_eligible"]!="True":
            matched[plant]={"status":"WCVP_NON_SPECIES_OR_UNACCEPTED","speciesKey":""}
        elif len(keys[plant])>1:
            matched[plant]={"status":"CONFLICTED_GBIF_PLANT_KEY","speciesKey":""}
        elif keys[plant]:
            matched[plant]={"status":"VERIFIED_PREVIOUS_GBIF_SPECIES_KEY","speciesKey":next(iter(keys[plant]))}
        else:
            obj=gbif.key_from_gbif(tax["wcvp_accepted_name"])
            matched[plant]={"status":obj["status"],"speciesKey":obj.get("key") or "",
                            "detail":obj.get("detail","")}
        time.sleep(.13)
    counts=Counter(item["status"] for item in matched.values())
    print("GBIF matched taxonomy",dict(counts),flush=True)
    results=[]
    missing=[]
    for r in pairs:
        plant,code=r["plant_id"],r["island_code"]
        ck=(plant,code)
        if ck in cached:
            base=dict(cached[ck])
            base["status"]="DATED_OCCURRENCE_FOUND" if base["year"] is not None else gbif.STATUS_NO_EVIDENCE
        elif geos[code]["status"]!="GEOMETRY_READY":
            base={"status":geos[code]["status"],"year":None,"source":"NEVER_QUERIED_COMPLEX_REGION"}
        elif not matched[plant]["speciesKey"]:
            base={"status":matched[plant]["status"],"year":None,"source":"NEVER_QUERIED_UNRESOLVED_TAXON"}
        else:
            missing.append((plant,code,matched[plant]["speciesKey"],geos[code]["wkt"]))
            continue
        results.append({**r,**{"wcvp_plant_name":names[plant]["wcvp_accepted_name"],
            "source_date_status":base["status"],
            "gbif_species_key":matched[plant]["speciesKey"],
            "first_gbif_indexed_plant_year":base["year"] if base["year"] is not None else "",
            "source_reuse":base["source"],"error":""}})
    print("GBIF fresh queries required",len(missing),flush=True)
    for k,(plant,code,taxon,wkt) in enumerate(missing,1):
        error=""
        try:
            obj=gbif.facet_year(taxon,wkt)
            status=obj["status"]
            year=obj["earliest_indexed_year"]
        except Exception as e:
            status="GBIF_QUERY_FAILED"
            year=""
            error=str(e)[:200]
        results.append({"plant_id":plant,"island_code":code,
          "matched_case_source_species":next(r["matched_case_source_species"] for r in pairs if r["plant_id"]==plant and r["island_code"]==code),
          "selected_case_or_control":next(r["selected_case_or_control"] for r in pairs if r["plant_id"]==plant and r["island_code"]==code),
          "matched_on":next(r["matched_on"] for r in pairs if r["plant_id"]==plant and r["island_code"]==code),
          "wcvp_plant_name":names[plant]["wcvp_accepted_name"],
          "source_date_status":status,"gbif_species_key":taxon,
          "first_gbif_indexed_plant_year":year,
          "source_reuse":"NEW_DOI_MATCHED_QUERY","error":error})
        if k%12==0:print("GBIF matching source dates",k,len(missing),dict(Counter(x["source_date_status"] for x in results)),flush=True)
        # Sequential, conservative public API use. Original Q4 had 429 errors
        # under six workers; do not create another concurrent burst.
        time.sleep(.75)
    if len(results)!=100:raise RuntimeError("Not all 100 plant-island matching statuses accounted")
    results.sort(key=lambda x:(x["island_code"],x["plant_id"]))
    csv_save(a.output_dir/"matched_host_plant_gbif_firstyears.csv",results)
    byhost={(r["plant_id"],r["island_code"]):r for r in results}
    years={(r["species"],r["area_code_l3"]):exact_year(r["first_indexed_year_in_DOI_subset"]) for r in pos}
    cellrows=[]
    for row in cases:
        sp,code=row["species"],row["island_code"]
        index_case=row["index_case"]
        # Controls inherit the original indexed case's first DOI year as a
        # *comparison cutoff*. This is not a predicted arrival year.
        case_cutoff=years.get((sp,index_case))
        if not case_cutoff:
            raise RuntimeError(f"No original case record year for matched block {sp}:{index_case}")
        hosts=pair_membership.get((sp,code),set())
        if not hosts:raise RuntimeError(f"Missing Q2 hosts for {sp} {code}")
        audited=[byhost[(h,code)] for h in hosts if (h,code) in byhost]
        dated=[int(x["first_gbif_indexed_plant_year"]) for x in audited
               if x["source_date_status"]=="DATED_OCCURRENCE_FOUND"]
        before=[y for y in dated if y<case_cutoff]
        n_undated=sum(x["source_date_status"]!="DATED_OCCURRENCE_FOUND" for x in audited)
        cellrows.append({
          **row,"case_reference_first_indexed_butterfly_year":case_cutoff,
          "butterfly_resource_host_taxa_in_Q2":len(hosts),
          "host_pairs_in_exact_matched_retrieval":len(audited),
          "host_pairs_with_GBIF_year":len(dated),
          "host_pairs_without_GBIF_year_or_failed":n_undated,
          "earliest_GBIF_indexed_host_year":min(dated) if dated else "",
          "host_first_index_year_before_case_butterfly_index":int(bool(before)),
          "earliest_GBIF_host_year_before_case_DOI_index":min(before) if before else "",
          "document_index_year_gap_at_most":case_cutoff-min(before) if before else "",
          "all_Q2_candidate_hosts_queried":int(len(audited)==len(hosts)),
        })
    if sum(x["recorded"]=="1" for x in cellrows)!=13:
        raise RuntimeError("DOI positives not retained")
    csv_save(a.output_dir/"matched_species_island_index_calendar.csv",cellrows)
    audit={
      "schema":"chocho_q4_matched_case_control_plant_year_retrieval_v0.1",
      "status":"INDEXED_HOST_SOURCE_DATES_NOT_COLONIZATION",
      "scope":{"matched_butterfly_region_cells":52,"indexed_butterfly_cases":13,
        "nonindexed_DOI_controls_NOT_absences":39,"host_plant_island_pairs":100,
        "source_plant_taxa":44,
        "GBIF_new_fresh_plant_region_queries":len(missing),
        "GBIF_source_cached_pairs":sum((x["plant_id"],x["island_code"]) in cached for x in pairs)},
      "taxon":dict(counts),
      "geometry":dict(Counter(x["status"] for x in geos.values())),
      "outcomes":dict(Counter(x["source_date_status"] for x in results)),
      "source_reuse":dict(Counter(x["source_reuse"] for x in results)),
      "matched_index_cases_with_any_dated_host":
        sum(x["recorded"]=="1" and x["host_pairs_with_GBIF_year"]>0 for x in cellrows),
      "matched_index_cases_with_host_first_index_before_butterfly_DOI":
        sum(x["recorded"]=="1" and x["host_first_index_year_before_case_butterfly_index"] for x in cellrows),
      "matched_nonindexed_control_cells_with_prior_host_documentation":
        sum(x["recorded"]=="0" and x["host_first_index_year_before_case_butterfly_index"] for x in cellrows),
      "matched_species_island_cells_with_all_source_host_pairs_queried":
        sum(x["all_Q2_candidate_hosts_queried"] for x in cellrows),
      "not_full_year_riskset":True,
      "notes":[
         "The source sample is 13 DOI observed and 39 DOI nonindexed cells selected by butterfly outcome, not independent arrivals or absences.",
         "GBIF first year is an upper bound on arrival time, subject to incomplete historical digitization and taxonomic mismatch.",
         "The 100 pairs are union of potential host resources across selected butterfly cells. A downstream join must use the full corrected Q2 edges; do not associate all plants on an island with every butterfly.",
         "Unresolved plant taxon / complex island polygon / GBIF query failure is missing, never botanical zero.",
         "Regional WGSRPD3 host record can be on a different physical island from the butterfly DOI index.",
         "Port exposure has already entered matching; no port coefficient or shipping-effect estimate is identified from this design."
      ]
    }
    (a.output_dir/"matched_casecontrol_plant_year_audit.json").write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RESULT",json.dumps(audit,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
