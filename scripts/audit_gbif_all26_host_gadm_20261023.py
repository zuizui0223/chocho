#!/usr/bin/env python3
"""GBIF GADM1 plant-record coverage for all 26 added-host butterfly-positive country cells.

All 44 original source WCVP accepted plant-ID × country queries and 24 exact
botanical taxon names are fixed before GBIF observation outcomes. No
individual record bodies, local feeding, recruitment or abundance inference.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

SOURCE_SCHEMA="chocho_all26_positive_country_source_hosts_v1"
OUTPUT_SCHEMA="chocho_all26_gbif_plant_gadm1_feasibility_20261023_v1"
EXPECTED_COUNTS=(48,26,46,24,44)
ORIGINAL_FIVE_SCHEMA="chocho_strict_plant_taxon_gadm_cocoverage_five_juvenile_country_cases_20261021_v1"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def load_inputs(source_path, butterfly_path, original_five_path):
    raw=Path(source_path).read_bytes()
    source=json.loads(raw)
    bf=json.loads(Path(butterfly_path).read_bytes())
    old=json.loads(Path(original_five_path).read_bytes())
    if source.get("schema")!=SOURCE_SCHEMA or source.get("status")!="SOURCE_ONLY_ALL26_NO_GBIF_PLANT_RECORDS":
        raise ValueError("Full 26-cell source-only botanical cohort not fixed")
    if (source.get("frozen_48_positive_species_country_cells"),source.get("frozen_added_host_cells"),
        source.get("candidate_butterfly_country_plant_links"),source.get("distinct_source_plant_ids"),
        source.get("distinct_plant_country_source_queries"))!=EXPECTED_COUNTS:
        raise ValueError("Full 26 botanical source cohort numerical margins changed")
    if bf.get("status")!="COMPLETE_48" or len(bf.get("observations",[]))!=48:
        raise ValueError("All 48 GBIF butterfly GADM1 facets must remain complete")
    if old.get("schema")!=ORIGINAL_FIVE_SCHEMA or old.get("status")!="ALL_FIVE_CASES_SOURCE_ASSESSED":
        raise ValueError("Previous 20 original botanical GADM1 source facets not complete")
    all48={(r["species"],r["country"]):r for r in bf["observations"]}
    if len(all48)!=48:
        raise ValueError("Original GADM1 source positive identities duplicate")
    hosts=source["all_48_positive_country_cells"]
    positive={(r["butterfly"],r["country"]):r for r in hosts
              if r["added_host_resource_country"]}
    if len(positive)!=26 or len(hosts)!=48:
        raise ValueError("Missing added host source cases")
    ids={}
    for r in source["unique_plant_country_query_basis"]:
        if r["source_name_gate"]!="EXACT_SOURCE_BINOMIAL_NEEDS_GBIF_TAXON_RESOLUTION" or len(r["accepted_names"])!=1:
            raise ValueError("Unresolved WCVP source botanical name: no nonexact GBIF query")
        pid=r["wcvp_accepted_id"];name=r["accepted_names"][0]
        if pid in ids and ids[pid]!=name:
            raise ValueError("One WCVP accepted ID has conflicting names")
        ids[pid]=name
    if len(ids)!=24:
        raise ValueError("Twenty-four unique accepted botanical species required")
    pairs={}
    total_links=0
    for (sp,ct),row in positive.items():
        if (sp,ct) not in all48 or all48[sp,ct]["added_resource_country"]!=1:
            raise ValueError("Source butterfly positive added flag differs")
        plants=row["candidate_accepted_plants"]
        if len(plants)!=row["candidate_source_plant_count"] or not plants:
            raise ValueError("Candidate plant-country source not complete")
        for p in plants:
            pid=p["wcvp_accepted_id"]
            if ids.get(pid) not in p["accepted_names"]:
                raise ValueError("Candidate plant taxon is not source-resolved")
            key=(pid,ct)
            pairs.setdefault(key,{"plant_name":ids[pid],"country":ct,"wcvp_accepted_id":pid,
                                  "butterflies_linked":set()})["butterflies_linked"].add(sp)
            total_links+=1
    if len(pairs)!=44 or total_links!=46:
        raise ValueError("Source links/country plant deduplication changed")
    source_query_set={(p["wcvp_accepted_id"],p["country"]) for p in source["unique_plant_country_query_basis"]}
    if set(pairs)!=source_query_set:
        raise ValueError("Source-only 44 country/plant queries do not match cases")
    prior={(r["plant_name"],r["country"]):r for r in old["per_source_host_country"]}
    if len(prior)!=20:
        raise ValueError("Previous pilot 20 plant/country comparison missing")
    return source,all48,old,ids,pairs,prior,{
        "source_inventory_sha256":sha(raw),
        "source_positive48_sha256":sha(Path(butterfly_path).read_bytes()),
        "previous_20_plant_facets_sha256":sha(Path(original_five_path).read_bytes()),
        "source_original_zip_sha256":source["provenance_sha256"]["original_archive_inner_zip"]}


def audit(source,all48,old,ids,pairs,prior,taxonomy,results,metadata):
    if set(taxonomy)!=set(ids):
        raise ValueError("Missing botanical WCVP ID taxonomy attempts")
    if set(results)!=set(pairs):
        raise ValueError("Missing fixed source plant/country queries")
    output=[]
    drift=[]
    for key in sorted(pairs):
        spec=pairs[key]
        result=results[key]
        if result["status"]=="COMPLETE_INDEXED_PLANT_FACET":
            n=result["indexed_country_plant_records"]
            g=result["gadm1_counts"]
            if n<0 or sum(g.values())+result["gadm1_unmapped_records"]!=n:
                raise ValueError("GADM1 botanical facet margins not complete")
        if (spec["plant_name"],spec["country"]) in prior:
            original=prior[spec["plant_name"],spec["country"]]
            if result["status"]=="COMPLETE_INDEXED_PLANT_FACET":
                delta=result["indexed_country_plant_records"]-original["indexed_country_plant_records"]
                shift={k:result["gadm1_counts"].get(k,0)-original["gadm1_counts"].get(k,0)
                       for k in sorted(set(result["gadm1_counts"])|set(original["gadm1_counts"]))
                       if result["gadm1_counts"].get(k,0)!=original["gadm1_counts"].get(k,0)}
                drift.append({"plant":spec["plant_name"],"country":spec["country"],
                              "old_record_count":original["indexed_country_plant_records"],
                              "new_record_count":result["indexed_country_plant_records"],
                              "count_delta":delta,"gadm1_nonzero_changes":shift})
        output.append({
            "wcvp_accepted_id":key[0],"country":key[1],
            "plant_name":spec["plant_name"],
            "linked_butterflies":sorted(spec["butterflies_linked"]),
            "gbif_taxonomy":taxonomy[key[0]], **result
        })
    bycase=[]
    for row in source["all_48_positive_country_cells"]:
        if not row["added_host_resource_country"]:
            continue
        sp,ct=row["butterfly"],row["country"]
        sub=[v for v in output if v["country"]==ct and sp in v["linked_butterflies"]]
        if len(sub)!=row["candidate_source_plant_count"]:
            raise ValueError("Source candidate plant selection drifted")
        good=[p for p in sub if p["status"]=="COMPLETE_INDEXED_PLANT_FACET"]
        all_recorded=len(good)==len(sub)
        occupied=set(g for p in good for g,n in p["gadm1_counts"].items() if n>0)
        gf=all48[sp,ct]["facets"]["gadm1"]["categories"]
        N=sum(gf.values())
        exposed_count=sum(n for g,n in gf.items() if g in occupied)
        co=len(set(gf)&occupied)
        bycase.append({
            "butterfly":sp,"country":ct,
            "source_candidate_host_species":len(sub),
            "complete_botanical_query_count":len(good),
            "index_plant_positive_source_links":sum(x["indexed_country_plant_records"]>0 for x in good),
            "total_botanical_indexed_records":sum(x["indexed_country_plant_records"] for x in good),
            "candidate_plants_detected_in_country":sum(x["indexed_country_plant_records"]>0 for x in good),
            "all_butterfly_mapped_gadm1_records":N,
            "butterfly_positive_gadm1_units":len(gf),
            "candidate_plant_positive_gadm1_units":len(occupied),
            "butterfly_gadm1_units_with_candidate_plant_records":co if all_recorded else None,
            "butterfly_records_in_candidate_plant_positive_gadm1":exposed_count if all_recorded else None,
            "source_all_plant_facets_complete":all_recorded,
            "interpretation":"Administrative co-detection only, not verified larval host use",
        })
    completed=[r for r in output if r["status"]=="COMPLETE_INDEXED_PLANT_FACET"]
    fully=[r for r in bycase if r["source_all_plant_facets_complete"]]
    complete=len(completed)==44 and all(r["status"]=="EXACT_ACCEPTED_GBIF_PLANT" for r in taxonomy.values())
    return {
        "schema":OUTPUT_SCHEMA,"status":"ALL26_BOTANICAL_GADM_QUERIES_COMPLETE" if complete else "PARTIAL_GBIF_OR_TAXONOMY_NOT_PLANT_ABSENCE",
        "source":metadata,"executed_utc":now(),
        "taxa_exact_accepted_in_GBIF":sum(r["status"]=="EXACT_ACCEPTED_GBIF_PLANT" for r in taxonomy.values()),
        "taxa_attempted":24,"plant_country_source_pairs":44,
        "plant_country_GADM_facets_complete":len(completed),
        "positive_plant_country_pairs":sum(r["indexed_country_plant_records"]>0 for r in completed),
        "indexed_plant_source_records":sum(r["indexed_country_plant_records"] for r in completed),
        "original_five_case_plant_pairs_compared":len(drift),
        "original_five_case_pair_count_drifts":sum(x["count_delta"]!=0 for x in drift),
        "original_five_case_gadm_bin_drifts":sum(bool(x["gadm1_nonzero_changes"]) for x in drift),
        "old_source_drift_receipts":drift,
        "added_butterfly_country_cases":len(bycase),
        "complete_butterfly_country_cases":len(fully),
        "cases_with_any_GBIF_botanical_source_host_country_record":sum(r["candidate_plants_detected_in_country"]>0 for r in fully),
        "cases_with_butterfly_GADM1_overlap_with_indexed_host":sum((r["butterfly_gadm1_units_with_candidate_plant_records"] or 0)>0 for r in fully),
        "source_only_host_geography_not_local_host_use":True,
        "taxonomy":taxonomy,
        "plant_country_results":output,
        "butterfly_country_results":bycase
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for f in ("--source-inventory","--positive48","--previous-five-plant-source","--output"):
        p.add_argument(f,required=True,type=Path)
    args=p.parse_args()
    from audit_gbif_five_host_gadm_cocoverage_20261021 import (
        get,match_url,strict_match,occurrence_url,parse_plant_gadm
    )
    source,all48,old,ids,pairs,prior,meta=load_inputs(
        args.source_inventory,args.positive48,args.previous_five_plant_source)
    taxonomy={}
    for pid,name in sorted(ids.items()):
        taxonomy[pid]=strict_match(name,get(match_url(name)))
        time.sleep(1.05)
    results={}
    for key,pair in sorted(pairs.items()):
        info=taxonomy[key[0]]
        if info["status"]!="EXACT_ACCEPTED_GBIF_PLANT":
            results[key]={"status":"SOURCE_PLANT_TAXON_UNRESOLVED_NOT_ABSENCE"}
            continue
        try:
            v=parse_plant_gadm(get(occurrence_url(pair,info["accepted_species_key"])))
        except ValueError as err:
            v={"status":"INVALID_GBIF_PLANT_GADM_NOT_ABSENCE","error":str(err)[:200]}
        results[key]=v
        time.sleep(1.05)
    r=audit(source,all48,old,ids,pairs,prior,taxonomy,results,meta)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    summary={k:v for k,v in r.items() if k not in ("taxonomy","plant_country_results","butterfly_country_results","old_source_drift_receipts")}
    print("GBIF_ALL26_GADM_INDEX_SUMMARY="+json.dumps(summary,separators=(",",":")),flush=True)
    print("GBIF_ALL26_GADM_INDEX_FULL="+json.dumps(r,separators=(",",":")),flush=True)
    if r["status"]!="ALL26_BOTANICAL_GADM_QUERIES_COMPLETE":
        raise SystemExit("Incomplete original plant/GADM source: not interpretable as plant absence")


if __name__=="__main__":
    main()
