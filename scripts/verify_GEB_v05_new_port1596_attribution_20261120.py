#!/usr/bin/env python3
"""Verify permanently archived NEW GBIF all1596-port datasetKey source index.

No API calls, no source modifications. The historical 2026 WKT-search GBIF
publisher attribution stays unresolved. This independently re-computes
counts and source identifiers from the JSON archived by Actions.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/geb_v05_gbif_source_attribution_20261120"
DATA=BASE/"NEW_GBIF_ORIGINAL1596_PORTS_DATASETKEY_INDEX_20261120.json"
SUMMARY=BASE/"VERIFIED_NEW1596_GBIF_SOURCE_SUMMARY_20261120.json"

def run():
    src=DATA.read_bytes()
    summary=json.loads(SUMMARY.read_text(encoding="utf-8"))
    q=json.loads(src)
    s=summary["new_2026_11_20_GBIF_source"]
    c=summary["citation_and_scientific_guards"]
    if hashlib.sha256(src).hexdigest()!=s["source_file_sha256"]:
        raise AssertionError("Original 1596-source new snapshot SHA256 changed")
    if hashlib.sha1(b"blob "+str(len(src)).encode()+b"\0"+src).hexdigest()!=s["source_file_git_blob"]:
        raise AssertionError("Original Git blob differs from expected complete current-source file")
    if q["schema"]!="chocho_all1596_original_WPI2019_ports_NEW_current_GBIF_datasetKey_facets_20261120_v2":
        raise AssertionError("Source schema changed")
    rows=q["new_source_port_index_facets"]
    ids=[int(x["original_WPI2019_OBJECTID"]) for x in rows]
    if ids!=sorted(ids) or len(ids)!=1596 or len(set(ids))!=1596:
        raise AssertionError("Original source 1596 port roster missing, duplicated or not sorted")
    sums=[]
    keys={}
    for row in rows:
        p=row["original_WPI2019_OBJECTID"]
        providers=row["current_GBIF_record_counts_by_datasetKey"]
        names={a["datasetKey"] for a in providers}
        if len(names)!=len(providers) or len(providers)!=row["new_datasetKey_count"]:
            raise AssertionError(f"Source port {p} has duplicated or missing publisher datasetKeys")
        total=sum(int(a["count"]) for a in providers)
        if total!=row["new_2026_GBIF_index_count"]:
            raise AssertionError(f"Port {p} summed publisher count mismatched")
        for a in providers:
            key=str(a["datasetKey"])
            keys[key]=keys.get(key,0)+int(a["count"])
        sums.append(total)
    counts=q["new_source_provider_counts"]
    indexed={d["datasetKey"]:int(d["sum_port_buffer_index_hits_NOT_UNIQUE_OCCURRENCES"]) for d in counts}
    if len(indexed)!=len(counts) or indexed!=keys:
        raise AssertionError("Provider aggregate receipts disagree with original per-port source")
    if len(keys)!=s["distinct_current_GBIF_datasetKeys"] or len(keys)!=1482:
        raise AssertionError("Complete NEW GBIF datasetKey number drifted")
    if sum(sums)!=q["new_source_sum_port_buffer_indexed_hits_WITH_OVERLAP"] or sum(sums)!=s["sum_GBIF_port_buffer_record_hits_including_overlap"] or sum(sums)!=3138675:
        raise AssertionError("All-port counted public source GBIF hits changed")
    positive=sum(n>0 for n in sums)
    if positive!=s["GBIF_positive_port_search_buffers"] or positive!=1461 or 1596-positive!=s["GBIF_zero_index_port_search_buffers"]:
        raise AssertionError("Positive vs empty indexing source panel drift")
    top=sorted(indexed.items(),key=lambda t:-t[1])
    if top[0]!=(s["top1_source_datasetKey"],s["top1_source_provider_port_buffer_record_hits"]):
        raise AssertionError("Current source top provider changed")
    top1=top[0][1]/sum(sums)
    top10=sum(x[1] for x in top[:10])/sum(sums)
    if not math.isclose(top1,s["top1_aggregate_share_of_port_buffer_hits"],abs_tol=1e-10) or not math.isclose(top10,s["top10_provider_aggregate_share_of_port_buffer_hits"],abs_tol=1e-10):
        raise AssertionError("Source provider concentration metrics changed")
    byport={int(z["original_WPI2019_OBJECTID"]):z["new_2026_GBIF_index_count"] for z in rows}
    if {str(i):byport[int(i)] for i in s["seven_port_previous_GBIF_pilot_index_counts_reproduced"]}!=s["seven_port_previous_GBIF_pilot_index_counts_reproduced"]:
        raise AssertionError("Seven selected original port pilot counts no longer match")
    shards=q["new_source_shard_receipts"]
    if len(shards)!=4 or sorted(t["part"] for t in shards)!=[0,1,2,3] or any(t["original_port_ids"]!=399 or t["source_workflow_run"]!=44742493164 for t in shards):
        raise AssertionError("Original four 399-port source-shard identity missing")
    if q["source_original_2026_GBIF_contributing_datasetKeys_recovered"] or q["source_original_2026_plant_GBIF_derived_DOI_registered"] or c["historic_2026_original_GBIF_faceted_1596_datasetKey_attribution_completed"] or c["retrospectively_minted_GBIF_derived_DOI"] or c["source_WPI_original_5km_WKT_buffer_is_exactly_current_GBIF_geoDistance_5km"]:
        raise AssertionError("Source snapshot mistakenly relabelled historical DOI or geometry equivalence")
    print("ARCHIVED_GBIF_1596_NEW_SOURCE_PORTS_AND_DATASETKEYS_INDEPENDENT_QA_PASS")
    print(json.dumps({"ports":len(ids),"positive":positive,"zero":len(ids)-positive,
                      "new_dataset_keys":len(keys),
                      "sum_port_buffer_indexed_hits_with_overlap":sum(sums),
                      "top1_provider_hit_share":top1,
                      "top10_provider_hit_share":top10,
                      "historical_source_GBIF_datasetkey_claim":"NOT_CLAIMED"},sort_keys=True))
if __name__=="__main__":
    run()
