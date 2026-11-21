#!/usr/bin/env python3
"""Consolidate four ALREADY SUCCESSFUL 2026-11-20 GBIF Plantae provider-facet shards.

This is a pure stdlib, no-network, source-exact reconstruction from original
source-complete WPI2019 port IDs and four immutable Actions artifacts.
The current GBIF datasetKey index DOES NOT reconstruct the publisher
contributions for the historical original speciesKey-based 1596-port query,
and DOES NOT license a retrospectively assigned GBIF DOI.

No original private WPI GeoJSON, GBIF occurrence ID, or author identity
is included in the compact long-term source receipt.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ORIGINAL_PORTS = Path("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv")
ORIGINAL_SOURCE_BLOB = "fdf73b11ddcb1644f46730123350ae52cc61cae6"
ORIGINAL_GIT_COMMIT = "3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250"
SOURCE_WORKFLOW_RUN = 44742493164
EXPECTED_SCHEMA = "chocho_new_20261120_all1596_WPI_plants_GBIF_datasetkey_facets_shard_v1"

def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def git_blob(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def iso_utc(value: str) -> datetime:
    dt = datetime.fromisoformat(value.replace("Z","+00:00"))
    if dt.tzinfo is None:
        raise AssertionError("Original GBIF query shard has naive retrieval timestamp")
    if dt.date().isoformat() != "2026-11-20":
        raise AssertionError("New source must be the explicitly dated 2026-11-20 query, not an unrelated snapshot")
    return dt

def audit(source_root: Path, shard_dir: Path) -> dict:
    data = (source_root / ORIGINAL_PORTS).read_bytes()
    if git_blob(data) != ORIGINAL_SOURCE_BLOB:
        raise AssertionError("Original source-locked 1596 port roster drift")
    source = list(csv.DictReader(data.decode("utf-8-sig").splitlines()))
    original = sorted(int(r["port_id"]) for r in source if r["query_complete"] == "True")
    if len(source) != 3630 or len(original) != 1596 or len(set(original)) != 1596:
        raise AssertionError("Original 3630/1596 source cohort changed")

    shards = sorted(shard_dir.rglob("shard_new_provider_*.json"))
    if len(shards) != 4:
        raise AssertionError("Expected four complete source shard JSONs: " + str([str(p) for p in shards]))
    allrows = []
    shard_receipts = []
    byprovider = defaultdict(int)
    byprovider_ports = defaultdict(set)
    shardparts = set()
    snapshot_start, snapshot_end = [], []
    for p in shards:
        raw=p.read_bytes()
        r=json.loads(raw)
        if r.get("schema") != EXPECTED_SCHEMA or r.get("status") != "SOURCE_COMPLETED_ONE_OF_FOUR_NEW_PUBLIC_GBIF_PROVIDER_SHARDS_NOT_ORIGINAL_2026_QUERY":
            raise AssertionError("Unexpected source shard schema/status")
        part=int(r["source_shard_index"])
        if part in shardparts or not 0<=part<4 or int(r["source_shard_total"])!=4:
            raise AssertionError("Duplicated/unrecognized original 4-way WPI shard")
        shardparts.add(part)
        if (r["original_git_source_commit"] != ORIGINAL_GIT_COMMIT or
            r["original_all_ports"] != 3630 or r["original_complete_ports"] != 1596 or
            not r["source_original_private_WPI2019_location_match_6dp"] or
            not r["new_provider_data_not_original_citation"] or
            not r["original_1596_cache_WKT_geometry_not_exactly_equivalent_to_new_geoDistance"] or
            r["historical_2026_GBIF_publisher_datasetKeys_recovered"] or
            r["new_GBIF_derived_DOI_minted"] or
            r["new_source_query_api"] != "https://api.gbif.org/v1/occurrence/search" or
            r["query_years"] != "2001,2024" or r["query_radius"] != "5km"):
            raise AssertionError("One current provider-facet shard has wrong original source/DOI/geometry scope")
        begin=iso_utc(r["new_index_timestamp_start_utc"])
        end=iso_utc(r["new_index_timestamp_end_utc"])
        if begin>end:
            raise AssertionError("Chronology of source acquisition is inconsistent")
        snapshot_start.append(begin)
        snapshot_end.append(end)
        rows = r["new_datasetKey_facets"]
        expect_ids=set(original[part::4])
        keys={int(x["original_WPI2019_OBJECTID"]) for x in rows}
        if len(rows)!=399 or len(keys)!=399 or keys!=expect_ids or r["n_port_queries"]!=399:
            raise AssertionError("Shard port roster is not complete and source-original")
        shard_total=0
        unique=set()
        for rec in rows:
            pid=int(rec["original_WPI2019_OBJECTID"])
            total=int(rec["new_2026_GBIF_index_count"])
            pages=int(rec["facet_pages"])
            providers=rec["current_GBIF_record_counts_by_datasetKey"]
            if total<0 or pages<1 or len(providers)!=int(rec["new_datasetKey_count"]):
                raise AssertionError("Invalid new source facet count")
            if total>0 and not providers:
                raise AssertionError("Indexed positive port without a publisher")
            if not str(rec["new_source_query_geoDistance"]).endswith(",5km"):
                raise AssertionError("One source query no longer used original 5km port centre")
            present=set()
            value=0
            for d in providers:
                key=str(d["datasetKey"])
                cnt=int(d["count"])
                if len(key)!=36 or key in present or cnt<0:
                    raise AssertionError("Duplicate/invalid provider name or nonnegative count")
                present.add(key)
                unique.add(key)
                value+=cnt
                byprovider[key]+=cnt
                byprovider_ports[key].add(pid)
            if total!=value:
                raise AssertionError("One port's publisher-source indexed counts do not reconcile")
            shard_total+=total
            allrows.append({
                "original_WPI2019_OBJECTID":pid,
                "new_2026_GBIF_index_count":total,
                "new_datasetKey_count":len(providers),
                "current_GBIF_record_counts_by_datasetKey":providers,
                "facet_pages":pages
            })
        if shard_total!=int(r["new_indexed_port_buffers_sum_not_unique_occurrence_records"]) or len(unique)!=int(r["new_unique_datasetKey_in_shard"]):
            raise AssertionError("Original independently audited shard total/unique source dataset counts changed")
        shard_receipts.append({
            "part":part,
            "source_workflow_run":SOURCE_WORKFLOW_RUN,
            "sha256":sha(raw),
            "original_port_ids":399,
            "new_current_GBIF_port_indexed_hit_sum_WITH_OVERLAP":shard_total,
            "new_unique_datasetKeys_in_shard":len(unique),
            "new_index_start_utc":r["new_index_timestamp_start_utc"],
            "new_index_end_utc":r["new_index_timestamp_end_utc"]
        })
    ids = [x["original_WPI2019_OBJECTID"] for x in allrows]
    if shardparts!={0,1,2,3} or len(ids)!=1596 or len(set(ids))!=1596 or set(ids)!=set(original):
        raise AssertionError("Consolidated 1596 original complete WPI port source population absent/duplicated")
    allrows.sort(key=lambda x:x["original_WPI2019_OBJECTID"])
    total=sum(row["new_2026_GBIF_index_count"] for row in allrows)
    positive=sum(row["new_2026_GBIF_index_count"]>0 for row in allrows)
    if sum(byprovider.values())!=total:
        raise AssertionError("Consolidated provider-source sum differs from indexed port hit sum")
    providers=[{
        "datasetKey":key,
        "sum_port_buffer_index_hits_NOT_UNIQUE_OCCURRENCES":n,
        "number_source_port_buffers_indexing_dataset":len(byprovider_ports[key])
    } for key,n in sorted(byprovider.items())]
    if len(providers)<1 or total<=0:
        raise AssertionError("No positive GBIF source providers")
    report={
       "schema":"chocho_all1596_original_WPI2019_ports_NEW_current_GBIF_datasetKey_facets_20261120_v2",
       "status":"COMPLETE_1596_NEW_CURRENT_SOURCE_PORT_BUFFER_PROVIDERS_NOT_ORIGINAL_HISTORICAL_2026_GBIF",
       "source_original_Git_commit":ORIGINAL_GIT_COMMIT,
       "source_original_workflow_run":SOURCE_WORKFLOW_RUN,
       "source_original_complete_port_IDs_Git_blob":ORIGINAL_SOURCE_BLOB,
       "source_original_WPI_candidates":3630,
       "source_original_WPI_plant_queries_complete":1596,
       "source_private_original_2019_WPI_centres_same_to_six_decimals":True,
       "source_cached_1596_original_GBIF_query_used_WKT_buffer_not_exactly_same_as_current_geoDistance":True,
       "source_original_2026_GBIF_contributing_datasetKeys_recovered":False,
       "source_original_2026_plant_GBIF_derived_DOI_registered":False,
       "new_source_query_years":"2001-2024",
       "new_source_query_central_WPI_geoDistance":"5km",
       "new_source_query":"GBIF occurrence/search; Plantae taxonKey6; PRESENT; hasCoordinate=true; hasGeospatialIssue=false; facet=datasetKey",
       "new_source_starts_utc":min(snapshot_start).isoformat(),
       "new_source_ends_utc":max(snapshot_end).isoformat(),
       "new_source_shard_receipts":sorted(shard_receipts,key=lambda x:x["part"]),
       "new_source_original_WPI_portIDs_queried":len(ids),
       "new_source_port_buffers_with_positive_GBIF_index":positive,
       "new_source_port_buffers_with_zero_GBIF_index":1596-positive,
       "new_source_datasetKey_contributors":len(providers),
       "new_source_sum_port_buffer_indexed_hits_WITH_OVERLAP":total,
       "new_source_provider_index_counts_sum_WITH_OVERLAP":sum(byprovider.values()),
       "new_source_provider_counts":providers,
       "new_source_port_index_facets":allrows,
       "audit_no_missing_duplicate_port_ID":True,
       "audit_individual_port_and_four_shard_provider_count_exact":True,
       "review_state":"AUDITED_NEW_SOURCE_SNAPSHOT_NOT_ANONYMIZED_OR_JOURNAL_DATA_DEPOSIT",
       "limits":[
        "The original 2026 port speciesKey faceted query stored no occurrence gbifID or publisher datasetKey and cannot be retrospectively attributed with the new snapshot.",
        "NEW current WPI centre geoDistance 5km is not necessarily the same boundary as original WKT5km buffered search; port source is a partial geographically ordered 1596/3630 candidate selection.",
        "Four queries ran at separate UTC times and GBIF indexing is mutable, so this is not a simultaneous original 2026 GBIF snapshot.",
        "Adjacent port search buffers can overlap; summing port-buffer source counts duplicates occurrences; total is NOT unique plant observations.",
        "Counts measure public GBIF indexed botanical observations, NOT abundance, true introduction timing, observed larval host feeding, vessels, or butterfly colonisation.",
        "This source does not inherit 10.15468/dl.umqmsf, the original 155/224 butterfly download DOI; no new derived DOI has been registered.",
        "Public source port coordinates and datasetKey counts must still be checked for archive rights and review anonymity before publication."
       ]
    }
    return report

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--original-root",required=True,type=Path)
    arg.add_argument("--shards-root",required=True,type=Path)
    arg.add_argument("--out",required=True,type=Path)
    args=arg.parse_args()
    r=audit(args.original_root,args.shards_root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("NEW_GBIF_COMPLETE_1596_PORTS_CONSOLIDATED_SOURCE_EXACT",
          json.dumps({"all_ports":r["new_source_original_WPI_portIDs_queried"],
              "positive_ports":r["new_source_port_buffers_with_positive_GBIF_index"],
              "unique_datasetKeys":r["new_source_datasetKey_contributors"],
              "overlap_count_sum":r["new_source_sum_port_buffer_indexed_hits_WITH_OVERLAP"],
              "historical_DOI_gate":r["source_original_2026_GBIF_contributing_datasetKeys_recovered"]}))
if __name__=="__main__":
    main()
