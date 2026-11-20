#!/usr/bin/env python3
"""Re-audit original source-snapshot 225 butterfly × 51 island GBIF facet index.

NO GBIF requery is performed and NO individual record/Darwin Core record or
contributing datasetKey is silently inferred from aggregate facet counts.

Reads an immutable checkout of commit 365526cc70e57b8a6634ee29494bd678f95961ff,
which contains all six original 2026-11-09 source shards, query receipts, and
the 11,475-row merged index. Outputs a citation-readiness report, not a DOI.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

PINNED_SOURCE_COMMIT = "365526cc70e57b8a6634ee29494bd678f95961ff"
BASE = Path("data/exploratory/public51_butterfly_index_20261109")
EXPECTED_MERGED_SHA256 = "d8fbb3fbf936bfdbedc6c91f1e19f11deda882d3125ccb6600a9c8d800a57887"
EXPECTED_SOURCE_ISLAND_BLOB = "18cecde5bf3066f276380167f451f1fa5cc1caa1"
EXPECTED_SOURCE_SPECIES_BLOB = "39fe84981cab2b0e47063399f0a5c0c9731e800d"
EXPECTED_KEYS = (
    "species", "speciesKey", "landmass_id",
    "pre2001_2010_indexed_records", "future2011_2017_indexed_records",
    "future2018_2022_indexed_records", "eligible_2011_2017",
    "first_indexed_in_2011_2017", "eligible_2018_2022",
    "first_indexed_in_2018_2022", "source_not_future_selected",
    "index_not_true_absence_or_colonization"
)
PERIODS = {"pre2001_2010","first2011_2017","holdout2018_2022"}
BLOCKED_RECORD_FIELDS = {"gbifid","occurrenceid","datasetkey","datasetdoi"}
SHARD_SCHEMA = "chocho_public_GBIF_unselected_225_butterfly_x_original51_island_window_v1"

def load_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rd = csv.DictReader(stream)
        columns = tuple(rd.fieldnames or ())
        return columns, list(rd)

def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def find_record_fields(value, loc="root"):
    """Search actual keys, not cautionary prose mentioning datasetKey."""
    hits=[]
    if isinstance(value, dict):
        for k,v in value.items():
            if k.lower().replace("_","") in BLOCKED_RECORD_FIELDS:
                hits.append(loc + "/" + k)
            hits.extend(find_record_fields(v, loc + "/" + k))
    elif isinstance(value, list):
        for i,v in enumerate(value):
            hits.extend(find_record_fields(v, loc + "/" + str(i)))
    return hits

def audit(source_root: Path):
    root=source_root / BASE
    if not root.is_dir():
        raise FileNotFoundError("Original public51 source branch not available: " + str(root))
    merged=root/"full51_225butterfly_period_index.csv"
    doc=root/"full51_source_validation.json"
    fields, rows=load_csv(merged)
    if fields!=EXPECTED_KEYS or sha(merged)!=EXPECTED_MERGED_SHA256:
        raise AssertionError("Merged source population/schema/hash changed")
    if len(rows)!=11475 or len({(r["speciesKey"],r["landmass_id"]) for r in rows})!=11475:
        raise AssertionError("Missing or duplicate fixed 225×51 source opportunity cells")
    sp={r["speciesKey"]:r["species"] for r in rows}
    isl={r["landmass_id"] for r in rows}
    if len(sp)!=225 or len(isl)!=51:
        raise AssertionError("Frozen 225-butterfly/51-island source sample drift")
    if any(not r["source_not_future_selected"].lower()=="true" or
           not r["index_not_true_absence_or_colonization"].lower()=="true" for r in rows):
        raise AssertionError("An index source has been relabelled as biological arrival")
    if any(int(r[k])<0 for r in rows for k in fields[3:6]):
        raise AssertionError("Negative original indexed-record counts")
    obs={
       "eligible_2011_2017":sum(int(r["eligible_2011_2017"]) for r in rows),
       "first_indexed_in_2011_2017":sum(int(r["first_indexed_in_2011_2017"]) for r in rows),
       "eligible_2018_2022":sum(int(r["eligible_2018_2022"]) for r in rows),
       "first_indexed_in_2018_2022":sum(int(r["first_indexed_in_2018_2022"]) for r in rows)
    }
    expected_totals={"eligible_2011_2017":10739,"first_indexed_in_2011_2017":221,
                     "eligible_2018_2022":10518,"first_indexed_in_2018_2022":96}
    if obs != expected_totals:
        raise AssertionError("Original 2026 indexed event window totals drift")
    receipt=json.loads(doc.read_text(encoding="utf-8"))
    if (receipt["cells"],receipt["source_islands"],receipt["source_butterflies"])!=(11475,51,225):
        raise AssertionError("Original merged source receipt changed")
    if receipt["source_island_git_blob"]!=EXPECTED_SOURCE_ISLAND_BLOB or \
       receipt["source_butterflies_git_blob"]!=EXPECTED_SOURCE_SPECIES_BLOB or \
       receipt["merged_table_sha256"]!=EXPECTED_MERGED_SHA256 or \
       receipt["totals"]!=obs:
        raise AssertionError("Original merged query receipt not source-frozen")
    shards=receipt["shard_audits"]
    if len(shards)!=6:
        raise AssertionError("Expected six independent historic queries")
    seen_islands=set()
    values=[]
    union_rows={}
    query_count=0
    snapshot_times=[]
    filters=set()
    for s in shards:
        prefix=s["prefix"]
        csvpath=root/"source_shards"/(prefix+"_butterfly_island_225_atrisk.csv")
        jsonpath=root/"source_shards"/(prefix+"_query_receipt.json")
        if sha(csvpath)!=s["csv_sha256"] or sha(jsonpath)!=s["receipt_sha256"]:
            raise AssertionError("One of six immutable query-shard SHA256s changed")
        columns, shardrows=load_csv(csvpath)
        if columns!=fields or len(shardrows)!=s["risk_cells"]:
            raise AssertionError("Source shard unexpected schema/row population")
        j=json.loads(jsonpath.read_text(encoding="utf-8"))
        if (j["schema"]!=SHARD_SCHEMA or
            j["data_source"]!="https://api.gbif.org/v1/occurrence/search" or
            j["risk_cells"]!=s["risk_cells"] or
            j["species_candidates_frozen_before_future_225"]!=225):
            raise AssertionError("Source no longer original unselected GBIF facet query")
        if j.get("snapshot_timestamp_utc")!=s["index_timestamp_utc"]:
            raise AssertionError("Original GBIF query timestamp mismatch")
        ids=set(map(str,j["islands_tested"]))
        if seen_islands & ids or ids!=set(map(str,s["ids"])):
            raise AssertionError("Source island present in multiple query shards")
        seen_islands |= ids
        if len(j["GBIF_index_query_receipts"])!=len(ids)*3:
            raise AssertionError("Query receipt missing source island×period facets")
        qkeys={(str(q["island_id"]),str(q["period"])) for q in j["GBIF_index_query_receipts"]}
        if qkeys != {(i,p) for i in ids for p in PERIODS}:
            raise AssertionError("Missing one of the original three-index periods")
        if j["first_indexed_2011_2017"]!=s["indexed_new_2011_2017"] or \
           j["first_indexed_2018_2022"]!=s["indexed_new_2018_2022"]:
            raise AssertionError("Original query facet positive count mismatch")
        if j["butterfly_taxon_filter"]["method"] not in (
           "GBIF_SUPERFAMILY_TAXONKEY","REPEATED_GBIF_FAMILY_TAXONKEY_OR_FILTER"):
            raise AssertionError("Unknown source taxon filter")
        filters.add(json.dumps(j["butterfly_taxon_filter"],sort_keys=True))
        for row in shardrows:
            key=(row["speciesKey"],row["landmass_id"])
            if row["landmass_id"] not in ids or key in union_rows:
                raise AssertionError("Shard contains duplicated or foreign island")
            union_rows[key]=row
        query_count+=len(j["GBIF_index_query_receipts"])
        snapshot_times.append(j["snapshot_timestamp_utc"])
        values.append({"batch":prefix,"source_islands":len(ids),
                       "cells":len(shardrows),"source_query_receipts":len(j["GBIF_index_query_receipts"]),
                       "snapshot_utc":j["snapshot_timestamp_utc"],
                       "source_csv_sha256":s["csv_sha256"],"source_receipt_sha256":s["receipt_sha256"]})
    if seen_islands!=isl or len(union_rows)!=len(rows) or len(filters)!=1:
        raise AssertionError("Not identical source population or taxon filter across shards")
    if any(union_rows[(r["speciesKey"],r["landmass_id"])]!=r for r in rows):
        raise AssertionError("Independent original shard reconstruction does not match merged source")
    if find_record_fields(receipt) or any(k.lower().replace("_","") in BLOCKED_RECORD_FIELDS for k in fields):
        raise AssertionError("Unexpected individual-record identifiers in a claimed facet-only source")
    if len(snapshot_times)!=6 or not receipt["six_shards_not_simultaneous_snapshot"]:
        raise AssertionError("Missing original index time heterogeneity warning")
    drift=receipt["original17_reindexed_snapshot_drift_audit"]
    if drift["comparable_original_species_island_cells"]!=3825 or \
       drift["changed_binary_first_indexing"]!={"first_indexed_in_2011_2017":1,
                                               "first_indexed_in_2018_2022":0}:
        raise AssertionError("Original 2026 GBIF index drift audit changed")
    return {
        "schema":"chocho_public51_archive_GBIF_facet_attribution_audit_20261120_v1",
        "status":"SIX_ORIGINAL_FACET_SHARDS_SOURCE_EXACT_NO_RECORD_PUBLISHERS",
        "pinned_git_commit":PINNED_SOURCE_COMMIT,
        "source_endpoint":"GBIF occurrence/search facet=speciesKey",
        "source_date":"2026-11-09",
        "merged_sha256":EXPECTED_MERGED_SHA256,
        "source_population":{"cells":len(rows),"butterflies":len(sp),"islands":len(isl),
                             "query_periods_per_island":3,"total_source_query_receipts":query_count},
        "source_window_count":obs,
        "shards":values,
        "distinct_source_query_times":len(set(snapshot_times)),
        "original17_index_snapshot_drift":drift,
        "verified":"Six source CSVs and six query receipts independently equal the merged table and the source-pinned original recorded checksums",
        "cannot_recover_from_facets":["original occurrence gbifID","original datasetKey/DOI",
           "original publisher and provider-specific counts","original per-record precision"],
        "DOI_gate":"NOT_MET: GBIF source publisher identities missing; no occurrence download or derived DOI minted",
        "scope_limit":"This 225 candidate×51 island 2026 re-index is not the distinct 224-butterfly historical private Parquet or observed-community 224-butterfly published GBIF DOI",
        "scientific_limit":"SpeciesKey facet counts are indexed occurrences, not true first arrival, breeding, or real absences"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    j=audit(args.source_root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(j,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PUBLIC51_GBIF_FACET_SIX_SHARDS_REPLAY_PASS")
    print(json.dumps({k:j[k] for k in ("status","source_population","DOI_gate")},ensure_ascii=False))

if __name__=="__main__":
    main()
