#!/usr/bin/env python3
"""Conditional GADM1 facets of source-locked juvenile stage labels: 15 queries.

Not verified feeding, host use, wild reproduction, recruitment or colonization.
The queries contain no occurrence record body, coordinates or observer IDs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.gbif.org/v1/occurrence/search"
SOURCE_SHA = "190ff7efabeed547f39d214263a0bde85fe8a3126d191359b188ec44de43f582"
SOURCE_SCHEMA = "chocho_subcountry_48_positive_GBIF_GADM1_lifestage_facet_v1"
SCHEMA = "chocho_juvenile_conditioned_GADM1_source_facet_20261021_v1"
STAGES = frozenset(("Egg", "Larva", "Caterpillar", "Pupa"))
BASIS = ("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
         "MACHINE_OBSERVATION", "OCCURRENCE")
SEED_QUERY_LIMIT = 300


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def cohort(path):
    raw = Path(path).read_bytes()
    # The byte SHA is from the earlier replay artifact. The committed
    # human-readable JSON may differ only in serialization/whitespace;
    # the original source fields must match the locked contract below.
    source = json.loads(raw)
    if source.get("schema") != SOURCE_SCHEMA or source.get("status") != "COMPLETE_48":
        raise ValueError("Not the completed October 20 48-cell GBIF stage facet source")
    if source.get("source_count_exact_cells") != 48 or source.get("source_count_drifted_cells") != 0:
        raise ValueError("Original GBIF source cohort has incomplete or drifting counts")
    if (source.get("indexed_records_across_completed_queries"),
        source.get("life_stage_metadata_records"),
        source.get("establishment_metadata_records")) != (756, 500, 0):
        raise ValueError("Frozen source 756/500/0 records changed")
    records = source["observations"]
    if len(records) != 48:
        raise ValueError("Source positive cell count changed")
    pairs, queries = {}, []
    for r in records:
        key = (r["species"], r["country"])
        if key in pairs or r["status"] != "COMPLETE_GADM1" or r["source_pooled_delta"] != 0:
            raise ValueError("Duplicate or drifting positive cell")
        pairs[key] = r
        stage = r["facets"].get("life_stage", {}).get("categories", {})
        for label, count in sorted(stage.items()):
            if label not in STAGES or not count:
                continue
            if type(count) is not int or count < 0:
                raise ValueError("Original stage counts invalid")
            queries.append({
                "species": r["species"], "country": r["country"],
                "taxon_key": r["taxon_key"], "life_stage": label,
                "expected_stage_count": count,
                "added_host_country": r["added_resource_country"],
            })
    if (len(pairs), len(queries),
        len({(q["species"], q["country"]) for q in queries}),
        sum(q["expected_stage_count"] for q in queries)) != (48, 15, 11, 39):
        raise ValueError("Locked juvenile 48/15/11/39 source totals changed")
    if len({q["species"] for q in queries}) != 3:
        raise ValueError("Source juvenile taxonomic composition changed")
    return queries, pairs, {"source_file_sha256": sha(raw),
                            "archived_original_replay_sha256": SOURCE_SHA,
                            "source_schema": SOURCE_SCHEMA,
                            "frozen_48_cells": 48, "juvenile_pairs": 11,
                            "conditional_queries": 15,
                            "expected_juvenile_label_records": 39}


def make_url(row):
    params = [
        ("taxonKey", str(row["taxon_key"])), ("country", row["country"]),
        ("year", "2010,2025"), ("lifeStage", row["life_stage"]),
        ("hasCoordinate", "true"), ("hasGeospatialIssue", "false"),
        ("occurrenceStatus", "PRESENT"),
        ("limit", "0"), ("facetLimit", str(SEED_QUERY_LIMIT)),
        ("facet", "gadmLevel1Gid")]
    params.extend(("basisOfRecord", x) for x in BASIS)
    return API + "?" + urlencode(params)


def parse_gadm(raw):
    try:
        doc = json.loads(raw)
    except (ValueError, UnicodeDecodeError) as e:
        raise ValueError("GBIF returned invalid JSON") from e
    if not isinstance(doc, dict) or type(doc.get("count")) is not int or doc["count"] < 0:
        raise ValueError("GBIF total is not nonnegative integer")
    if doc.get("limit") != 0 or doc.get("results") != []:
        raise ValueError("GBIF occurrence records were returned unexpectedly")
    facets = doc.get("facets")
    if not isinstance(facets, list) or len(facets) != 1:
        raise ValueError("GBIF GADM facet missing or duplicated")
    f = facets[0]
    if str(f.get("field")).upper().replace("_", "") != "GADMLEVEL1GID":
        raise ValueError("Unexpected GBIF facet name")
    entries = f.get("counts")
    if not isinstance(entries, list) or len(entries) >= SEED_QUERY_LIMIT:
        raise ValueError("Truncated or missing GADM1 counts")
    out = {}
    for r in entries:
        if not isinstance(r, dict):
            raise ValueError("Invalid GBIF GADM1 entry")
        gid, n = r.get("name"), r.get("count")
        if not isinstance(gid, str) or "." not in gid or type(n) is not int or n < 0 or gid in out:
            raise ValueError("Invalid, non-GADM1 or duplicate GBIF GADM1 category")
        out[gid] = n
    if sum(out.values()) > doc["count"]:
        raise ValueError("GBIF GADM1 mapped counts exceed matching stage total")
    return {"record_count": doc["count"], "gadm1": dict(sorted(out.items())),
            "gadm1_mapped": sum(out.values()),
            "gadm1_unmapped": doc["count"] - sum(out.values())}


def fetch(row, retries=3):
    url = make_url(row)
    history = []
    for trial in range(retries):
        try:
            req = Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "chocho-juvenile-gadm-source-conditional/1.0"})
            with urlopen(req, timeout=35) as response:
                status = response.status
                body = response.read(1_000_000)
            if status != 200 or len(body) >= 1_000_000:
                raise ValueError("GBIF HTTP status or facet size invalid")
            x = parse_gadm(body)
            return {**row, **x, "status": (
                "COMPLETE_MATCHED_STAGE" if x["record_count"] == row["expected_stage_count"]
                else "LIFESTAGE_SOURCE_MISMATCH_UNRESOLVED"),
                "http": status, "url": url, "sha256": sha(body),
                "retrieved_utc": utc(), "attempts": trial + 1}
        except (HTTPError, URLError, TimeoutError, OSError, ValueError) as e:
            history.append({"attempt": trial + 1, "error": type(e).__name__,
                            "message": str(e)[:200]})
            if trial < retries - 1:
                time.sleep(2 + trial)
    return {**row, "status": "API_UNAVAILABLE_NOT_ZERO", "url": url,
            "attempts": retries, "errors": history, "retrieved_utc": utc()}


def summarize(items, cells, meta):
    if len(items) != 15 or len({(r["species"], r["country"], r["life_stage"]) for r in items}) != 15:
        raise ValueError("Missing conditional source query")
    accepted = 0
    for r in items:
        if r["status"] != "COMPLETE_MATCHED_STAGE":
            continue
        src = cells[r["species"], r["country"]]
        if r["record_count"] != r["expected_stage_count"]:
            raise ValueError("Accepted result violates frozen stage margin")
        bins = src["facets"]["gadm1"]["categories"]
        for name, n in r["gadm1"].items():
            if n > bins.get(name, 0):
                r["status"] = "GADM_SUBSET_NOT_NESTED_UNRESOLVED"
                break
        if r["status"] == "COMPLETE_MATCHED_STAGE":
            accepted += 1
    valid = [r for r in items if r["status"] == "COMPLETE_MATCHED_STAGE"]
    grouped = defaultdict(list)
    for r in valid:
        grouped[r["species"], r["country"]].append(r)
    cell_summary = []
    for key, qs in sorted(grouped.items()):
        source = cells[key]
        stage_total = sum(r["record_count"] for r in qs)
        mapped = sum(r["gadm1_mapped"] for r in qs)
        # Different stage labels are mutually exclusive source facet values.
        bins = Counter()
        for q in qs:
            bins.update(q["gadm1"])
        mismatch = {k: v for k, v in bins.items() if v > source["facets"]["gadm1"]["categories"].get(k, 0)}
        cell_summary.append({
            "species": key[0], "country": key[1],
            "added_host_country": source["added_resource_country"],
            "juvenile_stage_labeled": stage_total,
            "mapped_juvenile_records": mapped,
            "juvenile_gadm1_units": len(bins),
            "all_life_stage_gadm1_units": source["mapped_gadm_units"],
            "juvenile_gadm1_counts": dict(sorted(bins.items())),
            "cross_category_gadm_mismatch": mismatch,
            "status": "VALID_GADM_STAGE_SUBSET" if not mismatch else "MULTISTAGE_NOT_NESTED"
        })
    full_complete = accepted == 15 and len(cell_summary) == 11 and all(
        c["status"] == "VALID_GADM_STAGE_SUBSET" for c in cell_summary)
    crosstab = {}
    for added in (0, 1):
        rows = [x for x in cell_summary if x["added_host_country"] == added]
        crosstab[str(added)] = {
            "juvenile_cells": len(rows),
            "stage_labeled_juvenile_records": sum(x["juvenile_stage_labeled"] for x in rows),
            "mapped_juvenile_records": sum(x["mapped_juvenile_records"] for x in rows)
        }
    status = "COMPLETE_CONDITIONAL_JUVENILE_GADM1" if full_complete else "INCOMPLETE_STAGE_GADM_DATA_NOT_ZERO"
    return {"schema": SCHEMA, "status": status, "executed_utc": utc(),
            "source": meta, "attempts": items, "matched_conditionals": accepted,
            "juvenile_stage_cells_resolved": len(cell_summary),
            "all_frozen_juvenile_stage_records": 39,
            "resolved_juvenile_records": sum(x["juvenile_stage_labeled"] for x in cell_summary),
            "gadm1_mapped_juvenile_records": sum(x["mapped_juvenile_records"] for x in cell_summary),
            "cell_gadm1_juvenile_summary": cell_summary,
            "by_source_country_added_resource_flag": crosstab,
            "no_verified_plant_host_or_breeding": True}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--frozen-positive48", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    queries, pairs, meta = cohort(args.frozen_positive48)
    output = []
    for row in queries:
        output.append(fetch(row))
        time.sleep(1.0)
    result = summarize(output, pairs, meta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("GBIF_JUVENILE_GADM_CONDITIONAL="+json.dumps(
        {k:v for k,v in result.items() if k!="attempts"}, separators=(",", ":"), flush=True))
    if result["status"] != "COMPLETE_CONDITIONAL_JUVENILE_GADM1":
        raise SystemExit("Juvenile facet query conditional mismatch, preserve but do not interpret")


if __name__ == "__main__":
    main()
