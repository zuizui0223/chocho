#!/usr/bin/env python3
"""Small, source-locked provenance facet audit for 7 GBIF country-detection cells.

Does NOT establish butterfly colonization, wild reproduction or local host use.
Only reads 3 group-by facets per frozen focal taxon/country, no occurrence bodies.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.gbif.org/v1/occurrence/search"
BASIS = ("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
         "MACHINE_OBSERVATION", "OCCURRENCE")
FACETS = ("datasetKey", "year", "basisOfRecord")
ORIGINAL_SHA = "18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2"
COUNTRY_ARCHIVE_SHA = "d9f5207cd611e737279b7b579fe8f2347132d68bab36fc78a56e727f498c0359"
FACET_LIMIT = 300


def sha(data):
    return hashlib.sha256(data).hexdigest()


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def load_candidates(temporal, original, country_counts):
    original_bytes = Path(original).read_bytes()
    if sha(original_bytes) != ORIGINAL_SHA or sha(Path(country_counts).read_bytes()) != COUNTRY_ARCHIVE_SHA:
        raise ValueError("Frozen original names or country facet CSV changed")
    names = json.loads(original_bytes)
    resolved = {row["species"]: int(row["taxon"]["key"])
                for row in names["species_receipts"] if row["status"] != "TAXON_UNRESOLVED"}
    if len(resolved) != 231:
        raise ValueError("Expected 231 original resolved names")
    raw = json.loads(Path(temporal).read_text(encoding="utf-8"))
    if raw["source"]["country_counts_sha256"] != COUNTRY_ARCHIVE_SHA:
        raise ValueError("Temporal source has a different GBIF country source")
    if raw["source"]["original_manifest_sha256"] != ORIGINAL_SHA:
        raise ValueError("Temporal taxonomy receipt changed")
    depth = raw["record_depth_posthoc_diagnostics"]["2"]
    result = []
    for row in depth["positive_species_country_details"]:
        sp, ct = row["species"], row["country"]
        if row["early_record_count"] != 0 or row["late_record_count"] < 2:
            raise ValueError("Not an early-zero / late >=2 pair")
        if sp not in resolved:
            raise ValueError(f"Unresolved or changed source name: {sp}")
        if len(ct) != 2 or not ct.isupper():
            raise ValueError("Non-country ISO code")
        result.append({"species": sp, "country": ct, "taxon_key": resolved[sp],
                       "frozen_late_records": int(row["late_record_count"]),
                       "frozen_early_records": int(row["early_record_count"]),
                       "added_host_country": int(row["added_host_country"]),
                       "other_taxa_early": int(row["other_taxa_early"]),
                       "other_taxa_late": int(row["other_taxa_late"])})
    if len(result) != 7 or len({(x["species"], x["country"]) for x in result}) != 7:
        raise ValueError("Expected seven unique paired-effort late-only >=2 cells")
    if sum(x["added_host_country"] for x in result) != 3:
        raise ValueError("Frozen added-host exposure cohort changed")
    if sorted(x["frozen_late_records"] for x in result if x["frozen_late_records"] >= 5) != [9, 86]:
        raise ValueError("Frozen high-count cohort changed")
    return sorted(result, key=lambda r: (r["species"], r["country"]))


def query_url(pair):
    params = [
        ("taxonKey", str(pair["taxon_key"])), ("country", pair["country"]),
        ("year", "2018,2025"), ("hasCoordinate", "true"),
        ("hasGeospatialIssue", "false"), ("occurrenceStatus", "PRESENT"),
        ("limit", "0"), ("facetLimit", str(FACET_LIMIT))
    ]
    params += [("basisOfRecord", basis) for basis in BASIS]
    params += [("facet", facet) for facet in FACETS]
    return API + "?" + urlencode(params)


def normalized_field(raw):
    return str(raw).upper().replace("_", "")


def validate_response(raw):
    try:
        value = json.loads(raw)
    except (ValueError, UnicodeDecodeError) as exc:
        raise ValueError("Invalid/truncated GBIF JSON body") from exc
    count = value.get("count")
    if type(count) is not int or count < 0 or value.get("limit") != 0 or value.get("results") != []:
        raise ValueError("Expected a nonnegative count and limit=0 facets without record bodies")
    facets = value.get("facets")
    if not isinstance(facets, list):
        raise ValueError("Missing GBIF facets")
    fields = {}
    required = {normalized_field(field) for field in FACETS}
    for row in facets:
        if not isinstance(row, dict):
            raise ValueError("Malformed facet")
        field = normalized_field(row.get("field"))
        if field not in required or field in fields:
            raise ValueError("Unknown or duplicate GBIF facet: " + field)
        entries = row.get("counts")
        if not isinstance(entries, list) or len(entries) >= FACET_LIMIT:
            raise ValueError("Missing/truncated facet entries: " + field)
        values = {}
        for entry in entries:
            if not isinstance(entry, dict):
                raise ValueError("Malformed facet count")
            name, n = entry.get("name"), entry.get("count")
            if not isinstance(name, str) or type(n) is not int or n < 0 or name in values:
                raise ValueError("Invalid or duplicate facet category")
            if field == "YEAR" and (not name.isdigit() or not 2018 <= int(name) <= 2025):
                raise ValueError("Year outside source window")
            values[name] = n
        if sum(values.values()) != count:
            raise ValueError(f"Incomplete or mismatched {field} facet totals")
        fields[field] = values
    if set(fields) != required:
        raise ValueError("Missing one or more requested GBIF facets")
    return {"record_count": count, "facet_counts": fields,
            "dataset_count": len(fields["DATASETKEY"]),
            "years_recorded": sorted(fields["YEAR"]),
            "basis_types_recorded": sorted(fields["BASISOFRECORD"]),
            "largest_dataset_fraction": max(fields["DATASETKEY"].values()) / count if count else None,
            "largest_year_fraction": max(fields["YEAR"].values()) / count if count else None}


def fetch_one(pair, *, max_attempts=3):
    url = query_url(pair)
    errors = []
    for attempt in range(max_attempts):
        request = Request(url, headers={
            "User-Agent": "chocho-gbif-country-provenance-audit/1.0",
            "Accept": "application/json"})
        try:
            with urlopen(request, timeout=35) as response:
                status = response.status
                body = response.read(1_000_000)
            if status != 200 or len(body) >= 1_000_000:
                raise ValueError(f"GBIF HTTP status/response length invalid: {status}/{len(body)}")
            parsed = validate_response(body)
            return {**pair, "status": "COMPLETE_LIVE_FACET_QC",
                    "retrieved_utc": utcnow(), "source_query_url": url,
                    "raw_response_sha256": sha(body),
                    "source_drift_records": parsed["record_count"] - pair["frozen_late_records"],
                    **parsed, "attempts": attempt + 1}
        except (HTTPError, URLError, TimeoutError, ValueError, OSError) as exc:
            errors.append({"attempt": attempt + 1, "error_type": type(exc).__name__,
                           "message": str(exc)[:400]})
            if attempt < max_attempts - 1:
                time.sleep(3 * (attempt + 1))
    return {**pair, "status": "LIVE_FACET_QC_UNAVAILABLE_NOT_ZERO",
            "source_query_url": url, "retrieved_utc": utcnow(),
            "attempts": max_attempts, "errors": errors}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--temporal-summary", type=Path, required=True)
    p.add_argument("--original-manifest", type=Path, required=True)
    p.add_argument("--frozen-country-counts", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    cohort = load_candidates(args.temporal_summary, args.original_manifest,
                             args.frozen_country_counts)
    print("LOCKED_CANDIDATES=" + json.dumps(cohort, separators=(",", ":")), flush=True)
    results = []
    for pair in cohort:
        results.append(fetch_one(pair))
        time.sleep(1.0)
    doc = {
        "schema": "gbif_country_record_provenance_facet_qc_v1",
        "source_protocol": "provenance/reviewer_defenses/protocols/GBIF_COUNTRY_SOURCE_INDEPENDENCE_GATE_20261020.md",
        "executed_utc": utcnow(), "frozen_candidates": len(cohort),
        "quality_audit_only_not_establishment": True,
        "live_facets_complete": sum(x["status"] == "COMPLETE_LIVE_FACET_QC" for x in results),
        "results": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("GBIF_COUNTRY_FACET_QC_SUMMARY=" + json.dumps(doc, separators=(",", ":")), flush=True)
    if doc["live_facets_complete"] != 7:
        raise SystemExit("LIVE_FACET_QC_INCOMPLETE; inspect receipts, do not infer zeros")


if __name__ == "__main__":
    main()
