#!/usr/bin/env python3
"""Cross-tab all 11 dataset/year series for 7 GBIF country-level butterfly detections.

No raw occurrence IDs, coordinates, observer names, host-use or colonization
claims. This is a *posthoc source-provenance check* of a frozen ecological panel.
See GBIF_COUNTRY_DATASET_YEAR_CONFOUND_20261020.md before interpreting output.
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
PROTOCOL = "provenance/reviewer_defenses/protocols/GBIF_COUNTRY_DATASET_YEAR_CONFOUND_20261020.md"
YEAR_RANGE = tuple(str(y) for y in range(2018, 2026))
MATRIX_SIZE = 11
PAIR_SIZE = 7
SCHEMA = "gbif_country_seven_candidate_dataset_by_year_qc_v1"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def load_fixed_pairs(source_qc: Path, registry: Path):
    source_bytes = source_qc.read_bytes()
    src = json.loads(source_bytes)
    reg = json.loads(registry.read_bytes())
    if src.get("schema") != "gbif_country_record_provenance_facet_qc_v1":
        raise ValueError("Unexpected frozen source provenance schema")
    if src.get("live_facets_complete") != PAIR_SIZE or len(src.get("results", [])) != PAIR_SIZE:
        raise ValueError("The prior seven-cell source audit is incomplete")
    if reg.get("schema") != "gbif_country_dataset_registry_followup_posthoc_v1":
        raise ValueError("Unexpected GBIF registry schema")
    if reg.get("source_qc_sha256") != sha(source_bytes):
        raise ValueError("Dataset registry is tied to a different source-QC receipt")
    registry_by_key = {r["dataset_key"]: r for r in reg["dataset_registry"]}
    if len(registry_by_key) != 6 or reg.get("n_complete") != 6:
        raise ValueError("Six dataset keys must have complete registry metadata")
    pairs = []
    source_keys = set()
    for row in src["results"]:
        if row["status"] != "COMPLETE_LIVE_FACET_QC" or row["source_drift_records"] != 0:
            raise ValueError("Frozen seven-cell source audit includes a missing or drifted result")
        if row["frozen_early_records"] != 0 or row["frozen_late_records"] < 2:
            raise ValueError("Not the original non-singleton early-zero cohort")
        dataset_counts = row["facet_counts"]["DATASETKEY"]
        year_counts = row["facet_counts"]["YEAR"]
        if sum(dataset_counts.values()) != row["frozen_late_records"] or sum(year_counts.values()) != row["frozen_late_records"]:
            raise ValueError("Frozen dataset/year marginals changed")
        pair_id = row["species"], row["country"]
        if pair_id in source_keys:
            raise ValueError("Repeated focal country cell")
        source_keys.add(pair_id)
        for key, count in dataset_counts.items():
            if key not in registry_by_key or count <= 0:
                raise ValueError("Unknown or non-positive dataset key in candidate cell")
            pairs.append({
                "species": row["species"], "country": row["country"],
                "taxon_key": row["taxon_key"], "dataset_key": key,
                "source_dataset_count": int(count),
                "source_country_total": int(row["frozen_late_records"]),
                "source_country_years": dict(year_counts),
                "added_host_country": int(row["added_host_country"]),
                "dataset_title": registry_by_key[key]["title"],
                "dataset_doi": registry_by_key[key]["doi"],
                "dataset_created_utc": registry_by_key[key]["created"],
                "publishing_organization_key": registry_by_key[key]["publishing_organization_key"],
            })
    if len(pairs) != MATRIX_SIZE or len(source_keys) != PAIR_SIZE:
        raise ValueError("Original 7 country-cells must contain exactly 11 dataset-specific combinations")
    expected_high = {
        ("Aguna asander", "DO"): {"50c9509d-22c7-4a22-a47d-8c48425ef4a7": 9},
        ("Papilio demoleus", "IL"): {
            "50c9509d-22c7-4a22-a47d-8c48425ef4a7": 43,
            "6c8d87be-829d-4364-b2bc-578258761cc2": 42,
            "cf3bdc30-370c-48d3-8fff-b587a39d72d6": 1
        }
    }
    for cell, expected in expected_high.items():
        actual = {p["dataset_key"]: p["source_dataset_count"] for p in pairs
                  if (p["species"], p["country"]) == cell}
        if actual != expected:
            raise ValueError("Two strongest source candidates changed unexpectedly")
    return sorted(pairs, key=lambda r: (r["species"], r["country"], r["dataset_key"])), {
        "source_qc_sha256": sha(source_bytes),
        "registry_sha256": sha(registry.read_bytes()),
        "seven_candidate_cells": len(source_keys),
        "dataset_country_cells": len(pairs),
        "dataset_registry_keys": len(registry_by_key),
        "source_qc_executed_utc": src["executed_utc"],
    }


def query_url(pair):
    params = [
        ("taxonKey", str(pair["taxon_key"])), ("country", pair["country"]),
        ("datasetKey", pair["dataset_key"]), ("year", "2018,2025"),
        ("hasCoordinate", "true"), ("hasGeospatialIssue", "false"),
        ("occurrenceStatus", "PRESENT"), ("limit", "0"),
        ("facet", "year"), ("facetLimit", "300")
    ]
    params += [("basisOfRecord", basis) for basis in BASIS]
    return API + "?" + urlencode(params)


def validate_year_facet(raw):
    try:
        payload = json.loads(raw)
    except (ValueError, UnicodeDecodeError) as exc:
        raise ValueError("Truncated or invalid occurrence facet JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError("Expected JSON object")
    count = payload.get("count")
    if type(count) is not int or count < 0 or payload.get("limit") != 0 or payload.get("results") != []:
        raise ValueError("Expected limit-zero facet and nonnegative integer total")
    facets = payload.get("facets")
    if not isinstance(facets, list) or len(facets) != 1:
        raise ValueError("Missing or duplicated year facet")
    facet = facets[0]
    if not isinstance(facet, dict) or str(facet.get("field")).upper().replace("_", "") != "YEAR":
        raise ValueError("Unexpected facet name")
    entries = facet.get("counts")
    if not isinstance(entries, list) or len(entries) >= 300:
        raise ValueError("Missing/truncated year entries")
    out = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Invalid year entry")
        year, n = entry.get("name"), entry.get("count")
        if not isinstance(year, str) or year not in YEAR_RANGE or type(n) is not int or n <= 0 or year in out:
            raise ValueError("Duplicated, missing or out-of-range year category")
        out[year] = n
    if sum(out.values()) != count:
        raise ValueError("Year counts do not sum to the GBIF query total")
    return {"record_count": count, "year_counts": dict(sorted(out.items()))}


def fetch_dataset_year(pair, *, retries=3):
    url = query_url(pair)
    failures = []
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "chocho-butterfly-dataset-year-provenance/1.0",
                                        "Accept": "application/json"})
            with urlopen(req, timeout=35) as response:
                http = response.status
                raw = response.read(1000000)
            if http != 200 or len(raw) >= 1000000:
                raise ValueError("GBIF HTTP response incomplete")
            value = validate_year_facet(raw)
            return {
                **pair, **value,
                "status": ("COMPLETE_MATCHED_SOURCE_COUNT" if value["record_count"] == pair["source_dataset_count"]
                           else "INDEX_DRIFT_NOT_USED_FOR_SOURCE_INFERENCE"),
                "retrieved_utc": utc_now(), "query_url": url,
                "raw_response_sha256": sha(raw), "http_status": http, "attempts": attempt + 1,
            }
        except (HTTPError, URLError, OSError, TimeoutError, ValueError) as err:
            failures.append({"attempt": attempt + 1, "exception": type(err).__name__,
                             "message": str(err)[:300]})
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
    return {**pair, "status": "UNAVAILABLE_NOT_ZERO", "query_url": url,
            "retrieved_utc": utc_now(), "attempts": retries, "failures": failures}


def assemble(pairs, fetched, meta):
    lookup = {(p["species"], p["country"], p["dataset_key"]): p for p in pairs}
    entries = {(p["species"], p["country"], p["dataset_key"]): p for p in fetched}
    if len(entries) != len(fetched) or set(entries) != set(lookup):
        raise ValueError("Missing or duplicate dataset × country × species source query")
    results = []
    for sp, ct in sorted({(p["species"], p["country"]) for p in pairs}):
        source = [lookup[k] for k in lookup if k[0] == sp and k[1] == ct]
        observed = [entries[k] for k in lookup if k[0] == sp and k[1] == ct]
        if any(r["status"] != "COMPLETE_MATCHED_SOURCE_COUNT" for r in observed):
            results.append({"species": sp, "country": ct,
                            "status": "SOURCE_DRIFT_OR_UNAVAILABLE_NO_INFERENCE",
                            "partial_sources": observed})
            continue
        year_sums = Counter()
        matrix = []
        for p in observed:
            year_sums.update(p["year_counts"])
            matrix.append(p)
        expected_years = source[0]["source_country_years"]
        expected_total = source[0]["source_country_total"]
        if dict(year_sums) != expected_years or sum(year_sums.values()) != expected_total:
            results.append({"species": sp, "country": ct, "status": "YEAR_MARGINAL_DRIFT_NO_INFERENCE",
                            "source_country_years": expected_years,
                            "current_country_years": dict(sorted(year_sums.items())),
                            "datasets": matrix})
            continue
        results.append({"species": sp, "country": ct, "status": "COMPLETE_RECONCILED",
                        "record_count": expected_total,
                        "country_years": dict(sorted(year_sums.items())),
                        "added_host_country": source[0]["added_host_country"],
                        "datasets": matrix})
    all_complete = all(row["status"] == "COMPLETE_RECONCILED" for row in results)
    report = {
        "schema": SCHEMA, "source_protocol": PROTOCOL, "completed_utc": utc_now(),
        "source": meta, "expected_dataset_queries": MATRIX_SIZE,
        "complete_dataset_queries": sum(p["status"] == "COMPLETE_MATCHED_SOURCE_COUNT" for p in fetched),
        "complete_country_cells": sum(r["status"] == "COMPLETE_RECONCILED" for r in results),
        "status": "ALL_SEVEN_SOURCE_CELLS_RECONCILED" if all_complete else "INCOMPLETE_NO_INFERENCE",
        "outcomes_cannot_establish_butterfly_colonization": True,
        "country_cells": results,
    }
    if all_complete:
        il = next(r for r in results if (r["species"], r["country"]) == ("Papilio demoleus", "IL"))
        do = next(r for r in results if (r["species"], r["country"]) == ("Aguna asander", "DO"))
        il_2025 = il["country_years"].get("2025", 0)
        by_dset = {p["dataset_title"]: p for p in il["datasets"]}
        inpa = by_dset["INPA - Observations"]
        other = [p for p in il["datasets"] if p["dataset_title"] != "INPA - Observations"]
        by_year_other = Counter()
        for p in other:
            by_year_other.update(p["year_counts"])
        report["primary_two_case_readout"] = {
            "israel_papilio": {
                "record_count": il["record_count"], "2025_total": il_2025,
                "2024_total": il["country_years"].get("2024", 0),
                "2023_total": il["country_years"].get("2023", 0),
                "2025_by_dataset_title": {p["dataset_title"]: p["year_counts"].get("2025", 0)
                                          for p in il["datasets"]},
                "inpa_2025_contribution": inpa["year_counts"].get("2025", 0),
                "inpa_dataset_created": inpa["dataset_created_utc"],
                "inpa_fraction_of_2025_records": (inpa["year_counts"].get("2025", 0) / il_2025 if il_2025 else None),
                "non_inpa_year_marginals": dict(sorted(by_year_other.items())),
                "non_inpa_total": sum(by_year_other.values()),
                "added_host_country": il["added_host_country"],
            },
            "dominican_aguna": {
                "record_count": do["record_count"],
                "country_years": do["country_years"],
                "dataset_count": len(do["datasets"]),
                "added_host_country": do["added_host_country"]
            }
        }
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source-qc", type=Path, required=True)
    p.add_argument("--dataset-registry", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    pairs, meta = load_fixed_pairs(args.source_qc, args.dataset_registry)
    print("DATASET_YEAR_LOCKED_PAIRS=" + json.dumps([
        {"species": p["species"], "country": p["country"], "dataset_key": p["dataset_key"],
         "expected_records": p["source_dataset_count"]} for p in pairs]), flush=True)
    fetched = []
    for pair in pairs:
        fetched.append(fetch_dataset_year(pair))
        time.sleep(1)
    result = assemble(pairs, fetched, meta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("DATASET_YEAR_AUDIT_RESULT=" + json.dumps(result, separators=(",", ":")), flush=True)
    if result["status"] != "ALL_SEVEN_SOURCE_CELLS_RECONCILED":
        raise SystemExit("INCOMPLETE_OR_SOURCE_DRIFT: preserve receipts; do not infer absence or growth")


if __name__ == "__main__":
    main()
