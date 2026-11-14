#!/usr/bin/env python3
"""Bounded original-GSHHG island search for larva/egg GBIF records in 15 cases.

This is an **indexed immature-stage evidence audit**, NOT a test of successful
colonization or establishment. Geometry is borrowed from the exact 225 x 51
source reconstruction. An API nonindex is UNKNOWN breeding status, NOT absence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from build_public51_butterfly_risk_20261109 import (
    GBIF, GSHHG_SHA256, get_json, load_polygons,
)

PERIODS = {"first_index_2011_2017": "2011,2017",
           "reindex_2018_2022": "2018,2022"}
STAGES = ("Larva", "Egg")


def write_csv(path, rows):
    fields = list(rows[0]) if rows else ["speciesKey", "island", "status"]
    with path.open("w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def search(params):
    data = get_json(GBIF + "/occurrence/search", params)
    # Count can be nonzero even if a species is not part of focal original 15.
    n = int(data.get("count", -1))
    if n < 0:
        raise ValueError("GBIF response omitted usable count")
    facets = {}
    for facet in data.get("facets", []):
        if str(facet.get("field", "")).upper() == "SPECIES_KEY":
            for entry in facet.get("counts", []):
                k, count = str(entry["name"]), int(entry["count"])
                if k in facets:
                    raise ValueError("Repeated GBIF speciesKey in facet")
                facets[k] = count
    if n > 0 and not facets:
        raise ValueError("Nonzero search without SPECIES_KEY facet")
    return n, facets, hashlib.sha256(
        json.dumps(data, sort_keys=True, ensure_ascii=False).encode("utf8")
    ).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gshhg-zip", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--cases", type=Path,
                   default=Path("results/port_host_reindex_20261114/audit.json"))
    p.add_argument("--out-dir", type=Path,
                   default=Path("results/port_host_immature_pilot_20261114"))
    p.add_argument("--max-wkt-chars", type=int, default=2300)
    p.add_argument("--max-facets", type=int, default=190)
    args = p.parse_args()

    b = args.cases.read_bytes()
    saved = json.loads(b)
    cases = saved["introduced_case_index"]
    if len(cases) != 15 or len({(c["species"], c["island"]) for c in cases}) != 15:
        raise ValueError("Frozen main introduced-host first-index cases no longer 15 unique cells")
    by_island = defaultdict(set)
    for r in cases:
        by_island[str(r["island"])].add(str(r["species"]))
    # Candidate case index uses GBIF speciesKey under the key "species".
    # The accepted binomial is "name"; never use the latter as a taxonKey.
    assert len(by_island) == 11
    shape, shape_receipts, shape_digest = load_polygons(
        args.gshhg_zip, args.work, set(by_island), args.max_wkt_chars)
    if shape_digest != GSHHG_SHA256:
        raise ValueError("Original GSHHG source hash not verified")
    expected = sum(len(parts) * 2 * len(STAGES) for parts in shape.values())
    if expected > args.max_facets:
        raise ValueError(
            f"Source-geometry target requires {expected} facet queries > "
            f"predeclared cap {args.max_facets}; no partial search or false zeros")
    print("PREDECLARED_FACET_QUERY_COUNT", expected, flush=True)

    cases_index = {(str(z["island"]), str(z["species"])): z for z in cases}
    cumulative = defaultdict(int)
    details, receipts, errors = [], [], []
    total_queries = 0
    for island in sorted(by_island, key=int):
        species = sorted(by_island[island], key=int)
        for index, polygon in enumerate(shape[island]):
            for era, interval in PERIODS.items():
                for stage in STAGES:
                    params = {
                        "geometry": polygon, "year": interval,
                        "taxonKey": species,
                        "hasCoordinate": "true", "hasGeospatialIssue": "false",
                        "occurrenceStatus": "PRESENT", "lifeStage": stage,
                        "facet": "speciesKey", "facetLimit": 1000, "limit": 0,
                    }
                    total_queries += 1
                    try:
                        n, found, digest = search(params)
                        if any(key not in species for key in found):
                            raise ValueError(
                                f"GBIF speciesKey facet escaped requested species: {sorted(found)}")
                        if sum(found.values()) > n:
                            raise ValueError("Facet species counts exceed all indexed records")
                        for key, count in found.items():
                            cumulative[(island, key, era, stage)] += count
                        receipts.append({
                            "island": island, "fragment_number": index+1,
                            "era": era, "stage": stage, "indexed_total": n,
                            "indexed_focal_species_records": sum(found.values()),
                            "response_sha256": digest,
                            "source_WKT_sha256": hashlib.sha256(polygon.encode()).hexdigest(),
                            "query_status": "OK",
                        })
                    except Exception as error:
                        errors.append({
                            "island": island, "fragment_number": index + 1,
                            "era": era, "stage": stage, "error": repr(error)
                        })
                        receipts.append({
                            "island": island, "fragment_number": index+1,
                            "era": era, "stage": stage, "indexed_total": "",
                            "indexed_focal_species_records": "",
                            "response_sha256": "",
                            "source_WKT_sha256": hashlib.sha256(polygon.encode()).hexdigest(),
                            "query_status": "FAILED_UNRESOLVED",
                        })
                    print("FACET", island, index+1, era, stage,
                          receipts[-1]["query_status"],
                          receipts[-1]["indexed_focal_species_records"], flush=True)

    for case in cases:
        island, key = str(case["island"]), str(case["species"])
        for era in PERIODS:
            for stage in STAGES:
                matching = [r for r in receipts
                            if r["island"] == island and r["era"] == era
                            and r["stage"] == stage]
                all_covered = all(r["query_status"] == "OK" for r in matching)
                details.append({
                    "species": case["name"],
                    "speciesKey": key, "island": island, "era": era, "stage": stage,
                    "status": "COMPLETE_FACET" if all_covered else "INCOMPLETE_UNKNOWN",
                    "indexed_stage_records": cumulative[(island, key, era, stage)]
                        if all_covered else "",
                    "GBIF_nonindex_NOT_true_breeding_absence": True,
                })
    assert len(details) == 15 * len(PERIODS) * len(STAGES)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "case_immature_stage_facets.csv", details)
    write_csv(args.out_dir / "api_fragment_receipts.csv", receipts)
    audit = {
        "analysis": "GBIF_IMMATURE_SOURCE_INDEX_NOT_ESTABLISHMENT",
        "utc": datetime.now(timezone.utc).isoformat(),
        "case_input_sha256": hashlib.sha256(b).hexdigest(),
        "gshhg_zip_sha256": shape_digest,
        "geometry_audit": shape_receipts,
        "cases": len(cases),
        "island_groups": len(by_island),
        "facet_queries_planned": expected, "facet_queries_attempted": total_queries,
        "facet_query_errors": errors,
        "case_stage_rows": len(details),
        "complete_case_stage_rows": sum(r["status"] == "COMPLETE_FACET" for r in details),
        "complete_positive_case_stage_rows": [
            r for r in details
            if r["status"] == "COMPLETE_FACET" and int(r["indexed_stage_records"]) > 0
        ],
        "zero_is_not_breeding_absence": True,
        "limitations": [
            "lifeStage=Larva/Egg is incompletely populated across GBIF publishers",
            "original case set selected after observational host indexing and butterfly records",
            "GBIF facet matches are stage-indexed locality observations, not proved multi-season breeding",
            "origin statuses are modern regional labels, not local historic dates",
            "speciesKey and isotope/photographic verification would need original record review",
            "this pilot is not a globally comprehensive breeding comparison",
        ],
    }
    (args.out_dir / "source_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    print(json.dumps({
        "status": "PARTIAL_SOURCE_ERRORS" if errors else "COMPLETE",
        "planned": expected, "attempted": total_queries,
        "errors": len(errors), "positive_stage_cells": len(audit["complete_positive_case_stage_rows"]),
        "output": str(args.out_dir),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
