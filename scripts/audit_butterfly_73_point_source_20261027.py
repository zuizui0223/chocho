#!/usr/bin/env python3
"""Source-complete live GBIF points for all 73 pre-frozen PR72 positive cases.

Only metadata/aggregated source receipts leave memory. The country assignment
describes host-resource opportunity; NOT observed butterfly point climate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from analyze_butterfly_host_point_proximity_20261026 import (
    BASIS, valid_point,
)

SOURCE_SCHEMA = "chocho_pr72_exact73_point_followup_roster_20261027_v01"
OUT_SCHEMA = "chocho_source_verified_73_positive_butterfly_points_20261027_v01"
MIN_INTERVAL = 1.05
MAX_QUERY_ROWS = 5000
MAX_PANEL_ROWS = 40000
PAGE = 300


def source_roster(path):
    source = json.loads(Path(path).read_text(encoding="utf-8"))
    if source.get("schema") != SOURCE_SCHEMA or source.get("status") != \
            "SOURCE_IDENTITY_FROZEN_BEFORE_NEW_POINT_OBSERVATION":
        raise ValueError("The original PR72 cohort was not frozen")
    s = source["summary"]
    if (s["positive_country_species_cases"], s["novel_resource_country_cases"],
            s["analog_resource_country_cases"], s["distinct_butterfly_species"],
            s["original_GBIF_indexed_record_total_in_positive_cases"]) != (
            73, 32, 41, 20, 6076):
        raise ValueError("Source cohort counters differ from PR72")
    roster = source["roster"]
    if len(roster) != 73 or len({(x["species"], x["country"]) for x in roster}) != 73:
        raise ValueError("Missing or duplicated PR72 original positive countries")
    if not all(x["resource_class"] in {"novel", "analogue"} and
               isinstance(x["gbif_accepted_species_key"], int) and
               x["original_gbif_record_count"] >= 1
               for x in roster):
        raise ValueError("Unverified original accepted taxonKey / classification")
    if sum(x["original_gbif_record_count"] for x in roster) != 6076:
        raise ValueError("Original source occurrence index count changed")
    return source


def query_url(case, offset=0):
    if offset < 0 or offset % PAGE != 0 or offset + PAGE > 100000:
        raise ValueError("Invalid GBIF bounded page offset")
    if case["resource_class"] not in ("novel", "analogue"):
        raise ValueError("Unknown climate-host opportunity class")
    p = [
        ("taxonKey", str(case["gbif_accepted_species_key"])),
        ("country", case["country"]),
        ("year", "2010,2025"),
        ("hasCoordinate", "true"),
        ("hasGeospatialIssue", "false"),
        ("occurrenceStatus", "PRESENT"),
        ("limit", str(PAGE)), ("offset", str(offset)),
    ] + [("basisOfRecord", b) for b in sorted(BASIS)]
    return "https://api.gbif.org/v1/occurrence/search?" + urlencode(p)


class Reader:
    def __init__(self):
        self.at = None
        self.calls = 0
        self.retries = 0

    def get(self, url):
        for attempt in range(4):
            delay = 0 if self.at is None else MIN_INTERVAL - (
                time.monotonic() - self.at)
            if delay > 0:
                time.sleep(delay)
            self.at = time.monotonic()
            self.calls += 1
            try:
                req = Request(url, headers={
                    "Accept": "application/json",
                    "User-Agent": "chocho-source-climate-niche-audit/20261027 (+https://github.com/zuizui0223/chocho)",
                })
                with urlopen(req, timeout=40) as response:
                    if response.status != 200:
                        raise ValueError("GBIF returned non-200")
                    raw = response.read()
                data = json.loads(raw)
                if not isinstance(data, dict) or "results" not in data or \
                        "count" not in data:
                    raise ValueError("Unrecognized GBIF occurrence JSON")
                return data, hashlib.sha256(raw).hexdigest()
            except (HTTPError, URLError, TimeoutError) as e:
                if isinstance(e, HTTPError) and e.code not in (
                    429, 500, 502, 503, 504,
                ):
                    raise
                self.retries += 1
                if attempt == 3:
                    raise RuntimeError("Repeated transient source failure") from e
                time.sleep(min(2 ** (attempt + 1), 12))
        raise RuntimeError("Source request retry exhaustion")


def valid_metadata(record, key, country):
    """Validate accepted source identity independently from GBIF point validity."""
    point, status = valid_point(record, key, country)
    if point is None:
        return None, status
    stage = str(record.get("lifeStage") or "").lower().strip()
    basis = record.get("basisOfRecord")
    if basis not in BASIS:
        return None, "unaccepted_basis"
    if stage in ("adult", "imago"):
        stage_class = "adult"
    elif any(s in stage for s in ("egg", "larva", "caterpillar", "pupa", "nymph")):
        stage_class = "juvenile_source_label"
    elif not stage:
        stage_class = "missing_stage"
    else:
        stage_class = "other_stage_label"
    source_wild = str(record.get("degreeOfEstablishment") or
                      record.get("establishmentMeans") or "").strip()
    coord = (round(point["lat"], 2), round(point["lon"], 2))
    return {
        "year": point["year"], "stage": stage_class,
        "datasetKey": point["dataset"], "strict": point["strict"],
        "coordinate_cell": coord,
        "reported_establishment_annotation": bool(source_wild),
    }, "valid"


def summarize_one(case, obs, pages, count):
    exclusion = Counter()
    rows = []
    ids = set()
    for rec in obs:
        index = rec.get("key")
        if index is None or index in ids:
            raise ValueError("Missing or duplicate source GBIF index key")
        ids.add(index)
        point, status = valid_metadata(
            rec, case["gbif_accepted_species_key"], case["country"])
        exclusion[status] += 1
        if point is not None:
            rows.append(point)
    if len(obs) != count or len(ids) != count:
        raise ValueError("GBIF index responses not exhaustive or duplicated")
    year = Counter("2010_2017" if x["year"] <= 2017 else "2018_2025"
                   for x in rows)
    stages = Counter(x["stage"] for x in rows)
    basis_tally = Counter(str(rec.get("basisOfRecord")) for rec in obs)
    dataset = Counter(str(x["datasetKey"]) for x in rows)
    cells = Counter(x["coordinate_cell"] for x in rows)
    return {
        "butterfly": case["species"], "country": case["country"],
        "resource_country_class": case["resource_class"],
        "original_GBIF_indexed_count": case["original_gbif_record_count"],
        "new_live_indexed_count": count,
        "original_to_live_count_delta": count - case["original_gbif_record_count"],
        "full_source_status": "COMPLETE_NEW_LIVE_GBIF_POINTS",
        "accepted_source_gbif_species_key": case["gbif_accepted_species_key"],
        "source_pages": pages,
        "source_record_validation": dict(sorted(exclusion.items())),
        "validated_point_indices": len(rows),
        "reported_butterfly_coord_uncertainty_le1000m": sum(
            x["strict"] for x in rows),
        "point_source_reporting_year_counts": dict(sorted(year.items())),
        "point_source_stage_labels": dict(sorted(stages.items())),
        "point_source_basis_of_record": dict(sorted(basis_tally.items())),
        "point_source_establishment_metadata_indices": sum(
            x["reported_establishment_annotation"] for x in rows),
        "point_0p01degree_reporting_cells": len(cells),
        "max_indices_in_one_0p01degree_cell": max(cells.values(), default=0),
        "distinct_index_dataset_keys": len(dataset),
        "largest_single_dataset_index_records": max(dataset.values(), default=0),
        "climate_annotation_status": "NOT_MEASURED_PENDING_COHERENT_POINT_AND_BASELINE_RASTER",
    }


def retrieve_case(case, reader, remaining):
    receipt = []
    try:
        first, first_sha = reader.get(query_url(case, 0))
        n = int(first["count"])
        if n < 0 or n > 100000:
            raise ValueError("Invalid source index count")
        if n > MAX_QUERY_ROWS or n > remaining:
            return {
                "butterfly": case["species"], "country": case["country"],
                "resource_country_class": case["resource_class"],
                "full_source_status": "NOT_EVALUABLE_PREDECLARED_RECORD_CAP",
                "new_live_indexed_count": n,
                "original_GBIF_indexed_count": case["original_gbif_record_count"]
            }
        pages = list(range(0, n, PAGE)) or [0]
        seen = []
        for offset in pages:
            data, digest = (first, first_sha) if offset == 0 else \
                reader.get(query_url(case, offset))
            if int(data["count"]) != n:
                raise ValueError("Source indexed count changed inside pagination")
            rec = data["results"]
            expected = min(PAGE, max(0, n - offset))
            if not isinstance(rec, list) or len(rec) != expected:
                raise ValueError("Source page truncated or mismatched")
            receipt.append({"offset": offset, "records": len(rec),
                            "source_body_sha256": digest})
            seen.extend(rec)
        return summarize_one(case, seen, receipt, n)
    except (HTTPError, URLError, ValueError, RuntimeError,
            TimeoutError) as exc:
        return {
            "butterfly": case["species"], "country": case["country"],
            "resource_country_class": case["resource_class"],
            "full_source_status": "INCOMPLETE_GBIF_SOURCE_QUERY",
            "original_GBIF_indexed_count": case["original_gbif_record_count"],
            "source_partial_page_hashes": receipt,
            "reason": str(exc)[:180],
        }


def make_report(source, rows, reader):
    complete = [r for r in rows
                if r["full_source_status"] == "COMPLETE_NEW_LIVE_GBIF_POINTS"]
    all_complete = len(complete) == 73
    class_report = {}
    for cls in ("novel", "analogue"):
        part = [r for r in complete if r["resource_country_class"] == cls]
        n = sum(r["validated_point_indices"] for r in part)
        class_report[cls] = {
            "source_positive_cases_total": 32 if cls == "novel" else 41,
            "complete_cases": len(part),
            "sum_new_live_indexed_count": sum(r["new_live_indexed_count"]
                                              for r in part),
            "validated_point_indices": n,
            "point_coord_uncertainty_le1000m": sum(
                r["reported_butterfly_coord_uncertainty_le1000m"] for r in part),
            "year_counts": dict(sum((Counter(r["point_source_reporting_year_counts"])
                                     for r in part), Counter())),
            "life_stage_counts": dict(sum((Counter(r["point_source_stage_labels"])
                                         for r in part), Counter())),
            "median_records_per_complete_case": sorted(
                (r["new_live_indexed_count"] for r in part))[
                    len(part)//2] if part else None,
        }
    complete_by_cases = sorted(
        complete, key=lambda r: r["validated_point_indices"], reverse=True)
    all_points = sum(r["validated_point_indices"] for r in complete)
    top2 = sum(r["validated_point_indices"] for r in complete_by_cases[:2])
    return {
        "schema": OUT_SCHEMA,
        "status": "COMPLETE_ALL_73_NEW_LIVE_GBIF_QUERY_SNAPSHOT" if all_complete
                  else "INCOMPLETE_73_CASE_QUERY_SNAPSHOT",
        "source_roster_schema": source["schema"],
        "original_73_cohort_source_positive_class_counts": {
            "novel": 32, "analogue": 41,
        },
        "source_original_positive_index_total": 6076,
        "network_requests": reader.calls,
        "network_retry_count": reader.retries,
        "cases_complete": len(complete),
        "cases_incomplete_or_above_cap": 73 - len(complete),
        "class_summary": class_report,
        "source_all_validated_point_indices": all_points,
        "point_record_share_from_two_largest_country_species_cases": (
            top2 / all_points if all_points else None),
        "top_two_case_labels": [
            {"butterfly": r["butterfly"], "country": r["country"],
             "validated_indices": r["validated_point_indices"]}
            for r in complete_by_cases[:2]
        ],
        "cases": rows,
        "source_bounded_query_contract": {
            "records_per_case_max": MAX_QUERY_ROWS,
            "panel_record_max": MAX_PANEL_ROWS,
            "GBIF_period": "2010–2025",
            "query_limit": PAGE, "rate_limit_seconds": MIN_INTERVAL,
            "same_five_original_basis_groups": sorted(BASIS),
        },
        "point_climate_not_measured": True,
        "biological_inference": (
            "These are GBIF index recording points only. Country resource class is "
            "NOT the point's measured climate, habitat, larval host use, "
            "population establishment or a dated niche shift."),
    }


def run(args):
    source = source_roster(args.roster)
    reader = Reader()
    rows = []
    remaining = MAX_PANEL_ROWS
    for case in sorted(source["roster"],
                       key=lambda r: (r["species"], r["country"])):
        result = retrieve_case(case, reader, remaining)
        rows.append(result)
        if result["full_source_status"] == "COMPLETE_NEW_LIVE_GBIF_POINTS":
            remaining -= result["new_live_indexed_count"]
    report = make_report(source, rows, reader)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = {k: report[k] for k in (
        "status", "cases_complete", "cases_incomplete_or_above_cap",
        "source_all_validated_point_indices", "network_requests",
        "network_retry_count", "point_record_share_from_two_largest_country_species_cases"
    )}
    summary["class_summary"] = report["class_summary"]
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--roster", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    run(p.parse_args())
