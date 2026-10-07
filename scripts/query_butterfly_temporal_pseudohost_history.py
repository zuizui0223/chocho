#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

from audit_butterfly_temporal_hindcast_full_history import (
    first_exact_record_year,
    load_geometry,
    resolve_species,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matches-csv", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--start-year", type=int, default=1800)
    ap.add_argument("--end-year", type=int, default=2025)
    ap.add_argument("--request-spacing-seconds", type=float, default=0.25)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    geom = load_geometry(args.level3_geojson)
    rows = list(csv.DictReader(args.matches_csv.open(newline="", encoding="utf-8")))
    if not rows:
        raise RuntimeError("empty pseudo-host match table")

    host_names = sorted({str(r["candidate_name"]).strip() for r in rows})
    butterflies = sorted({str(r["species"]).strip() for r in rows})
    resolved = {}
    resolution_errors = {}

    for name in butterflies + host_names:
        if name in resolved or name in resolution_errors:
            continue
        try:
            resolved[name] = resolve_species(name)
        except Exception as exc:
            resolution_errors[name] = str(exc)
        time.sleep(args.request_spacing_seconds)

    first_cache = {}
    query_errors = {}

    def first_for(name: str, region: str):
        key = (name, region)
        if key in first_cache:
            return first_cache[key]
        if name in resolution_errors:
            result = {"first_year": None, "first_record_count": 0, "first_basis": "", "first_dataset_count": 0}
            first_cache[key] = result
            return result
        if region not in geom:
            query_errors[key] = "missing WGSRPD3 geometry"
            result = {"first_year": None, "first_record_count": 0, "first_basis": "", "first_dataset_count": 0}
            first_cache[key] = result
            return result
        try:
            taxon_key, _ = resolved[name]
            result = first_exact_record_year(
                taxon_key,
                geom[region]["geometry"],
                args.start_year,
                args.end_year,
            )
        except Exception as exc:
            query_errors[key] = str(exc)
            result = {"first_year": None, "first_record_count": 0, "first_basis": "", "first_dataset_count": 0}
        first_cache[key] = result
        time.sleep(args.request_spacing_seconds)
        return result

    butterfly_first = {}
    for sp in butterflies:
        for region in sorted({r["wgsrpd3_code"] for r in rows if r["species"] == sp}):
            butterfly_first[(sp, region)] = first_for(sp, region)

    out = []
    for n, row in enumerate(rows, start=1):
        sp = str(row["species"]).strip()
        region = str(row["wgsrpd3_code"]).strip()
        host = str(row["candidate_name"]).strip()
        hf = first_for(host, region)
        bf = butterfly_first[(sp, region)]
        hy = hf["first_year"]
        by = bf["first_year"]
        out.append({
            **row,
            "local_host_first_record_year": "" if hy is None else hy,
            "local_host_first_basis": hf["first_basis"],
            "local_host_first_dataset_count": hf["first_dataset_count"],
            "butterfly_first_record_year_full": "" if by is None else by,
            "butterfly_first_basis": bf["first_basis"],
            "butterfly_first_dataset_count": bf["first_dataset_count"],
            "local_record_lag_years": "" if hy is None or by is None else int(by) - int(hy),
            "host_query_error": query_errors.get((host, region), ""),
            "butterfly_query_error": query_errors.get((sp, region), ""),
        })
        if n % 25 == 0:
            print(json.dumps({"completed_assignments": n, "total": len(rows)}), flush=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    actual = [r for r in out if r["assignment_type"] == "actual"]
    pseudo = [r for r in out if r["assignment_type"] == "pseudo"]
    payload = {
        "schema": "chocho_butterfly_temporal_pseudohost_local_history_v0.1",
        "status": "POSTHOC_PSEUDOHOST_LOCAL_RECORD_HISTORY",
        "record_window": [args.start_year, args.end_year],
        "assignment_rows": len(out),
        "actual_rows": len(actual),
        "pseudo_rows": len(pseudo),
        "unique_host_region_queries": sum(1 for (name, region) in first_cache if name in host_names),
        "actual_rows_with_local_record": sum(bool(str(r["local_host_first_record_year"]).strip()) for r in actual),
        "pseudo_rows_with_local_record": sum(bool(str(r["local_host_first_record_year"]).strip()) for r in pseudo),
        "taxon_resolution_errors": resolution_errors,
        "local_query_error_count": len(query_errors),
        "local_query_errors": [
            {"name": k[0], "region": k[1], "error": v}
            for k, v in sorted(query_errors.items())
        ],
        "claim_boundary": [
            "First GBIF record is a detection/digitization date, not an establishment date.",
            "The same 1800-2025 window and exact WGSRPD3 geometry are used for actual and pseudo plants.",
            "This is a conditional post-hoc identity-specific timing null; event cells were selected after apparent butterfly outcomes."
        ],
    }
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
