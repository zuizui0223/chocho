#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from audit_butterfly_temporal_hindcast_full_history import (
    build_snapshot_candidates,
    first_exact_record_year,
    load_geometry,
    load_pairs,
    load_units,
    null_summary,
    resolve_species,
)

EXPECTED = {
    "candidate_cells": 1477,
    "treatment_cells": 442,
    "control_cells": 1035,
    "treatment_snapshot_hits": 9,
    "control_snapshot_hits": 11,
}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def prepare(args) -> None:
    geom = load_geometry(args.level3_geojson)
    pairs = load_pairs(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    cont = load_units(args.contemporary_distribution_csv)
    rows = build_snapshot_candidates(
        args.occurrence_mapped_csv, pairs, native, cont, geom, args.effort_threshold
    )
    snapshot = {
        "candidate_cells": len(rows),
        "treatment_cells": sum(int(r["treatment"]) for r in rows),
        "control_cells": sum(1 - int(r["treatment"]) for r in rows),
        "treatment_snapshot_hits": sum(
            int(r["treatment"]) and int(r["snapshot_outcome"]) for r in rows
        ),
        "control_snapshot_hits": sum(
            (1 - int(r["treatment"])) and int(r["snapshot_outcome"]) for r in rows
        ),
    }
    if snapshot != EXPECTED:
        raise RuntimeError(f"frozen candidate reconstruction drift: {snapshot} != {EXPECTED}")
    for i, row in enumerate(rows):
        row["candidate_index"] = i
    write_csv(args.output_csv, rows)
    print(json.dumps(snapshot, sort_keys=True))


def query(args) -> None:
    geom = load_geometry(args.level3_geojson)
    with args.candidates_csv.open(newline="", encoding="utf-8") as f:
        all_rows = list(csv.DictReader(f))
    selected = [
        row for row in all_rows
        if int(row["candidate_index"]) % args.shard_count == args.shard_index
    ]
    resolved = {s: resolve_species(s) for s in sorted({r["species"] for r in selected})}

    def task(row):
        key, _ = resolved[row["species"]]
        result = first_exact_record_year(
            key,
            geom[row["wgsrpd3_code"]]["geometry"],
            args.history_start,
            args.history_end,
        )
        out = dict(row)
        out["historical_first_record_year"] = (
            "" if result["first_year"] is None else result["first_year"]
        )
        out["historical_first_basis"] = result["first_basis"]
        out["historical_first_dataset_count"] = result["first_dataset_count"]
        return out

    out = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(task, row) for row in selected]
        for n, fut in enumerate(as_completed(futures), start=1):
            out.append(fut.result())
            if n % 25 == 0:
                print(json.dumps({
                    "shard": args.shard_index,
                    "completed": n,
                    "total": len(selected),
                }), flush=True)
    out.sort(key=lambda r: int(r["candidate_index"]))
    write_csv(args.output_csv, out)
    print(json.dumps({
        "shard": args.shard_index,
        "shard_count": args.shard_count,
        "rows": len(out),
    }, sort_keys=True))


def aggregate(args) -> None:
    rows = []
    for path in sorted(args.shard_dir.glob("history_shard_*.csv")):
        with path.open(newline="", encoding="utf-8") as f:
            rows.extend(csv.DictReader(f))
    rows.sort(key=lambda r: int(r["candidate_index"]))
    if len(rows) != EXPECTED["candidate_cells"]:
        raise RuntimeError(f"expected 1477 shard rows, got {len(rows)}")
    keys = {(r["species"], r["wgsrpd3_code"]) for r in rows}
    if len(keys) != len(rows):
        raise RuntimeError("duplicate candidate rows after shard merge")

    for r in rows:
        raw = str(r["historical_first_record_year"]).strip()
        first = None if not raw else int(float(raw))
        r["historical_record_pre2018"] = int(first is not None and first <= 2017)
        r["full_history_eligible"] = int(first is None or first >= 2018)
        r["full_history_outcome"] = int(first is not None and 2018 <= first <= 2025)

    clean = [r for r in rows if int(r["full_history_eligible"])]
    treat = [r for r in clean if int(r["treatment"])]
    control = [r for r in clean if not int(r["treatment"])]
    th = sum(int(r["full_history_outcome"]) for r in treat)
    ch = sum(int(r["full_history_outcome"]) for r in control)
    tr = th / len(treat) if treat else None
    cr = ch / len(control) if control else None

    null_species = null_summary(clean, ["species"], args.permutations, args.seed)
    null_level1 = null_summary(
        clean, ["species", "level1"], args.permutations, args.seed
    )
    null_full = null_summary(
        clean,
        [
            "species", "level1", "baseline_effort_bin",
            "test_effort_bin", "distance_bin"
        ],
        args.permutations,
        args.seed,
    )

    payload = {
        "schema": "chocho_butterfly_temporal_resource_hindcast_full_history_parallel_v0.1",
        "status": "POSTHOC_FULL_HISTORY_CORRECTED_TEMPORAL_HINDCAST",
        "snapshot_candidate_reconstruction": EXPECTED,
        "history_window": [1800, 2025],
        "history_audit": {
            "cells_with_pre2018_record_removed": sum(
                int(r["historical_record_pre2018"]) for r in rows
            ),
            "treatment_pre2018_removed": sum(
                int(r["historical_record_pre2018"]) and int(r["treatment"])
                for r in rows
            ),
            "control_pre2018_removed": sum(
                int(r["historical_record_pre2018"]) and not int(r["treatment"])
                for r in rows
            ),
            "clean_cells": len(clean),
            "clean_treatment_cells": len(treat),
            "clean_control_cells": len(control),
        },
        "observed": {
            "treatment_first_records_2018_2025": th,
            "treatment_rate": tr,
            "control_first_records_2018_2025": ch,
            "control_rate": cr,
            "raw_risk_difference": None if tr is None or cr is None else tr - cr,
            "raw_risk_ratio": None if tr is None or not cr else tr / cr,
        },
        "null_ladder": {
            "species_only": null_species,
            "species_plus_level1": null_level1,
            "full_protocol_match": null_full,
        },
        "decision": {
            "placement_signal_beyond_level1": bool(
                null_level1.get("p_high", 1) <= 0.05
            ),
            "placement_signal_beyond_full_match": bool(
                null_full.get("p_high", 1) <= 0.05
            ),
        },
        "claim_boundary": [
            "The capped 53,434-record snapshot defines the outcome-blind candidate and effort strata only; full GBIF history is queried for every candidate before temporal eligibility and outcome are defined.",
            "A first GBIF record is a first documented detection, not a colonization or establishment date.",
            "Contemporary WCVP introduced ranges remain time-invariant and do not establish host presence by 2017.",
        ],
    }
    write_csv(args.output_csv, rows)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prepare")
    p.add_argument("--occurrence-mapped-csv", type=Path, required=True)
    p.add_argument("--insect-host-csv", type=Path, required=True)
    p.add_argument("--native-distribution-csv", type=Path, required=True)
    p.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    p.add_argument("--level3-geojson", type=Path, required=True)
    p.add_argument("--effort-threshold", type=int, default=10)
    p.add_argument("--output-csv", type=Path, required=True)
    p.set_defaults(func=prepare)

    q = sub.add_parser("query")
    q.add_argument("--candidates-csv", type=Path, required=True)
    q.add_argument("--level3-geojson", type=Path, required=True)
    q.add_argument("--shard-index", type=int, required=True)
    q.add_argument("--shard-count", type=int, required=True)
    q.add_argument("--history-start", type=int, default=1800)
    q.add_argument("--history-end", type=int, default=2025)
    q.add_argument("--workers", type=int, default=4)
    q.add_argument("--output-csv", type=Path, required=True)
    q.set_defaults(func=query)

    a = sub.add_parser("aggregate")
    a.add_argument("--shard-dir", type=Path, required=True)
    a.add_argument("--permutations", type=int, default=99999)
    a.add_argument("--seed", type=int, default=20261007)
    a.add_argument("--output-csv", type=Path, required=True)
    a.add_argument("--output-json", type=Path, required=True)
    a.set_defaults(func=aggregate)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
