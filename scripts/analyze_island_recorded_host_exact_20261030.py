#!/usr/bin/env python3
"""Exact, exploratory co-recorded known-host comparison on the frozen island fiber.

This is NOT a test of biological host absence or butterfly establishment.
The exposed indicator denotes a taxonomically conservative known-host species
with an indexed pre-2018 botanical record in the exact target island.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]
HOST = "data/exploratory/island_dated_host_20261030/dated_host_primary_risk_cells.csv"
PANEL = "data/exploratory/island_expansion_20261029/expanded_island_panel_selection.csv"
LATITUDE = "data/exploratory/island_latitude_null_20261030/latitude_null_diagnostic.json"
EXPECTED_FIBER = 147_117_600
EXPECTED_CELLS = 2385
EXPECTED_POSITIVES = 66


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def zone_for(latitude: float) -> str:
    if latitude < -23.5:
        return "south_extratropical"
    if latitude > 23.5:
        return "north_extratropical"
    return "tropical"


def convolve(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    out = Counter()
    for h, n in a.items():
        for k, m in b.items():
            out[h + k] += n * m
    return dict(sorted(out.items()))


def summary_for_distribution(dist: dict[int, int], observed: int) -> dict:
    n = sum(dist.values())
    assert n > 0
    tail = sum(freq for h, freq in dist.items() if h >= observed)
    return {
        "observed_recorded_host_cells": observed,
        "null_expected_recorded_host_cells": sum(h * freq for h, freq in dist.items()) / n,
        "exact_fiber_count": n,
        "inclusive_upper_tail_count": tail,
        "inclusive_upper_tail_probability": tail / n,
        "score_support_min": min(dist),
        "score_support_max": max(dist),
    }


def block_distribution(rows: list[tuple[int, tuple[tuple[int, int, float], ...]]],
                       caps: tuple[int, ...]) -> dict[int, int]:
    """Dynamic program over row choices and remaining fixed island capacities.

    A row entry is (number of positives, ((column index, host flag, distance), ...)).
    Its output values are exact integer *numbers of binary tables*.
    """
    @lru_cache(None)
    def dp(row_index: int, remaining: tuple[int, ...]):
        if row_index == len(rows):
            return ((0, 1),) if all(v == 0 for v in remaining) else ()
        positives, cells = rows[row_index]
        out = Counter()
        for chosen in combinations(cells, positives):
            remaining_next = list(remaining)
            host_increment = 0
            for col, host, _distance in chosen:
                remaining_next[col] -= 1
                host_increment += host
            if min(remaining_next, default=0) < 0:
                continue
            for score, frequency in dp(row_index + 1, tuple(remaining_next)):
                out[score + host_increment] += frequency
        return tuple(sorted(out.items()))
    return dict(dp(0, caps))


def block_moments(rows: list[tuple[int, tuple[tuple[int, int, float], ...]]],
                  caps: tuple[int, ...]) -> dict:
    """Independent direct enumeration of complete tables and H/D moments."""
    count = 0
    sum_h = sum_h2 = sum_d = sum_d2 = sum_hd = 0.0
    def walk(index: int, remaining: tuple[int, ...], h: int, d: float):
        nonlocal count, sum_h, sum_h2, sum_d, sum_d2, sum_hd
        if index == len(rows):
            if any(remaining):
                return
            count += 1
            sum_h += h
            sum_h2 += h * h
            sum_d += d
            sum_d2 += d * d
            sum_hd += h * d
            return
        margin, cells = rows[index]
        for chosen in combinations(cells, margin):
            next_caps = list(remaining)
            for col, _host, _distance in chosen:
                next_caps[col] -= 1
            if min(next_caps, default=0) < 0:
                continue
            walk(index + 1, tuple(next_caps),
                 h + sum(c[1] for c in chosen),
                 d + sum(c[2] for c in chosen))
    walk(0, caps, 0, 0.0)
    assert count > 0
    mh, md = sum_h / count, sum_d / count
    return {
        "count": count,
        "expected_H": mh,
        "variance_H": max(0.0, sum_h2 / count - mh * mh),
        "expected_D": md,
        "variance_D": max(0.0, sum_d2 / count - md * md),
        "covariance_H_D": sum_hd / count - mh * md,
    }


def run(root: Path = ROOT) -> dict:
    inputs = {name: root / rel for name, rel in
              (("host", HOST), ("panel", PANEL), ("latitude", LATITUDE))}
    hashes = {key: hashlib.sha256(path.read_bytes()).hexdigest()
              for key, path in inputs.items()}
    host_cells = read_csv(inputs["host"])
    panel = {row["island_key"]: row for row in read_csv(inputs["panel"])}
    diag = json.loads(inputs["latitude"].read_text(encoding="utf-8"))["analyses"]["primary"]
    expected_blocks = {b["longitude_block"]: b for b in diag["blocks"]}
    if len(expected_blocks) != 12:
        raise ValueError("Expected the twelve previously fixed latitude blocks")

    if len(host_cells) != EXPECTED_CELLS:
        raise ValueError("Frozen eligible risk set changed")
    seen_keys = set()
    grouped = {k: [] for k in expected_blocks}
    observed_h = 0
    observed_d = 0.0
    total_positive = 0
    for c in host_cells:
        key = (c["current_butterfly_key"], c["target_island"])
        if key in seen_keys or key[1] not in panel:
            raise ValueError(f"Bad or missing primary key: {key}")
        seen_keys.add(key)
        if c["late_positive"] not in ("0", "1") or c["confirmed_host_record"] not in ("True", "False"):
            raise ValueError(f"Unexpected outcome or host state: {key}")
        expected_flag = c["host_record_status"] == "RECORDED_KNOWN_HOST"
        if expected_flag != (c["confirmed_host_record"] == "True"):
            raise ValueError(f"Source host flags disagree: {key}")
        positive = int(c["late_positive"])
        host = int(c["confirmed_host_record"] == "True")
        dist = math.log1p(float(c["nearest_source_upper_km"]))
        block = f'{c["longitude_block"]}:{zone_for(float(panel[key[1]]["representative_lat"]))}'
        if block not in grouped:
            raise ValueError(f"Unrecognized frozen latitude block {block}")
        grouped[block].append((key[0], key[1], positive, host, dist))
        total_positive += positive
        observed_h += positive * host
        observed_d += positive * dist

    if total_positive != EXPECTED_POSITIVES or observed_h != 20:
        raise ValueError("Previously recorded observed outcomes changed")
    if not math.isclose(observed_d / total_positive, float(diag["observed_T"]), abs_tol=1e-9):
        raise ValueError("Distance-score audit does not reproduce prior observed statistic")

    blocks = []
    full = {0: 1}
    global_expected_h = global_expected_d = global_variance_h = global_variance_d = global_covariance = 0.0
    for name, expected in expected_blocks.items():
        cells = grouped[name]
        per_row, per_col = defaultdict(list), Counter()
        observed_block_h = 0
        block_positives = 0
        for species, island, positive, host, dist in cells:
            per_row[species].append((island, positive, host, dist))
            per_col[island] += positive
            observed_block_h += positive * host
            block_positives += positive
        col_names = sorted(k for k, n in per_col.items() if n > 0)
        col_position = {name: i for i, name in enumerate(col_names)}
        rows = []
        active_cell_count = 0
        for species in sorted(per_row):
            values = per_row[species]
            margin = sum(x[1] for x in values)
            if margin == 0:
                continue
            allowed = tuple((col_position[x[0]], x[2], x[3])
                            for x in values if x[0] in col_position)
            active_cell_count += len(allowed)
            if margin > len(allowed):
                raise ValueError("A fixed row exceeds allowed eligible columns")
            rows.append((margin, allowed))
        capacities = tuple(per_col[k] for k in col_names)
        if (len(cells) != int(expected["allowed_cells"])
                or block_positives != int(expected["positive_cells"])
                or len(rows) != int(expected["active_rows"])
                or len(col_names) != int(expected["active_columns"])
                or active_cell_count != int(expected["active_allowed_cells"])):
            raise ValueError(f"Frozen structural margins changed: {name}")

        hist = block_distribution(rows, capacities)
        moments = block_moments(rows, capacities)
        n = sum(hist.values())
        if n != moments["count"] or str(n) != expected["fiber_cardinality_decimal"]:
            raise ValueError(f"Mismatch with old exactly counted block {name}")
        full = convolve(full, hist)
        global_expected_h += moments["expected_H"]
        global_expected_d += moments["expected_D"]
        global_variance_h += moments["variance_H"]
        global_variance_d += moments["variance_D"]
        global_covariance += moments["covariance_H_D"]
        blocks.append({
            "block": name, "allowed_cells": len(cells),
            "observed_recorded_host_cells": observed_block_h,
            "distribution": {str(h): n for h, n in hist.items()},
            **moments,
        })
    if sum(full.values()) != EXPECTED_FIBER or str(sum(full.values())) != diag["fiber_cardinality_decimal"]:
        raise ValueError("Full table count does not match the already verified fiber")

    main = summary_for_distribution(full, observed_h)
    if abs(main["null_expected_recorded_host_cells"] - global_expected_h) > 1e-10:
        raise ValueError("Independent direct enumeration and score DP disagree")

    nonconstant = [b for b in blocks if len(b["distribution"]) > 1]
    sensitivity = {}
    for removed in nonconstant:
        subset = {0: 1}
        obs = 0
        for b in blocks:
            if b["block"] == removed["block"]:
                continue
            subset = convolve(subset, {int(h): n for h, n in b["distribution"].items()})
            obs += b["observed_recorded_host_cells"]
        sensitivity[removed["block"]] = summary_for_distribution(subset, obs)

    score_correlation = (global_covariance / math.sqrt(global_variance_h * global_variance_d)
                         if global_variance_h > 0 and global_variance_d > 0 else None)
    return {
        "schema": "chocho_exploratory_recorded_host_exact_fiber_20261030_v1",
        "not_biological_host_absence": True,
        "not_independent_of_distance_or_recording": True,
        "not_confirmatory": True,
        "input_sha256": hashes,
        "original_candidate_cells": len(host_cells),
        "original_later_positive_cells": total_positive,
        "host_recording_statistic": main,
        "score_probability_denominator": EXPECTED_FIBER,
        "score_probability_numerator": main["inclusive_upper_tail_count"],
        "score_distribution": {str(h): n for h, n in full.items()},
        "blocks": blocks,
        "informative_host_blocks": len(nonconstant),
        "leave_one_informative_block_out": sensitivity,
        "distance_score_companion": {
            "observed_total_log1p_gap": observed_d,
            "null_expected_total_log1p_gap": global_expected_d,
            "null_H_D_covariance": global_covariance,
            "null_H_D_correlation": score_correlation,
            "direction": "negative correlation means closer-source tables tend to have more host-recorded cells",
            "not_distance_adjusted_host_test": True,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, help="Write JSON result, otherwise print to stdout")
    args = parser.parse_args()
    result = run(args.root)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else args.root / args.output
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(rendered, encoding="utf-8")
        print(f"Wrote {target}")
        print(json.dumps(result["host_recording_statistic"], sort_keys=True))
    else:
        print(rendered)


if __name__ == "__main__":
    main()
