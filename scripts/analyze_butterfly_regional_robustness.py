#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np


REGION_NAMES = {
    "1": "Europe",
    "2": "Africa",
    "3": "Temperate Asia",
    "4": "Tropical Asia",
    "5": "Australasia",
    "6": "Pacific",
    "7": "Northern America",
    "8": "Southern America",
    "9": "Antarctic",
}


def average_ranks(values):
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def spearman(a, b):
    if len(a) < 3:
        return None
    x, y = average_ranks(a), average_ranks(b)
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def load_level1(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = {}
    for feature in payload.get("features", []):
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        region = str(props.get("LEVEL1_COD") or "").strip()
        if code and region:
            result[code] = region
    return result


def summarize(rows):
    return {
        "species": len(rows),
        "rho_host_family_vs_log_expansion": spearman(
            [row["host_family_count"] for row in rows],
            [row["log_resource_expansion"] for row in rows],
        ),
        "median_log_expansion": (
            None
            if not rows
            else float(np.median([row["log_resource_expansion"] for row in rows]))
        ),
        "expanded_fraction": (
            None
            if not rows
            else sum(row["introduced_added_units"] > 0 for row in rows) / len(rows)
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-footprints-json", type=Path, required=True)
    parser.add_argument("--anthropogenic-csv", type=Path, required=True)
    parser.add_argument("--level3-geojson", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    footprints = json.loads(
        args.native_footprints_json.read_text(encoding="utf-8")
    )["species"]
    level1 = load_level1(args.level3_geojson)

    rows = []
    with args.anthropogenic_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            units = [code for code in footprints[row["species"]] if code in level1]
            counts = Counter(level1[code] for code in units)
            if not counts:
                continue
            dominant_region, dominant_count = sorted(
                counts.items(), key=lambda item: (-item[1], item[0])
            )[0]
            rows.append(
                {
                    "species": row["species"],
                    "host_family_count": float(row["host_family_count"]),
                    "log_resource_expansion": float(row["log_resource_expansion"]),
                    "introduced_added_units": int(row["introduced_added_units"]),
                    "dominant_region": dominant_region,
                    "dominant_region_share": dominant_count / len(units),
                }
            )

    by_region = {}
    for region in sorted(set(row["dominant_region"] for row in rows)):
        subset = [row for row in rows if row["dominant_region"] == region]
        strict = [
            row for row in subset if row["dominant_region_share"] >= 0.70
        ]
        by_region[region] = {
            "name": REGION_NAMES.get(region, region),
            "dominant_region": summarize(subset),
            "dominant_region_share_ge_0_70": summarize(strict),
        }

    payload = {
        "schema": "chocho_butterfly_regional_robustness_v0.1",
        "status": "POSTHOC_GEOGRAPHIC_BIAS_SENSITIVITY",
        "species": len(rows),
        "overall": summarize(rows),
        "dominant_native_resource_region": by_region,
        "interpretation": (
            "Regional stratification asks whether the near-zero global association "
            "between host-family breadth and proportional anthropogenic resource "
            "expansion is driven by the geographic composition of the LepTraits/HOSTS panel."
        ),
        "claim_boundary": (
            "Dominant native host-resource region is a resource-geographic stratifier, "
            "not the butterfly's complete native range and not a correction for uneven "
            "data collection intensity within regions."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
