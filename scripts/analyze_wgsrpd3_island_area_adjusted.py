#!/usr/bin/env python3
"""Area-adjusted descriptive island–mainland resource opportunity checks.

Uses frozen Q1 region data and pinned WGSRPD3 geometry. This is not a causal
effect estimate or an independent confirmation of historical colonization.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyproj import Geod
from shapely.geometry import shape


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def area_by_code(path: Path):
    obj = json.loads(path.read_text(encoding="utf-8"))
    geod = Geod(ellps="WGS84")
    areas = defaultdict(float)
    for feat in obj["features"]:
        code = str((feat.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        geo = feat.get("geometry")
        if not code or not geo:
            continue
        polygon = shape(geo)
        area, _ = geod.geometry_area_perimeter(polygon)
        areas[code] += abs(area) / 1e6
    return dict(areas)


def mean(seq):
    return sum(seq) / len(seq)


def match_block_diff(data, label_names, blocks):
    by_block = defaultdict(lambda: {"island": [], "mainland": []})
    for r, label, block in zip(data, label_names, blocks):
        by_block[block][label].append(r["share"])
    contrasts = []
    for name, groups in sorted(by_block.items()):
        if groups["island"] and groups["mainland"]:
            contrasts.append({
                "stratum": name,
                "island_n": len(groups["island"]),
                "mainland_n": len(groups["mainland"]),
                "difference": mean(groups["island"]) - mean(groups["mainland"]),
                "weight": min(len(groups["island"]), len(groups["mainland"])),
            })
    weight_total = sum(x["weight"] for x in contrasts)
    return (None if weight_total == 0 else
            sum(x["weight"] * x["difference"] for x in contrasts) / weight_total), contrasts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--q1-regions", required=True, type=Path)
    p.add_argument("--q1-summary", required=True, type=Path)
    p.add_argument("--wgsrpd3-geojson", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    summaries = json.loads(a.q1_summary.read_text(encoding="utf-8"))
    region_rows = list(rows(a.q1_regions))
    if len(region_rows) != 369 or summaries["panel_total"]["n_panel_mismatches"] != 0:
        raise ValueError("Q1 region counts or source audit mismatch")
    area = area_by_code(a.wgsrpd3_geojson)
    missing = sorted(set(r["area_code_l3"] for r in region_rows) - area.keys())
    zero_areas = sorted(k for k, v in area.items() if v <= 0)
    if missing:
        raise ValueError(f"Missing official WGSRPD3 areas: {missing}")
    if zero_areas:
        raise ValueError(f"Non-positive WGSRPD3 areas: {zero_areas}")

    frame = []
    for r in region_rows:
        if r["class"] == "mixed" or int(r["contemporary"]) == 0:
            continue
        code = r["area_code_l3"]
        frame.append({
            "code": code, "level1": r["level1"],
            "class": r["class"], "area_km2": area[code],
            "ln_area": math.log(area[code]),
            "share": int(r["added"]) / int(r["contemporary"]),
            "native": int(r["native"]),
            "contemporary": int(r["contemporary"]),
            "added": int(r["added"]),
        })
    labels = [r["class"] for r in frame]
    y = np.array([r["share"] for r in frame])
    island = np.array([1 if r["class"] == "island" else 0 for r in frame])
    log_area = np.array([r["ln_area"] for r in frame])
    level1 = sorted({r["level1"] for r in frame})
    level1_indicators = [np.array([1.0 if r["level1"] == level else 0.0 for r in frame])
                         for level in level1[1:]]
    X = np.column_stack([np.ones(len(frame)), island, log_area, *level1_indicators])
    coefficients = np.linalg.lstsq(X, y, rcond=None)[0]
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise RuntimeError("Non-full-rank design matrix")
    # Within same level1, match areas by fixed quartile bins of global log-area.
    q = np.quantile(log_area, [.25, .50, .75])
    blocks = [f"{r['level1']}|{int(np.searchsorted(q, r['ln_area'], side='right'))}"
              for r in frame]
    observed, strata = match_block_diff(frame, labels, blocks)
    if observed is None:
        raise RuntimeError("No island–mainland area-overlap strata")
    inds_by_block = defaultdict(list)
    for i, block in enumerate(blocks):
        inds_by_block[block].append(i)
    rng = random.Random(20261102)
    draws = []
    for _ in range(4999):
        perm = labels.copy()
        for inds in inds_by_block.values():
            vals = [perm[i] for i in inds]
            rng.shuffle(vals)
            for j, val in zip(inds, vals):
                perm[j] = val
        draw, _ = match_block_diff(frame, perm, blocks)
        if draw is not None:
            draws.append(draw)

    p_two = (1 + sum(abs(x) >= abs(observed) for x in draws)) / (1 + len(draws))
    result = {
        "schema": "chocho_wgsrpd3_area_adjusted_island_resource_v0.1",
        "source": "Pinned WGSRPD3 Level 3 geometry + exact audited Q1 output",
        "area_metric": "WGS84 geodesic polygon area km2, entire WGSRPD3 Level3 region",
        "active_regions": len(frame),
        "active_islands": sum(r["class"] == "island" for r in frame),
        "active_mainland": sum(r["class"] == "mainland" for r in frame),
        "ols_model": {
            "response": "introduced_added / contemporary opportunity",
            "predictors": ["intercept", "island", "ln(region area km2)", *[f"level1=={x}" for x in level1[1:]]],
            "coefficients": [float(x) for x in coefficients],
            "island_coefficient": float(coefficients[1]),
            "interpretation": "Descriptive conditional association, not independence-adjusted inference",
        },
        "stratum_test": {
            "definition": "Level1 biogeographic region x pooled global log-area quartile",
            "area_quartile_cutoffs_log_km2": [float(x) for x in q],
            "strata_with_both_classes": strata,
            "weighted_island_minus_mainland_share": observed,
            "within_stratum_permutation_n": len(draws),
            "within_stratum_descriptive_p_two_sided": p_two,
            "note": "Only overlap strata contribute; cannot support causal or spatial-independent inference",
        },
        "limits": [
            "WGSRPD3 polygons include coarse region boundaries; natural island area differs from precise ecological land area.",
            "Out-of-range island areas have no matched mainland analogue and are excluded from strict matched strata.",
            "Area and broad biogeographic region do not adjust climate, plant introduction effort, HOSTS observation bias or human land use.",
            "Repeated butterflies across regions create geographic and interaction-network dependence; permutation p values are descriptive."
        ],
    }
    (a.out_dir / "area_adjusted_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (a.out_dir / "region_area_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frame[0].keys()))
        writer.writeheader()
        writer.writerows(frame)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
