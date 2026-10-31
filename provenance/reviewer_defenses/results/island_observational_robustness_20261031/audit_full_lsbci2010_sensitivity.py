#!/usr/bin/env python3
"""Sensitivity of island-pair composition to a *country-pair* shipping index.

Uses already source-verified full_lsbci2010_island_pairs.csv (49 islands,
1,130 eligible island pairs, 325 repeated country-pair LSBCI scores). This is
NOT observed port-to-port voyages and does not identify shipping mediation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
SOURCE = BASE / "full_lsbci2010_island_pairs.csv"
RECEIPT = BASE / "full_lsbci2010_audit.json"


def records(source=SOURCE):
    with Path(source).open(encoding="utf8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1275
    selected = [r for r in rows if r["status"] == "INDEX_2010_OBSERVED"]
    assert len(selected) == 1130
    for r in selected:
        r["log_port_distance"] = math.log1p(float(r["port_distance_km"]))
    return selected


def rank(values):
    ordered = sorted(range(len(values)), key=lambda i: values[i])
    out = np.empty(len(values), dtype=float)
    k = 0
    while k < len(values):
        j = k + 1
        while j < len(values) and values[ordered[j]] == values[ordered[k]]:
            j += 1
        out[ordered[k:j]] = (k + 1 + j) / 2
        k = j
    return out


def partial(rows, left, right, controls=()):
    v = {key: rank([float(x[key]) for x in rows])
         for key in set([left, right, *controls])}
    x, y = v[left], v[right]
    if controls:
        design = np.column_stack([np.ones(len(rows))] + [v[k] for k in controls])
        x = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
        y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return float(np.corrcoef(x, y)[0, 1])


def summary(rows, shipping):
    basics = ("log_port_distance", "earlier_butterfly_jaccard")
    return {
        "ship_host_raw_rho": partial(rows, shipping, "host_jaccard"),
        "ship_host_adjusted": partial(rows, shipping, "host_jaccard", basics),
        "ship_butterfly_raw_rho": partial(rows, shipping, "new_butterfly_jaccard"),
        "ship_butterfly_adjusted_host": partial(
            rows, shipping, "new_butterfly_jaccard", (*basics, "host_jaccard")),
        "host_butterfly_adjusted_shipping": partial(
            rows, "host_jaccard", "new_butterfly_jaccard", (*basics, shipping)),
    }


def evaluate():
    rows = records()
    original = json.loads(RECEIPT.read_text(encoding="utf8"))
    assert original["source_sha256"] == (
        "edab41de33b9a6b4fc568bf454544dedfcbf32117e2bad7ccc77203e1ee6c33e")
    assert len(rows) == original["sample_size"]["island_pairs"]
    assert all(r["lsbci_2006_2010_mean"] != "" for r in rows)

    countries = sorted({c for r in rows for c in (r["country_a"], r["country_b"])})
    full2010 = summary(rows, "lsbci_2010")
    fullmean = summary(rows, "lsbci_2006_2010_mean")
    country_drop = []
    for country in countries:
        sub = [r for r in rows if country not in (r["country_a"], r["country_b"])]
        stats = summary(sub, "lsbci_2010")
        country_drop.append({
            "removed_country": country, "remaining_island_pairs": len(sub),
            "ship_butterfly_adjusted_host": stats["ship_butterfly_adjusted_host"],
            "host_butterfly_adjusted_shipping": stats["host_butterfly_adjusted_shipping"],
        })

    shipping_vals = [z["ship_butterfly_adjusted_host"] for z in country_drop]
    host_vals = [z["host_butterfly_adjusted_shipping"] for z in country_drop]
    return {
        "status": "OBSERVATIONAL_GLOBAL_COUNTRY_BILATERAL_SENSITIVITY_NOT_PORT_VOYAGES",
        "source_audit": "full_lsbci2010_audit.json",
        "source_sha256": original["source_sha256"],
        "island_pairs": len(rows),
        "unique_islands": len({x for r in rows for x in (r["island_a"], r["island_b"])}),
        "countries": len(countries),
        "country_pairs": len({tuple(sorted((r["country_a"], r["country_b"]))) for r in rows}),
        "main_2010": full2010,
        "mean_2006_2010": fullmean,
        "leave_one_country_out": {
            "n_countries": len(countries),
            "ship_adjusted_min_max": [min(shipping_vals), max(shipping_vals)],
            "host_adjusted_min_max": [min(host_vals), max(host_vals)],
            "ship_adjusted_negative_folds": sum(v < 0 for v in shipping_vals),
            "host_adjusted_positive_folds": sum(v > 0 for v in host_vals),
            "per_country": country_drop,
        },
        "interpretation_limits": [
            "Cross-country scheduled liner connectivity is not an observed island-to-island maritime route.",
            "2010 country-pair LSBCI repeats across island pairs and cannot identify independent port effects.",
            "All 1,130 pairs reuse 49 islands and only 26 country nodes; do not use pair-independent p-values.",
            "Compared 2006-2010 and 2010 measures are both source-vintage 2026 archived history.",
            "Host plants indexed 2001-2024 can postdate butterfly first records in 2011-2022.",
            "Historical community, sea barriers, survey effort, regional climate and general flora remain potential confounders.",
            "Observational rank associations cannot distinguish host shipping, passenger transport, self-dispersal or plant-mediated arrival.",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    d = evaluate()
    raw = json.dumps(d, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(raw, encoding="utf8")
    else:
        print(raw)


if __name__ == "__main__":
    main()
