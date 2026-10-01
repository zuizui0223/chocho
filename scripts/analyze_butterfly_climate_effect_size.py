#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np


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


def partial_spearman(response, predictor, control):
    ry = average_ranks(response)
    rx = average_ranks(predictor)
    rc = average_ranks(control)
    design = np.column_stack([np.ones(len(rc)), rc])
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    if np.std(ey) == 0 or np.std(ex) == 0:
        return None
    return float(np.corrcoef(ey, ex)[0, 1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-result-json", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=30000)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.primary_result_json.read_text(encoding="utf-8"))
    rows = payload["species"]
    y = np.asarray([row["climate_filtering_score"] for row in rows], dtype=float)
    x = np.asarray([row["host_family_count"] for row in rows], dtype=float)
    control = np.log1p(
        np.asarray(
            [row["contemporary_host_resource_units"] for row in rows], dtype=float
        )
    )

    observed = partial_spearman(y, x, control)
    rng = np.random.default_rng(20260928)
    bootstrap = []
    while len(bootstrap) < args.bootstrap:
        indices = rng.integers(0, len(rows), len(rows))
        value = partial_spearman(y[indices], x[indices], control[indices])
        if value is not None and math.isfinite(value):
            bootstrap.append(value)
    bootstrap = np.asarray(bootstrap, dtype=float)

    alpha_z = NormalDist().inv_cdf(0.95)
    power_z = NormalDist().inv_cdf(0.80)
    fisher_se_denominator = len(rows) - 4
    mde = math.tanh(
        (alpha_z + power_z) / math.sqrt(fisher_se_denominator)
    )

    result = {
        "schema": "chocho_butterfly_climate_effect_size_v0.1",
        "status": "POSTHOC_PRECISION_DIAGNOSTIC",
        "species": len(rows),
        "observed_partial_spearman": observed,
        "bootstrap_replicates": len(bootstrap),
        "bootstrap_percentile_95_ci": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "bootstrap_median": float(np.median(bootstrap)),
        "approximate_minimum_detectable_absolute_partial_correlation": {
            "power": 0.80,
            "one_sided_alpha": 0.05,
            "fisher_z_approximation": mde,
            "effective_denominator": fisher_se_denominator,
        },
        "interpretation": (
            "The independent panel rules out neither a moderate negative association "
            "nor a modest positive one. With n=24, the test is primarily powered for "
            "large partial correlations, so the unsupported directional prediction "
            "should be reported as imprecise rather than evidence of equivalence."
        ),
        "claim_boundary": (
            "The confidence interval is a species bootstrap percentile interval and "
            "the minimum detectable effect is a Fisher-z approximation; neither was "
            "part of the frozen confirmatory decision rule."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
