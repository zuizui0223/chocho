#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

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


def correlation(a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def residualize_exact_count(values, counts):
    values = np.asarray(values, dtype=float)
    counts = np.asarray(counts)
    residuals = values.copy()
    for count in np.unique(counts):
        idx = np.flatnonzero(counts == count)
        residuals[idx] = values[idx] - np.mean(values[idx])
    return residuals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=9999)
    args = parser.parse_args()

    rows = []
    with args.metrics_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["host_taxonomy_lower_bound_adequate"] != "True":
                continue
            if int(row["introduced_added_units"]) <= 0:
                continue
            rows.append(row)

    host_family = np.asarray([float(r["host_family_count"]) for r in rows])
    host_count = np.asarray([int(r["resolved_host_species"]) for r in rows])
    x_rank = average_ranks(host_family)
    x_residual = residualize_exact_count(x_rank, host_count)
    groups = [np.flatnonzero(host_count == c) for c in np.unique(host_count)]
    metric_names = (
        "effective_contributor_number",
        "maximum_single_host_fractional_share",
    )

    y_residuals = {}
    observed = {}
    for name in metric_names:
        y = np.asarray([float(r[name]) for r in rows])
        marginal = correlation(average_ranks(host_family), average_ranks(y))
        residual = residualize_exact_count(average_ranks(y), host_count)
        y_residuals[name] = residual
        observed[name] = {
            "marginal_spearman": marginal,
            "exact_host_count_adjusted_rank_correlation": correlation(x_residual, residual),
        }

    rng = np.random.default_rng(20260928)
    null = {name: [] for name in metric_names}
    for _ in range(args.permutations):
        permuted = x_rank.copy()
        for idx in groups:
            if len(idx) > 1:
                permuted[idx] = rng.permutation(permuted[idx])
        x_perm = residualize_exact_count(permuted, host_count)
        for name in metric_names:
            null[name].append(correlation(x_perm, y_residuals[name]))

    for name in metric_names:
        values = np.asarray(null[name], dtype=float)
        obs = observed[name]["exact_host_count_adjusted_rank_correlation"]
        observed[name].update(
            {
                "permutation_p_two_sided": float(
                    (1 + np.sum(np.abs(values) >= abs(obs))) / (args.permutations + 1)
                ),
                "null_q025": float(np.quantile(values, 0.025)),
                "null_median": float(np.median(values)),
                "null_q975": float(np.quantile(values, 0.975)),
            }
        )

    payload = {
        "schema": "chocho_butterfly_architecture_hostcount_adjusted_v0.1",
        "status": "POSTHOC_STRUCTURAL_DIAGNOSTIC",
        "species": len(rows),
        "permutations": args.permutations,
        "adjustment": "Exact resolved host-species count fixed effects on rank scale; host-family ranks permuted only within exact host-count strata.",
        "metrics": observed,
        "interpretation": "The marginal family-breadth architecture gradients largely disappear once exact resolved host-species richness is held constant.",
        "claim_boundary": "This diagnoses structural dependence on host count; it does not test whether the identities or geographic redistribution histories of observed hosts differ from matched alternative hosts.",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
