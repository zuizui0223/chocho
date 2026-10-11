#!/usr/bin/env python3
"""Post-hoc association between resource gain and geographic host redundancy.

Uses ONLY the frozen per-butterfly figure CSVs. Fixed-host-count permutation
controls taxonomic recording portfolio size, not plant geography or climate.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from analyze_butterfly_resource_redundancy_hostcount import (
    Mulberry32,
    correlation,
    csv_index,
    prepare,
    quantile,
    ranks,
    residualize_by_groups,
)


def analyze(expansion_path: Path, contributions_path: Path, permutations=9999, seed=20261011):
    source = csv_index(expansion_path)
    rows = prepare(source, csv_index(contributions_path))
    if len(rows) != 191:
        raise RuntimeError("Exact frozen quality filter drift")
    groups_dict = defaultdict(list)
    for i, r in enumerate(rows):
        groups_dict[r["resolved_host_species"]].append(i)
    groups = [groups_dict[k] for k in sorted(groups_dict)]
    y = ranks([r["multiple_fraction"] for r in rows])
    ry = residualize_by_groups(y, groups)
    rng = Mulberry32(seed)
    endpoints = {}
    for key, source_key in [
        ("relative_log_expansion", "log_resource_expansion"),
        ("absolute_added_region_units", "introduced_added_units"),
    ]:
        x = ranks([float(source[r["species"]][source_key]) for r in rows])
        rx = residualize_by_groups(x, groups)
        observed = correlation(rx, ry)
        random_correlations = []
        for _ in range(permutations):
            values = rx.copy()
            for g in groups:
                for j in range(len(g) - 1, 0, -1):
                    k = int(rng.random() * (j + 1))
                    values[g[j]], values[g[k]] = values[g[k]], values[g[j]]
            random_correlations.append(correlation(values, ry))
        n_extreme = sum(
            abs(v) + 1e-12 >= abs(observed) for v in random_correlations
        )
        endpoints[key] = {
            "marginal_spearman": correlation(x, y),
            "exact_host_species_count_adjusted_rank_correlation": observed,
            "null_q025": quantile(random_correlations, 0.025),
            "null_q975": quantile(random_correlations, 0.975),
            "permutation_p_two_sided": (n_extreme + 1) / (permutations + 1),
        }
    return {
        "schema": "chocho_resource_gain_vs_redundancy_exact_host_count_v0.1",
        "status": "EXPLORATORY_POSTHOC_STRUCTURAL_ASSOCIATION",
        "species": len(rows),
        "permutations": permutations,
        "seed": seed,
        "endpoints": endpoints,
        "interpretation": (
            "Among butterflies with the same number of resolved hosts, greater "
            "introduced-host geographic resource gain is associated with more "
            "added regions supported by multiple recorded hosts."
        ),
        "claim_boundary": (
            "Both axes are derived from the same underlying plant range unions. "
            "The association can be a mathematical and spatial-coverage effect. "
            "It is not a causal tradeoff, demographic insurance result or an "
            "independently validated ecological process. Evaluate against the "
            "plant-family x continent geography-constrained introduction null "
            "before making strong claims."
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expansion-csv", type=Path, default=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv"))
    ap.add_argument("--contributions-csv", type=Path, default=Path("data/frozen/figure_sources/host_contribution_metrics.csv"))
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=9999)
    ap.add_argument("--seed", type=int, default=20261011)
    args = ap.parse_args()
    result = analyze(args.expansion_csv, args.contributions_csv, args.permutations, args.seed)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
