#!/usr/bin/env python3
"""Post-hoc geographic redundancy diagnostic from the frozen chocho figure CSVs.

Run:
python scripts/analyze_butterfly_resource_redundancy_hostcount.py \
  --expansion-csv data/frozen/figure_sources/anthropogenic_species_metrics.csv \
  --contributions-csv data/frozen/figure_sources/host_contribution_metrics.csv \
  --output-json provenance/reviewer_defenses/results/butterfly_resource_redundancy_hostcount_20261011.json

This is a structural sensitivity, not evidence of demographic dependence.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def csv_index(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"species", "introduced_added_units"} <= set(reader.fieldnames or []):
            # The contribution CSV has a compatible introduced_added_units column.
            raise ValueError(f"Unexpected CSV columns: {path}")
        out = {}
        for row in reader:
            key = row["species"].strip()
            if key in out:
                raise ValueError(f"Duplicate species {key}")
            out[key] = row
        return out


def ranks(values: list[float]) -> list[float]:
    sorted_index = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    lo = 0
    while lo < len(values):
        hi = lo + 1
        while hi < len(values) and values[sorted_index[hi]] == values[sorted_index[lo]]:
            hi += 1
        rank = (lo + hi + 1) / 2.0  # average of 1-based tied ranks
        for k in range(lo, hi):
            out[sorted_index[k]] = rank
        lo = hi
    return out


def correlation(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or not x:
        raise ValueError("Invalid correlation vectors")
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sx = sum((a - mx) ** 2 for a in x)
    sy = sum((b - my) ** 2 for b in y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy) ** 0.5


def residualize_by_groups(values: list[float], groups: list[list[int]]) -> list[float]:
    out = values.copy()
    for indices in groups:
        mean = sum(values[i] for i in indices) / len(indices)
        for i in indices:
            out[i] -= mean
    return out


class Mulberry32:
    """Cross-runtime deterministic PRNG for the reported exact-count permutations."""

    def __init__(self, seed: int):
        self.state = seed & 0xFFFFFFFF

    def random(self) -> float:
        self.state = (self.state + 0x6D2B79F5) & 0xFFFFFFFF
        x = self.state
        x = ((x ^ (x >> 15)) * (x | 1)) & 0xFFFFFFFF
        x = (x ^ ((x + (((x ^ (x >> 7)) * (x | 61)) & 0xFFFFFFFF)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        return ((x ^ (x >> 14)) & 0xFFFFFFFF) / (2 ** 32)


def quantile(values: list[float], p: float) -> float:
    vals = sorted(values)
    point = p * (len(vals) - 1)
    lo, hi = int(point), min(int(point) + 1, len(vals) - 1)
    return vals[lo] + (vals[hi] - vals[lo]) * (point - lo)


def prepare(expansion: dict[str, dict[str, str]], contributions: dict[str, dict[str, str]]):
    if set(expansion) != set(contributions):
        raise ValueError("Species IDs differ across frozen inputs")
    rows = []
    for species, row in expansion.items():
        if not int(row["introduced_added_units"]):
            continue
        if row["host_taxonomy_lower_bound_adequate"] != "True":
            continue
        other = contributions[species]
        if int(row["introduced_added_units"]) != int(other["introduced_added_units"]):
            raise ValueError(f"Added-unit count differs for {species}")
        rows.append({
            "species": species,
            "host_family_count": float(row["host_family_count"]),
            "resolved_host_species": int(row["resolved_host_species"]),
            "added_units": int(row["introduced_added_units"]),
            "multiple_fraction": float(other["fraction_added_units_with_multiple_contributing_hosts"]),
        })
    return rows


def analyze(rows, permutations: int, seed: int, bootstraps: int):
    families = ranks([r["host_family_count"] for r in rows])
    redundancies = ranks([r["multiple_fraction"] for r in rows])
    strata = defaultdict(list)
    for i, row in enumerate(rows):
        strata[row["resolved_host_species"]].append(i)
    groups = [strata[n] for n in sorted(strata)]
    x = residualize_by_groups(families, groups)
    y = residualize_by_groups(redundancies, groups)
    observed = correlation(x, y)
    rng = Mulberry32(seed)
    null = []
    for _ in range(permutations):
        shuffled = x.copy()
        for indices in groups:
            for j in range(len(indices) - 1, 0, -1):
                k = int(rng.random() * (j + 1))
                shuffled[indices[j]], shuffled[indices[k]] = (
                    shuffled[indices[k]], shuffled[indices[j]]
                )
        null.append(correlation(shuffled, y))
    exceed = sum(abs(value) + 1e-12 >= abs(observed) for value in null)

    bins = []
    for lo, hi, label in [(1, 2, "1-2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10**10, "10+")]:
        subset = [r for r in rows if lo <= r["resolved_host_species"] <= hi]
        total = sum(r["added_units"] for r in subset)
        single = sum(r["added_units"] * (1 - r["multiple_fraction"]) for r in subset)
        bins.append({
            "host_species_group": label,
            "species": len(subset),
            "added_units": total,
            "single_host_added_units": round(single),
            "single_host_added_share": single / total,
            "species_with_ge90pct_single_host_added_units": sum(
                r["multiple_fraction"] <= 0.1 for r in subset
            ),
        })
    rich = [r for r in rows if r["resolved_host_species"] >= 10]
    bootstrap_rng = Mulberry32(seed + 1)
    bootstrap = []
    for _ in range(bootstraps):
        sampled = [rich[int(bootstrap_rng.random() * len(rich))] for _ in rich]
        num = sum(r["added_units"] * (1 - r["multiple_fraction"]) for r in sampled)
        den = sum(r["added_units"] for r in sampled)
        bootstrap.append(num / den)
    return {
        "status": "EXPLORATORY_POSTHOC_DIAGNOSTIC_NOT_MANUSCRIPT_PRIMARY",
        "n_species": len(rows),
        "filter": "introduced_added_units > 0 and host_taxonomy_lower_bound_adequate == True",
        "marginal_spearman_family_breadth_vs_multiple_host_fraction": correlation(
            families, redundancies
        ),
        "exact_host_species_count_adjusted_rank_correlation": observed,
        "exact_host_species_count_group_permutations": permutations,
        "permutation_seed": seed,
        "permutation_p_two_sided": (exceed + 1) / (permutations + 1),
        "permutation_null_q025": quantile(null, 0.025),
        "permutation_null_q975": quantile(null, 0.975),
        "host_count_bins": bins,
        "rich_host_ge10_species_bootstrap_replicates": bootstraps,
        "rich_host_ge10_species_bootstrap_q025": quantile(bootstrap, 0.025),
        "rich_host_ge10_species_bootstrap_q975": quantile(bootstrap, 0.975),
        "n_informative_species": sum(
            len(indices) for indices in groups
            if len({rows[i]["host_family_count"] for i in indices}) > 1
        ),
        "interpretation": "Local added-resource redundancy is largely related to recorded host-species count. Host-family breadth adds little detectable independent rank-scale association at exactly fixed recorded host-species count. Even butterflies with >=10 recorded hosts have substantial single-recorded-host added geography.",
        "claim_boundary": "The geography is WGSRPD3 plant-range reconstruction, not observed site-level co-occurrence, host suitability, abundance, realized larval use, population sensitivity or removal outcome. The global HOSTS panel is incomplete. Both tests are post-hoc, lack butterfly phylogenetic correction and should not imply a universal zero effect.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expansion-csv", type=Path, default=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv"))
    ap.add_argument("--contributions-csv", type=Path, default=Path("data/frozen/figure_sources/host_contribution_metrics.csv"))
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=9999)
    ap.add_argument("--bootstraps", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261011)
    args = ap.parse_args()
    if args.permutations <= 0 or args.bootstraps <= 0:
        raise ValueError("Resampling counts must be positive")
    rows = prepare(csv_index(args.expansion_csv), csv_index(args.contributions_csv))
    if len(rows) != 191:
        raise ValueError(f"Frozen filter drift: expected 191; got {len(rows)}")
    result = analyze(rows, args.permutations, args.seed, args.bootstraps)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
