#!/usr/bin/env python3
"""Reanalyse fixed 51-island pair observations without claiming port causality.

Uses only the two public, frozen 2026-10-30 source CSV files. The 1,275 pairs
are NOT independent replicates. Geographically restricted row/column
permutations are post-hoc sensitivity diagnostics, not calibrated proof of
causality or species-level establishment.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "provenance/reviewer_defenses/results/island_host_network_20261030"
PAIRS = SOURCE / "island_host_network_pairs.csv"
ISLANDS = SOURCE / "island_host_network_islands.csv"


def read_csv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def rank(values):
    """One-based midranks with exact (frozen-CSV) ties."""
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(order)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        mid = (i + j + 1) / 2.0
        for k in range(i, j):
            out[order[k]] = mid
        i = j
    return out


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def project(values, controls):
    """Residualize against intercept and each control via Gram-Schmidt QR."""
    mean = sum(values) / len(values)
    residual = [float(v - mean) for v in values]
    basis = []
    for raw in controls:
        center = sum(raw) / len(raw)
        v = [float(x - center) for x in raw]
        for unit in basis:
            alpha = dot(v, unit)
            v = [x - alpha * y for x, y in zip(v, unit)]
        norm = math.sqrt(dot(v, v))
        if norm > 1e-8:
            basis.append([x / norm for x in v])
    for unit in basis:
        alpha = dot(residual, unit)
        residual = [x - alpha * y for x, y in zip(residual, unit)]
    return residual


def corr(x, y):
    denominator = math.sqrt(dot(x, x) * dot(y, y))
    return dot(x, y) / denominator if denominator else None


def spearman(rows, x, y):
    return corr(project(rank([r[x] for r in rows]), []),
                project(rank([r[y] for r in rows]), []))


def load_observations():
    eligible = {r["landmass_id"]: r for r in read_csv(ISLANDS)
                if r["eligible"] == "True"}
    ids = sorted(eligible)
    index = {name: i for i, name in enumerate(ids)}
    dist = [[0.0 if i == j else math.inf for j in range(len(ids))]
            for i in range(len(ids))]
    rows = []
    for r in read_csv(PAIRS):
        a, b = r["island_a"], r["island_b"]
        if a not in eligible or b not in eligible:
            raise ValueError("Pair references a non-eligible island")
        ia, ib = index[a], index[b]
        if ia == ib:
            raise ValueError("Self-pair in observed dataset")
        distance = float(r["port_distance_km"])
        dist[ia][ib] = dist[ib][ia] = distance
        ai, bi = eligible[a], eligible[b]
        rows.append({
            "i": ia, "j": ib,
            "h": float(r["host_jaccard"]),
            "y": float(r["new_butterfly_jaccard"]),
            "old": float(r["earlier_butterfly_jaccard"]),
            "dist": distance,
            "ha": int(ai["alien_host_species"]), "hb": int(bi["alien_host_species"]),
            "na": int(ai["new_butterfly_species"]), "nb": int(bi["new_butterfly_species"]),
            "ea": int(ai["earlier_butterfly_species"]), "eb": int(bi["earlier_butterfly_species"]),
            "pa": int(ai["ports"]), "pb": int(bi["ports"]),
        })
    n = len(ids)
    assert n == 51 and len(rows) == n * (n - 1) // 2 == 1275
    assert all(math.isfinite(dist[i][j]) for i in range(n)
               for j in range(n))
    assert len({tuple(sorted((r["i"], r["j"]))) for r in rows}) == len(rows)
    return ids, rows, dist


def complete_link_groups(dist, threshold_km):
    """Geographic clusters fixed using ONLY distances, never response values."""
    groups = [[i] for i in range(len(dist))]
    while True:
        best = (math.inf, -1, -1)
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                furthest = max(dist[i][j] for i in groups[a] for j in groups[b])
                if furthest < best[0]:
                    best = (furthest, a, b)
        if best[1] < 0 or best[0] > threshold_km:
            return groups
        groups[best[1]].extend(groups[best[2]])
        del groups[best[2]]


class Xorshift32:
    """Explicit RNG so JS and Python source audits use the same permutations."""

    def __init__(self, seed):
        if not seed:
            raise ValueError("seed must be nonzero")
        self.state = seed & 0xffffffff

    def randbelow(self, maximum):
        x = self.state
        x ^= (x << 13) & 0xffffffff
        x ^= x >> 17
        x ^= (x << 5) & 0xffffffff
        self.state = x & 0xffffffff
        return self.state % maximum


def permutation_audit(rows, n, xr, yr, groups, repetitions, seed):
    """Permute row+column island labels of residualized butterfly matrix.

    This is an exploratory residual-label QAP-like reference, NOT an exact
    conditional permutation p-value. Spatial block exchangeability and the
    residual's network structure cannot be guaranteed from these data.
    """
    mat = [[0.0] * n for _ in range(n)]
    for k, row in enumerate(rows):
        a, b = row["i"], row["j"]
        mat[a][b] = mat[b][a] = yr[k]
    observed = corr(xr, yr)
    denominator = math.sqrt(dot(xr, xr) * dot(yr, yr))
    rng = Xorshift32(seed)
    labels = list(range(n))
    values = []
    ge = both = 0
    for _ in range(repetitions):
        for group in groups:
            shuffled = group.copy()
            for k in range(len(shuffled) - 1, 0, -1):
                j = rng.randbelow(k + 1)
                shuffled[k], shuffled[j] = shuffled[j], shuffled[k]
            for a, b in zip(group, shuffled):
                labels[a] = b
        t = sum(xr[k] * mat[labels[row["i"]]][labels[row["j"]]]
                for k, row in enumerate(rows)) / denominator
        values.append(t)
        if t >= observed - 1e-12:
            ge += 1
        if abs(t) >= abs(observed) - 1e-12:
            both += 1
    values.sort()
    return {
        "observed": observed, "permutations": repetitions, "seed": seed,
        "one_sided_positive_p": (ge + 1) / (repetitions + 1),
        "two_sided_abs_p": (both + 1) / (repetitions + 1),
        "null_mean": sum(values) / len(values),
        "null_q025": values[int(.025 * (repetitions - 1))],
        "null_q975": values[int(.975 * (repetitions - 1))],
        "movable_islands": sum(len(g) for g in groups if len(g) > 1),
        "groups": len(groups), "largest_group": max(map(len, groups)),
    }


def analyse(repetitions=4999):
    ids, rows, dist = load_observations()
    n = len(ids)
    host_rank = rank([r["h"] for r in rows])
    newer_rank = rank([r["y"] for r in rows])
    controls = [rank([math.log1p(r["dist"]) for r in rows]),
                rank([r["old"] for r in rows])]
    richness = [
        rank([min(r["ha"], r["hb"]) for r in rows]),
        rank([max(r["ha"], r["hb"]) for r in rows]),
        rank([min(r["na"], r["nb"]) for r in rows]),
        rank([max(r["na"], r["nb"]) for r in rows]),
        rank([min(r["ea"], r["eb"]) for r in rows]),
        rank([max(r["ea"], r["eb"]) for r in rows]),
        rank([max(r["pa"], r["pb"]) for r in rows]),
    ]
    xr = project(host_rank, controls)
    yr = project(newer_rank, controls)
    rr = project(host_rank, controls + richness)
    sr = project(newer_rank, controls + richness)
    distances = sorted(r["dist"] for r in rows)
    bounds = [distances[int(len(rows) * q)] for q in (.25, .5, .75)]
    by_distance = []
    for k in range(4):
        sample = [r for r in rows if (r["dist"] <= bounds[0] if k == 0
                  else r["dist"] > bounds[2] if k == 3
                  else bounds[k - 1] < r["dist"] <= bounds[k])]
        by_distance.append({
            "group": k + 1, "n": len(sample), "rho": spearman(sample, "h", "y"),
            "host_pos": sum(r["h"] > 0 for r in sample),
            "new_pos": sum(r["y"] > 0 for r in sample),
        })
    subsets = {
        "host_positive": [r for r in rows if r["h"] > 0],
        "both_positive": [r for r in rows if r["h"] > 0 and r["y"] > 0],
        "new_positive": [r for r in rows if r["y"] > 0],
    }
    global_group = [list(range(n))]
    geo1000 = complete_link_groups(dist, 1000)
    geo2000 = complete_link_groups(dist, 2000)
    geo3000 = complete_link_groups(dist, 3000)
    return {
        "status": "EXPLORATORY_OBSERVATIONAL_ISLAND_NETWORK_SENSITIVITY",
        "source": "frozen 51-island pair and eligible-island CSV, 2026-10-30",
        "eligible_islands": n, "island_pairs_not_independent": len(rows),
        "source_git_blobs": {
            "island_host_network_pairs.csv": "5663d2ffe42419f5c4496b92c9f3bbf76175907b",
            "island_host_network_islands.csv": "18cecde5bf3066f276380167f451f1fa5cc1caa1",
        },
        "observed": {
            "host_new_rho": spearman(rows, "h", "y"),
            "host_prior_rho": spearman(rows, "h", "old"),
            "partial_distance_prior": corr(xr, yr),
            "partial_distance_prior_richness": corr(rr, sr),
            "pairs_host_positive": sum(r["h"] > 0 for r in rows),
            "pairs_new_positive": sum(r["y"] > 0 for r in rows),
            "pairs_both_zero": sum(r["h"] == 0 and r["y"] == 0 for r in rows),
            "within_positive_subsets": {
                key: {"n": len(sample), "rho": spearman(sample, "h", "y")}
                for key, sample in subsets.items()
            },
            "distance_quartile_boundaries_km": bounds,
            "within_distance_quartiles": by_distance,
        },
        "permutation_sensitivity": {
            "distance_prior_global": permutation_audit(rows, n, xr, yr, global_group, repetitions, 31031),
            "distance_prior_geo1000": permutation_audit(rows, n, xr, yr, geo1000, repetitions, 31032),
            "distance_prior_geo2000": permutation_audit(rows, n, xr, yr, geo2000, repetitions, 31033),
            "distance_prior_geo3000": permutation_audit(rows, n, xr, yr, geo3000, repetitions, 31034),
            "richness_global": permutation_audit(rows, n, rr, sr, global_group, repetitions, 31035),
            "richness_geo2000": permutation_audit(rows, n, rr, sr, geo2000, repetitions, 31036),
        },
        "interpretation_limits": [
            "Ports index plants, but no matched nonport botanical host comparison is in this 51-island CSV.",
            "2001-2024 plant inventory overlaps and postdates 2011-2022 butterfly records.",
            "Earlier butterfly similarity correlates more strongly than new-record similarity.",
            "Pair count is not the independent sample size; 51 island nodes are reused.",
            "Distance groups and all post-hoc sensitivities are exploratory, not preregistered.",
            "Residual-label permutations assume unverified network/spatial exchangeability.",
            "Richness adjustments use observed catalogue coverage, not full observer effort.",
            "Neither true colonization, host use, shipping pathways nor relay is measured.",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=4999)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.permutations < 99:
        parser.error("At least 99 permutations for an interpretable exploratory audit")
    result = analyse(args.permutations)
    data = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
