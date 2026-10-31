#!/usr/bin/env python3
"""Retrospective 51-island leave-one-island-out Jaccard evaluation.

This tests the out-of-island *descriptive* association between port-indexed
host composition and butterfly compositional similarity. The candidate host
plant records cover 2001-2024, including years AFTER the butterfly outcome;
DO NOT describe this as a temporal forecast, colonization prediction, or
evidence of port route transport.

No new-butterfly richness or any other outcome-derived feature is used.
All predictors are standardized within each training fold.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parents[1] / "island_host_network_20261030"
PAIRS = BASE / "island_host_network_pairs.csv"
ISLANDS = BASE / "island_host_network_islands.csv"

BASIC = ("dist", "old")
EFFORT = ("dist", "old", "e0", "e1", "h0", "h1", "p0", "p1")
MODELS = {"M0": BASIC, "M1": BASIC + ("host",),
          "E0": EFFORT, "E1": EFFORT + ("host",)}


def load_csv(path):
    with path.open(encoding="utf8", newline="") as f:
        return list(csv.DictReader(f))


def load():
    eligible = {x["landmass_id"]: x for x in load_csv(ISLANDS)
                if x["eligible"] == "True"}
    ids = sorted(eligible)
    assert len(ids) == 51
    result = []
    for r in load_csv(PAIRS):
        a, b = eligible[r["island_a"]], eligible[r["island_b"]]
        minmax = lambda key: (math.log1p(min(int(a[key]), int(b[key]))),
                              math.log1p(max(int(a[key]), int(b[key]))))
        e0, e1 = minmax("earlier_butterfly_species")
        h0, h1 = minmax("alien_host_species")
        p0, p1 = minmax("ports")
        result.append({
            "a": r["island_a"], "b": r["island_b"],
            "y": float(r["new_butterfly_jaccard"]),
            "dist": math.log1p(float(r["port_distance_km"])),
            "old": float(r["earlier_butterfly_jaccard"]),
            "host": float(r["host_jaccard"]),
            "e0": e0, "e1": e1, "h0": h0, "h1": h1,
            "p0": p0, "p1": p1,
        })
    assert len(result) == 1275
    return ids, result


def solve(matrix, vector):
    """Partial-pivot Gaussian elimination for a small regularized system."""
    n = len(vector)
    a = [list(row) + [vector[i]] for i, row in enumerate(matrix)]
    for i in range(n):
        pivot = max(range(i, n), key=lambda j: abs(a[j][i]))
        a[i], a[pivot] = a[pivot], a[i]
        if abs(a[i][i]) < 1e-12:
            raise ValueError("singular least-squares system")
        f = a[i][i]
        for k in range(i, n + 1):
            a[i][k] /= f
        for j in range(n):
            if j == i:
                continue
            f = a[j][i]
            for k in range(i, n + 1):
                a[j][k] -= f * a[i][k]
    return [r[n] for r in a]


def fit_predict(train, test, features, lam):
    """Ridge linear regression (intercept unpenalized), clipped to [0, 1]."""
    n = len(train)
    means = [sum(r[f] for r in train) / n for f in features]
    scales = [math.sqrt(sum((r[f] - means[i]) ** 2 for r in train) / n)
              or 1.0 for i, f in enumerate(features)]

    def vector(row):
        return [1.0] + [(row[f] - means[i]) / scales[i]
                        for i, f in enumerate(features)]

    x = [vector(r) for r in train]
    p = len(features) + 1
    gram = [[0.] * p for _ in range(p)]
    target = [0.] * p
    for i, vec in enumerate(x):
        for k in range(p):
            target[k] += vec[k] * train[i]["y"]
            for j in range(p):
                gram[k][j] += vec[k] * vec[j]
    for i in range(1, p):
        gram[i][i] += lam
    beta = solve(gram, target)
    return [min(1., max(0., sum(beta[i] * z for i, z in enumerate(vector(r)))))
            for r in test]


def evaluate(lam=1.0):
    if lam < 0:
        raise ValueError("ridge lambda cannot be negative")
    ids, rows = load()
    fold = {key: [] for key in MODELS}
    for island in ids:
        train = [r for r in rows if r["a"] != island and r["b"] != island]
        test = [r for r in rows if r["a"] == island or r["b"] == island]
        assert len(train) == 1225 and len(test) == 50
        for model, feat in MODELS.items():
            preds = fit_predict(train, test, feat, lam)
            mse = sum((preds[i] - r["y"])**2 for i, r in enumerate(test)) / len(test)
            mae = sum(abs(preds[i] - r["y"]) for i, r in enumerate(test)) / len(test)
            fold[model].append({"island": island, "mse": mse, "mae": mae, "n": len(test)})
    model_stats = {}
    for m, values in fold.items():
        model_stats[m] = {
            "island_equal_mse": sum(x["mse"] for x in values) / len(values),
            "island_equal_mae": sum(x["mae"] for x in values) / len(values),
        }
    comparisons = {}
    for null, expanded in (("M0", "M1"), ("E0", "E1")):
        diffs = [a["mse"] - b["mse"]
                 for a, b in zip(fold[null], fold[expanded])]
        avg = sum(diffs) / len(diffs)
        comparisons[expanded + " vs " + null] = {
            "mse_reduction_pct": 100 * avg / model_stats[null]["island_equal_mse"],
            "mean_mse_reduction": avg,
            "islands_improved": sum(x > 0 for x in diffs),
            "islands_worsened": sum(x < 0 for x in diffs),
            "median_diff": sorted(diffs)[len(diffs)//2],
            "range": [min(diffs), max(diffs)],
        }
    return {
        "status": "RETROSPECTIVE_HOST_COMPOSITION_OUT_OF_ISLAND_DESCRIPTIVE_EVALUATION",
        "data": "2001-2024 port plants; 2011-2022 newly indexed butterflies; no temporal precedence",
        "islands": 51, "pairs": 1275, "folds": 51,
        "training_pairs_per_fold": 1225, "evaluation_pairs_per_fold": 50,
        "ridge_lambda": lam,
        "outcome": "Jaccard similarity of butterfly species first indexed 2011-2022",
        "scoring": "island-equal mean of held-out pairwise squared error (and MAE)",
        "feature_groups": {key: list(value) for key, value in MODELS.items()},
        "model_stats": model_stats, "comparisons": comparisons,
        "limitations": [
            "51 islands are the sampling units; held-out pairs are reused across folds.",
            "Host plant indexing through 2024 can postdate 2011-2022 butterfly observations.",
            "All general non-host flora and actual maritime origin-destination links are absent.",
            "Catalogued host counts are not direct botanical observation-effort measurement.",
            "Predictive improvement in pairwise assemblage similarity does not imply host use, arrival, reproduction or transport.",
            "The richer model excludes outcome-derived new butterfly richness to avoid feature leakage.",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ridge-lambda", type=float, default=1.)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    out = evaluate(args.ridge_lambda)
    encoded = json.dumps(out, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
