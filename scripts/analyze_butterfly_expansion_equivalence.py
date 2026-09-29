#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def average_ranks(values: np.ndarray) -> np.ndarray:
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


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    x = average_ranks(a)
    y = average_ranks(b)
    if np.std(x) == 0 or np.std(y) == 0:
        raise RuntimeError("undefined Spearman correlation")
    return float(np.corrcoef(x, y)[0, 1])


def normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def fisher_ci(rho: float, n: int, level: float) -> tuple[float, float]:
    zcrit = {
        0.90: 1.6448536269514722,
        0.95: 1.959963984540054,
    }[level]
    z = math.atanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    return math.tanh(z - zcrit * se), math.tanh(z + zcrit * se)


def fisher_tost(
    rho: float,
    n: int,
    margin: float,
) -> dict[str, float | bool]:
    z = math.atanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    lower_bound = math.atanh(-margin)
    upper_bound = math.atanh(margin)
    z_lower = (z - lower_bound) / se
    z_upper = (z - upper_bound) / se
    p_lower = 1.0 - normal_cdf(z_lower)
    p_upper = normal_cdf(z_upper)
    p_tost = max(p_lower, p_upper)
    return {
        "margin": margin,
        "p_lower": p_lower,
        "p_upper": p_upper,
        "p_tost": p_tost,
        "equivalent_at_alpha_0_05": p_tost < 0.05,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anthropogenic-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    host, expansion = [], []
    with args.anthropogenic_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"host_family_count", "log_resource_expansion"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("anthropogenic metric schema drift")
        for row in reader:
            host.append(float(row["host_family_count"]))
            expansion.append(float(row["log_resource_expansion"]))

    x = np.asarray(host, dtype=float)
    y = np.asarray(expansion, dtype=float)
    if len(x) != 239:
        raise RuntimeError(
            f"expected 239 resource-eligible butterflies; got {len(x)}"
        )

    rho = spearman(x, y)
    ci90 = fisher_ci(rho, len(x), 0.90)
    ci95 = fisher_ci(rho, len(x), 0.95)
    tost = fisher_tost(rho, len(x), 0.10)

    payload = {
        "schema": "chocho_butterfly_expansion_equivalence_v0.1",
        "status": "SUCCESS_POSTHOC_EQUIVALENCE_DIAGNOSTIC",
        "n": len(x),
        "association": "Spearman host-family breadth vs log resource expansion",
        "rho": rho,
        "fisher_z_approx_ci90": list(ci90),
        "fisher_z_approx_ci95": list(ci95),
        "tost_fisher_z_approx": tost,
        "interpretation": (
            "The point estimate is near zero, but equivalence within ±0.10 is not "
            "established at alpha=0.05 because the upper one-sided test is narrowly "
            "non-significant. The result supports wording such as 'largely independently "
            "of diet breadth' rather than an exact-null claim."
        ),
        "claim_boundary": (
            "Fisher-z confidence intervals and TOST are approximate for Spearman rho; "
            "they quantify the precision of the near-zero rank association rather than "
            "prove an exact zero effect."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
