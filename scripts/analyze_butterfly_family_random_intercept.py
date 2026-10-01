#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import statsmodels.api as sm


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


def zscore(values):
    x = np.asarray(values, dtype=float)
    return (x - np.mean(x)) / np.std(x)


def load_family(path):
    result = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            result[str(row["Species"]).strip()] = str(row["Family"]).strip()
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--anthropogenic-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    family = load_family(args.leptraits_csv)
    rows = []
    with args.anthropogenic_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            species = row["species"]
            fam = family.get(species)
            if not fam:
                raise RuntimeError(f"missing butterfly Family: {species}")
            rows.append(
                {
                    "species": species,
                    "butterfly_family": fam,
                    "host_family_count": float(row["host_family_count"]),
                    "log_resource_expansion": float(row["log_resource_expansion"]),
                }
            )

    x = zscore(average_ranks([r["host_family_count"] for r in rows]))
    y = zscore(average_ranks([r["log_resource_expansion"] for r in rows]))
    X = sm.add_constant(x)
    groups = np.asarray([r["butterfly_family"] for r in rows], dtype=object)

    model = sm.MixedLM(y, X, groups=groups)
    fit = model.fit(reml=True, method="lbfgs", maxiter=1000, disp=False)

    beta = float(fit.params[1])
    se = float(fit.bse[1])
    pvalue = float(fit.pvalues[1])
    ci = [float(v) for v in fit.conf_int()[1]]
    random_intercept_variance = float(np.asarray(fit.cov_re)[0, 0])
    residual_variance = float(fit.scale)
    icc = random_intercept_variance / (
        random_intercept_variance + residual_variance
    )

    payload = {
        "schema": "chocho_butterfly_family_random_intercept_sensitivity_v0.1",
        "status": "POSTHOC_TAXONOMIC_NONINDEPENDENCE_SENSITIVITY",
        "species": len(rows),
        "butterfly_families": len(set(groups)),
        "model": (
            "Standardized rank(log resource expansion) ~ standardized "
            "rank(host-family breadth) + (1 | butterfly Family)"
        ),
        "fixed_host_breadth": {
            "coefficient": beta,
            "standard_error": se,
            "wald_p_value": pvalue,
            "wald_95_ci": ci,
        },
        "random_intercept": {
            "variance": random_intercept_variance,
            "residual_variance": residual_variance,
            "icc": icc,
        },
        "converged": bool(fit.converged),
        "interpretation": (
            "The host-family-breadth association with proportional resource expansion "
            "remains near zero after including butterfly Family as a random intercept."
        ),
        "claim_boundary": (
            "This satisfies a coarse taxonomic random-effect sensitivity but is not "
            "equivalent to species-level PGLS or a dated phylogenetic covariance model."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
