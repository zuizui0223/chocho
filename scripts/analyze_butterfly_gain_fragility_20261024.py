#!/usr/bin/env python3
"""Exploratory resource extent versus host support redundancy, 239 butterfly panel.

Reconstruct source-locked descriptive (not causal) association exactly.
Protocol: provenance/reviewer_defenses/protocols/BUTTERFLY_GEOGRAPHIC_GAIN_FRAGILITY_20261024.md
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "data/frozen/figure_sources/anthropogenic_species_metrics.csv": "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
    "data/frozen/figure_sources/host_contribution_metrics.csv": "4cf4851499ca2233ff2adadf1c4ee7f0e863f23d",
    "data/frozen/figure_sources/s1_resource_descriptors.csv": "09541b5894529cc8629c7515ea9022ded9edead2",
}
SEED, N_BOOT = 20261024, 5000


def data_by_name(path: Path, sha: str) -> dict[str, dict]:
    raw = path.read_bytes()
    identity = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\x00" + raw).hexdigest()
    if identity != sha:
        raise ValueError(f"Unexpected frozen source identity: {path.name}")
    with path.open(newline="", encoding="utf-8") as fd:
        rows = {}
        for row in csv.DictReader(fd):
            sp = row["species"]
            if sp in rows:
                raise ValueError(f"Duplicate accepted butterfly: {sp}")
            rows[sp] = row
        return rows


def panel() -> list[dict]:
    a, h, d = [data_by_name(ROOT / path, sha) for path, sha in SOURCES.items()]
    if (len(a), len(h), len(d)) != (239, 239, 339):
        raise ValueError("Input sample unexpectedly changed")
    if set(a) != set(h) or not set(a).issubset(d):
        raise ValueError("Scientific name join differs")
    joined = []
    for name, row in a.items():
        mix = {**row, **h[name], **d[name]}
        n = int(mix["native_resource_units"])
        added = int(mix["introduced_added_units"])
        now = int(mix["contemporary_resource_units"])
        if n <= 0 or added < 0 or n + added != now:
            raise ValueError(f"Resource union identity failed for {name}")
        if added == 0:
            continue
        multi = float(mix["fraction_added_units_with_multiple_contributing_hosts"])
        peak = float(mix["maximum_single_host_fractional_share"])
        if not (0 <= multi <= 1 and 0 <= peak <= 1):
            raise ValueError(f"Invalid contributor fractions for {name}")
        joined.append({
            "species": name,
            "genus": name.split()[0],
            "adequate": mix["host_taxonomy_lower_bound_adequate"] == "True",
            "gain": math.log(now / n),
            "single_share": 1.0 - multi,
            "max_contributor_share": peak,
            "log_native": math.log(n),
            "log_hosts": math.log1p(int(mix["resolved_host_species"])),
            "families": float(mix["host_family_count"]),
        })
    if len(joined) != 206 or len([v for v in joined if v["adequate"]]) != 191:
        raise ValueError("Expanded butterfly eligible panel differs")
    if sum(int(z["introduced_added_units"]) for z in a.values()) != 14553:
        raise ValueError("Original new resource total unexpectedly changed")
    return joined


class Xorshift32:
    def __init__(self, seed):
        self.s = seed & 0xffffffff

    def sample_index(self, n):
        v = self.s
        v ^= (v << 13) & 0xffffffff
        v ^= v >> 17
        v ^= (v << 5) & 0xffffffff
        self.s = v & 0xffffffff
        return (self.s * n) // 4294967296


def ranks(arr: np.ndarray) -> np.ndarray:
    order = np.argsort(arr, kind="stable")
    ordered = arr[order]
    out = np.empty(len(arr), dtype=float)
    i = 0
    while i < len(arr):
        j = i + 1
        while j < len(arr) and ordered[j] == ordered[i]:
            j += 1
        out[order[i:j]] = (i + j + 1) / 2
        i = j
    return out


def fit(X: np.ndarray, y: np.ndarray):
    if len(y) <= X.shape[1]:
        return None
    beta, residuals, rank, singular = np.linalg.lstsq(X, y, rcond=None)
    if rank < X.shape[1] or not np.isfinite(beta).all():
        return None
    return beta


def arrays(rows, secondary=False):
    X = np.array([[1.0, r["gain"], r["log_native"], r["log_hosts"], r["families"]] for r in rows])
    y = np.array([r["max_contributor_share"] if secondary else r["single_share"] for r in rows])
    return X, y


def examine(rows, secondary=False):
    X, y = arrays(rows, secondary)
    beta = fit(X, y)
    if beta is None:
        raise ValueError("Observed OLS design singular")
    genera = defaultdict(list)
    for i, r in enumerate(rows):
        genera[r["genus"]].append(i)
    groups = [np.array(idx, dtype=int) for idx in genera.values()]
    sample = Xorshift32(SEED)
    bootstrap, failed = [], 0
    for _ in range(N_BOOT):
        idx = np.concatenate([groups[sample.sample_index(len(groups))] for __ in groups])
        coef = fit(X[idx], y[idx])
        if coef is None:
            failed += 1
        else:
            bootstrap.append(float(coef[1]))
    if not bootstrap:
        raise ValueError("All bootstrap fits failed")
    ci = np.quantile(bootstrap, [.025, .975], method="linear").tolist()
    leave_one_genus = []
    for g in genera.values():
        keep = np.ones(len(rows), dtype=bool)
        keep[np.array(g, dtype=int)] = False
        coef = fit(X[keep], y[keep])
        if coef is not None:
            leave_one_genus.append(float(coef[1]))

    control = X[:, [0, 2, 3, 4]]
    residual_gain = X[:, 1] - control @ np.linalg.lstsq(control, X[:, 1], rcond=None)[0]
    residual_y = y - control @ np.linalg.lstsq(control, y, rcond=None)[0]
    terc1, terc2 = np.quantile(X[:, 1], [1/3, 2/3], method="linear")
    low = X[:, 1] <= terc1
    high = X[:, 1] >= terc2
    return {
        "species": len(rows),
        "genera": len(genera),
        "outcome": ("max_single_host_fractional_credit" if secondary else
                    "fraction_added_cells_single_recorded_host"),
        "model": "equal_species_OLS_y~1+log_gain+log_native_units+log1p_recorded_hosts+host_family_count",
        "coefficient_log_gain": float(beta[1]),
        "all_coefficients": beta.tolist(),
        "raw_spearman_gain_vs_outcome": float(np.corrcoef(ranks(X[:, 1]), ranks(y))[0, 1]),
        "partial_residual_pearson": float(np.corrcoef(residual_gain, residual_y)[0, 1]),
        "genus_bootstrap_95_percentile": ci,
        "genus_bootstrap_draws": N_BOOT,
        "genus_bootstrap_successful": len(bootstrap),
        "genus_bootstrap_failed": failed,
        "positive_coefficient_bootstrap_fraction": float(np.mean(np.array(bootstrap) > 0)),
        "leave_one_genus_out_coefficient_range": [min(leave_one_genus), max(leave_one_genus)],
        "median_resource_ratio": float(np.exp(np.median(X[:, 1]))),
        "median_outcome": float(np.median(y)),
        "exploratory_terciles": {
            "low_n": int(np.sum(low)), "high_n": int(np.sum(high)),
            "low_mean_outcome": float(np.mean(y[low])),
            "high_mean_outcome": float(np.mean(y[high])),
        },
    }


def result():
    rows = panel()
    return {
        "schema": "chocho_geographic_gain_vs_host_support_fragility_macro_20261024",
        "status": "PRESPECIFIED_POSITIVE_FRAGILITY_HYPOTHESIS_FALSIFIED_EXPLORATORY_NEGATIVE_CATALOGUE_ASSOCIATION",
        "source_git_blob_sha": SOURCES,
        "full_butterfly_panel": 239,
        "added_resource_regions_total": 14553,
        "expanded_species": len(rows),
        "analyses": {
            "primary_all": examine(rows),
            "primary_adequate_catalogue": examine([x for x in rows if x["adequate"]]),
            "secondary_max_host_share": examine(rows, secondary=True),
        },
        "poaceae_status": "NOT_EVALUABLE_NO_FROZEN_PER_SPECIES_POACEAE_HOST_FLAG_IN_THREE_INPUT_CSVS",
        "interpretive_limits": [
            "Extent and redundancy are derived from the same union-of-known-hosts resource map",
            "Evidence is about recorded introduced hosts supporting potential regions, not larval performance",
            "Taxonomic catalogue incompleteness can overstate single-host support",
            "Butterfly genus clusters are not a phylogenetic random effect or independent study replication",
            "No time of plant arrival, butterfly abundance, adult emergence or conservation intervention measured",
        ],
    }


def verify(a, b):
    for panel, actual in a["analyses"].items():
        expected = b["analyses"][panel]
        for key in ("species", "genera", "genus_bootstrap_draws",
                    "genus_bootstrap_successful", "genus_bootstrap_failed"):
            if actual[key] != expected[key]:
                raise AssertionError(f"{panel}.{key}: count mismatch")
        for key in ("coefficient_log_gain", "raw_spearman_gain_vs_outcome",
                    "partial_residual_pearson", "positive_coefficient_bootstrap_fraction",
                    "median_resource_ratio", "median_outcome"):
            if not math.isclose(actual[key], expected[key], abs_tol=1e-7):
                raise AssertionError(f"{panel}.{key}: {actual[key]} vs {expected[key]}")
        for key in ("genus_bootstrap_95_percentile", "leave_one_genus_out_coefficient_range"):
            if not np.allclose(actual[key], expected[key], atol=1e-7, rtol=0):
                raise AssertionError(f"{panel}.{key}: range mismatch")
    print("PASS: source hashes, accounting, regressions, 15000 genus bootstrap fits verified")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    p.add_argument("--verify", type=Path)
    args = p.parse_args()
    r = result()
    if args.verify:
        verify(r, json.loads(args.verify.read_text(encoding="utf-8")))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(r, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    else:
        print(json.dumps(r, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
