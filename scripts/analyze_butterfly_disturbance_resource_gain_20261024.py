#!/usr/bin/env python3
"""Source-locked, exploratory disturbance-affinity butterfly macro analysis.

Protocol: provenance/reviewer_defenses/protocols/DISTURBANCE_TRAIT_RESOURCE_GAIN_20261024.md
Only frozen three species CSVs; no GBIF query or author-only external source needed.
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
A = "Disturbance-associated (strong)"
V = "Disturbance-avoidant (strong)"
FILES = (
    "data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "data/frozen/figure_sources/s1_resource_descriptors.csv",
    "data/frozen/figure_sources/host_contribution_metrics.csv",
)
EXPECTED_GIT_BLOBS = {
    FILES[0]: "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
    FILES[1]: "09541b5894529cc8629c7515ea9022ded9edead2",
    FILES[2]: "4cf4851499ca2233ff2adadf1c4ee7f0e863f23d",
}
BOOTSTRAPS = 5000
SEED = 20261024


class XorShift32:
    """Simple identical cross-language 32-bit bootstrap sampler."""

    def __init__(self, seed: int):
        self.state = seed & 0xFFFFFFFF
        if not self.state:
            raise ValueError("Seed must be nonzero")

    def next_uint(self) -> int:
        x = self.state
        x ^= (x << 13) & 0xFFFFFFFF
        x ^= x >> 17
        x ^= (x << 5) & 0xFFFFFFFF
        self.state = x & 0xFFFFFFFF
        return self.state

    def index(self, n: int) -> int:
        return (self.next_uint() * n) // (2 ** 32)


def read_source(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        out = {}
        for row in csv.DictReader(fh):
            name = row["species"]
            if name in out:
                raise ValueError(f"Duplicate species: {name} in {path}")
            out[name] = row
        return out


def make_panel() -> list[dict[str, str]]:
    for relpath, expected in EXPECTED_GIT_BLOBS.items():
        data = (ROOT / relpath).read_bytes()
        git_sha = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()
        if git_sha != expected:
            raise ValueError(f"Source blob changed: {relpath}: {git_sha} != {expected}")
    sources = [read_source(ROOT / p) for p in FILES]
    a, traits, host = sources
    if len(a) != 239 or len(host) != 239 or len(traits) != 339:
        raise ValueError("Unexpected frozen data panel sizes")
    if not set(a).issubset(traits) or set(a) != set(host):
        raise ValueError("Species join mismatch")
    rows = []
    for name in a:
        merged = dict(a[name])
        merged.update(traits[name])
        merged.update(host[name])
        if int(merged["native_resource_units"]) <= 0:
            raise ValueError(f"Nonpositive native resource: {name}")
        if (int(merged["native_resource_units"]) +
            int(merged["introduced_added_units"]) !=
            int(merged["contemporary_resource_units"])):
            raise ValueError(f"Resource accounting mismatch: {name}")
        rows.append(merged)
    if sum(int(r["introduced_added_units"]) for r in rows) != 14553:
        raise ValueError("Frozen resource total mismatch")
    return rows


def make_xy(rows: list[dict[str, str]], *, reduced=False, secondary=False):
    X, y, genus, group = [], [], [], []
    for r in rows:
        if r["disturbance_affinity"] not in {A, V}:
            continue
        n = int(r["native_resource_units"])
        added = int(r["introduced_added_units"])
        if secondary and added == 0:
            continue
        covariates = [
            1.0,
            float(r["disturbance_affinity"] == A),
            math.log(n),
            math.log1p(int(r["resolved_host_species"])),
            float(r["host_family_count"]),
        ]
        if not reduced:
            covariates.extend([
                float(r["wing_size_proxy"]),
                float(r["voltinism"] == "B"),
                float(r["voltinism"] == "U"),
            ])
        if secondary:
            covariates.append(math.log1p(added))
            outcome = float(r["fraction_added_units_with_multiple_contributing_hosts"])
        else:
            outcome = math.log(int(r["contemporary_resource_units"]) / n)
        if not all(math.isfinite(x) for x in covariates) or not math.isfinite(outcome):
            continue
        X.append(covariates)
        y.append(outcome)
        genus.append(r["species"].split()[0])
        group.append(r["disturbance_affinity"])
    return np.asarray(X, dtype=float), np.asarray(y, dtype=float), genus, group


def fit_coef(X: np.ndarray, y: np.ndarray):
    if X.shape[0] <= X.shape[1]:
        return None
    coef, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1] or not np.all(np.isfinite(coef)):
        return None
    return float(coef[1])


def analyze(rows, *, reduced=False, secondary=False):
    X, y, genera, labels = make_xy(rows, reduced=reduced, secondary=secondary)
    if len(y) <= X.shape[1]:
        raise ValueError("Insufficient complete-case data")
    groups = defaultdict(list)
    for i, name in enumerate(genera):
        groups[name].append(i)
    cluster_ids = list(groups)
    cluster_indices = [np.asarray(groups[g], dtype=int) for g in cluster_ids]
    n_clusters = len(cluster_indices)
    obs = fit_coef(X, y)
    if obs is None:
        raise ValueError("Singular observed design")
    rng = XorShift32(SEED)
    betas = []
    failures = 0
    for _ in range(BOOTSTRAPS):
        ix = np.concatenate(
            [cluster_indices[rng.index(n_clusters)] for _ in range(n_clusters)]
        )
        coef = fit_coef(X[ix], y[ix])
        if coef is None:
            failures += 1
        else:
            betas.append(coef)
    betas = np.asarray(betas)
    if not len(betas):
        raise ValueError("All bootstrap models singular")
    lo, hi = np.quantile(betas, [0.025, 0.975], method="linear").tolist()
    loo = []
    for cluster in cluster_ids:
        ix = np.asarray([i for i, g in enumerate(genera) if g != cluster], dtype=int)
        b = fit_coef(X[ix], y[ix])
        if b is not None:
            loo.append(b)
    which_a = [i for i, label in enumerate(labels) if label == A]
    which_v = [i for i, label in enumerate(labels) if label == V]
    return {
        "n": len(y),
        "n_genera": n_clusters,
        "n_associated": len(which_a),
        "n_avoidant": len(which_v),
        "raw_mean_difference": float(np.mean(y[which_a]) - np.mean(y[which_v])),
        "adjusted_coefficient": obs,
        "genus_bootstrap_95_interval": [lo, hi],
        "bootstrap_positive_fraction": float(np.mean(betas > 0)),
        "bootstrap_draws": BOOTSTRAPS,
        "bootstrap_successes": len(betas),
        "bootstrap_singular": failures,
        "leave_one_genus_out_range": [min(loo), max(loo)],
        "response": (
            "fraction_added_units_with_multiple_contributing_hosts"
            if secondary else "log(contemporary_resource_units/native_resource_units)"
        ),
        "control_set": "reduced" if reduced else "full",
    }


def group_readout(rows, label):
    r = [x for x in rows if x["disturbance_affinity"] == label]
    y = np.asarray([
        math.log(int(x["contemporary_resource_units"]) / int(x["native_resource_units"]))
        for x in r
    ], dtype=float)
    return {
        "n": len(r),
        "mean_log_ratio": float(np.mean(y)),
        "median_log_ratio": float(np.median(y)),
        "geometric_mean_ratio": float(math.exp(np.mean(y))),
        "median_ratio": float(math.exp(np.median(y))),
        "mean_added_regions": float(np.mean([int(x["introduced_added_units"]) for x in r])),
        "median_added_regions": float(np.median([int(x["introduced_added_units"]) for x in r])),
        "zero_gain_species": sum(int(x["introduced_added_units"]) == 0 for x in r),
    }


def result():
    all_rows = make_panel()
    pairs = [x for x in all_rows if x["disturbance_affinity"] in {A, V}]
    adequate = [x for x in pairs if x["host_taxonomy_lower_bound_adequate"] == "True"]
    return {
        "status": "EXPLORATORY_TRAIT_SORTING_PRIMARY_SUGGESTIVE_ADEQUATE_HOST_SENSITIVITY_UNRESOLVED",
        "source_paths": list(FILES),
        "source_panel_species": len(all_rows),
        "source_added_resource_units": sum(int(r["introduced_added_units"]) for r in all_rows),
        "prespecified_groups": {A: group_readout(pairs, A), V: group_readout(pairs, V)},
        "analyses": {
            "primary_full": analyze(pairs),
            "primary_reduced": analyze(pairs, reduced=True),
            "primary_taxonomy": analyze(adequate),
            "secondary_full": analyze(pairs, secondary=True),
        },
        "limitations": [
            "Potential resource geography only, not recruitment or actual host use",
            "One trait contrast was prespecified on an already examined source panel",
            "Host inventory completeness and correlated taxa limit inference",
            "Bootstrap genus clusters are not a phylogenetic comparative model",
            "Posthoc observational association cannot identify causal mechanism",
        ],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=None, help="Write full result JSON")
    ap.add_argument("--compare-json", type=Path, default=None, help="Verify frozen receipt")
    args = ap.parse_args()
    output = result()
    if args.compare_json:
        reference = json.loads(args.compare_json.read_text(encoding="utf-8"))
        for key in output["analyses"]:
            actual, expected = output["analyses"][key], reference["analyses"][key]
            for field in ("n", "n_genera", "n_associated", "n_avoidant",
                          "bootstrap_successes", "bootstrap_singular"):
                if actual[field] != expected[field]:
                    raise AssertionError(f"{key}.{field}: {actual[field]} != {expected[field]}")
            for field in ("adjusted_coefficient", "raw_mean_difference",
                          "bootstrap_positive_fraction"):
                if not math.isclose(actual[field], expected[field], abs_tol=1e-6):
                    raise AssertionError(f"{key}.{field}: {actual[field]} != {expected[field]}")
            for field in ("genus_bootstrap_95_interval", "leave_one_genus_out_range"):
                if not np.allclose(actual[field], expected[field], rtol=0, atol=1e-6):
                    raise AssertionError(f"{key}.{field}: {actual[field]} != {expected[field]}")
        print("VERIFIED: all four source-locked result receipts match")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    else:
        print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
