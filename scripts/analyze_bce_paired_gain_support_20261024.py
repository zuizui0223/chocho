#!/usr/bin/env python3
"""Posthoc BCE/Clarke host-catalogue source sensitivity of potential resource gain.

Do not conflate this mean support intensity with the per-cell single-host share.
Requires previously source-verified Actions artifact 12084636746.
Protocol: provenance/reviewer_defenses/protocols/BUTTERFLY_BCE_GAIN_SUPPORT_SOURCE_SENSITIVITY_20261024.md
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
SPECIES_CSV = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
EXPECTED_GIT_SHA = "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"
RUN, ARTIFACT_ID = 38988757959, 12084636746
SEED, N_BOOT = 20261024, 5000


class Xorshift32:
    def __init__(self, seed):
        self.s = seed & 0xffffffff
    def next_index(self, n):
        x = self.s
        x ^= (x << 13) & 0xffffffff
        x ^= x >> 17
        x ^= (x << 5) & 0xffffffff
        self.s = x & 0xffffffff
        return (self.s * n) // 4294967296


def read_pinned_species() -> dict[str, dict]:
    contents = SPECIES_CSV.read_bytes()
    sha = hashlib.sha1(b"blob " + str(len(contents)).encode() + b"\x00" + contents).hexdigest()
    if sha != EXPECTED_GIT_SHA:
        raise ValueError(f"Source species SHA changed unexpectedly: {sha}")
    with SPECIES_CSV.open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    if len(rows) != 239 or len({r["species"] for r in rows}) != 239:
        raise ValueError("Original 239-species names changed")
    return {r["species"]: r for r in rows}


def read_archived_source(artifact_dir: Path) -> tuple[dict, Path]:
    options = []
    for path in artifact_dir.rglob("*.json"):
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except (UnicodeError, json.JSONDecodeError):
            continue
        if isinstance(obj, dict) and obj.get("schema") == "chocho_butterfly_bce_host_inventory_resource_novelty_v0.1":
            options.append((obj, path))
    if len(options) != 1:
        raise ValueError(f"Expected exactly one original BCE source JSON, found {len(options)}")
    obj, path = options[0]
    if obj.get("status") != "SUCCESS_POSTHOC_EUROPE_BIASED_INVENTORY_SENSITIVITY_NOT_SPATIAL_NULL":
        raise ValueError("Unexpected BCE source receipt status")
    if len(obj.get("affected_butterflies_sorted", [])) != 81:
        raise ValueError("BCE changed-species cohort not 81")
    if obj["source"]["candidate_additional_host_links"] != 1027:
        raise ValueError("BCE expected accepted host links changed")
    if obj["base"]["focal_butterflies"] != 239 or obj["BCE_augmented_partial"]["focal_butterflies"] != 239:
        raise ValueError("Source full panel changed")
    for scenario, expected in (
        ("base", {"candidate_incidences": 75346, "native_absorbed_incidences": 46387,
                  "distinct_added_cells": 14553, "redundant_new_cell_incidences": 14406}),
        ("BCE_augmented_partial", {"candidate_incidences": 103349,
                                   "native_absorbed_incidences": 67740,
                                   "distinct_added_cells": 14967,
                                   "redundant_new_cell_incidences": 20642}),
    ):
        for key, value in expected.items():
            if obj[scenario][key] != value:
                raise ValueError(f"Unexpected BCE {scenario} {key}")
    return obj, path


def build_paired(original: dict, metrics: dict[str, dict]):
    out = []
    for z in original["affected_butterflies_sorted"]:
        name = z["species"]
        if name not in metrics:
            raise ValueError(f"Unexpected BCE focal name: {name}")
        source = metrics[name]
        extra = int(z["BCE_candidate_host_links_added"])
        if extra <= 0:
            raise ValueError("An affected row has no newly recognized accepted hosts")
        baseline_hosts = int(source["resolved_host_species"])
        family_count = float(source["host_family_count"])
        genus = name.split()[0]
        scenarios = {}
        for tag in ("base", "augmented"):
            native = int(z[f"{tag}_native"])
            gained = int(z[f"{tag}_added"])
            H, A, R = (int(z[f"{tag}_{k}"]) for k in ("H", "A", "R"))
            if native <= 0 or gained < 0 or H - A != gained + R:
                raise ValueError(f"Host-incidence source identity failed: {name},{tag}")
            host_count = baseline_hosts + (extra if tag == "augmented" else 0)
            if gained > 0:
                scenarios[tag] = {
                    "gain": math.log((native + gained) / native),
                    "log_native": math.log(native),
                    "log_hosts": math.log1p(host_count),
                    "families": family_count,
                    "support": 1.0 + R / gained,
                    "log_support": math.log1p(R / gained),
                    "U": gained,
                    "native": native,
                    "hosts": host_count,
                }
        out.append({
            "species": name,
            "genus": genus,
            "base": scenarios.get("base"),
            "augmented": scenarios.get("augmented"),
        })
    if len(out) != 81 or len({x["species"] for x in out}) != 81:
        raise ValueError("BCE affected sample differs")
    return out


def fit(X, y):
    if len(y) <= X.shape[1]:
        return None
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1] or not np.isfinite(beta).all():
        return None
    return beta


def estimate(rows, which):
    xs, ys, genus = [], [], []
    for row in rows:
        r = row[which]
        if r is None:
            raise ValueError("Paired panel contains a scenario with U=0")
        xs.append([1.0, r["gain"], r["log_native"], r["log_hosts"], r["families"]])
        ys.append(r["log_support"])
        genus.append(row["genus"])
    X, y = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    beta = fit(X, y)
    if beta is None:
        raise ValueError("BCE regression singular")
    groups = defaultdict(list)
    for i, g in enumerate(genus):
        groups[g].append(i)
    idx_groups = [np.asarray(v, dtype=int) for v in groups.values()]
    rnd = Xorshift32(SEED)
    draws = []
    failed = 0
    for _ in range(N_BOOT):
        ind = np.concatenate([idx_groups[rnd.next_index(len(idx_groups))]
                              for __ in idx_groups])
        est = fit(X[ind], y[ind])
        if est is None:
            failed += 1
        else:
            draws.append(float(est[1]))
    if not draws:
        raise ValueError("All BCE bootstrap samples singular")
    loo = []
    for g, v in groups.items():
        mask = np.ones(len(y), dtype=bool)
        mask[np.asarray(v, dtype=int)] = False
        co = fit(X[mask], y[mask])
        if co is not None:
            loo.append(float(co[1]))
    return {
        "species": len(rows),
        "genera": len(groups),
        "slope_log_support_intensity_on_log_gain": float(beta[1]),
        "genus_bootstrap_95_percentile": np.quantile(draws, [.025, .975], method="linear").tolist(),
        "genus_bootstrap_successes": len(draws),
        "genus_bootstrap_failed": failed,
        "positive_bootstrap_fraction": float(np.mean(np.asarray(draws) > 0)),
        "leave_one_genus_coefficient_range": [min(loo), max(loo)],
        "median_introduced_support_incidence_per_new_cell": float(np.median(np.exp(y))),
        "median_geographic_resource_expansion_ratio": float(np.exp(np.median(X[:, 1]))),
        "model": "log(1+R/U)~1+log_growth+log_native+log1p_resolved_hosts+host_family_breadth",
    }


def make_result(artifact_dir):
    source, path = read_archived_source(artifact_dir)
    rows = build_paired(source, read_pinned_species())
    paired = [z for z in rows if z["base"] is not None and z["augmented"] is not None]
    before, after = estimate(paired, "base"), estimate(paired, "augmented")
    delta_E = np.asarray([r["augmented"]["gain"] - r["base"]["gain"] for r in paired])
    delta_S = np.asarray([r["augmented"]["log_support"] - r["base"]["log_support"] for r in paired])
    full_totals = {
        "base_81": source["subgroups"]["81_BCE_candidate_augmented_butterflies"]["base"],
        "augmented_81": source["subgroups"]["81_BCE_candidate_augmented_butterflies"]["BCE_augmented"],
    }
    return {
        "schema": "chocho_bce_paired_gain_support_intensity_diagnostic_20261024",
        "status": "POSTHOC_SOURCE_SENSITIVITY_NOT_DIRECT_SINGLE_HOST_FRACTION_REPLICATION",
        "source_run": RUN, "source_artifact_id": ARTIFACT_ID,
        "source_artifact_name": path.name,
        "source_uncompressed_json_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "original_239_record_git_blob": EXPECTED_GIT_SHA,
        "all_augmented_species": len(rows),
        "positive_gain_both_catalogues": len(paired),
        "nonpaired_due_to_zero_U": len(rows) - len(paired),
        "analyses": {"before_extra_links": before, "after_extra_links": after},
        "paired_source_change": {
            "median_change_log_geographic_gain": float(np.median(delta_E)),
            "median_change_log_support_intensity": float(np.median(delta_S)),
            "mean_change_log_geographic_gain": float(np.mean(delta_E)),
            "mean_change_log_support_intensity": float(np.mean(delta_S)),
            "fraction_species_with_higher_intensity_after": float(np.mean(delta_S > 0)),
            "correlation_within_species_changes": float(np.corrcoef(delta_E, delta_S)[0, 1])
                if np.std(delta_E) > 0 and np.std(delta_S) > 0 else None,
        },
        "source_full_81_aggregates": full_totals,
        "scientific_limits": [
            "BCE catalogue Europe-biased, not independent global cohort",
            "Mean introduced supports per cell is not the fraction of cells with ≥2 supports",
            "H, A, R and U are mathematically linked components of same map",
            "A changed host list can move regions from newly added to natively supported",
            "No field verification of local plants, larval survival or butterfly recruitment",
            "Posthoc exploratory source sensitivity, not causal host buffering",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--artifact-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = make_result(args.artifact_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "paired_species": result["positive_gain_both_catalogues"],
        "before": result["analyses"]["before_extra_links"]["slope_log_support_intensity_on_log_gain"],
        "after": result["analyses"]["after_extra_links"]["slope_log_support_intensity_on_log_gain"],
    }))


if __name__ == "__main__":
    main()
