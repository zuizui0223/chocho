#!/usr/bin/env python3
"""Offline Figure 2A-C model replay on projected original inputs (NumPy only)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median

import numpy as np

from audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 import evaluate as one_held
from audit_WCVP126_native_alien_double_unseen_islands_20261117 import evaluate as both_unseen
from audit_WCVP126_uniform_circuit_double_unseen_20261118 import BOOT_N, FOLD_N, SEED, contrasts, matrix
from verify_island_port_native_alien_matched_20261118 import fit

ROOT = Path(__file__).resolve().parent
PORT_FIELDS = ("code", "level1", "area_km2", "log_area", "native", "introduced",
               "wpi_port_count", "log2ports", "resource_active")
FLORA_FIELDS = ("source_region", "destination_region", "native_shared", "native_union",
                "introduced_shared", "introduced_union")
MARINE_FIELDS = ("island1", "island2", "n1", "n2", "a1", "a2", "native", "alien", "log_geo",
                 "sea_detour", "circuit", "port_pair_proxy", "same_level1", "abs_lat")


def read_table(root, name, fields):
    with (root / "data" / name).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames) != fields:
            raise ValueError(f"Unexpected projection schema for {name}")
        return list(reader)


def compare(actual, expected, path="root"):
    """Assert all original analytical values and scope fields, including folds."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise AssertionError(f"Result fields differ at {path}")
        return sum(compare(actual[k], v, f"{path}/{k}") for k, v in expected.items())
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise AssertionError(f"Result length differs at {path}")
        return sum(compare(a, e, f"{path}/{i}") for i, (a, e) in enumerate(zip(actual, expected)))
    if type(expected) is float:
        if type(actual) not in (int, float) or not math.isclose(
            actual, expected, rel_tol=1e-8, abs_tol=1e-11
        ):
            raise AssertionError(f"{path}: {actual!r} differs from {expected!r}")
        return 1
    if type(actual) is not type(expected) or actual != expected:
        raise AssertionError(f"{path}: {actual!r} differs from {expected!r}")
    return 1


def port_models(rows):
    def fit_origins(region_rows):
        return {name: fit(region_rows, name) for name in ("native", "introduced", "difference")}

    active = [row for row in rows if row["resource_active"] == "1"]
    positive = [row for row in rows if int(row["wpi_port_count"]) > 0]
    levels = sorted({int(row["level1"]) for row in rows})
    if (len(rows), len({r["code"] for r in rows}), len(active), len(positive), len(levels)) != (
        126, 126, 123, 91, 9
    ):
        raise ValueError("Figure 2A source cohort changed")
    deletions = []
    for level in levels:
        retained = [row for row in rows if int(row["level1"]) != level]
        deletions.append({"excluded_TDWG_Level1": level, "sample_n": len(retained),
                          "models": fit_origins(retained)})
    ranges = {}
    for outcome in ("native", "introduced", "difference"):
        slopes = [d["models"][outcome]["beta"] for d in deletions]
        ranges[outcome] = [min(slopes), max(slopes)]
    return {
        "model_scopes": {"all_126": fit_origins(rows), "resource_active_123": fit_origins(active)},
        "geographic_sensitivity": {
            "leave_one_TDWG_Level1_out": deletions,
            "range_by_outcome": ranges,
            "port_positive_only_91": fit_origins(positive),
        },
    }


def floristic_overlap(rows, island_codes):
    keys = [tuple(sorted((row["source_region"], row["destination_region"]))) for row in rows]
    if len(keys) != 7875 or len(set(keys)) != 7875:
        raise ValueError("Figure 2B requires 7,875 distinct island pairs")
    if {code for pair in keys for code in pair} != island_codes or any(a == b for a, b in keys):
        raise ValueError("Figure 2B island roster differs from Figure 2A")
    results = {}
    undefined = {}
    for origin in ("native", "introduced"):
        values = []
        undefined[origin] = []
        for key, row in zip(keys, rows):
            shared = int(row[f"{origin}_shared"])
            union = int(row[f"{origin}_union"])
            if not 0 <= shared <= union:
                raise ValueError(f"Invalid {origin} shared/union counts for {key}")
            if union == 0:
                undefined[origin].append(key)
            else:
                values.append(shared / union)
        results[origin] = {"defined": len(values), "median": median(values)}
    if undefined != {"native": [], "introduced": [("BOU", "SSA")]}:
        raise ValueError(f"Unexpected undefined Jaccard: {undefined}")
    return {
        "source_island_pairs": len(keys),
        "defined_pairs": {key: value["defined"] for key, value in results.items()},
        "undefined_introduced_pair": list(undefined["introduced"][0]),
        "median_Jaccard": {key: value["median"] for key, value in results.items()},
    }


def circuit_models(rows, areas):
    # These are the exact pre-fit scalars from the original join, in source row order.
    pairs = [dict(island1=row["island1"], island2=row["island2"],
                  **{name: float(row[name]) for name in MARINE_FIELDS[2:]}) for row in rows]
    keys = {tuple(sorted((entry["island1"], entry["island2"]))) for entry in pairs}
    islands = {code for pair in keys for code in pair}
    if len(pairs) != 5355 or len(keys) != 5355 or len(islands) != 104 or not islands <= areas.keys():
        raise ValueError("Figure 2C requires 5,355 defined source dyads on 104 islands")
    if ("BOU", "SSA") in keys or any(a == b for a, b in keys):
        raise ValueError("Undefined introduced or invalid self-pair entered marine fits")
    shortest = {name: both_unseen(pairs, name, areas) for name in ("native", "alien")}
    for entry in pairs:
        entry["sea_detour"] = entry.pop("circuit")
    native_one = one_held(pairs, "native", areas)
    introduced_one = one_held(pairs, "alien", areas)
    native_both = both_unseen(pairs, "native", areas)
    introduced_both = both_unseen(pairs, "alien", areas)
    for name, fit_result in (("native", native_both), ("alien", introduced_both)):
        reference = shortest[name]
        if not fit_result["test_islands_absent_from_train_endpoints"]:
            raise ValueError("Held island endpoint entered the source training split")
        for before, after in zip(reference["folds"], fit_result["folds"]):
            for key in ("held_island_codes", "n_train_pairs", "n_test_double_unseen_pairs",
                        "MSE_B", "MSE_P"):
                compare(after[key], before[key], f"{name} fixed circuit/shortest fold {key}")
    for native_fold, introduced_fold in zip(native_both["folds"], introduced_both["folds"]):
        for key in ("held_island_codes", "n_train_pairs", "n_test_double_unseen_pairs"):
            compare(native_fold[key], introduced_fold[key], f"paired strict fold {key}")
    nmat = matrix(native_both)
    imat = matrix(introduced_both)
    native_gain = contrasts(nmat)["circuit_beyond_geography_pct"]
    introduced_gain = contrasts(imat)["circuit_beyond_geography_pct"]
    rng = np.random.default_rng(SEED)
    native_draws, introduced_draws = [], []
    for _ in range(BOOT_N):
        sample = rng.integers(0, FOLD_N, FOLD_N)
        native_draws.append(contrasts(nmat[sample])["circuit_beyond_geography_pct"])
        introduced_draws.append(contrasts(imat[sample])["circuit_beyond_geography_pct"])

    def interval(draws):
        return [float(np.quantile(draws, .025)), float(np.quantile(draws, .975))]

    return {
        "one_held": {"native_flora": native_one, "introduced_flora": introduced_one},
        "strict": {
            "source_common_regions": len(islands),
            "source_common_pairs": len(pairs),
            "source_n_test_double_unseen_pairs": native_both["n_double_unseen_test_pairs"],
            "n_source_folds": FOLD_N,
            "original_shortest_double_unseen_verified_before_electrical_fit": {
                name: value["macro_equal_fold_relative_mse_improvement_percent"]
                for name, value in shortest.items()
            },
            "native_uniform_circuit": native_both,
            "alien_uniform_circuit": introduced_both,
            "origin_paired_comparison": {
                "same_folds_same_pairs": True,
                "native_circuit_gain_pct": native_gain,
                "alien_circuit_gain_pct": introduced_gain,
                "alien_minus_native_relative_mse_gain_points": introduced_gain - native_gain,
                "descriptive_shared_5fold_resample_4999_interval": interval(
                    np.asarray(introduced_draws) - np.asarray(native_draws)
                ),
                "bootstrap_seed": SEED,
                "bootstrap_draws": BOOT_N,
                "native_MSE_gain_descriptive_4999_interval": interval(native_draws),
                "alien_MSE_gain_descriptive_4999_interval": interval(introduced_draws),
            },
        },
    }


def replay(root):
    ports = read_table(root, "port_flora_regions.csv", PORT_FIELDS)
    flora = read_table(root, "flora_jaccard_counts.csv", FLORA_FIELDS)
    marine = read_table(root, "marine_fit_inputs.csv", MARINE_FIELDS)
    areas = {row["code"]: float(row["area_km2"]) for row in ports}
    return {
        "A": port_models(ports),
        "B": floristic_overlap(flora, set(areas)),
        "C": circuit_models(marine, areas),
    }


def summary(results, n_checked):
    one = results["C"]["one_held"]
    strict = results["C"]["strict"]
    contrast = strict["origin_paired_comparison"]
    return {
        "status": "FIG2ABC_DERIVED_REPLAY_PASS_INTERNAL_ONLY",
        "release_ready": False,
        "full_SI_replayed": False,
        "checked_original_result_fields": n_checked,
        "A": results["A"]["model_scopes"]["all_126"],
        "B": results["B"],
        "C": {
            "common_islands": strict["source_common_regions"],
            "common_pairs": strict["source_common_pairs"],
            "one_held_folds": one["native_flora"]["held_island_exactly_once_folds"],
            "strict_pairs": strict["source_n_test_double_unseen_pairs"],
            "strict_folds": strict["n_source_folds"],
            "one_held_gain_percent": {
                name: one[f"{name}_flora"]["relative_delta_MSE_percent_positive_is_improvement"]["G_vs_B_pct"]
                for name in ("native", "introduced")
            },
            "one_held_intervals": {
                name: one[f"{name}_flora"]["descriptive_island_bootstrap_95pct_MSE_delta"]["G_vs_B_pct"]
                for name in ("native", "introduced")
            },
            "strict_gain_percent": {"native": contrast["native_circuit_gain_pct"],
                                    "introduced": contrast["alien_circuit_gain_pct"]},
            "strict_intervals": {
                "native": contrast["native_MSE_gain_descriptive_4999_interval"],
                "introduced": contrast["alien_MSE_gain_descriptive_4999_interval"],
            },
            "strict_paired_interval": contrast["descriptive_shared_5fold_resample_4999_interval"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    manifest = json.loads((ROOT / "MANIFEST.json").read_bytes())
    for name, expected_hash in manifest["archive_file_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected_hash:
            raise ValueError(f"SHA256 mismatch: {name}")
    expected = json.loads((ROOT / "EXPECTED_FIT_RESULTS.json").read_bytes())
    actual = replay(ROOT)
    checked = compare(actual, expected)
    print(json.dumps(summary(actual, checked), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
