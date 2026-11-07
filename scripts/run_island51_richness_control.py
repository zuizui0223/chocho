#!/usr/bin/env python3
"""Post-hoc host composition versus marginal host-richness sensitivity for 51 islands."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from exploratory_island51_circuit_network_holdout import source, train_predict
from run_island51_host_circuit_joint import (
    INPUT_PATHS, TARGETS, load_original_port_counts, load_resistance_matrix,
    make_joint_features, validate_input_lock,
)
from run_island51_observed_ocean_circuits import sha256

PARENT = ROOT / "results/island51_host_circuit_joint_20261106"
ATLAS_LOCK = ROOT / "results/island51_host_circuit_atlas_20261106/input_lock.json"
ROSTER = ROOT / INPUT_PATHS["islands"]
PAIRS = ROOT / INPUT_PATHS["pairs"]
PARENT_FOLDS = PARENT / "model_fold_losses.csv"
MODELS = ("C", "H", "HC", "CR", "HR", "HCR")
COMPARISONS = (
    ("C_to_HC", "C", "HC", "parent_host_increment"),
    ("H_to_HC", "H", "HC", "parent_circuit_increment"),
    ("C_to_CR", "C", "CR", "candidate_host_richness_alone"),
    ("CR_to_HCR", "CR", "HCR", "primary_host_identity_beyond_richness"),
    ("H_to_HR", "H", "HR", "richness_beyond_host_identity"),
    ("HR_to_HCR", "HR", "HCR", "secondary_circuit_beyond_hosts_and_richness"),
    ("HC_to_HCR", "HC", "HCR", "richness_beyond_joint_model"),
)


def host_counts(path, ids):
    frame = pd.read_csv(path, dtype={"landmass_id": str})
    eligible = frame.loc[frame.eligible.astype(str).str.lower().eq("true")]
    if len(eligible) != 51 or eligible.landmass_id.nunique() != 51 or set(eligible.landmass_id) != set(ids):
        raise ValueError("Host counts do not cover the fixed 51 islands")
    numbers = {}
    for row in eligible.itertuples(index=False):
        x = float(row.alien_host_species)
        if not np.isfinite(x) or x < 0 or x != int(x):
            raise ValueError("Nonnegative integer host counts required")
        numbers[str(row.landmass_id)] = int(x)
    return numbers


def richness_features(q, ids, counts):
    if set(counts) != set(ids):
        raise ValueError("Host richness island identity mismatch")
    log = np.log1p(np.array([counts[k] for k in ids], dtype=float))
    l, r = log[q.i.to_numpy(dtype=int)], log[q.j.to_numpy(dtype=int)]
    values = {
        "sum_log1p_candidate_host_taxa": l + r,
        "abs_difference_log1p_candidate_host_taxa": np.abs(l - r),
    }
    if any(v.shape != (len(q),) or not np.isfinite(v).all() for v in values.values()):
        raise ValueError("Unaligned richness predictors")
    return values


def compare_folds(folds):
    if len(folds) != len(TARGETS) * len(MODELS) * 51 or folds.duplicated(["target", "model", "heldout_island_id"]).any():
        raise ValueError("Incomplete 51-island model-fold support")
    rows = []
    for target in TARGETS:
        f = folds.loc[folds.target.eq(target)]
        for name, ref, aug, role in COMPARISONS:
            a = f.loc[f.model.eq(ref)].set_index("heldout_island_id")
            b = f.loc[f.model.eq(aug)].set_index("heldout_island_id")
            if len(a) != 51 or len(b) != 51 or set(a.index) != set(b.index):
                raise ValueError("Unmatched fold identities")
            b = b.reindex(a.index)
            g = a.MSE.to_numpy() - b.MSE.to_numpy()
            g_mae = a.MAE.to_numpy() - b.MAE.to_numpy()
            x, y = float(a.MSE.mean()), float(b.MSE.mean())
            xa, ya = float(a.MAE.mean()), float(b.MAE.mean())
            rows.append({
                "target": target, "comparison": name, "role": role,
                "reference_model": ref, "augmented_model": aug,
                "islands": 51, "heldout_endpoint_predictions": 2550,
                "MSE_reference": x, "MSE_augmented": y,
                "MSE_reduction_percent": 100 * (x-y) / x,
                "islands_mse_improved": int(np.sum(g > 0)),
                "islands_mse_worse": int(np.sum(g < 0)),
                "islands_mse_tied": int(np.sum(g == 0)),
                "MAE_reference": xa, "MAE_augmented": ya,
                "MAE_reduction_percent": 100 * (xa-ya) / xa,
                "islands_mae_improved": int(np.sum(g_mae > 0)),
                "islands_mae_worse": int(np.sum(g_mae < 0)),
            })
    return pd.DataFrame(rows)


def run(out_dir):
    out = Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT / p).resolve()) for p in
           ("data", "provenance", "results/island51_host_circuit_joint_20261106",
            "results/island51_host_circuit_atlas_20261106")):
        raise ValueError("Cannot overwrite source, parent, or atlas artifacts")
    validated = validate_input_lock(PARENT / "input_lock.json")
    atlas = json.loads(ATLAS_LOCK.read_text())
    pin = {item["path"]: item["sha256"] for item in atlas["files"]}
    relative_fold = str(PARENT_FOLDS.relative_to(ROOT))
    if pin.get(relative_fold) != sha256(PARENT_FOLDS):
        raise ValueError("Parent fold checksum changed from atlas-pinned version")
    q, ids, _ = source(PAIRS, ROSTER)
    if len(ids) != 51 or len(q) != 1275:
        raise ValueError("Original island/pair support changed")
    prep = json.loads((ROOT / INPUT_PATHS["preparation"]).read_text())
    ports = load_original_port_counts(ROSTER, ids, prep)
    counts = host_counts(ROSTER, ids)
    richness = richness_features(q, ids, counts)
    matrices = {key: load_resistance_matrix(ROOT / INPUT_PATHS[key], ids) for key in ("uniform", "ais", "wind")}
    parent = pd.read_csv(PARENT_FOLDS, dtype={"heldout_island_id": str}).set_index(
        ["target", "model", "heldout_island_id"], verify_integrity=True)
    folds, predictions, features_by_target = [], [], {}
    max_mse = max_mae = 0.0
    ii, jj = q.i.to_numpy(dtype=int), q.j.to_numpy(dtype=int)
    for target in TARGETS:
        old = make_joint_features(q, ids, ports, matrices, target)
        features = {
            "C": old["C"], "H": old["H"], "HC": old["HC"],
            "CR": dict(old["C"], **richness),
            "HR": dict(old["H"], **richness),
            "HCR": dict(old["HC"], **richness),
        }
        features_by_target[target] = {m: list(features[m]) for m in MODELS}
        design = {m: np.column_stack(list(features[m].values())) for m in MODELS}
        y = q[target].to_numpy(dtype=float)
        for i, island in enumerate(ids):
            test = (ii == i) | (jj == i)
            train = ~test
            if (int(test.sum()), int(train.sum())) != (50, 1225):
                raise ValueError("Invalid leave-one-island-out fold")
            for m in MODELS:
                estimate = train_predict(design[m], y, train, test)
                obs = y[test]
                if len(estimate) != 50 or not np.isfinite(estimate).all():
                    raise ValueError("Invalid prediction")
                mse = float(np.mean((obs-estimate)**2))
                mae = float(np.mean(np.abs(obs-estimate)))
                if m in ("C", "H", "HC"):
                    p = parent.loc[(target, m, island)]
                    dx, dy = abs(mse-float(p.MSE)), abs(mae-float(p.MAE))
                    if dx > 1e-12 or dy > 1e-12:
                        raise ValueError("Failed parent fold reproduction: " + target + "/" + m + "/" + island)
                    max_mse, max_mae = max(max_mse, dx), max(max_mae, dy)
                folds.append({"target": target, "model": m, "heldout_island_id": island,
                              "train_pairs": 1225, "test_pairs": 50, "MSE": mse, "MAE": mae})
                for pidx, value, actual in zip(np.flatnonzero(test), estimate, obs):
                    predictions.append({"target": target, "model": m, "heldout_island_id": island,
                                        "source_pair_index": int(pidx), "observed": float(actual),
                                        "predicted": float(value),
                                        "squared_error": float((actual-value)**2),
                                        "absolute_error": float(abs(actual-value))})
    if len(predictions) != 30600 or len(folds) != 612:
        raise ValueError("Missing per-endpoint predictions or model folds")
    fold_frame = pd.DataFrame(folds)
    summary = compare_folds(fold_frame)
    if len(summary) != 14:
        raise ValueError("Missing fixed model comparisons")
    out.mkdir(parents=True, exist_ok=True)
    fold_frame.to_csv(out/"model_fold_losses.csv", index=False, float_format="%.17g")
    summary.to_csv(out/"comparison_summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(predictions).to_csv(out/"heldout_predictions.csv", index=False, float_format="%.17g")
    receipt = {
        "schema": "chocho_island51_host_richness_sensitivity_v1",
        "status": "ALL_51_ISLANDS_AND_SIX_FIXED_MODELS_EVALUATED",
        "status_role": "exploratory_posthoc_sensitivity_not_independent_validation",
        "source_lock_status": validated["status"], "parent_folds_sha256": sha256(PARENT_FOLDS),
        "max_parent_mse_difference": max_mse, "max_parent_mae_difference": max_mae,
        "island_count": 51, "unique_island_pairs": 1275, "original_port_count": sum(ports.values()),
        "candidate_host_counts": {k: counts[k] for k in ids},
        "new_marginal_terms": list(richness), "feature_names": features_by_target,
        "models": list(MODELS), "fixed_comparisons": [list(x) for x in COMPARISONS],
        "fold_rows": 612, "endpoint_prediction_rows": 30600,
        "estimator": "train-only standardized OLS with intercept, predictions clipped to [0,1]",
        "folds_independent": False,
        "plant_inventory_pooled_2001_2024": True,
        "causal_colonization_or_transport_identified": False,
        "output_sha256": {p.name: sha256(p) for p in sorted(out.glob("*.csv"))},
    }
    (out/"result_receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2)+"\n")
    print("RICHNESS_RESULT_JSON="+summary.to_json(orient="records", double_precision=12), flush=True)
    print("RICHNESS_CHECK_JSON="+json.dumps({
        "parent_max_mse_error": max_mse, "parent_max_mae_error": max_mae,
        "model_fold_rows": len(folds), "predictions": len(predictions),
    }, sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, required=True)
    run(ap.parse_args().out_dir)
