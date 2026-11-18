#!/usr/bin/env python3
"""Fixed exploratory joint host/circuit prediction on the original51 islands.

Reuses the frozen distributed-port matrices and existing train-only OLS
routine. Each model is fitted once per island/response; every held-out
prediction is retained before any loss aggregation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy

from exploratory_island51_circuit_network_holdout import source, train_predict
from island51_ocean_circuit import resistance_feature
from run_island51_observed_ocean_circuits import make_base_features, sha256, validate_source_bundle


REPO = Path(__file__).resolve().parents[1]
RESULT_DIRECTORY = "results/island51_host_circuit_joint_20261106"
PRIOR_DIRECTORY = "results/island51_observed_ocean_circuit_20261106"
NATIVE_DIRECTORY = PRIOR_DIRECTORY + "/native_source_sensitivity"
SOURCE_DIRECTORY = "provenance/reviewer_defenses/results/island_host_network_20261030"
DEFAULT_LOCK = REPO / RESULT_DIRECTORY / "input_lock.json"
TARGETS = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
MODELS = ("G", "H", "C", "HC", "I")
COMPARISONS = (
    ("G_to_H", "G", "H", "contextual"),
    ("G_to_C", "G", "C", "contextual"),
    ("C_to_HC", "C", "HC", "primary"),
    ("H_to_HC", "H", "HC", "primary"),
    ("HC_to_I", "HC", "I", "secondary_interaction"),
)
INPUT_PATHS = {
    "pairs": SOURCE_DIRECTORY + "/island_host_network_pairs.csv",
    "islands": SOURCE_DIRECTORY + "/island_host_network_islands.csv",
    "uniform": PRIOR_DIRECTORY + "/uniform_distributed_port_resistance.npz",
    "wind": PRIOR_DIRECTORY + "/wind_distributed_port_resistance.npz",
    "ais": NATIVE_DIRECTORY + "/native_AIS_distributed_port_resistance.npz",
    "preparation": PRIOR_DIRECTORY + "/source_preparation.json",
    "prior_summary": NATIVE_DIRECTORY + "/model_comparison_summary.csv",
    "prior_folds": NATIVE_DIRECTORY + "/all_island_holdout_losses.csv",
    "official_primary": PRIOR_DIRECTORY + "/official_circuitscape_validation.json",
    "official_native": NATIVE_DIRECTORY + "/native_AIS_official_validation.json",
    "source_bundle": "data/exploratory/island51_observed_ocean_20261106/source_bundle_manifest.json",
    "prior_protocol": PRIOR_DIRECTORY + "/PROTOCOL.md",
}


def validate_input_lock(lock_path, *, repo=REPO):
    root = Path(repo).resolve()
    lock = json.loads(Path(lock_path).read_text())
    if (lock.get("schema") != "chocho_island51_joint_host_circuit_input_lock_v1" or
            lock.get("status") != "FIXED_BEFORE_NEW_JOINT_MODEL_FITS" or
            lock.get("prior_circuit_and_host_results_already_seen") is not True or
            lock.get("new_joint_model_results_already_seen") is not False or
            (lock.get("fixed_islands"), lock.get("fixed_pairs")) != (51, 1275) or
            lock.get("primary_comparisons") != ["C_to_HC", "H_to_HC"]):
        raise ValueError("Joint source lock does not preserve the fixed exploratory design")
    files = lock.get("files", [])
    if len(files) != 12 or {entry["path"] for entry in files} != set(INPUT_PATHS.values()):
        raise ValueError("All12 declared source and prior-result files must be locked")
    for entry in files:
        path = (root / entry["path"]).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Locked input is outside the repository")
        if (not path.is_file() or path.stat().st_size != entry["bytes"] or
                sha256(path) != entry["sha256"]):
            raise ValueError("Locked input size or checksum changed: " + entry["path"])
    protocol = (root / lock["protocol_path"]).resolve()
    if not protocol.is_relative_to(root) or sha256(protocol) != lock.get("protocol_sha256"):
        raise ValueError("Frozen joint protocol checksum changed")
    implementation = protocol.with_name("IMPLEMENTATION_PLAN.md")
    if sha256(implementation) != lock.get("implementation_plan_sha256"):
        raise ValueError("Frozen joint implementation plan checksum changed")
    return {"status": "ALL_LOCKED_JOINT_INPUTS_AND_PROTOCOL_VERIFIED", "locked_files_verified": len(files),
            "input_lock_sha256": sha256(lock_path), "protocol_sha256": lock["protocol_sha256"],
            "implementation_plan_sha256": lock["implementation_plan_sha256"],
            "frozen_at_utc": lock["frozen_at_utc"], "base_commit": lock["base_commit"],
            "prior_results_already_seen_when_frozen": True,
            "new_joint_results_already_seen_when_frozen": False, "files": files}


def load_resistance_matrix(path, ids):
    with np.load(path, allow_pickle=False) as data:
        matrix_ids = np.asarray(data["island_ids"])
        matrix = np.asarray(data["resistance"], dtype=float)
    if matrix_ids.ndim != 1 or matrix_ids.tolist() != list(ids):
        raise ValueError("Circuit matrix island ID order differs from canonical source order")
    n = len(ids)
    if (matrix.shape != (n, n) or not np.isfinite(matrix).all() or
            not np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-12) or
            not np.allclose(np.diag(matrix), 0., rtol=0, atol=1e-12) or
            np.any(matrix[np.triu_indices(n, 1)] <= 0)):
        raise ValueError("Circuit matrix must be finite, symmetric, zero-diagonal and positive between distinct islands")
    return matrix


def load_original_port_counts(path, ids, preparation):
    roster = pd.read_csv(path, dtype={"landmass_id": str})
    eligible = roster[roster.eligible.astype(str).str.lower().eq("true")]
    if (len(eligible) != 51 or eligible.landmass_id.nunique() != 51 or
            set(eligible.landmass_id) != set(ids) or
            preparation.get("original_ports_audited") != 194 or
            preparation.get("marine_ports_included") != 182 or preparation.get("islands_retained") != 51):
        raise ValueError("The original port-count roster or source scope changed")
    counts = {}
    for row in eligible.itertuples(index=False):
        value = float(row.ports)
        if not np.isfinite(value) or value <= 0 or value != int(value):
            raise ValueError("Invalid original port count")
        counts[row.landmass_id] = int(value)
    recorded = preparation["by_island"]
    if (set(recorded) != set(ids) or sum(counts.values()) != 194 or
            any(counts[key] != recorded[key]["original_ports"] for key in ids) or
            sum(recorded[key]["marine_ports"] for key in ids) != 182):
        raise ValueError("Archived island and source-preparation original port counts disagree")
    return counts


def make_joint_features(q, ids, port_counts, matrices, target):
    if target not in TARGETS:
        raise ValueError("Only the two fixed butterfly targets may be fitted")
    if set(matrices) != {"uniform", "ais", "wind"}:
        raise ValueError("The three fixed circuit matrices are required")
    circuit = {key: resistance_feature(matrix, q.i.to_numpy(), q.j.to_numpy())
               for key, matrix in matrices.items()}
    host = q.host_jaccard.to_numpy(dtype=float)
    if not np.isfinite(host).all() or np.any((host < 0) | (host > 1)):
        raise ValueError("Original host Jaccard is outside its recorded range")
    base = make_base_features(q, ids, port_counts, target)
    geographic = dict(base, log1p_normalized_uniform_circuit=circuit["uniform"])
    host_model = dict(geographic, host_jaccard=host)
    circuit_model = dict(geographic, log1p_normalized_ais_circuit=circuit["ais"],
                         log1p_normalized_wind_circuit=circuit["wind"])
    joint = dict(host_model, log1p_normalized_ais_circuit=circuit["ais"],
                 log1p_normalized_wind_circuit=circuit["wind"])
    interaction = dict(joint, host_times_AIS_connectivity=host*(-circuit["ais"]),
                       host_times_wind_connectivity=host*(-circuit["wind"]))
    result = {"G": geographic, "H": host_model, "C": circuit_model, "HC": joint, "I": interaction}
    if any(np.asarray(value).shape != (len(q),) or not np.isfinite(value).all()
           for model in result.values() for value in model.values()):
        raise ValueError("A declared joint feature is not aligned or finite")
    return result


def fit_joint_models(q, ids, port_counts, matrices):
    if len(ids) != 51 or len(set(ids)) != 51 or len(q) != 1275:
        raise ValueError("Joint prediction requires the exact51-island/1275-pair support")
    left, right = q.i.to_numpy(dtype=int), q.j.to_numpy(dtype=int)
    if (np.any(left < 0) or np.any(right < 0) or np.any(left >= 51) or np.any(right >= 51) or
            np.any(left == right) or len({tuple(sorted(pair)) for pair in zip(left, right)}) != 1275):
        raise ValueError("Joint prediction pair endpoints are incomplete or duplicated")
    pair_ids = np.array([ids[min(i, j)]+":"+ids[max(i, j)] for i, j in zip(left, right)])
    source_rows = np.arange(len(q))
    observed_a, observed_b = q.island_a.to_numpy(), q.island_b.to_numpy()
    predictions, folds, feature_sets = [], [], {}
    for target in TARGETS:
        features = make_joint_features(q, ids, port_counts, matrices, target)
        feature_sets[target] = {model: list(values) for model, values in features.items()}
        designs = {model: np.column_stack(list(values.values())) for model, values in features.items()}
        response = q[target].to_numpy(dtype=float)
        if not np.isfinite(response).all() or np.any((response < 0) | (response > 1)):
            raise ValueError("A fixed butterfly response is invalid")
        for index, island in enumerate(ids):
            test = (left == index) | (right == index)
            train = ~test
            if (int(train.sum()), int(test.sum())) != (1225, 50):
                raise ValueError("Every joint fold must retain1225 training and50 test pairs")
            for model in MODELS:
                # Exactly one call per model/fold/target; comparisons share it.
                predicted = train_predict(designs[model], response, train, test)
                if predicted.shape != (50,) or not np.isfinite(predicted).all():
                    raise ValueError("Joint model produced nonfinite or misaligned predictions")
                observed = response[test]
                squared = (observed-predicted)**2
                absolute = np.abs(observed-predicted)
                predictions.append(pd.DataFrame({
                    "target": target, "model": model, "heldout_island_id": island,
                    "source_pair_index": source_rows[test], "pair_id": pair_ids[test],
                    "island_a": observed_a[test], "island_b": observed_b[test],
                    "observed": observed, "predicted": predicted,
                    "squared_error": squared, "absolute_error": absolute,
                }))
                folds.append({"target": target, "model": model, "heldout_island_id": island,
                              "train_pairs": 1225, "test_pairs": 50,
                              "MSE": float(np.mean(squared)), "MAE": float(np.mean(absolute))})
    return pd.concat(predictions, ignore_index=True), pd.DataFrame(folds), feature_sets


def summarize_results(model_folds):
    keys = ["target", "model", "heldout_island_id"]
    if (len(model_folds) != 510 or model_folds.duplicated(keys).any() or
            set(model_folds.target) != set(TARGETS) or set(model_folds.model) != set(MODELS) or
            not model_folds.groupby(["target", "model"]).size().eq(51).all()):
        raise ValueError("All510 unique fixed model-fold losses are required")
    model_rows, summary_rows, comparison_rows = [], [], []
    for target in TARGETS:
        selected = model_folds[model_folds.target.eq(target)]
        for model in MODELS:
            values = selected[selected.model.eq(model)]
            model_rows.append({"target": target, "model": model, "MSE": float(values.MSE.mean()),
                               "MAE": float(values.MAE.mean()), "n_folds": 51, "n_predictions": 2550})
        for comparison, reference, augmented, role in COMPARISONS:
            a = selected[selected.model.eq(reference)]
            b = selected[selected.model.eq(augmented)]
            joined = a.merge(b, on="heldout_island_id", suffixes=("_ref", "_aug"), validate="one_to_one")
            if len(joined) != 51:
                raise ValueError("Compared models have different island fold identities")
            current = []
            for row in joined.itertuples(index=False):
                if (row.train_pairs_ref, row.train_pairs_aug, row.test_pairs_ref, row.test_pairs_aug) != (1225, 1225, 50, 50):
                    raise ValueError("Compared model folds use different source support")
                current.append({"target": target, "comparison": comparison, "reference_model": reference,
                    "augmented_model": augmented, "heldout_island_id": row.heldout_island_id,
                    "train_pairs": 1225, "test_pairs": 50, "reference_mse": row.MSE_ref,
                    "augmented_mse": row.MSE_aug, "delta_mse_augmented_minus_reference": row.MSE_aug-row.MSE_ref,
                    "reference_mae": row.MAE_ref, "augmented_mae": row.MAE_aug,
                    "improved": row.MSE_aug < row.MSE_ref, "worse": row.MSE_aug > row.MSE_ref,
                    "tied": row.MSE_aug == row.MSE_ref})
            mse_a = float(np.mean([row["reference_mse"] for row in current]))
            mse_b = float(np.mean([row["augmented_mse"] for row in current]))
            summary_rows.append({"comparison": comparison, "target": target, "comparison_role": role,
                "reference_model": reference, "augmented_model": augmented,
                "MSE_reference": mse_a, "MSE_augmented": mse_b,
                "relative_MSE_reduction_percent": 100*(mse_a-mse_b)/mse_a if mse_a > 0 else None,
                "MAE_reference": float(np.mean([row["reference_mae"] for row in current])),
                "MAE_augmented": float(np.mean([row["augmented_mae"] for row in current])),
                "absolute_MSE_gain": mse_a-mse_b,
                "islands_improved": sum(row["improved"] for row in current),
                "islands_worse": sum(row["worse"] for row in current),
                "islands_tied": sum(row["tied"] for row in current)})
            comparison_rows.extend(current)
    return pd.DataFrame(model_rows), pd.DataFrame(summary_rows), pd.DataFrame(comparison_rows)


def verify_prior_circuit_results(model_folds, summary, prior_folds, prior_summary):
    previous = prior_folds[(prior_folds.comparison == "uniform_vs_ais_and_wind") & prior_folds.target.isin(TARGETS)]
    if (len(previous) != 102 or previous.duplicated(["target", "heldout_island_id"]).any() or
            not previous.groupby("target").size().eq(51).all()):
        raise ValueError("The prior native circuit reference must have102 unique butterfly fold rows")
    current = model_folds.set_index(["target", "model", "heldout_island_id"], verify_integrity=True)
    differences = []

    def compare(observed, expected, context):
        difference = abs(float(observed)-float(expected))
        if not np.isfinite(difference) or not np.isclose(observed, expected, atol=1e-12, rtol=1e-12):
            raise ValueError("New G/C result differs from prior native circuit " + context)
        differences.append(difference)

    for old in previous.itertuples(index=False):
        try:
            geographic = current.loc[(old.target, "G", str(old.heldout_island_id))]
            circuit = current.loc[(old.target, "C", str(old.heldout_island_id))]
        except KeyError as exc:
            raise ValueError("A prior circuit fold identity is absent") from exc
        for new, mse, mae in [(geographic, old.baseline_mse, old.MAE_base),
                               (circuit, old.circuit_augmented_mse, old.MAE_circuit)]:
            compare(new.MSE, mse, "fold MSE")
            compare(new.MAE, mae, "fold MAE")
            if (new.train_pairs, new.test_pairs, old.train_pairs, old.test_pairs) != (1225, 50, 1225, 50):
                raise ValueError("A prior circuit fold has different train/test support")
        compare(circuit.MSE-geographic.MSE, old.delta_mse_circuit_minus_baseline, "fold difference")
        if bool(circuit.MSE < geographic.MSE) != bool(old.circuit_improves_this_island):
            raise ValueError("A prior circuit fold improvement flag changed")
    earlier_summary = prior_summary[(prior_summary.comparison == "uniform_vs_ais_and_wind") & prior_summary.target.isin(TARGETS)]
    if len(earlier_summary) != 2 or earlier_summary.target.nunique() != 2:
        raise ValueError("The prior native circuit must have both butterfly summaries")
    for target in TARGETS:
        new = summary[(summary.comparison == "G_to_C") & summary.target.eq(target)]
        old = earlier_summary[earlier_summary.target.eq(target)]
        if len(new) != 1 or len(old) != 1:
            raise ValueError("A prior circuit summary is missing or duplicated")
        for key in ["MSE_reference", "MSE_augmented", "MAE_reference", "MAE_augmented", "relative_MSE_reduction_percent", "islands_improved"]:
            compare(new.iloc[0][key], old.iloc[0][key], "summary " + key)
    return {"status": "PRIOR_NATIVE_G_AND_C_FOLDS_AND_SUMMARIES_REPRODUCED",
            "prior_comparison_fold_rows_checked": 102, "prior_model_fold_records_checked": 204,
            "prior_summary_rows_checked": 2, "maximum_absolute_numeric_difference": max(differences),
            "relative_tolerance": 1e-12, "absolute_tolerance": 1e-12}


def write_predictions(path, predictions):
    with Path(path).open("wb") as binary:
        with gzip.GzipFile(filename="", fileobj=binary, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                predictions.to_csv(stream, index=False, float_format="%.17g")


def run(out_dir, input_lock=DEFAULT_LOCK):
    started = time.monotonic()
    out = Path(out_dir).resolve()
    for protected in [REPO / "data", REPO / "provenance", REPO / PRIOR_DIRECTORY]:
        if out.is_relative_to(protected.resolve()):
            raise ValueError("Joint outputs must be separate from sources and prior results")
    validation = validate_input_lock(input_lock)
    paths = {key: REPO / value for key, value in INPUT_PATHS.items()}
    bundle = validate_source_bundle(paths["source_bundle"].parent, paths["prior_protocol"])
    if not bundle.get("manifest_present"):
        raise ValueError("The prior frozen source bundle is required")
    q, ids, _ = source(paths["pairs"], paths["islands"])
    preparation = json.loads(paths["preparation"].read_text())
    counts = load_original_port_counts(paths["islands"], ids, preparation)
    matrices = {key: load_resistance_matrix(paths[key], ids) for key in ("uniform", "ais", "wind")}
    prior_folds = pd.read_csv(paths["prior_folds"], dtype={"heldout_island_id": str})
    prior_summary = pd.read_csv(paths["prior_summary"])
    validation.update({"completed_before_joint_fit_at_utc": datetime.now(timezone.utc).isoformat(),
        "prior_frozen_source_bundle": bundle, "canonical_island_ids": ids,
        "original_islands": len(ids), "original_pairs": len(q), "original_ports": sum(counts.values()),
        "marine_circuit_ports": preparation["marine_ports_included"],
        "original_host_inventory_and_zero_similarities_retained": True,
        "host_zero_similarity_pairs": int(q.host_jaccard.eq(0).sum()),
        "original_port_counts": counts, "matrix_ids_aligned": True,
        "matrices_finite_symmetric_zero_diagonal_positive_off_diagonal": True,
        "native_AIS_matrix_used": INPUT_PATHS["ais"], "circuits_rebuilt": False})
    out.mkdir(parents=True, exist_ok=True)
    (out / "source_validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    print("JOINT_SOURCE_GATES_PASS", json.dumps({"locked_files": 12, "islands": 51, "pairs": 1275,
          "original_ports": 194, "marine_circuit_ports": 182}), flush=True)
    predictions, model_folds, feature_sets = fit_joint_models(q, ids, counts, matrices)
    model_summary, summary, comparison_folds = summarize_results(model_folds)
    reproduction = verify_prior_circuit_results(model_folds, summary, prior_folds, prior_summary)
    if (len(predictions), len(model_folds), len(comparison_folds), len(summary)) != (25500, 510, 510, 10):
        raise ValueError("Fixed joint analysis output support changed")
    repeated = predictions.groupby(["target", "model", "source_pair_index"]).size()
    if len(repeated) != 12750 or not repeated.eq(2).all():
        raise ValueError("Every fixed model/target pair must have two endpoint-held-out predictions")
    if validate_input_lock(input_lock)["input_lock_sha256"] != validation["input_lock_sha256"]:
        raise ValueError("Frozen joint lock changed during fitting")
    validate_source_bundle(paths["source_bundle"].parent, paths["prior_protocol"])
    write_predictions(out / "heldout_predictions.csv.gz", predictions)
    outputs = {"model_fold_losses.csv": model_folds, "comparison_fold_losses.csv": comparison_folds,
               "model_summary.csv": model_summary, "comparison_summary.csv": summary}
    for filename, frame in outputs.items():
        frame.to_csv(out / filename, index=False, float_format="%.17g")
    (out / "prior_circuit_reproduction.json").write_text(json.dumps(reproduction, indent=2)+"\n")
    output_names = ["heldout_predictions.csv.gz", *outputs, "source_validation.json", "prior_circuit_reproduction.json"]
    receipt = {"schema": "chocho_fixed_island51_joint_host_circuit_run_v1",
        "status": "FIVE_FIXED_JOINT_MODELS_AND_ALL_HELDOUT_PREDICTIONS_COMPLETE",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(), "elapsed_seconds": time.monotonic()-started,
        "prior_results_already_seen_when_design_frozen": True,
        "new_joint_results_known_when_design_frozen": False,
        "independent_validation_population": False, "models_or_transforms_selected_by_new_scores": False,
        "source_validation": validation, "prior_circuit_reproduction": reproduction,
        "targets": list(TARGETS), "models": list(MODELS), "feature_sets": feature_sets,
        "model_fits": 510, "individual_prediction_rows": len(predictions),
        "model_fold_rows": len(model_folds), "comparison_fold_rows": len(comparison_folds),
        "comparison_summary_rows": len(summary),
        "loss_aggregation": "Mean of50 individual losses within each island fold, then mean over51 folds; paired predictions are not averaged before loss.",
        "improved_worse_tied_definition": "Strict comparisons of fold MSE: augmented<reference, augmented>reference, augmented==reference.",
        "estimation": "Existing train_predict: training-only standardization, intercept, ordinary least squares, clipping to[0,1].",
        "interaction_definition": "Raw host_jaccard times negative normalized log1p AIS resistance; same for wind; both enter together with all main effects.",
        "comparisons": summary.to_dict(orient="records"), "model_summary": model_summary.to_dict(orient="records"),
        "python_environment": {"python": platform.python_version(), "numpy": np.__version__,
                               "pandas": pd.__version__, "scipy": scipy.__version__},
        "implementation_sha256": {str(Path(__file__).relative_to(REPO)): sha256(__file__),
            "scripts/exploratory_island51_circuit_network_holdout.py": sha256(REPO / "scripts/exploratory_island51_circuit_network_holdout.py"),
            "scripts/run_island51_observed_ocean_circuits.py": sha256(REPO / "scripts/run_island51_observed_ocean_circuits.py"),
            "scripts/island51_ocean_circuit.py": sha256(REPO / "scripts/island51_ocean_circuit.py")},
        "outputs_sha256": {name: sha256(out / name) for name in output_names},
        "interpretation_limits": [
            "Exploratory extension informed by known host/circuit results, fixed before these new joint fits.",
            "Host inventory pools2001-2024 and can postdate2011-2022 butterfly records; no prospective or temporal mediation claim.",
            "Candidate host composition and undirected marine fields do not establish plant carriage, butterfly movement or breeding.",
            "Pair and fold dependence precludes independent-pair or independent-fold p-values and confidence intervals.",
            "One held-out island is paired with otherwise represented islands; two unseen endpoints are not evaluated.",
            "Conditional increments do not erase prior marginal host associations or identify a causal mechanism.",
        ]}
    (out / "joint_run_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False)+"\n")
    print("JOINT_PRIOR_G_C_REPRODUCED", json.dumps(reproduction), flush=True)
    print("JOINT_FIXED_COMPARISONS_COMPLETE", json.dumps(summary.to_dict(orient="records")), flush=True)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--input-lock", type=Path, default=DEFAULT_LOCK)
    args = parser.parse_args()
    run(args.out_dir, args.input_lock)


if __name__ == "__main__":
    main()
