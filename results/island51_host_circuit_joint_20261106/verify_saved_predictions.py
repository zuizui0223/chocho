#!/usr/bin/env python3
"""Independently audit joint-model prediction files without fitting any model.

Uses original source observations, csv parsing and math.fsum loss aggregation.
Does not import the joint runner or any of its fitting/aggregation helpers.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[2]
SPEC_DIRECTORY = Path(__file__).resolve().parent
ORIGINAL = Path("provenance/reviewer_defenses/results/island_host_network_20261030")
CIRCUIT = Path("results/island51_observed_ocean_circuit_20261106")
NATIVE = CIRCUIT / "native_source_sensitivity"
TARGETS = ("new_butterfly_jaccard", "earlier_butterfly_jaccard")
MODELS = ("G", "H", "C", "HC", "I")
COMPARISONS = {
    "G_to_H": ("G", "H", "contextual"),
    "G_to_C": ("G", "C", "contextual"),
    "C_to_HC": ("C", "HC", "primary"),
    "H_to_HC": ("H", "HC", "primary"),
    "HC_to_I": ("HC", "I", "secondary_interaction"),
}
TABLE_KEYS = {
    "model_fold_losses": ("target", "model", "heldout_island_id"),
    "comparison_fold_losses": ("target", "comparison", "heldout_island_id"),
    "model_summary": ("target", "model"),
    "comparison_summary": ("target", "comparison"),
}
ATOL = RTOL = 1e-12


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing or duplicate CSV headers: {path}")
        rows = list(reader)
    require(all(None not in row and None not in row.values() for row in rows),
            f"Malformed CSV row: {path}")
    return rows


def finite(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric {label}: {value!r}") from exc
    require(math.isfinite(number), f"Nonfinite {label}")
    return number


def boolean(value, label):
    if isinstance(value, bool):
        return value
    text = str(value).lower()
    require(text in {"true", "false"}, f"Invalid boolean {label}: {value!r}")
    return text == "true"


def close(value, expected, label, maxima=None, *, atol=ATOL, rtol=RTOL):
    actual = finite(value, label)
    error = abs(actual - expected)
    if maxima is not None:
        maxima[label] = max(maxima.get(label, 0.0), error)
    require(math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol),
            f"Mismatch {label}: saved={actual!r}, reconstructed={expected!r}")


def reconstruct(predictions, source_pairs, island_ids):
    """Score each endpoint-held-out prediction first, then aggregate losses."""
    n = len(island_ids)
    require(len(set(island_ids)) == n, "Duplicate source island ID")
    require(len(source_pairs) == n * (n - 1) // 2, "Source pair coverage changed")
    required_rows = len(TARGETS) * len(MODELS) * n * (n - 1)
    require(len(predictions) == required_rows, "Prediction row coverage changed")
    source_seen = set()
    source = []
    for index, row in enumerate(source_pairs):
        a, b = str(row["island_a"]), str(row["island_b"])
        require(a in island_ids and b in island_ids and a != b,
                "Invalid original source endpoints")
        pair_id = ":".join(sorted((a, b), key=int))
        require(pair_id not in source_seen, "Duplicate original source pair")
        source_seen.add(pair_id)
        observations = {target: finite(row[target], f"source.{target}") for target in TARGETS}
        require(all(0 <= value <= 1 for value in observations.values()),
                "Original source observation outside Jaccard range")
        source.append((a, b, pair_id, observations))

    groups = defaultdict(list)
    seen = set()
    endpoint_coverage = defaultdict(set)
    maxima = {}
    for row in predictions:
        target, model = str(row["target"]), str(row["model"])
        heldout = str(row["heldout_island_id"])
        require(target in TARGETS and model in MODELS, "Unexpected prediction target/model")
        value = finite(row["source_pair_index"], "source_pair_index")
        require(value.is_integer() and 0 <= value < len(source), "Invalid source_pair_index")
        index = int(value)
        a, b, pair_id, observations = source[index]
        require((str(row["island_a"]), str(row["island_b"])) == (a, b),
                "Prediction endpoints disagree with original source pair index/orientation")
        require(str(row["pair_id"]) == pair_id, "Prediction pair_id is not the canonical source pair")
        require(heldout in (a, b), "Held-out island is not a source pair endpoint")
        key = (target, model, heldout, index)
        require(key not in seen, "duplicate prediction replaces an endpoint observation")
        seen.add(key)
        endpoint_coverage[target, model, index].add(heldout)
        observed = finite(row["observed"], "observed")
        close(observed, observations[target], "observed_vs_source", maxima,
              atol=2e-15, rtol=2e-15)
        predicted = finite(row["predicted"], "predicted")
        require(0 <= predicted <= 1, "Prediction outside the clipped Jaccard range")
        error = observed - predicted
        squared, absolute = error * error, abs(error)
        close(row["squared_error"], squared, "individual.squared_error", maxima)
        close(row["absolute_error"], absolute, "individual.absolute_error", maxima)
        groups[target, model, heldout].append((squared, absolute))

    for target in TARGETS:
        for model in MODELS:
            for index, (a, b, _, _) in enumerate(source):
                require(endpoint_coverage[target, model, index] == {a, b},
                        "Both held-out endpoints must represent every source pair")
    result = {name: [] for name in TABLE_KEYS}
    lookup = {}
    for target in TARGETS:
        for model in MODELS:
            losses = []
            for heldout in island_ids:
                rows = groups[target, model, heldout]
                require(len(rows) == n - 1, "Incomplete model-fold pair coverage")
                mse = math.fsum(x[0] for x in rows) / (n - 1)
                mae = math.fsum(x[1] for x in rows) / (n - 1)
                row = {"target": target, "model": model, "heldout_island_id": heldout,
                       "train_pairs": (n - 1) * (n - 2) // 2, "test_pairs": n - 1,
                       "MSE": mse, "MAE": mae}
                result["model_fold_losses"].append(row)
                lookup[target, model, heldout] = row
                losses.append((mse, mae))
            result["model_summary"].append({
                "target": target, "model": model,
                "MSE": math.fsum(x[0] for x in losses) / n,
                "MAE": math.fsum(x[1] for x in losses) / n,
                "n_folds": n, "n_predictions": n * (n - 1),
            })
        for label, (reference, augmented, role) in COMPARISONS.items():
            comparison_rows = []
            for heldout in island_ids:
                ref, aug = lookup[target, reference, heldout], lookup[target, augmented, heldout]
                row = {"target": target, "comparison": label,
                       "reference_model": reference, "augmented_model": augmented,
                       "heldout_island_id": heldout, "train_pairs": ref["train_pairs"],
                       "test_pairs": ref["test_pairs"],
                       "reference_mse": ref["MSE"], "augmented_mse": aug["MSE"],
                       "delta_mse_augmented_minus_reference": aug["MSE"] - ref["MSE"],
                       "reference_mae": ref["MAE"], "augmented_mae": aug["MAE"],
                       "improved": aug["MSE"] < ref["MSE"],
                       "worse": aug["MSE"] > ref["MSE"],
                       "tied": aug["MSE"] == ref["MSE"]}
                comparison_rows.append(row)
                result["comparison_fold_losses"].append(row)
            rm = math.fsum(row["reference_mse"] for row in comparison_rows) / n
            am = math.fsum(row["augmented_mse"] for row in comparison_rows) / n
            require(rm > 0, "Zero reference MSE makes the fixed percentage reduction undefined")
            result["comparison_summary"].append({
                "comparison": label, "target": target, "comparison_role": role,
                "reference_model": reference, "augmented_model": augmented,
                "MSE_reference": rm, "MSE_augmented": am,
                "relative_MSE_reduction_percent": 100 * (rm - am) / rm,
                "MAE_reference": math.fsum(row["reference_mae"] for row in comparison_rows) / n,
                "MAE_augmented": math.fsum(row["augmented_mae"] for row in comparison_rows) / n,
                "absolute_MSE_gain": rm - am,
                "islands_improved": sum(row["improved"] for row in comparison_rows),
                "islands_worse": sum(row["worse"] for row in comparison_rows),
                "islands_tied": sum(row["tied"] for row in comparison_rows),
            })
    result["prediction_checks"] = {
        "rows": len(predictions), "unique_source_pairs": len(source), "islands": n,
        "models": list(MODELS), "targets": list(TARGETS),
        "each_pair_predicted_once_per_endpoint_per_model_and_target": True,
        "maximum_absolute_numeric_errors": maxima,
        "aggregation": "individual prediction loss, then mean within held-out island, then unweighted mean across islands",
        "endpoint_predictions_averaged_before_scoring": False,
    }
    return result


def indexed(rows, keys, label):
    result = {}
    for row in rows:
        key = tuple(str(row[col]) for col in keys)
        require(key not in result, f"duplicate key in {label}: {key}")
        result[key] = row
    return result


def compare_tables(reconstructed, saved):
    checks = {}
    for name, keys in TABLE_KEYS.items():
        expected = indexed(reconstructed[name], keys, name)
        actual = indexed(saved[name], keys, name)
        require(expected.keys() == actual.keys(), f"{name} coverage differs from saved predictions")
        maxima = {}
        for key, row in expected.items():
            for column, want in row.items():
                require(column in actual[key], f"Missing {name}.{column}")
                value = actual[key][column]
                if isinstance(want, bool):
                    require(boolean(value, f"{name}.{column}") == want,
                            f"Mismatch {name}.{column} at {key}")
                elif isinstance(want, int):
                    require(finite(value, f"{name}.{column}") == want,
                            f"Mismatch {name}.{column} at {key}")
                elif isinstance(want, float):
                    close(value, want, column, maxima)
                else:
                    require(str(value) == want, f"Mismatch {name}.{column} at {key}")
        checks[name] = {"checked_rows": len(expected), "all_values_match": True,
                        "maximum_absolute_numeric_errors": maxima}
    return checks


def verify_sources(repo, spec_directory):
    """Independently check the frozen lock, source bundle and circuit labels."""
    import numpy as np

    lock_path = spec_directory / "input_lock.json"
    lock = json.loads(lock_path.read_text())
    require(lock["status"] == "FIXED_BEFORE_NEW_JOINT_MODEL_FITS" and
            lock["new_joint_model_results_already_seen"] is False,
            "Joint specification is not recorded as frozen before the new fits")
    require((lock["fixed_islands"], lock["fixed_pairs"]) == (51, 1275), "Locked support changed")
    require(lock["primary_comparisons"] == ["C_to_HC", "H_to_HC"], "Primary comparisons changed")
    require(digest(repo / lock["protocol_path"]) == lock["protocol_sha256"], "Locked protocol hash changed")
    require(digest(spec_directory / "IMPLEMENTATION_PLAN.md") == lock["implementation_plan_sha256"],
            "Locked implementation plan changed")
    require(len(lock["files"]) == 12, "Expected all 12 frozen inputs")
    input_hashes = {}
    for entry in lock["files"]:
        path = repo / entry["path"]
        sha = digest(path)
        require(path.stat().st_size == entry["bytes"] and sha == entry["sha256"],
                f"Locked input changed: {entry['path']}")
        input_hashes[entry["path"]] = sha
    require(len(input_hashes) == 12, "Duplicate path in frozen source lock")

    manifest_path = repo / "data/exploratory/island51_observed_ocean_20261106/source_bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    require(digest(repo / manifest["protocol_path"]) == manifest["protocol_sha256_before_outcome_fit"],
            "Prior frozen source protocol changed")
    for entry in manifest["files"]:
        path = manifest_path.parent / entry["path"]
        require(path.stat().st_size == entry["bytes"] and digest(path) == entry["sha256"],
                f"Prior frozen bundle input changed: {entry['path']}")

    source_pairs = read_csv(repo / ORIGINAL / "island_host_network_pairs.csv")
    island_rows = read_csv(repo / ORIGINAL / "island_host_network_islands.csv")
    eligible = {row["landmass_id"]: row for row in island_rows if boolean(row["eligible"], "eligible")}
    ids = sorted(eligible, key=int)
    require(len(ids) == 51 and len(source_pairs) == 1275, "Original source support changed")
    prep = json.loads((repo / CIRCUIT / "source_preparation.json").read_text())
    require(set(prep["by_island"]) == set(ids), "Source-preparation island IDs changed")
    for island in ids:
        require(int(eligible[island]["ports"]) == prep["by_island"][island]["original_ports"],
                f"Original port count mismatch at island {island}")
    matrices = {
        "uniform": CIRCUIT / "uniform_distributed_port_resistance.npz",
        "wind": CIRCUIT / "wind_distributed_port_resistance.npz",
        "native_ais": NATIVE / "native_AIS_distributed_port_resistance.npz",
    }
    matrix_checks = {}
    for label, path in matrices.items():
        with np.load(repo / path, allow_pickle=False) as archive:
            names = archive["island_ids"].astype(str).tolist()
            values = archive["resistance"]
            require(names == ids, f"{label} circuit island order mismatch")
            require(values.shape == (51, 51) and np.isfinite(values).all(), f"Invalid {label} matrix")
            require(np.allclose(values, values.T, rtol=RTOL, atol=ATOL), f"Asymmetric {label} matrix")
            require(np.array_equal(np.diag(values), np.zeros(51)), f"Nonzero {label} diagonal")
            require((values[~np.eye(51, dtype=bool)] > 0).all(), f"Nonpositive {label} off-diagonal")
        matrix_checks[label] = {"path": str(path), "sha256": digest(repo / path),
                                "canonical_island_order_verified": True}
    return source_pairs, ids, {
        "input_lock_sha256": digest(lock_path), "protocol_sha256": lock["protocol_sha256"],
        "implementation_plan_sha256": lock["implementation_plan_sha256"],
        "frozen_at_utc": lock["frozen_at_utc"], "locked_inputs_verified": input_hashes,
        "prior_bundle_files_verified": len(manifest["files"]),
        "matrix_checks": matrix_checks,
        "original_port_counts_verified_for_all_islands": True,
    }


def verify_prior(reconstructed, prior_folds, prior_summaries):
    """Compare independently reconstructed G/C losses with frozen native results."""
    prior = indexed([row for row in prior_folds
                     if row["comparison"] == "uniform_vs_ais_and_wind" and row["target"] in TARGETS],
                    ("target", "heldout_island_id"), "prior_native_folds")
    summaries = indexed([row for row in prior_summaries
                         if row["comparison"] == "uniform_vs_ais_and_wind" and row["target"] in TARGETS],
                        ("target",), "prior_native_summaries")
    maxima = {}
    checked = 0
    for row in reconstructed["model_fold_losses"]:
        if row["model"] not in ("G", "C"):
            continue
        key = (row["target"], row["heldout_island_id"])
        require(key in prior, f"Missing prior native fold {key}")
        old = prior[key]
        require(int(old["train_pairs"]) == row["train_pairs"] and
                int(old["test_pairs"]) == row["test_pairs"], "Prior fold support mismatch")
        old_mse = "baseline_mse" if row["model"] == "G" else "circuit_augmented_mse"
        old_mae = "MAE_base" if row["model"] == "G" else "MAE_circuit"
        close(old[old_mse], row["MSE"], f"{row['model']}_fold_MSE", maxima)
        close(old[old_mae], row["MAE"], f"{row['model']}_fold_MAE", maxima)
        checked += 1
    require(checked == 2 * len(prior), "Prior native fold coverage differs")
    for row in reconstructed["comparison_summary"]:
        if row["comparison"] != "G_to_C":
            continue
        key = (row["target"],)
        require(key in summaries, "Missing prior native summary")
        for col in ("MSE_reference", "MSE_augmented", "relative_MSE_reduction_percent",
                    "MAE_reference", "MAE_augmented"):
            close(summaries[key][col], row[col], f"prior_summary.{col}", maxima)
        require(int(summaries[key]["islands_improved"]) == row["islands_improved"],
                "Prior native improved-island count differs")
    require(len(summaries) == len(TARGETS), "Prior summary target coverage differs")
    return {"verdict": "PASS", "model_folds_checked": checked,
            "comparison_summaries_checked": len(summaries),
            "absolute_tolerance": ATOL, "relative_tolerance": RTOL,
            "maximum_absolute_numeric_errors": maxima}


def verify_receipt_links(results_dir, source_checks, *, repo=REPO):
    """Bind saved outputs to the recorded source gate and implementation."""
    run = json.loads((results_dir / "joint_run_receipt.json").read_text())
    source = json.loads((results_dir / "source_validation.json").read_text())
    prior = json.loads((results_dir / "prior_circuit_reproduction.json").read_text())
    require(run["schema"] == "chocho_fixed_island51_joint_host_circuit_run_v1" and
            run["status"] == "FIVE_FIXED_JOINT_MODELS_AND_ALL_HELDOUT_PREDICTIONS_COMPLETE",
            "Joint run receipt does not declare a complete fixed comparison")
    output_names = {"heldout_predictions.csv.gz", "source_validation.json",
                    "prior_circuit_reproduction.json", *[name + ".csv" for name in TABLE_KEYS]}
    require(set(run["outputs_sha256"]) == output_names, "Run receipt output hash coverage changed")
    for name, sha in run["outputs_sha256"].items():
        require(digest(results_dir / name) == sha, f"Run receipt output hash mismatch: {name}")
    code_paths = {
        "scripts/run_island51_host_circuit_joint.py",
        "scripts/exploratory_island51_circuit_network_holdout.py",
        "scripts/run_island51_observed_ocean_circuits.py",
        "scripts/island51_ocean_circuit.py",
    }
    require(set(run["implementation_sha256"]) == code_paths, "Run receipt implementation hash coverage changed")
    for path, sha in run["implementation_sha256"].items():
        require(digest(repo / path) == sha, f"Run receipt implementation hash mismatch: {path}")
    require(run["source_validation"] == source and run["prior_circuit_reproduction"] == prior,
            "Embedded run receipts differ from standalone receipts")
    for key in ("input_lock_sha256", "protocol_sha256", "implementation_plan_sha256", "frozen_at_utc"):
        require(source[key] == source_checks[key], f"Run receipt source gate disagrees on {key}")
    freeze = datetime.fromisoformat(source["frozen_at_utc"])
    preflight = datetime.fromisoformat(source["completed_before_joint_fit_at_utc"])
    completion = datetime.fromisoformat(run["completed_at_utc"])
    require(all(value.tzinfo is not None for value in (freeze, preflight, completion)),
            "Run receipt timestamps must state their UTC offsets")
    require(freeze <= preflight <= completion, "Declared freeze/preflight/completion order is invalid")
    return {"output_hashes_verified": len(output_names), "implementation_hashes_verified": len(code_paths),
            "embedded_and_standalone_receipts_identical": True,
            "freeze_precedes_preflight_and_completion": True,
            "declared_source_preflight_at_utc": source["completed_before_joint_fit_at_utc"],
            "declared_run_completed_at_utc": run["completed_at_utc"]}


def audit(results_dir):
    sources, ids, source_checks = verify_sources(REPO, SPEC_DIRECTORY)
    receipt_links = verify_receipt_links(results_dir, source_checks)
    predictions_path = results_dir / "heldout_predictions.csv.gz"
    predictions = read_csv(predictions_path)
    reconstructed = reconstruct(predictions, sources, ids)
    tables = {name: read_csv(results_dir / f"{name}.csv") for name in TABLE_KEYS}
    table_checks = compare_tables(reconstructed, tables)
    prior = verify_prior(reconstructed,
                         read_csv(REPO / NATIVE / "all_island_holdout_losses.csv"),
                         read_csv(REPO / NATIVE / "model_comparison_summary.csv"))
    outputs = [predictions_path] + [results_dir / f"{name}.csv" for name in TABLE_KEYS]
    receipt_metadata = {}
    for name in ("source_validation.json", "prior_circuit_reproduction.json", "joint_run_receipt.json"):
        path = results_dir / name
        payload = json.loads(path.read_text())
        require(isinstance(payload, dict), f"Invalid runner receipt: {name}")
        outputs.append(path)
        receipt_metadata[name] = {"schema": payload.get("schema"), "status": payload.get("status"),
                                   "sha256": digest(path)}
    return {
        "schema": "chocho_island51_joint_independent_prediction_audit_v1",
        "status": "PASS_SAVED_PREDICTIONS_AND_LOCKED_COMPARISONS_INDEPENDENTLY_RECONSTRUCTED",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "auditor_sha256": digest(Path(__file__)),
        "specification_compliance": {
            "verdict": "PASS_FOR_AUDITED_OUTPUT_CONTRACT_AND_FROZEN_INPUTS",
            "checked": ["frozen protocol, plan and all 12 inputs", "prior frozen source bundle",
                        "native AIS plus unchanged uniform/wind matrix identity and island order",
                        "original 51-island/1275-pair support and original port counts",
                        "two fixed targets and five fixed model labels",
                        "all endpoint-held-out test rows and their complementary training-set sizes",
                        "original observed target values and clipped predictions",
                        "all five fixed comparison directions, roles and strict tie definitions",
                        "recorded output and implementation hashes plus declared source-gate timing"],
            "scope_limit": "Prediction files cannot alone establish the fitted feature equations, actual training exclusion or train-only scaling. Those require the runner implementation and leakage tests; this auditor does not refit any model.",
        },
        "numerical_audit": {"verdict": "PASS", "prediction_checks": reconstructed["prediction_checks"],
                            "saved_table_checks": table_checks, "prior_G_C_reproduction": prior},
        "source_checks": source_checks,
        "runner_receipt_link_checks": receipt_links,
        "audited_output_sha256": {path.name: digest(path) for path in outputs},
        "runner_receipts_recorded_not_used_as_numeric_truth": receipt_metadata,
        "reconstructed_model_fold_losses": reconstructed["model_fold_losses"],
        "reconstructed_comparison_fold_losses": reconstructed["comparison_fold_losses"],
        "reconstructed_model_summaries": reconstructed["model_summary"],
        "reconstructed_comparison_summaries": reconstructed["comparison_summary"],
        "independence": "Original source observations plus independent CSV parsing and math.fsum. No joint-runner imports, no model fitting, no averaging of endpoint predictions before scoring.",
        "interpretation_limit": "This is an arithmetic, support and provenance audit, not independent ecological replication or a causal test. Pair and fold dependence is unchanged.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=SPEC_DIRECTORY)
    args = parser.parse_args()
    receipt_path = args.results_dir / "independent_prediction_audit.json"
    try:
        receipt = audit(args.results_dir)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        receipt = {"schema": "chocho_island51_joint_independent_prediction_audit_v1",
                   "status": "FAIL", "audited_at_utc": datetime.now(timezone.utc).isoformat(),
                   "error_type": type(exc).__name__, "error": str(exc)}
        if args.results_dir.is_dir():
            receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        print(json.dumps(receipt), file=sys.stderr)
        return 1
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": receipt["status"], "receipt": str(receipt_path),
                      "prediction_rows": receipt["numerical_audit"]["prediction_checks"]["rows"],
                      "model_folds": len(receipt["reconstructed_model_fold_losses"]),
                      "comparison_folds": len(receipt["reconstructed_comparison_fold_losses"]),
                      "model_summaries": len(receipt["reconstructed_model_summaries"]),
                      "comparison_summaries": len(receipt["reconstructed_comparison_summaries"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
