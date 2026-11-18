#!/usr/bin/env python3
"""One post hoc source-quality sensitivity: native AIS means, fixed design.

The primary frozen bundle and results are read-only. Only AIS changes; all
source geometry, terminal weights, transformations, models and folds are reused.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from numbers import Real
from pathlib import Path
import tempfile
import time

import numpy as np
import pandas as pd
from scipy import sparse

import run_island51_observed_ocean_circuits as primary


DEFAULT_PRIMARY = primary.REPO / "results/island51_observed_ocean_circuit_20261106"


def validate_plan(plan_path, lock_path):
    lock = json.loads(Path(lock_path).read_text())
    if primary.sha256(plan_path) != lock.get("plan_sha256"):
        raise ValueError("Native sensitivity plan checksum differs from its frozen lock")
    if (lock.get("status") != "POSTHOC_SOURCE_QUALITY_PLAN_FIXED_BEFORE_NATIVE_AIS_OUTCOME_FIT" or
            lock.get("primary_outcomes_already_seen") is not True or
            lock.get("native_outcomes_already_seen") is not False):
        raise ValueError("Native sensitivity must retain its honest post hoc plan state")
    return dict(lock, plan_lock_sha256=primary.sha256(lock_path))


def verify_native_sources(grid_path, receipt_path, audit_path):
    receipt = json.loads(Path(receipt_path).read_text())
    if receipt.get("schema") != "chocho_native_AIS_0p005_to_0p5_source_quality_v1":
        raise ValueError("Native source receipt has an unrecognized derivation schema")
    if receipt.get("source_original_SHA256") != primary.AIS_SOURCE_SHA:
        raise ValueError("Native AIS archive identity differs from the frozen primary archive")
    if primary.sha256(grid_path) != receipt.get("native_grid_sha256"):
        raise ValueError("The native grid checksum differs from its derivation receipt")
    if primary.sha256(audit_path) != receipt.get("independent_overlap_audit_sha256"):
        raise ValueError("Independent native overlap audit checksum differs from the derivation receipt")
    audit = json.loads(Path(audit_path).read_text())
    metadata = receipt["source_metadata"]
    transform = np.asarray(metadata["native_transform"], dtype=float)
    if metadata.get("native_overview_factors") != []:
        raise ValueError("Direct native source unexpectedly includes an overview")
    if (metadata.get("native_dtype") != "int32" or
            metadata.get("native_nodata") != 2147483647. or
            metadata.get("native_scales") != [1.] or metadata.get("native_offsets") != [0.] or
            metadata.get("native_width") != 72006 or metadata.get("native_height") != 33998 or
            transform.shape != (6,) or
            not np.allclose(transform, [.005, 0., -180.015311275, 0., -.005, 85.002647937],
                            rtol=0, atol=1e-10)):
        raise ValueError("Native source identity, scale or original0.005-degree geometry changed")
    member = audit["archive_members"]["ShipDensity_Commercial1.tif"]
    if (int(metadata["native_crc32_verified"], 16) != int(member["CRC32"], 16) or
            metadata.get("native_bytes") != member["uncompressed_bytes"]):
        raise ValueError("Verified native source CRC or length differs from the audited archive member")
    cases = audit.get("cases", [])
    if (audit.get("schema") != "chocho_published_AIS_native_overview_coarse_source_quality_v1" or
            audit.get("source_ZIP_sha256") != primary.AIS_SOURCE_SHA or
            audit.get("diagnostic_locations") != 39 or len(cases) != 39):
        raise ValueError("The complete independent39-cell native source audit is required")
    values, valid = primary.circuit.load_ais(grid_path)
    if values.shape != (340, 720):
        raise ValueError("Native quality sensitivity must use the unchanged0.5-degree grid")
    seen, checked = set(), []
    for case in cases:
        row, col = int(case["coarse_row"]), int(case["coarse_col"])
        if (row, col) in seen or not (0 <= row < 340 and 0 <= col < 720):
            raise ValueError("Independent native audit has duplicated or invalid cell coordinates")
        seen.add((row, col))
        expected = case["native_direct_coarse_overlap"]["mean"]
        if expected is None or not np.isfinite(expected) or expected < 0 or not valid[row, col]:
            raise ValueError("An independent native overlap audit cell lacks an observed mean")
        observed = float(values[row, col])
        error = abs(observed-float(expected)) / max(abs(float(expected)), 1.)
        if error > 3e-7:
            raise ValueError(f"New grid fails independent native pixel overlap at cell({row},{col})")
        checked.append({"row": row, "col": col, "independent_native_mean": expected,
                        "derived_mean": observed, "relative_error": error})
    return values, valid, {
        "native_grid_sha256": primary.sha256(grid_path),
        "native_receipt_sha256": primary.sha256(receipt_path),
        "independent_audit_sha256": primary.sha256(audit_path),
        "independent_native_cells_verified": len(checked),
        "maximum_relative_overlap_error": max(row["relative_error"] for row in checked),
        "overlap_relative_tolerance": 3e-7,
        "all_39_checks": checked, "native_source_receipt": receipt,
    }


def assert_fixed_comparisons(original, sensitivity):
    """Verify every control metric, fold identity and fold loss, not just MSE."""
    differences = []

    def same(a, b, path):
        if isinstance(a, dict):
            if not isinstance(b, dict) or set(a) != set(b):
                raise ValueError("Unchanged comparison fields differ: " + path)
            for key in a:
                same(a[key], b[key], path+"/"+str(key))
        elif isinstance(a, list):
            if not isinstance(b, list) or len(a) != len(b):
                raise ValueError("Unchanged comparison folds differ: " + path)
            for index, (left, right) in enumerate(zip(a, b)):
                same(left, right, path+"/"+str(index))
        elif isinstance(a, Real) and not isinstance(a, bool):
            if not isinstance(b, Real) or not np.isclose(a, b, rtol=1e-10, atol=1e-12):
                raise ValueError("Unchanged comparison numerical result differs: " + path)
            differences.append(abs(float(a)-float(b)))
        elif a != b:
            raise ValueError("Unchanged comparison identity differs: " + path)

    labels = ["base_vs_uniform", "uniform_vs_wind"]
    for label in labels:
        same(original[label], sensitivity[label], label)
    return {"comparisons_verified": labels,
            "all_targets_and_individual_fold_losses_checked": True,
            "maximum_absolute_numeric_difference": max(differences, default=0.),
            "relative_tolerance": 1e-10, "absolute_tolerance": 1e-12}


def run(args):
    began = time.monotonic()
    out, reference = Path(args.out_dir), Path(args.primary_results_dir)
    if out.resolve() == reference.resolve():
        raise ValueError("Sensitivity output must not replace the primary result directory")
    plan = validate_plan(args.plan, args.plan_lock)
    native, valid, source_check = verify_native_sources(args.native_grid, args.native_receipt, args.native_audit)
    primary_files = [reference / "source_preparation.json", reference / "observed_ocean_circuit_holdout.json"]
    primary_files += [reference / f"{label}_distributed_port_resistance.npz" for label in ("uniform", "ais", "wind")]
    primary_files += [Path(args.input_root) / "source_bundle_manifest.json", primary.PROTOCOL]
    original_hashes = {str(path.resolve()): primary.sha256(path) for path in primary_files}
    original = json.loads((reference / "observed_ocean_circuit_holdout.json").read_text())
    old_preparation = json.loads((reference / "source_preparation.json").read_text())
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="chocho_native_sensitivity_primary_") as work:
        q, ids, counts, graph, old_edges, weights, preparation = primary.prepare_inputs(
            args.input_root, args.coast_dir or Path(args.input_root) / "coast_0p0625",
            Path(work) / "preparation", Path(work) / "primary_networks", marine_ports_only=True)
    if (preparation["exported_network"]["files_sha256"] != old_preparation["exported_network"]["files_sha256"] or
            preparation["frozen_source_bundle"]["manifest_sha256"] != old_preparation["frozen_source_bundle"]["manifest_sha256"]):
        raise ValueError("Primary source geometry, terminal weights or conductance exports changed")
    matrices = {}
    for label in ("uniform", "ais", "wind"):
        with np.load(reference / f"{label}_distributed_port_resistance.npz") as data:
            if data["island_ids"].tolist() != ids:
                raise ValueError("Original circuit matrix island order differs from the source roster")
            matrices[label] = np.asarray(data["resistance"], dtype=float)
    original_ais = matrices["ais"].copy()
    node_row, node_col = graph["node_block_row"], graph["node_block_col"]
    sea = np.zeros(native.shape, dtype=bool)
    sea[node_row, node_col] = True
    if not np.all(valid[sea]):
        raise ValueError("Native source missingness would change the fixed common marine domain")
    wind = primary.load_wind(Path(args.input_root) / "wind/NOAA_NCEP_R1_10m_wind2015to2021_1deg.npz", native.shape)
    _, fields = primary.circuit.conductance_fields(native, valid, wind, sea)
    native_edges = primary.finite_volume_edges(graph, fields["ais"][node_row, node_col])
    edges = dict(old_edges, ais=native_edges)
    component_nodes = np.flatnonzero(graph["node_component"] == preparation["marine_component"])
    exported = primary.export_networks(args.network_dir, graph, edges, weights, ids, component_nodes)
    for name in ["uniform_conductances.txt", "wind_conductances.txt", "island_node_weights.csv"]:
        if exported["files_sha256"][name] != preparation["exported_network"]["files_sha256"][name]:
            raise ValueError("Native sensitivity changed a fixed graph or terminal invariant: " + name)
    preflight = {
        "status": "NATIVE_SOURCE_AUDIT_AND_FIXED_GEOMETRY_VERIFIED_BEFORE_NATIVE_OUTCOME_FIT",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "post_hoc_source_quality_sensitivity": True, "plan": plan,
        "source_checks": source_check,
        "original_manifest_sha256": preparation["frozen_source_bundle"]["manifest_sha256"],
        "original_protocol_sha256": preparation["protocol_sha256"],
        "original_source_preparation_sha256": primary.sha256(reference / "source_preparation.json"),
        "original_sources_reverified_by_primary_prepare_inputs": True,
        "original_islands": len(ids), "original_pairs": len(q), "original_ports": 194, "marine_ports": 182,
        "unchanged_uniform_wind_networks_and_terminal_weights": True,
        "exported_network": exported,
    }
    (out / "native_source_preflight.json").write_text(json.dumps(preflight, indent=2)+"\n")
    print("NATIVE_SOURCE_SENSITIVITY_PREFLIGHT_PASS", json.dumps({
        "independent_native_cells": source_check["independent_native_cells_verified"],
        "maximum_native_mean_error": source_check["maximum_relative_overlap_error"],
        "islands": len(ids), "marine_ports": 182}), flush=True)
    i, j = graph["edge_i"], graph["edge_j"]
    g = native_edges["conductance"]
    n = weights.shape[1]
    adjacency = sparse.csr_matrix((np.r_[g, g], (np.r_[i, j], np.r_[j, i])), shape=(n, n))
    laplacian = sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency
    matrices["ais"], diagnostics = primary.circuit.distributed_effective_resistance(laplacian, weights)
    np.savez_compressed(out / "native_AIS_distributed_port_resistance.npz",
                        island_ids=np.array(ids), resistance=matrices["ais"])
    comparisons = primary.compare_models(q, ids, counts, matrices)
    invariants = assert_fixed_comparisons(original["comparisons"], comparisons)
    rows, folds, side_by_side = [], [], []
    for label, comparison in comparisons.items():
        for target, fitted in comparison["targets"].items():
            old = original["comparisons"][label]["targets"][target]
            rows.append({"comparison": label, "target": target,
                         "MSE_reference": fitted["MSE_baseline"], "MSE_augmented": fitted["MSE_augmented"],
                         "relative_MSE_reduction_percent": fitted["relative_MSE_reduction_percent"],
                         "MAE_reference": fitted["MAE_baseline"], "MAE_augmented": fitted["MAE_augmented"],
                         "islands_improved": fitted["islands_with_better_augmented_prediction"]})
            side_by_side.append({"comparison": label, "target": target,
                "primary_overview_MSE": old["MSE_augmented"], "native_AIS_MSE": fitted["MSE_augmented"],
                "primary_overview_MAE": old["MAE_augmented"], "native_AIS_MAE": fitted["MAE_augmented"],
                "primary_overview_MSE_reduction_percent": old["relative_MSE_reduction_percent"],
                "native_AIS_MSE_reduction_percent": fitted["relative_MSE_reduction_percent"],
                "MSE_reduction_change_percentage_points": fitted["relative_MSE_reduction_percent"]-old["relative_MSE_reduction_percent"],
                "primary_islands_improved": old["islands_with_better_augmented_prediction"],
                "native_islands_improved": fitted["islands_with_better_augmented_prediction"]})
            folds.extend(dict(comparison=label, target=target, **fold) for fold in fitted["per_island_LOIO"])
    if len(rows) != 12 or len(folds) != 612:
        raise ValueError("Native sensitivity did not retain all12 comparisons and612 folds")
    pair_output = q.drop(columns=["i", "j"]).copy()
    pi, pj = q.i.to_numpy(), q.j.to_numpy()
    for label, matrix in matrices.items():
        pair_output[f"{label}_distributed_port_resistance"] = matrix[pi, pj]
        pair_output[f"log1p_normalized_{label}_circuit"] = primary.circuit.resistance_feature(matrix, pi, pj)
    pair_output["primary_overview_ais_distributed_port_resistance"] = original_ais[pi, pj]
    pair_output.to_csv(out / "island51_observed_circuit_pairs.csv", index=False, float_format="%.17g")
    pd.DataFrame(rows).to_csv(out / "model_comparison_summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(folds).to_csv(out / "all_island_holdout_losses.csv", index=False, float_format="%.17g")
    pd.DataFrame(side_by_side).to_csv(out / "primary_vs_native_comparison.csv", index=False, float_format="%.17g")
    before, after = original_ais[pi, pj], matrices["ais"][pi, pj]
    matrix_comparison = {
        "pairs": len(before), "raw_resistance_Pearson_correlation": float(np.corrcoef(before, after)[0, 1]),
        "median_native_over_primary_raw_resistance": float(np.median(after/before)),
        "maximum_absolute_relative_resistance_change": float(np.max(np.abs(after-before)/before)),
        "native_matrix_sha256": primary.sha256(out / "native_AIS_distributed_port_resistance.npz"),
    }
    for path, expected_hash in original_hashes.items():
        if primary.sha256(path) != expected_hash:
            raise ValueError("A protected primary input/result changed during sensitivity: " + path)
    result = {
        "schema": "chocho_single_native_AIS_source_quality_sensitivity_v1",
        "status": "SINGLE_POSTHOC_NATIVE_AIS_SOURCE_SENSITIVITY_COMPLETE",
        "post_hoc_after_primary_outcomes": True, "native_plan_fixed_before_native_outcomes": True,
        "conductivity_tuning_or_model_sweep_performed": False,
        "preflight": preflight, "numerical_diagnostics": diagnostics, "control_invariants": invariants,
        "comparisons": comparisons, "primary_vs_native": side_by_side,
        "matrix_comparison": matrix_comparison,
        "protected_primary_inputs_and_results_unchanged": True,
        "protected_primary_files_sha256": original_hashes,
        "elapsed_seconds": time.monotonic()-began,
        "official_Circuitscape_native_check_performed_by_this_script": False,
        "interpretation": "One source-quality sensitivity chosen after seeing the primary outcomes; not independent validation and not a replacement of the frozen primary. Only AIS pixel averaging changed.",
    }
    (out / "native_source_sensitivity.json").write_text(json.dumps(result, indent=2)+"\n")
    print("NATIVE_SOURCE_SENSITIVITY_COMPLETE", json.dumps(side_by_side), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--coast-dir", type=Path)
    parser.add_argument("--primary-results-dir", type=Path, default=DEFAULT_PRIMARY)
    parser.add_argument("--native-grid", required=True, type=Path)
    parser.add_argument("--native-receipt", required=True, type=Path)
    parser.add_argument("--native-audit", required=True, type=Path)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PRIMARY / "NATIVE_SOURCE_SENSITIVITY_PLAN.md")
    parser.add_argument("--plan-lock", type=Path, default=DEFAULT_PRIMARY / "native_source_sensitivity_plan_lock.json")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--network-dir", required=True, type=Path)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
