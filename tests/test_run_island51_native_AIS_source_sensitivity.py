"""The single native-source sensitivity must preserve its plan and controls."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("pandas")
pytest.importorskip("scipy")

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/run_island51_native_AIS_source_sensitivity.py"


def runner():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("native_AIS_sensitivity", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_changed_native_source_plan_is_rejected_before_analysis(tmp_path):
    m = runner()
    plan = tmp_path / "plan.md"
    plan.write_text("Only native AIS averaging changes; no outcome tuning.\n")
    lock = tmp_path / "plan_lock.json"
    lock.write_text(json.dumps({
        "status": "POSTHOC_SOURCE_QUALITY_PLAN_FIXED_BEFORE_NATIVE_AIS_OUTCOME_FIT",
        "plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
        "primary_outcomes_already_seen": True,
        "native_outcomes_already_seen": False,
        "frozen_at_utc": "2026-11-06T08:54:25.598496+00:00",
    }))
    checked = m.validate_plan(plan, lock)
    assert checked["primary_outcomes_already_seen"] is True
    plan.write_text("An outcome-selected new conductivity transform.\n")
    with pytest.raises(ValueError, match="plan.*checksum"):
        m.validate_plan(plan, lock)


def fixed_comparisons():
    values = {}
    for comparison in ["base_vs_uniform", "uniform_vs_wind", "uniform_vs_ais", "uniform_vs_ais_and_wind"]:
        values[comparison] = {"targets": {target: {
            "MSE_baseline": .1, "MSE_augmented": .09,
            "per_island_LOIO": [{"heldout_island_id": "9", "baseline_mse": .2, "circuit_augmented_mse": .19}],
        } for target in ["host_jaccard", "new_butterfly_jaccard", "earlier_butterfly_jaccard"]}}
    return values


def test_native_AIS_sensitivity_cannot_change_uniform_or_wind_only_folds():
    m = runner()
    primary = fixed_comparisons()
    changed = copy.deepcopy(primary)
    changed["uniform_vs_ais"]["targets"]["host_jaccard"]["MSE_augmented"] = .08
    checked = m.assert_fixed_comparisons(primary, changed)
    assert checked["maximum_absolute_numeric_difference"] == 0.
    changed["uniform_vs_wind"]["targets"]["host_jaccard"]["per_island_LOIO"][0]["circuit_augmented_mse"] += .01
    with pytest.raises(ValueError, match="uniform_vs_wind"):
        m.assert_fixed_comparisons(primary, changed)


def test_native_AIS_sensitivity_cannot_change_control_fold_identity():
    m = runner()
    primary = fixed_comparisons()
    changed = copy.deepcopy(primary)
    changed["base_vs_uniform"]["targets"]["host_jaccard"]["per_island_LOIO"][0]["heldout_island_id"] = "15"
    with pytest.raises(ValueError, match="base_vs_uniform"):
        m.assert_fixed_comparisons(primary, changed)


def native_source_fixture(tmp_path):
    grid, receipt, audit = (tmp_path / name for name in ("native.npz", "native_receipt.json", "audit.json"))
    values = np.ones((340, 720), dtype=np.float32)
    np.savez_compressed(grid, AIS_mean=values, source_valid_mask=np.ones_like(values, dtype=bool),
                        raster_transform=np.array([.5, 0., -180., 0., -.5, 85.]))
    archive_sha = "0cfdce41d9934982bbc8f16bd92ce7545cc088839327cbaa9b0ed1469dda31e8"
    transform = [.005, 0., -180.015311275, 0., -.005, 85.002647937]
    audited = {"schema": "chocho_published_AIS_native_overview_coarse_source_quality_v1",
               "source_ZIP_sha256": archive_sha, "diagnostic_locations": 39,
               "archive_members": {"ShipDensity_Commercial1.tif": {"CRC32": "0xca2fef6", "uncompressed_bytes": 9818136931}},
               "cases": [{"coarse_row": row, "coarse_col": 10,
                          "native_direct_coarse_overlap": {"mean": 1.}} for row in range(39)]}
    audit.write_text(json.dumps(audited))
    recorded = {"schema": "chocho_native_AIS_0p005_to_0p5_source_quality_v1",
                "source_original_SHA256": archive_sha,
                "native_grid_sha256": hashlib.sha256(grid.read_bytes()).hexdigest(),
                "independent_overlap_audit_sha256": hashlib.sha256(audit.read_bytes()).hexdigest(),
                "source_metadata": {"native_transform": transform, "native_dtype": "int32",
                    "native_nodata": 2147483647., "native_scales": [1.], "native_offsets": [0.],
                    "native_overview_factors": [], "native_width": 72006, "native_height": 33998,
                    "native_crc32_verified": "0ca2fef6", "native_bytes": 9818136931}}
    receipt.write_text(json.dumps(recorded))
    return grid, receipt, audit, recorded


def test_native_grid_is_checked_against_all_39_independent_overlap_means(tmp_path):
    m = runner()
    grid, receipt, audit, recorded = native_source_fixture(tmp_path)
    values, valid, verified = m.verify_native_sources(grid, receipt, audit)
    assert values.shape == valid.shape == (340, 720)
    assert verified["independent_native_cells_verified"] == 39
    assert verified["maximum_relative_overlap_error"] == 0.
    # Even a newly hashed grid must fail if it contradicts the independent
    # native-pixel overlap audit at one of the39 preselected locations.
    values[20, 10] = 2.
    np.savez_compressed(grid, AIS_mean=values, source_valid_mask=valid,
                        raster_transform=np.array([.5, 0., -180., 0., -.5, 85.]))
    recorded["native_grid_sha256"] = hashlib.sha256(grid.read_bytes()).hexdigest()
    receipt.write_text(json.dumps(recorded))
    with pytest.raises(ValueError, match="independent native.*overlap"):
        m.verify_native_sources(grid, receipt, audit)


def test_native_grid_checksum_and_true_native_source_are_required(tmp_path):
    m = runner()
    grid, receipt, audit, recorded = native_source_fixture(tmp_path)
    recorded["native_grid_sha256"] = "0"*64
    receipt.write_text(json.dumps(recorded))
    with pytest.raises(ValueError, match="native grid.*checksum"):
        m.verify_native_sources(grid, receipt, audit)
    recorded["native_grid_sha256"] = hashlib.sha256(grid.read_bytes()).hexdigest()
    recorded["source_metadata"]["native_overview_factors"] = [2]
    receipt.write_text(json.dumps(recorded))
    with pytest.raises(ValueError, match="native.*overview"):
        m.verify_native_sources(grid, receipt, audit)


def test_verified_native_CRC_must_equal_the_independent_archive_member(tmp_path):
    m = runner()
    grid, receipt, audit, recorded = native_source_fixture(tmp_path)
    recorded["source_metadata"]["native_crc32_verified"] = "deadbeef"
    receipt.write_text(json.dumps(recorded))
    with pytest.raises(ValueError, match="native.*CRC"):
        m.verify_native_sources(grid, receipt, audit)
