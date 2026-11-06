#!/usr/bin/env python3
"""Derive the separate, post hoc native-AIS source-quality sensitivity input.

This reads the verified native TIFF directly. It never edits the primary
overview-derived input, its frozen manifest, or biological results.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import zipfile
import zlib

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_island51_1degree_commercial_AIS_source import (  # noqa: E402
    ORIGINAL_SHA,
    SOURCE_URL,
    TIFF,
    average_overview_on_grid,
    global_grid,
    hash_file,
    validate_source_coverage,
)


REPO = Path(__file__).resolve().parents[2]
BUNDLE = Path("data/exploratory/island51_observed_ocean_20261106")
RESULTS = Path("results/island51_observed_ocean_circuit_20261106")
NATIVE_TRANSFORM = [0.005, 0.0, -180.015311275, 0.0, -0.005, 85.00264793700009]


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def repo_label(path):
    try:
        return str(path.resolve().relative_to(REPO))
    except ValueError:
        return str(path.resolve())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--native-tiff", type=Path, required=True)
    parser.add_argument("--primary-grid", type=Path,
                        default=REPO / BUNDLE / "ais/commercial_AIS_mean_0p5deg_2015_2021.npz")
    parser.add_argument("--independent-audit", type=Path,
                        default=REPO / RESULTS / "ais_source_quality/published_AIS_source_quality_receipt.json")
    parser.add_argument("--plan-lock", type=Path,
                        default=REPO / RESULTS / "native_source_sensitivity_plan_lock.json")
    parser.add_argument("--out-dir", type=Path,
                        default=REPO / BUNDLE / "ais_native_quality_check")
    args = parser.parse_args()
    start = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()

    lock = json.loads(args.plan_lock.read_text())
    plan = REPO / lock["plan_path"]
    if (hash_file(plan) != lock["plan_sha256"] or
            lock["native_outcomes_already_seen"] is not False):
        raise ValueError("Native source-quality sensitivity plan lock is invalid")
    audit = json.loads(args.independent_audit.read_text())
    if audit["source_ZIP_sha256"] != ORIGINAL_SHA or len(audit["cases"]) != 39:
        raise ValueError("Independent native overlap audit identity/count changed")
    primary_sha = hash_file(args.primary_grid)
    if primary_sha != audit["coarse_grid_sha256"]:
        raise ValueError("Primary grid differs from the independently audited grid")
    source_sha = hash_file(args.source_zip)
    if source_sha != ORIGINAL_SHA:
        raise ValueError("Official World Bank source archive checksum changed")
    with zipfile.ZipFile(args.source_zip) as archive:
        member = archive.getinfo(TIFF)
    if (member.file_size != 9818136931 or member.CRC != 0x0CA2FEF6 or
            args.native_tiff.stat().st_size != member.file_size):
        raise ValueError("Extracted native TIFF member identity/size changed")

    # One bounded-memory sequential pass validates the already-extracted member.
    crc = 0
    native_hash = hashlib.sha256()
    with args.native_tiff.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            crc = zlib.crc32(chunk, crc)
            native_hash.update(chunk)
    if crc != member.CRC:
        raise ValueError("Extracted native TIFF failed the official ZIP member CRC")

    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      GDAL_NUM_THREADS="2", GDAL_CACHEMAX=128 * 1024 * 1024):
        with rasterio.open(args.native_tiff) as source:
            metadata = {
                "native_width": source.width,
                "native_height": source.height,
                "native_dtype": source.dtypes[0],
                "native_nodata": source.nodata,
                "native_crs": str(source.crs),
                "native_transform": list(source.transform)[:6],
                "native_scales": list(source.scales),
                "native_offsets": list(source.offsets),
                "native_overview_factors": source.overviews(1),
                "native_crc32_verified": f"{crc:08x}",
                "native_bytes": member.file_size,
                "native_sha256": native_hash.hexdigest(),
            }
            if (source.count != 1 or source.width != 72006 or source.height != 33998 or
                    source.dtypes[0] != "int32" or source.nodata != 2147483647 or
                    str(source.crs) != "EPSG:4326" or
                    not np.allclose(list(source.transform)[:6], NATIVE_TRANSFORM,
                                    rtol=0, atol=1e-12) or
                    source.scales != (1.0,) or source.offsets != (0.0,) or
                    source.overviews(1)):
                raise ValueError("Native metadata changed or an overview is available")
            native_transform = source.transform
        print("VERIFIED_NATIVE_BASE_NO_OVERVIEWS", json.dumps(metadata), flush=True)
        # The existing helper reads through a base-band WarpedVRT. The source
        # reports no overviews in this environment, so no pyramid can be chosen.
        values, valid = average_overview_on_grid(
            args.native_tiff, native_transform, resolution_degrees=0.5)
    count, fraction = validate_source_coverage(valid)
    shape, transform, latitudes, longitudes = global_grid(0.5)
    if values.shape != shape:
        raise ValueError("Native target shape changed")

    comparisons = []
    for case in audit["cases"]:
        row, col = case["coarse_row"], case["coarse_col"]
        expected = float(case["native_direct_coarse_overlap"]["mean"])
        observed = float(values[row, col])
        error = abs(observed - expected)
        relative_error = error / max(abs(expected), 1.0)
        passed = bool(valid[row, col] and relative_error <= 3e-7)
        comparisons.append({
            "name": case["name"], "coarse_row": row, "coarse_col": col,
            "native_direct_coarse_overlap_mean": expected,
            "derived_native_coarse_mean": observed,
            "absolute_error": error, "relative_error": relative_error,
            "passed": passed,
        })
    if not all(row["passed"] for row in comparisons):
        raise ValueError("Native grid disagrees with independent overlap calculation: " +
                         json.dumps([row for row in comparisons if not row["passed"]]))

    with np.load(args.primary_grid) as primary:
        primary_values = np.asarray(primary["AIS_mean"], dtype=np.float64)
        primary_valid = primary["source_valid_mask"].astype(bool)
        if (primary_values.shape != shape or
                not np.array_equal(primary["raster_transform"], np.array(list(transform)[:6]))):
            raise ValueError("Primary and native grid alignment differs")
    common = valid & primary_valid
    native_common = values[common].astype(np.float64)
    overview_common = primary_values[common]
    difference = overview_common - native_common
    positive = native_common > 0
    relative = difference[positive] / native_common[positive]
    source_comparison = {
        "scope": "All common valid 0.5-degree cells, including valid land zeros; unweighted cell summaries, not global vessel totals",
        "difference_sign": "primary overview-derived minus direct-native",
        "common_valid_cells": int(common.sum()),
        "same_valid_mask": bool(np.array_equal(valid, primary_valid)),
        "native_valid_mean": float(values[valid].astype(np.float64).mean()),
        "primary_valid_mean": float(primary_values[primary_valid].mean()),
        "native_positive_cells": int((values[valid] > 0).sum()),
        "primary_positive_cells": int((primary_values[primary_valid] > 0).sum()),
        "native_observed_zero_cells": int((values[valid] == 0).sum()),
        "primary_observed_zero_cells": int((primary_values[primary_valid] == 0).sum()),
        "primary_zero_native_positive_cells": int(((overview_common == 0) & positive).sum()),
        "native_zero_primary_positive_cells": int(((native_common == 0) & (overview_common > 0)).sum()),
        "mean_signed_difference": float(difference.mean()),
        "mean_absolute_difference": float(np.abs(difference).mean()),
        "median_absolute_difference": float(np.median(np.abs(difference))),
        "positive_native_relative_difference_quantiles": {
            str(q): float(np.quantile(relative, q)) for q in (0, 0.25, 0.5, 0.75, 1)},
        "positive_native_median_absolute_relative_difference": float(np.median(np.abs(relative))),
        "native_maximum": float(values[valid].max()),
        "primary_maximum": float(primary_values[primary_valid].max()),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    grid = args.out_dir / "commercial_AIS_native_mean_0p5deg_2015_2021.npz"
    receipt_path = args.out_dir / "commercial_AIS_native_mean_0p5deg_2015_2021_receipt.json"
    check_path = args.out_dir / "native_overlap_validation_39_cells.json"
    np.savez_compressed(grid, AIS_mean=values, source_valid_mask=valid,
                        raster_transform=np.array(list(transform)[:6], dtype=float),
                        lat_centers=latitudes, lon_centers=longitudes)
    grid_sha = hash_file(grid)
    audit_sha = hash_file(args.independent_audit)
    dump(check_path, {
        "schema": "chocho_native_AIS_independent_overlap_check_v1",
        "independent_overlap_audit_path": repo_label(args.independent_audit),
        "independent_overlap_audit_sha256": audit_sha,
        "native_grid_sha256": grid_sha,
        "relative_tolerance": 3e-7,
        "denominator": "max(abs(independent native mean), 1)",
        "all_passed": True,
        "cases": comparisons,
    })
    receipt = {
        "schema": "chocho_native_AIS_0p005_to_0p5_source_quality_v1",
        "status": "DIRECT_NATIVE_SOURCE_QUALITY_SENSITIVITY_INPUT_VALIDATED_NO_BIOLOGICAL_FIT",
        "started_at_utc": started_at,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "source_original_ZIP_url": SOURCE_URL,
        "source_original_SHA256": source_sha,
        "source_metadata": metadata,
        "native_grid_path": repo_label(grid),
        "native_grid_sha256": grid_sha,
        "primary_grid_path": repo_label(args.primary_grid),
        "primary_grid_sha256": primary_sha,
        "independent_overlap_audit_path": repo_label(args.independent_audit),
        "independent_overlap_audit_sha256": audit_sha,
        "independent_overlap_validation_path": repo_label(check_path),
        "independent_overlap_validation_sha256": hash_file(check_path),
        "independent_overlap_cells": len(comparisons),
        "independent_overlap_all_passed": True,
        "independent_overlap_maximum_relative_error": max(row["relative_error"] for row in comparisons),
        "plan_lock_path": repo_label(args.plan_lock),
        "plan_lock_sha256": hash_file(args.plan_lock),
        "plan_sha256": lock["plan_sha256"],
        "plan_frozen_at_utc": lock["frozen_at_utc"],
        "derived_grid": {
            "shape": list(shape), "transform": list(transform)[:6], "crs": "EPSG:4326",
            "resampling": "GDAL geographic area-overlap average directly from every native 0.005-degree base-band pixel",
            "source_nodata_excluded": True, "valid_zero_included": True,
            "values_rescaled": False, "overview_used": False,
            "count_valid_pixels": count, "fraction_valid_pixels": fraction,
            "count_missing_pixels": int((~valid).sum()),
        },
        "memory_limits": {"gdal_cache_bytes": 128 * 1024 * 1024, "warp_memory_megabytes": 128},
        "native_to_primary_source_comparison": source_comparison,
        "reproducer_path": repo_label(Path(__file__)),
        "reproducer_sha256": hash_file(Path(__file__)),
        "units_note": "Published AIS density index retained without rescaling; independent physical calibration is unresolved.",
        "scope_note": "Post hoc source-quality sensitivity after primary outcomes; no primary input, manifest, protocol, result, or biological model changed by this derivation.",
    }
    dump(receipt_path, receipt)
    print("NATIVE_AIS_QUALITY_GRID_READY", json.dumps({
        "native_grid_path": repo_label(grid), "native_grid_sha256": grid_sha,
        "receipt_path": repo_label(receipt_path),
        "39_cell_maximum_relative_error": receipt["independent_overlap_maximum_relative_error"],
        "elapsed_seconds": receipt["elapsed_seconds"],
        "source_comparison": source_comparison,
    }), flush=True)


if __name__ == "__main__":
    main()
