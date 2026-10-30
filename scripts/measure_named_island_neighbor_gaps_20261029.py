#!/usr/bin/env python3
"""Compare eight named islands' direct-continental and nearest-other coast gaps.

Distances concern source geography, not butterfly sources, movement or actual
stepping stones. Only derived keys, areas and bounded measurements are exported.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from measure_island_port_geography_20261029 import (
    ROOT, R_KM, SHA, URL, CONTINENT_SEEDS, ISLAND_SEEDS, source_landmasses,
    shoreline_xyz, coast_gap_interval, derived_key, seed_indices, unit_xyz,
    chord_to_km, write_csv, digest,
)


def spherical_cap(points, step_km):
    """A cap enclosing sampled coast AND every retained continuous coast arc.

    Choose a sample as center; maximum sampled angular distance plus half the
    maximum arc subdivision step encloses every coast point by triangle
    inequality. A tiny outward numerical margin avoids pruning at roundoff ties.
    """
    center = points[0]
    angles = np.arctan2(np.linalg.norm(np.cross(points, center), axis=1),
                        np.clip(points @ center, -1, 1))
    radius = min(math.pi * R_KM, float(angles.max()) * R_KM + step_km / 2 + 1e-7)
    return center, radius


def nearest_other(query_index, samples, step_km, caps=None):
    """Exact candidate pruning, with a conservative interval for coast sampling."""
    if not 0 < step_km <= 5:
        raise ValueError("Sampling step must be in (0,5] km")
    caps = caps or {i: spherical_cap(p, step_km) for i, p in samples.items()}
    qa, qr = caps[query_index]
    candidates = []
    for index, (center, radius) in caps.items():
        if index == query_index:
            continue
        separation = float(np.arctan2(np.linalg.norm(np.cross(qa, center)),
                                     np.clip(qa @ center, -1, 1))) * R_KM
        lower_bound = max(0.0, separation - qr - radius - 1e-7)
        candidates.append((lower_bound, index))
    candidates.sort()
    best_upper, best_index = math.inf, None
    evaluated = []
    query_tree = cKDTree(samples[query_index])
    first_pruned_lower = None
    for cap_lower, index in candidates:
        if cap_lower >= best_upper:
            first_pruned_lower = cap_lower
            break
        # Query the smaller shoreline against a tree on the larger shoreline;
        # querying every continental sample against a small island is costly.
        if len(samples[index]) > len(samples[query_index]):
            lo, hi = coast_gap_interval(samples[query_index], cKDTree(samples[index]), step_km)
        else:
            lo, hi = coast_gap_interval(samples[index], query_tree, step_km)
        # Both bounds cover the continuous retained coast, so their maximum is
        # valid; using the minimum over candidates bounds the global minimum.
        evaluated.append((max(cap_lower, lo), hi, index))
        if hi < best_upper:
            best_upper, best_index = hi, index
    if best_index is None:
        raise ValueError("No other source landmass available")
    lower = min(x[0] for x in evaluated)
    possible = [i for lo, _, i in evaluated if lo <= best_upper]
    return {"neighbor_index": best_index, "lower_km": lower, "upper_km": best_upper,
            "candidates_evaluated": len(evaluated),
            "candidates_pruned": len(candidates) - len(evaluated),
            "first_pruned_cap_lower_km": first_pruned_lower,
            "potential_nearest_neighbor_indices": possible}


def run(cache, output, step=1.0):
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    land = source_landmasses(cache)
    named = seed_indices(land, ISLAND_SEEDS)
    continents = seed_indices(land, CONTINENT_SEEDS)
    candidates = land[land.area_km2 >= 5]
    if len(candidates) != 8265 or not set(continents.values()) <= set(candidates.index):
        raise ValueError("Candidate universe differs from 8260 islands + 5 continents")
    samples, caps = {}, {}
    for number, (ix, row) in enumerate(candidates.iterrows(), 1):
        points, _ = shoreline_xyz(row.geometry, step)
        samples[int(ix)] = points
        caps[int(ix)] = spherical_cap(points, step)
        if number % 2000 == 0:
            print(json.dumps({"coastlines_sampled": number}), flush=True)
    continental_tree = cKDTree(np.concatenate([samples[ix] for ix in continents.values()]))
    rows = []
    details = {}
    for name, ix in named.items():
        result = nearest_other(ix, samples, step, caps)
        other = result["neighbor_index"]
        direct_lo, direct_hi = coast_gap_interval(samples[ix], continental_tree, step)
        neighbor_names = [n for n, i in {**continents, **named}.items() if i == other]
        row = {"named_island": name,
               "island_key": derived_key(land.loc[ix, "parent_source_ID"]),
               "island_area_km2": float(land.loc[ix, "area_km2"]),
               "sampled_nearest_neighbor_key": derived_key(land.loc[other, "parent_source_ID"]),
               "sampled_nearest_neighbor_seed_label": "|".join(neighbor_names),
               "sampled_nearest_neighbor_area_km2": float(land.loc[other, "area_km2"]),
               "nearest_other_coast_gap_lower_km": result["lower_km"],
               "nearest_other_coast_gap_upper_km": result["upper_km"],
               "direct_continental_coast_gap_lower_km": direct_lo,
               "direct_continental_coast_gap_upper_km": direct_hi,
               "sampled_gap_difference_km": direct_hi - result["upper_km"],
               "neighbor_identity_ambiguous_under_sampling": int(len(result["potential_nearest_neighbor_indices"]) > 1),
               "potential_nearest_neighbor_units": len(result["potential_nearest_neighbor_indices"]),
               "candidates_evaluated": result["candidates_evaluated"],
               "candidates_pruned": result["candidates_pruned"]}
        rows.append(row)
        details[name] = {"first_pruned_cap_lower_km": result["first_pruned_cap_lower_km"],
                         "potential_nearest_neighbor_keys": [derived_key(land.loc[i, "parent_source_ID"])
                            for i in result["potential_nearest_neighbor_indices"]]}
        print(json.dumps(row), flush=True)
    csv_path = output / "named_island_neighbor_gaps.csv"
    write_csv(csv_path, rows)
    manifest = {"schema": "chocho_eight_named_island_nearest_other_coast_gap_20261029_v1",
        "status": "BOUNDED_SOURCE_GEOGRAPHY_DIAGNOSTIC_NO_BUTTERFLY_EFFECT",
        "gshhg_source": {"url": URL, "zip_sha256": SHA, "resolution": "h", "level": 1},
        "candidate_units": len(candidates), "candidate_min_area_km2": 5,
        "candidate_definition": "All 8265 dissolved GSHHG source L1 landmasses >=5km2, including 5 continental units; exclude the focal source unit itself. Antarctica absent from this source L1.",
        "method": {"sphere_radius_km": R_KM, "max_coast_segment_km": step,
            "distance_interval_width_at_most_km": step,
            "geometry": "Minimum sampled spherical coast-to-coast separation; retains all source coast samples and excludes artificial dateline cuts using the shared measurement code.",
            "pruning": "Each coast's spherical cap is centered on a coast sample, radius=max sampled great-circle radius + half-step + numerical inflation. Continuous coast gap lower bound=max(0, center distance - radii). Evaluate candidates by increasing cap lower bound; stop only when remaining lower bounds >= best sampled upper distance.",
            "global_interval": "Minimum across evaluated candidate lower intervals (max of cap lower and sampled upper minus step), and minimum sampled upper. Pruned candidates have continuous-coast lower bounds >= the chosen upper.",
            "neighbor_identity": "Reported neighbor is nearest among sampled source coasts; all units whose continuous-coast lower interval overlaps winning upper are listed separately. Unique identity depends on this sampling check.",
            "bounds_cover": "Coast densification only, not source positional accuracy, source age or spherical-versus-ellipsoidal distance."},
        "interpretation": "Direct-continental isolation and immediate nearest-other-landmass separation are distinct geographic quantities. Neighbors are not validated butterfly source populations, effective stepping stones or dispersal pathways. No port or biological model is fitted.",
        "measurements": rows, "candidate_pruning_qc": details,
        "csv_sha256": digest(csv_path)}
    (output / "named_island_neighbor_gaps.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/exploratory/island_port_20261029")
    parser.add_argument("--max-coast-segment-km", type=float, default=1.0)
    args = parser.parse_args()
    run(args.cache_dir, args.output_dir, args.max_coast_segment_km)
