#!/usr/bin/env python3
"""Measure two geographic neighbor definitions for the frozen island roster.

The 5 and 1000 km² candidate thresholds describe landmasses, not occupied
butterfly source areas or effective dispersal stepping stones.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from measure_island_port_geography_20261029 import (
    ROOT, R_KM, SHA, URL, source_landmasses, shoreline_xyz, coast_gap_interval,
    derived_key, unit_xyz, write_csv, digest,
)
from measure_named_island_neighbor_gaps_20261029 import spherical_cap


def candidate_indices(areas, focal, threshold):
    if threshold < 5:
        raise ValueError("Candidate landmasses must be >=5 km²")
    return sorted(i for i, area in areas.items() if i != focal and area >= threshold)


def measure_candidates(focal, candidates, samples, caps, step, trees):
    """Same conservative cap pruning as the named-island audit; shared trees."""
    if not candidates or focal in candidates or len(candidates) != len(set(candidates)):
        raise ValueError("Nonempty, unique, non-focal alternatives required")
    qa, qr = caps[focal]
    centers = np.array([caps[i][0] for i in candidates])
    radii = np.array([caps[i][1] for i in candidates])
    separation = np.arctan2(np.linalg.norm(np.cross(centers, qa), axis=1),
                            np.clip(centers @ qa, -1, 1))*R_KM
    bounds = np.maximum(0, separation-qr-radii-1e-7)
    order = sorted(zip(bounds.tolist(), candidates))
    best_upper, best_index = math.inf, None
    evaluated = []
    first_pruned = None
    for lower, ix in order:
        if lower >= best_upper:
            first_pruned = lower
            break
        larger, smaller = (ix, focal) if len(samples[ix]) > len(samples[focal]) else (focal, ix)
        if larger not in trees:
            trees[larger] = cKDTree(samples[larger])
        lo, hi = coast_gap_interval(samples[smaller], trees[larger], step)
        evaluated.append((max(lower, lo), hi, ix))
        if hi < best_upper:
            best_upper, best_index = hi, ix
    return {"neighbor_index": best_index, "lower_km": min(x[0] for x in evaluated),
            "upper_km": best_upper, "candidates_evaluated": len(evaluated),
            "candidates_pruned": len(candidates)-len(evaluated),
            "first_pruned_cap_lower_km": first_pruned,
            "possible_indices": [ix for lo, hi, ix in evaluated if lo <= best_upper]}


def run(source_cache, output):
    selection_path = output/"expanded_island_panel_selection.json"
    lock = json.loads(selection_path.read_bytes())
    selected = lock["selection"]["selected"]
    land = source_landmasses(source_cache)
    source = land[land.area_km2 >= 5]
    if len(source) != 8265:
        raise ValueError("Original >=5km² landmass inventory changed")
    bykey = {derived_key(row.parent_source_ID): int(ix) for ix, row in source.iterrows()}
    areas = {int(ix): float(row.area_km2) for ix, row in source.iterrows()}
    if any(r["island_key"] not in bykey for r in selected):
        raise ValueError("Frozen selected island is not in original landmass source")
    samples, caps, trees = {}, {}, {}
    for number, (ix, row) in enumerate(source.iterrows(), 1):
        points, _ = shoreline_xyz(row.geometry, 1)
        samples[int(ix)] = points
        caps[int(ix)] = spherical_cap(points, 1)
        if number % 2000 == 0:
            print(json.dumps({"coastlines_sampled": number}), flush=True)
    rows, details = [], {}
    for number, r in enumerate(selected, 1):
        focal = bykey[r["island_key"]]
        out = {"island_key": r["island_key"], "selection_stratum": r["selection_stratum"],
               "area_km2": float(r["area_km2"]),
               "continental_gap_lower_km": float(r["mainland_coast_gap_lower_km"]),
               "continental_gap_upper_km": float(r["mainland_coast_gap_upper_km"])}
        details[r["island_key"]] = {}
        for threshold in (5, 1000):
            candidates = candidate_indices(areas, focal, threshold)
            measured = measure_candidates(focal, candidates, samples, caps, 1, trees)
            other = measured["neighbor_index"]
            prefix = f"other_ge{threshold}km2_"
            out.update({prefix+"coast_gap_lower_km": measured["lower_km"],
                        prefix+"coast_gap_upper_km": measured["upper_km"],
                        prefix+"sampled_neighbor_key": derived_key(land.loc[other, "parent_source_ID"]),
                        prefix+"sampled_neighbor_area_km2": areas[other],
                        prefix+"identity_ambiguous": int(len(measured["possible_indices"]) > 1),
                        prefix+"candidate_count": len(candidates)})
            details[r["island_key"]][str(threshold)] = {
                "evaluated": measured["candidates_evaluated"], "pruned": measured["candidates_pruned"],
                "first_pruned_cap_lower_km": measured["first_pruned_cap_lower_km"],
                "possible_neighbor_keys": [derived_key(land.loc[i, "parent_source_ID"]) for i in measured["possible_indices"]]}
        if out["other_ge5km2_coast_gap_upper_km"] > out["other_ge1000km2_coast_gap_upper_km"]+1e-6:
            raise ValueError("Nested candidate universe produced a nonmonotone gap")
        if out["other_ge1000km2_coast_gap_upper_km"] > out["continental_gap_upper_km"]+1e-6:
            raise ValueError("Continents are a subset; nearest large-landmass gap cannot exceed continental gap")
        rows.append(out)
        if number % 16 == 0:
            print(json.dumps({"selected_islands_measured": number}), flush=True)
    csv_path = output/"expanded_island_access.csv"
    write_csv(csv_path, rows)
    remote = [r for r in rows if r["continental_gap_lower_km"] >= 500]
    summary = {"islands": len(rows), "ge5_candidate_landmasses": len(areas),
               "ge1000_candidate_landmasses": sum(a >= 1000 for a in areas.values()),
               "remote_ge500_islands": len(remote),
               "remote_with_other_ge5_certainly_within50_km": sum(r["other_ge5km2_coast_gap_upper_km"] <= 50 for r in remote),
               "remote_with_other_ge1000_certainly_within50_km": sum(r["other_ge1000km2_coast_gap_upper_km"] <= 50 for r in remote),
               "remote_with_other_ge1000_certainly_within100_km": sum(r["other_ge1000km2_coast_gap_upper_km"] <= 100 for r in remote)}
    manifest = {"schema": "chocho_expanded_island_access_20261029_v1", "summary": summary,
                "selection_file_sha256": digest(selection_path), "selection_sha256": lock["selection_sha256"],
                "gshhg": {"url": URL, "zip_sha256": SHA, "resolution": "h", "level": 1},
                "method": "Minimum spherical source coast-to-coast gaps; retained coastline arcs sampled <=1km; shared tested cap pruning; artificial dateline cuts removed. Focal unit excluded from each threshold candidate set. Bounds cover sampling only, not positional or sphere/ellipsoid error.",
                "interpretation": "Geographic candidate landmasses only; neither >=5 nor >=1000km² certifies an occupied butterfly source, host supply, habitat suitability or effective dispersal. No added biological regression chosen from these results.",
                "csv_sha256": digest(csv_path), "candidate_pruning": details}
    (output/"expanded_island_access.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n")
    print(json.dumps(summary, indent=2), flush=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"data/exploratory/island_expansion_20261029")
    args = parser.parse_args()
    run(args.source_cache, args.output_dir)
