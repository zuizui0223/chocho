#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from shapely.geometry import shape


def haversine_km(lat1, lon1, lat2, lon2):
    radius = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlon / 2) ** 2
    return 2 * radius * math.asin(min(1.0, math.sqrt(a)))


def average_ranks(values):
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def spearman(a, b):
    if len(a) < 3 or len(a) != len(b):
        return None
    x, y = average_ranks(a), average_ranks(b)
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def partial_spearman(y, x, control):
    if len(y) < 4:
        return None
    ry, rx, rc = average_ranks(y), average_ranks(x), average_ranks(control)
    design = np.column_stack([np.ones(len(rc)), rc])
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    if np.std(ey) == 0 or np.std(ex) == 0:
        return None
    return float(np.corrcoef(ey, ex)[0, 1])


def load_level3(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for feature in payload.get("features", []):
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        level1 = str(props.get("LEVEL1_COD") or "").strip()
        geometry = feature.get("geometry")
        if not code or not geometry:
            continue
        point = shape(geometry).representative_point()
        out[code] = (float(point.y), float(point.x), level1)
    if not out:
        raise RuntimeError("no WGSRPD3 geometries")
    if not all(v[2] for v in out.values()):
        raise RuntimeError("WGSRPD LEVEL1_COD missing; cannot run continental sensitivity")
    return out


def probability_greater(eval_rows, never_rows):
    if len(eval_rows) < 2 or len(never_rows) < 2:
        return None
    score = 0.0
    total = 0
    for observed in eval_rows:
        for unobserved in never_rows:
            score += 1.0 if unobserved["mismatch"] > observed["mismatch"] else 0.5 if unobserved["mismatch"] == observed["mismatch"] else 0.0
            total += 1
    return score / total


def greedy_distance_pairs(eval_rows, never_rows, caliper_km):
    candidates = []
    for i, observed in enumerate(eval_rows):
        for j, unobserved in enumerate(never_rows):
            gap = abs(observed["distance_to_train_km"] - unobserved["distance_to_train_km"])
            if gap <= caliper_km:
                candidates.append((gap, observed["code"], unobserved["code"], i, j))
    candidates.sort()
    used_eval, used_never, pairs = set(), set(), []
    for gap, _, __, i, j in candidates:
        if i in used_eval or j in used_never:
            continue
        used_eval.add(i)
        used_never.add(j)
        pairs.append((eval_rows[i], never_rows[j], gap))
    return pairs


def paired_score(pairs):
    if len(pairs) < 2:
        return None
    score = 0.0
    for observed, unobserved, _ in pairs:
        score += 1.0 if unobserved["mismatch"] > observed["mismatch"] else 0.5 if unobserved["mismatch"] == observed["mismatch"] else 0.0
    return score / len(pairs)


def summarize(rows, score_key, pair_key=None):
    chosen = [
        row for row in rows
        if row.get(score_key) is not None
        and (pair_key is None or int(row.get(pair_key, 0)) >= 2)
    ]
    if len(chosen) < 3:
        return {"species": len(chosen)}
    scores = [float(row[score_key]) for row in chosen]
    host_family = [float(row["host_family_count"]) for row in chosen]
    resource = [math.log1p(float(row["resource_units"])) for row in chosen]
    return {
        "species": len(chosen),
        "median_score": float(np.median(scores)),
        "species_above_0_5": sum(value > 0.5 for value in scores),
        "spearman_host_family": spearman(host_family, scores),
        "partial_spearman_host_family_controlling_log_resource_breadth": partial_spearman(scores, host_family, resource),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit-climate-csv", type=Path, required=True)
    parser.add_argument("--primary-result-json", type=Path, required=True)
    parser.add_argument("--level3-geojson", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    geometry = load_level3(args.level3_geojson)
    primary = json.loads(args.primary_result_json.read_text(encoding="utf-8"))
    meta = {row["species"]: row for row in primary["species"]}

    by_species = defaultdict(list)
    with args.unit_climate_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["wgsrpd3_code"] in geometry:
                by_species[row["species"]].append(row)

    output_rows = []
    for species, rows in by_species.items():
        if species not in meta:
            continue
        train_rows = [row for row in rows if row["observed_unit_split"] == "train_observed"]
        if not train_rows:
            continue
        train_points = [geometry[row["wgsrpd3_code"]] for row in train_rows]

        candidates = []
        for row in rows:
            split = row["observed_unit_split"]
            if split not in {"eval_observed", "never_observed"}:
                continue
            if str(row["effort_supported"]) != "1":
                continue
            value = str(row["climate_mismatch_to_train_niche"]).strip()
            if not value:
                continue
            lat, lon, level1 = geometry[row["wgsrpd3_code"]]
            distance = min(
                haversine_km(lat, lon, train_lat, train_lon)
                for train_lat, train_lon, _ in train_points
            )
            candidates.append(
                {
                    "code": row["wgsrpd3_code"],
                    "split": split,
                    "mismatch": float(value),
                    "distance_to_train_km": distance,
                    "level1": level1,
                }
            )

        eval_rows = [row for row in candidates if row["split"] == "eval_observed"]
        never_rows = [row for row in candidates if row["split"] == "never_observed"]
        observed_level1 = {geometry[row["wgsrpd3_code"]][2] for row in train_rows}
        observed_level1.update(row["level1"] for row in eval_rows)
        same_eval = [row for row in eval_rows if row["level1"] in observed_level1]
        same_never = [row for row in never_rows if row["level1"] in observed_level1]

        record = {
            "species": species,
            "host_family_count": meta[species]["host_family_count"],
            "resource_units": meta[species]["contemporary_host_resource_units"],
            "original_score_recomputed": probability_greater(eval_rows, never_rows),
            "same_level1_score": probability_greater(same_eval, same_never),
            "same_level1_eval_units": len(same_eval),
            "same_level1_never_units": len(same_never),
        }
        for caliper in (250, 500, 1000):
            pairs = greedy_distance_pairs(eval_rows, never_rows, caliper)
            record[f"distance_matched_score_{caliper}km"] = paired_score(pairs)
            record[f"distance_matched_pairs_{caliper}km"] = len(pairs)
            record[f"median_distance_difference_{caliper}km"] = (
                None if not pairs else float(np.median([pair[2] for pair in pairs]))
            )
        output_rows.append(record)

    payload = {
        "schema": "chocho_butterfly_climate_distance_sensitivity_v0.1",
        "status": "POSTHOC_DISTANCE_CONFOUNDING_SENSITIVITY",
        "same_wgsrpd_level1": summarize(output_rows, "same_level1_score"),
        "distance_matched_250km": summarize(output_rows, "distance_matched_score_250km", "distance_matched_pairs_250km"),
        "distance_matched_500km": summarize(output_rows, "distance_matched_score_500km", "distance_matched_pairs_500km"),
        "distance_matched_1000km": summarize(output_rows, "distance_matched_score_1000km", "distance_matched_pairs_1000km"),
        "claim_boundary": "Post-hoc spatial sensitivity. Matching reduces gross accessibility confounding but does not establish dispersal accessibility or causal climate filtering.",
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(sorted(output_rows, key=lambda row: row["species"]))
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
