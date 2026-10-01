#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def load_level3(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    level1 = {}
    for feature in payload.get("features", []):
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        region = str(props.get("LEVEL1_COD") or "").strip()
        if code and region:
            level1[code] = region
    if len(level1) < 300:
        raise RuntimeError("WGSRPD3 geometry/code drift")
    return level1


def hypergeom_upper_tail(N, K, n, x):
    if n < 0 or K < 0 or N < 0 or n > N or K > N:
        raise ValueError("invalid hypergeometric parameters")
    if x <= 0:
        return 1.0
    upper = min(K, n)
    if x > upper:
        return 0.0
    def logcomb(a, b):
        if b < 0 or b > a:
            return float("-inf")
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    logs = [
        logcomb(K, i) + logcomb(N - K, n - i) - logcomb(N, n)
        for i in range(x, upper + 1)
    ]
    peak = max(logs)
    return float(math.exp(peak) * sum(math.exp(value - peak) for value in logs))


def quantiles(values):
    x = np.asarray(values, dtype=float)
    return {
        "q025": float(np.quantile(x, 0.025)),
        "median": float(np.median(x)),
        "q975": float(np.quantile(x, 0.975)),
        "mean": float(np.mean(x)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-footprints-json", type=Path, required=True)
    parser.add_argument("--contemporary-footprints-json", type=Path, required=True)
    parser.add_argument("--unit-table-csv", type=Path, required=True)
    parser.add_argument("--level3-geojson", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=199999)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    native_payload = json.loads(args.native_footprints_json.read_text(encoding="utf-8"))
    contemporary_payload = json.loads(
        args.contemporary_footprints_json.read_text(encoding="utf-8")
    )
    native = {
        species: set(map(str, units))
        for species, units in native_payload["species"].items()
    }
    contemporary = {
        species: set(map(str, units))
        for species, units in contemporary_payload["species"].items()
    }
    level1 = load_level3(args.level3_geojson)
    universe = set(level1)

    observed_units = defaultdict(set)
    with args.unit_table_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if int(row["butterfly_observed"]):
                observed_units[row["species"]].add(row["wgsrpd3_code"])

    species_rows = []
    global_components = []
    region_components = []
    for species in sorted(contemporary):
        native_units = native.get(species, set()) & universe
        contemporary_units = contemporary[species] & universe
        added = contemporary_units - native_units
        outside = (observed_units.get(species, set()) & universe) - native_units
        overlap = added & outside
        if not outside:
            continue

        N = len(universe - native_units)
        K = len(outside)
        n = len(added)
        expected_global = 0.0 if N == 0 else K * n / N
        p_global = hypergeom_upper_tail(N, K, n, len(overlap))
        global_components.append((N, K, n))

        expected_region = 0.0
        region_terms = []
        for region in sorted(set(level1.values())):
            candidates = {
                code for code, value in level1.items()
                if value == region and code not in native_units
            }
            n_region = sum(level1[code] == region for code in added)
            K_region = sum(level1[code] == region for code in outside)
            if n_region == 0:
                continue
            N_region = len(candidates)
            if n_region > N_region:
                raise RuntimeError("added units exceed within-region candidate pool")
            expected_region += (
                0.0 if N_region == 0 else K_region * n_region / N_region
            )
            region_terms.append((N_region, K_region, n_region))
            region_components.append((N_region, K_region, n_region))

        species_rows.append(
            {
                "species": species,
                "native_units": len(native_units),
                "introduced_added_units": len(added),
                "outside_native_occurrence_units": len(outside),
                "observed_recovered_units": len(overlap),
                "fraction_recovered": len(overlap) / len(outside),
                "uniform_global_expected_overlap": expected_global,
                "uniform_global_p_high": p_global,
                "level1_composition_matched_expected_overlap": expected_region,
                "level1_component_count": len(region_terms),
            }
        )

    observed_total = sum(row["observed_recovered_units"] for row in species_rows)
    rng = np.random.default_rng(20260928)
    B = args.permutations
    global_total = np.zeros(B, dtype=np.int32)
    for N, K, n in global_components:
        global_total += rng.hypergeometric(K, N - K, n, size=B)
    region_total = np.zeros(B, dtype=np.int32)
    for N, K, n in region_components:
        region_total += rng.hypergeometric(K, N - K, n, size=B)

    global_p = float((1 + np.sum(global_total >= observed_total)) / (B + 1))
    region_p = float((1 + np.sum(region_total >= observed_total)) / (B + 1))

    payload = {
        "schema": "chocho_butterfly_occurrence_overlap_null_v0.1",
        "status": "POSTHOC_STRUCTURAL_AND_REGIONAL_OVERLAP_NULL",
        "species_with_outside_native_occurrences": len(species_rows),
        "observed_recovered_species_x_units": observed_total,
        "observed_fraction_recovered": observed_total / sum(
            row["outside_native_occurrence_units"] for row in species_rows
        ),
        "uniform_global_null": {
            **quantiles(global_total),
            "p_high": global_p,
            "definition": (
                "For each species, preserve native-envelope size, number of "
                "introduced-added units and number of outside-native butterfly "
                "occurrence units; place added units uniformly among non-native "
                "WGSRPD3 units."
            ),
        },
        "level1_composition_matched_null": {
            **quantiles(region_total),
            "p_high": region_p,
            "definition": (
                "For each species, preserve the number of introduced-added units "
                "within every WGSRPD level-1 region and sample them among non-native "
                "WGSRPD3 units in that same region."
            ),
        },
        "species_uniform_global_p_le_0_05": sum(
            row["uniform_global_p_high"] <= 0.05 for row in species_rows
        ),
        "species": species_rows,
        "claim_boundary": (
            "These nulls test whether occurrence recovery exceeds overlap expected "
            "from envelope enlargement and broad regional placement alone. They do "
            "not prove local larval use of introduced hosts or causal range expansion."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in payload.items() if k != "species"}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
