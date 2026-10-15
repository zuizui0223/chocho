#!/usr/bin/env python3
"""Exploratory climate-conditioned plant-range null for butterfly resource overlap.

Reuses the frozen CHELSA V2.1 (1981–2010) WGSRPD3 climate summaries from an
independent butterfly climate workflow. DO NOT claim climatic causal control,
realized butterfly occupancy, or larval fitness. The fixed native host graph,
the exact butterfly-host network and plant introduction counts are preserved.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import median

from analyze_butterfly_shared_host_structural_null import (
    attempt_range_swaps, build_model, percentiles, read_csv,
    site_contribution_stats,
)

VARIABLES = ("bio1", "bio7", "bio12", "bio15")
MODES = {
    # Prespecified primary climate filter: 2x2 temp and precipitation classes
    "temperature_precipitation_2x2": ("bio1", "bio12", 2),
    # Sensitivity: finer temp/precip co-stratification
    "temperature_precipitation_3x3": ("bio1", "bio12", 3),
    # Conservative four-axis stratum; may not be evaluable
    "four_bioclim_variables_2x2x2x2": ("all", "all", 2),
}
MIN_CLIMATE_EDGE_COVERAGE = 0.85
MIN_VIABLE_EDGE_FRACTION = 0.40
MIN_CHANGED_FRACTION = 0.05
MIN_ACCEPTED_PER_BURN = 100
EXPECTED = {
    "added_units": 14553,
    "introduced_host_butterfly_region_incidences": 75346,
    "introduced_incidences_in_native_supported_cells": 46387,
    "introduced_incidences_in_newly_supported_cells": 28959,
    "redundant_introduced_host_incidences_in_new_cells": 14406,
    "novel_exact_shared_pair_regions": 79412,
    "single_units": 8574,
}


def read_unit_climate(path: Path) -> tuple[dict[str, tuple[float, ...]], dict[str, object]]:
    """Deduplicate earlier species-level climate output by WGSRPD3 region.

    The climate raster samples for a region must be identical regardless of
    focal butterfly. Missing and nonfinite values are never imputed.
    """
    data: dict[str, tuple[float, ...]] = {}
    duplicates = 0
    total_rows = 0
    missing_codes: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        required = {"wgsrpd3_code", "climate_sample_points", *VARIABLES}
        if not required <= set(rows.fieldnames or ()):
            raise ValueError("climate artifact schema drift")
        for row in rows:
            total_rows += 1
            code = str(row["wgsrpd3_code"]).strip()
            if not code or any(str(row[k]).strip() == "" for k in VARIABLES):
                if code:
                    missing_codes.add(code)
                continue
            try:
                if int(row["climate_sample_points"]) < 1:
                    missing_codes.add(code)
                    continue
                values = tuple(float(row[k]) for k in VARIABLES)
            except (ValueError, TypeError) as exc:
                raise ValueError(f"invalid climate data for {code}") from exc
            if not all(math.isfinite(x) for x in values):
                missing_codes.add(code)
                continue
            if code in data:
                duplicates += 1
                if any(abs(a-b) > 1e-7 for a, b in zip(values, data[code])):
                    raise ValueError(f"nonidentical climate values across butterflies in {code}")
            else:
                data[code] = values
    if not data:
        raise ValueError("no valid WGSRPD3 climate summaries")
    if set(data) & missing_codes:
        # Some rows with a missing vector for a region are acceptable only
        # when all complete rows for that region match. Report these explicitly.
        missing_codes -= set(data)
    return data, {
        "source_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_rows": total_rows,
        "source_unique_valid_regions": len(data),
        "duplicate_rows_validated": duplicates,
        "regions_only_missing_climate": len(missing_codes),
        "data_source": (
            "CHELSA V2.1 1981-2010 median climate samples within WGSRPD3,"
            " frozen independent climate crossfit result,"
            " TTF run 36270581743 artifact 10919985332"
        ),
    }


def quantile_cutpoints(values: list[float], n: int) -> tuple[float, ...]:
    if n < 2 or not values:
        raise ValueError("invalid quantile inputs")
    ordered = sorted(values)
    cuts = []
    for k in range(1, n):
        t = k * (len(ordered)-1) / n
        lo = math.floor(t)
        hi = math.ceil(t)
        x = ordered[lo] if lo == hi else (
            ordered[lo] * (hi-t) + ordered[hi] * (t-lo)
        )
        cuts.append(x)
    return tuple(cuts)


def climate_strata(
    climate: dict[str, tuple[float, ...]],
    region_codes: list[str],
    level1: dict[int, str],
    mode: str,
) -> tuple[dict[int, tuple[int, ...]], dict[str, object]]:
    if mode not in MODES:
        raise ValueError(f"unknown mode: {mode}")
    keys = ("bio1", "bio12") if mode != "four_bioclim_variables_2x2x2x2" else VARIABLES
    n = MODES[mode][2]
    by_area: dict[str, list[tuple[int, tuple[float, ...]]]] = defaultdict(list)
    for i, code in enumerate(region_codes):
        if code in climate:
            by_area[level1[i]].append((i, climate[code]))
    cuts: dict[str, tuple[tuple[float, ...], ...]] = {}
    classes: dict[int, tuple[int, ...]] = {}
    for area, records in sorted(by_area.items()):
        thresholds = tuple(
            quantile_cutpoints(
                [vals[VARIABLES.index(key)] for _, vals in records], n
            )
            for key in keys
        )
        cuts[area] = thresholds
        for region, vec in records:
            classes[region] = tuple(
                bisect.bisect_left(thresholds[j], vec[VARIABLES.index(key)])
                for j, key in enumerate(keys)
            )
    return classes, {
        "mode": mode,
        "variables": list(keys),
        "n_quantile_classes_per_variable_within_level1": n,
        "cutpoints_by_level1": {
            a: {k: list(cs) for k, cs in zip(keys, values)}
            for a, values in cuts.items()
        },
        "climate_regions": len(classes),
        "bin_assignment_rule": (
            "within-Level1 equal-region-quantile cutpoints;"
            " exact threshold ties assigned to lower class"
        ),
    }


def grouped_introduced_edges(
    introduced: dict[str, set[int]],
    host_family: dict[str, str],
    level1: dict[int, str],
    climatic_classes: dict[int, tuple[int, ...]],
) -> tuple[list[list[tuple[str, int]]], list[int], dict[str, float | int]]:
    """Keep uncovered and non-swap-viable edges unchanged in null worlds."""
    by_group: dict[tuple[str, str, tuple[int, ...]], list[tuple[str, int]]] = defaultdict(list)
    all_edges = covered = 0
    for host, units in sorted(introduced.items()):
        family = host_family[host]
        for region in sorted(units):
            all_edges += 1
            if region not in climatic_classes:
                continue
            covered += 1
            by_group[(family, level1[region], climatic_classes[region])].append(
                (host, region)
            )

    viable = []
    stuck_groups = 0
    for group_key in sorted(by_group):
        edges = by_group[group_key]
        if len(edges) >= 2 and len({h for h, _ in edges}) >= 2 \
           and len({u for _, u in edges}) >= 2:
            viable.append(edges[:])
        else:
            stuck_groups += 1
    cumulative = []
    n = 0
    for group in viable:
        n += len(group)
        cumulative.append(n)
    summary = {
        "total_introduced_plant_region_edges": all_edges,
        "climate_covered_edges": covered,
        "climate_covered_fraction": covered/all_edges if all_edges else 0,
        "viable_edges": n,
        "viable_fraction_of_total": n/all_edges if all_edges else 0,
        "fixed_edges": all_edges-n,
        "viable_group_count": len(viable),
        "nonswappable_group_count": stuck_groups,
    }
    return viable, cumulative, summary


def evaluate(model: dict[str, object], introduced: dict[str, set[int]]) -> dict[str, float | int]:
    masks = {
        h: sum(1 << u for u in cells)
        for h, cells in introduced.items()
    }
    return site_contribution_stats(
        model["species_hosts"],
        model["native_bits"],
        masks,
        model["shared_pairs"],
        model["native_species"],
    )


def candidate_chain(
    model: dict[str, object],
    climate_classes: dict[int, tuple[int, ...]],
    reps: int,
    seed: int,
    burn_factor: int,
    spacing_factor: int,
) -> dict[str, object]:
    if reps < 99 or burn_factor < 1 or spacing_factor < 1:
        raise ValueError("insufficient requested shuffle replication")
    # Crucial: start EACH chain from the same original observed plant graph.
    introduced = {h: set(s) for h, s in model["introduced_sets"].items()}
    groups, cumulative, group_diag = grouped_introduced_edges(
        introduced, model["host_families"], model["level1"], climate_classes
    )
    observed = evaluate(model, introduced)
    for k, expected in EXPECTED.items():
        if observed[k] != expected:
            raise ValueError(f"frozen baseline drift for {k}: {observed[k]} != {expected}")
    if (
        group_diag["climate_covered_fraction"] < MIN_CLIMATE_EDGE_COVERAGE
        or group_diag["viable_fraction_of_total"] < MIN_VIABLE_EDGE_FRACTION
    ):
        return {
            "status": "NOT_EVALUABLE_CLIMATE_SUPPORT_OR_SWAP_GEOMETRY",
            "seed": seed,
            "observed": observed,
            "groups": group_diag,
            "thresholds": {
                "min_covered_fraction": MIN_CLIMATE_EDGE_COVERAGE,
                "min_viable_fraction": MIN_VIABLE_EDGE_FRACTION,
            },
        }

    rng = random.Random(seed)
    E = int(group_diag["total_introduced_plant_region_edges"])
    orig = {h: frozenset(v) for h, v in introduced.items()}
    burn = burn_factor*E
    burn_accepted = attempt_range_swaps(
        groups, cumulative, model["native_sets"], introduced, burn, rng
    )
    # Each swapped incidence is a pair of removed / added membership changes.
    moved_fraction = sum(
        len(orig[h] ^ introduced[h]) for h in introduced
    ) / (2*E)
    if moved_fraction < MIN_CHANGED_FRACTION or burn_accepted < MIN_ACCEPTED_PER_BURN:
        return {
            "status": "NOT_EVALUABLE_INADEQUATE_CLIMATE_NULL_MIXING",
            "seed": seed,
            "observed": observed,
            "groups": group_diag,
            "burn_acceptances": burn_accepted,
            "fraction_initial_introduced_edges_moved_after_burn": moved_fraction,
        }

    stats: dict[str, list[float]] = {k: [] for k in observed}
    accepted = []
    for _ in range(reps):
        accepted.append(attempt_range_swaps(
            groups, cumulative, model["native_sets"], introduced,
            spacing_factor*E, rng
        ))
        x = evaluate(model, introduced)
        for k, val in x.items():
            stats[k].append(float(val))
        for k in ("introduced_host_butterfly_region_incidences",):
            if x[k] != observed[k]:
                raise ValueError(f"fixed host × butterfly × region budget drift: {k}")
        if x["introduced_incidences_in_native_supported_cells"] \
           + x["introduced_incidences_in_newly_supported_cells"] \
           != x["introduced_host_butterfly_region_incidences"]:
            raise ValueError("native/new absorption conservation violated")
        if x["introduced_incidences_in_newly_supported_cells"] \
           - x["redundant_introduced_host_incidences_in_new_cells"] \
           != x["added_units"]:
            raise ValueError("unique addition accounting violated")

    outcomes = {}
    for key, null_values in stats.items():
        low, center, high = percentiles(null_values)
        obs = float(observed[key])
        mean_null = sum(null_values) / len(null_values)
        numerator = sum((a-mean_null)*(b-mean_null)
                        for a,b in zip(null_values[:-1],null_values[1:]))
        denominator = sum((a-mean_null)**2 for a in null_values)
        lag1 = numerator / denominator if denominator else 0.0
        outcomes[key] = {
            "observed": obs,
            "null_q025": low,
            "null_median": center,
            "null_q975": high,
            "observed_minus_null_median": obs-center,
            "nominal_one_sided_p_above": (
                1+sum(v >= obs for v in null_values)
            )/(reps+1),
            "nominal_one_sided_p_below": (
                1+sum(v <= obs for v in null_values)
            )/(reps+1),
            "lag1_null_draw_correlation": lag1,
        }
    return {
        "status": "EXPLORATORY_CLIMATE_CONDITIONED_NULL_EVALUABLE",
        "seed": seed,
        "null_draws": reps,
        "groups": group_diag,
        "burn_proposals": burn,
        "burn_accepted": burn_accepted,
        "changed_fraction_after_burn": moved_fraction,
        "spacing_proposals_per_draw": spacing_factor*E,
        "accepted_swaps_min": min(accepted),
        "accepted_swaps_median": median(accepted),
        "observed": observed,
        "outcomes": outcomes,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    for k in (
        "descriptors_csv", "insect_host_csv", "native_distribution_csv",
        "contemporary_distribution_csv", "level1_json", "unit_climate_csv",
        "output_json",
    ):
        p.add_argument("--"+k.replace("_", "-"), type=Path, required=True)
    p.add_argument("--draws", type=int, default=299)
    p.add_argument("--seeds", nargs="+", type=int, default=[20261015, 20261016])
    p.add_argument("--burn-factor", type=int, default=5)
    p.add_argument("--spacing-factor", type=int, default=2)
    args = p.parse_args()

    native_rows = read_csv(args.native_distribution_csv)
    contemporary_rows = read_csv(args.contemporary_distribution_csv)
    interaction_rows = read_csv(args.insect_host_csv)
    descriptors = read_csv(args.descriptors_csv)
    lookup = json.loads(args.level1_json.read_text(encoding="utf-8"))["units"]
    model = build_model(
        descriptors, interaction_rows, native_rows, contemporary_rows, lookup
    )
    focal_hosts = set(model["host_families"])
    region_codes = sorted({
        str(r["area_code_l3"])
        for r in native_rows+contemporary_rows
        if r["accepted_plant_name_id"] in focal_hosts
    })
    if set(region_codes) - set(lookup):
        raise ValueError("missing region in frozen WGSRPD hierarchy")
    if len(region_codes) != len(model["level1"]):
        raise ValueError("reconstructed region-index order mismatch")
    if any(lookup[c] != model["level1"][i] for i, c in enumerate(region_codes)):
        raise ValueError("hierarchy not aligned with original resource reconstruction")
    climate, prov = read_unit_climate(args.unit_climate_csv)
    observed = evaluate(model, model["introduced_sets"])
    if any(observed[k] != value for k, value in EXPECTED.items()):
        raise ValueError("source reconstruction differs from frozen null")
    output = {
        "schema": "chocho_climate_conditioned_plant_range_null_v0.1",
        "status": "CONDITIONAL_POSTHOC_CLIMATE_SENSITIVITY",
        "source_climate": prov,
        "study": {
            "focal_butterflies": len(model["species"]),
            "source_climate_regions": len(climate),
            "resource_regions": len(region_codes),
            "variables": list(VARIABLES),
            "source_climate_age": "1981–2010 climatology, not introduction dates",
            "source_geo": "WGSRPD3 are geographically coarse",
        },
        "baseline_observed": observed,
        "modes": {},
        "scientific_boundaries": [
            "Post-hoc climate-stratified sensitivity; both science questions and original null results already known.",
            "Plant-family × Level1 × coarse climate-bin swap preserves exact butterfly-host links and introduced counts per plant within each climate-bin; not species-level resource addition margins.",
            "Climatic classes are summaries from independent occurrence-panel host geography, not site-scale habitat or true plant climatic niches.",
            "Regions without a climate signature and nonviable strata are fixed rather than artificially imputed; results conditional on measured coverage.",
            "Climatic similarities within a bin are not controlled continuously, nor geography within continents, elevation, invasiveness, introduction age or recording intensity.",
            "Within-chain Monte Carlo randomizations are serially related, p values nominal and not post-hoc/multiplicity-corrected.",
            "No butterfly colonization, successful development or competition tested.",
        ],
    }
    for mode in MODES:
        classes, design = climate_strata(climate, region_codes, model["level1"], mode)
        chains = [
            candidate_chain(
                model, classes, args.draws, seed,
                args.burn_factor, args.spacing_factor
            )
            for seed in args.seeds
        ]
        output["modes"][mode] = {
            "climatic_blocks": design,
            "chains": chains,
            "status": (
                "EVALUABLE"
                if all(x["status"] == "EXPLORATORY_CLIMATE_CONDITIONED_NULL_EVALUABLE" for x in chains)
                else "NOT_EVALUABLE"
            ),
        }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "schema":output["schema"],
        "climate_source":prov,
        "study":output["study"],
        "baseline": {
            k:observed[k] for k in (
                "introduced_host_butterfly_region_incidences",
                "introduced_incidences_in_native_supported_cells",
                "added_units",
            )
        },
        "modes": {
            name:{
                "status": value["status"],
                "chains": [
                    {
                        "status": c["status"],
                        "seed": c["seed"],
                        "groups": c["groups"],
                        "burn_accepted": c.get("burn_accepted"),
                        "changed_fraction_after_burn": c.get("changed_fraction_after_burn"),
                        "absorption": c.get("outcomes",{}).get("introduced_incidences_in_native_supported_cells"),
                        "new_added_cells": c.get("outcomes",{}).get("added_units"),
                    } for c in value["chains"]
                ],
            } for name,value in output["modes"].items()
        },
        "output_json":str(args.output_json),
    },indent=2))


if __name__ == "__main__":
    main()
