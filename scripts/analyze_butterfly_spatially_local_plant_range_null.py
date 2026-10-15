#!/usr/bin/env python3
"""Geographically local climate-matched host-range null (exploratory only).

This spatial rearrangement is conditional on a frozen 239-butterfly HOSTS/WCVP
resource inventory. It never measures larval performance or colonization.
It preserves degrees and uses an original-footprint tether to prevent chains
from drifting outward through repeated locally accepted exchanges.
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

from shapely.geometry import shape

from analyze_butterfly_shared_host_structural_null import (
    build_model, percentiles, read_csv,
)
from analyze_butterfly_climate_matched_plant_range_null import (
    EXPECTED, climate_strata, evaluate, grouped_introduced_edges,
    read_unit_climate,
)

CALIPERS_KM = (1500, 750, 3000)
MAIN_METRIC = "introduced_incidences_in_native_supported_cells"
SECONDARY_METRIC = "added_units"
EARTH_KM = 6371.0088


def haversine_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    """Shortest great-circle distance between (lat, lon) coordinates."""
    lat_a, lon_a = (math.radians(x) for x in a)
    lat_b, lon_b = (math.radians(x) for x in b)
    x = math.sin((lat_b-lat_a)/2)**2 + math.cos(lat_a)*math.cos(lat_b) * (
        math.sin((lon_b-lon_a)/2)**2
    )
    return 2*EARTH_KM*math.asin(min(1., math.sqrt(max(0., x))))


def load_geography(
    path: Path, ordered_codes: list[str], level1: dict[int, str]
) -> tuple[dict[int, tuple[float, float]], dict[str, object]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    raw: dict[str, tuple[float, float]] = {}
    categories = {}
    for feature in data.get("features", []):
        prop = feature.get("properties") or {}
        code = str(prop.get("LEVEL3_COD") or "").strip()
        area = str(prop.get("LEVEL1_COD") or "").strip()
        geom = feature.get("geometry")
        if not code or not geom:
            continue
        if code in raw:
            raise ValueError(f"duplicate WGSRPD Level3 code: {code}")
        p = shape(geom).representative_point()
        lat, lon = float(p.y), float(p.x)
        if not (math.isfinite(lat) and math.isfinite(lon)
                and -90 <= lat <= 90 and -180 <= lon <= 180):
            raise ValueError(f"invalid representative-point geometry for {code}")
        raw[code] = (lat, lon)
        categories[code] = area

    missing = sorted(set(ordered_codes) - set(raw))
    if missing:
        raise ValueError("focal resource regions missing Level3 centroid: " + ",".join(missing[:10]))
    coords = {}
    for i, code in enumerate(ordered_codes):
        if categories[code] != level1[i]:
            raise ValueError(f"Level1 hierarchy mismatch for {code}")
        coords[i] = raw[code]
    return coords, {
        "source_geometry_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_geometry": "tdwg/wgsrpd geojson/level3.geojson at fixed commit 52da7828aba9d461dd133c27b3bd7a4407161f54",
        "mapped_focal_regions": len(coords),
        "available_geometry_codes": len(raw),
        "representative_point_note": "Within-polygon representative point, not population, geographic, or climatic centroid",
    }


def neighbour_index(
    coords: dict[int, tuple[float, float]], radius_km: int
) -> dict[int, frozenset[int]]:
    near = {}
    for i, a in coords.items():
        near[i] = frozenset(
            j for j, b in coords.items() if haversine_km(a, b) <= radius_km+1e-9
        )
        if i not in near[i]:
            raise ValueError("self geography missing")
    for i, neighbors in near.items():
        for j in neighbors:
            if i not in near[j]:
                raise ValueError("asymmetric distance graph")
    return near


def origin_tether(
    native_intro: dict[str, set[int]], near: dict[int, frozenset[int]]
) -> dict[str, frozenset[int]]:
    return {
        h: frozenset().union(*(near[u] for u in units))
        if units else frozenset()
        for h, units in native_intro.items()
    }


def admissible(
    left: tuple[str, int],
    right: tuple[str, int],
    introduced: dict[str, set[int]],
    native: dict[str, set[int]],
    neighbours: dict[int, frozenset[int]],
    tether: dict[str, frozenset[int]],
) -> bool:
    """Reversible cross-region two-edge swap with geographic original tether."""
    a, x = left
    b, y = right
    return (
        a != b and x != y
        and y in neighbours[x]
        and y in tether[a] and x in tether[b]
        and y not in introduced[a] and x not in introduced[b]
        and y not in native.get(a, set()) and x not in native.get(b, set())
    )


def local_swap(
    groups: list[list[tuple[str, int]]], cumulative: list[int],
    introduced: dict[str, set[int]], native: dict[str, set[int]],
    neighbours: dict[int, frozenset[int]],
    tether: dict[str, frozenset[int]],
    attempts: int, rng: random.Random,
) -> int:
    if not cumulative:
        return 0
    accepted = 0
    for _ in range(attempts):
        k = rng.randrange(cumulative[-1])
        gidx = bisect.bisect_right(cumulative, k)
        before = cumulative[gidx-1] if gidx else 0
        i = k-before
        group = groups[gidx]
        j = rng.randrange(len(group))
        if i == j:
            continue
        ai = group[i]
        bj = group[j]
        if not admissible(ai, bj, introduced, native, neighbours, tether):
            continue
        a, x = ai
        b, y = bj
        introduced[a].remove(x)
        introduced[a].add(y)
        introduced[b].remove(y)
        introduced[b].add(x)
        group[i] = (a, y)
        group[j] = (b, x)
        accepted += 1
    return accepted


def count_mobile_edges(
    groups: list[list[tuple[str, int]]],
    introduced: dict[str, set[int]], native: dict[str, set[int]],
    neighbours: dict[int, frozenset[int]],
    tether: dict[str, frozenset[int]],
) -> int:
    """Exact one-step mobility at the ORIGINAL arrangement, not a fit measure."""
    mobile = 0
    for group in groups:
        for i, a in enumerate(group):
            for j, b in enumerate(group):
                if i == j:
                    continue
                if admissible(a, b, introduced, native, neighbours, tether):
                    mobile += 1
                    break
    return mobile


def lag1(x: list[float]) -> float | None:
    n = len(x)
    if n < 3:
        return None
    m = sum(x)/n
    var = sum((v-m)**2 for v in x)
    if var == 0:
        return None
    return sum((x[k]-m)*(x[k+1]-m) for k in range(n-1))/var


def run_caliper(
    model: dict[str, object], categories: dict[int, tuple[int, ...]],
    coords: dict[int, tuple[float, float]],
    caliper_km: int, seeds: list[int], draws: int,
    burn_factor: int, spacing_factor: int,
) -> dict[str, object]:
    source = model["introduced_sets"]
    n_edges = sum(len(x) for x in source.values())
    near = neighbour_index(coords, caliper_km)
    tether = origin_tether(source, near)
    groups0, _, baseline = grouped_introduced_edges(
        source, model["host_families"], model["level1"], categories
    )
    init_mobile = count_mobile_edges(
        groups0, source, model["native_sets"], near, tether
    )
    geometry = {
        **baseline,
        "radius_km": caliper_km,
        "original_tethered_mobile_edges": init_mobile,
        "original_tethered_mobile_fraction": init_mobile/n_edges,
        "median_neighbour_regions": median([len(v) for v in near.values()]),
        "one_step_validity_definition": (
            "At least one reversible two-edge swap within family × Level1 "
            "× climate-bin, local to both involved region points and tethered "
            "to each host's originally introduced locations."
        ),
        "input_regions": len(coords),
    }
    gate = baseline["climate_covered_fraction"] >= 1.0 and (
        baseline["viable_fraction_of_total"] >= 0.4
    ) and init_mobile/n_edges >= 0.05
    if not gate:
        return {
            "status": "NOT_EVALUABLE_SPATIAL_GEOMETRY",
            "caliper_km": caliper_km,
            "geometry": geometry,
            "chains": [],
        }

    baseline_obs = evaluate(model, source)
    if any(baseline_obs[k] != v for k, v in EXPECTED.items()):
        raise ValueError("invariant observed host-resource numbers changed")

    all_chains = []
    for seed in seeds:
        # Every chain is independent: reset graph and group representation.
        introduced = {h: set(v) for h, v in source.items()}
        groups, cumulative, rep = grouped_introduced_edges(
            introduced, model["host_families"], model["level1"], categories
        )
        if rep != baseline:
            raise ValueError("original climate group reconstruction changed")
        rng = random.Random(seed)
        proposals = burn_factor*n_edges
        accepted = local_swap(
            groups, cumulative, introduced, model["native_sets"], near,
            tether, proposals, rng
        )
        changed = sum(
            len(introduced[h] ^ source[h]) for h in introduced
        )/(2*n_edges)
        if accepted < 100 or changed < 0.02:
            all_chains.append({
                "status": "NOT_EVALUABLE_SPATIAL_MIXING",
                "seed": seed,
                "burn_accepted": accepted,
                "fraction_original_edges_changed": changed,
            })
            continue
        stats: dict[str, list[float]] = {
            name: [] for name in baseline_obs
        }
        accepted_per_world = []
        for _ in range(draws):
            n = local_swap(
                groups, cumulative, introduced, model["native_sets"], near,
                tether, spacing_factor*n_edges, rng
            )
            accepted_per_world.append(n)
            result = evaluate(model, introduced)
            for key in stats:
                stats[key].append(float(result[key]))
            if result["introduced_host_butterfly_region_incidences"] != baseline_obs[
                "introduced_host_butterfly_region_incidences"
            ]:
                raise ValueError("fixed introduced host exposure budget drift")
            if result["introduced_incidences_in_native_supported_cells"] + (
                result["introduced_incidences_in_newly_supported_cells"]
            ) != result["introduced_host_butterfly_region_incidences"]:
                raise ValueError("native absorption conservation drift")
            if result["introduced_incidences_in_newly_supported_cells"] - (
                result["redundant_introduced_host_incidences_in_new_cells"]
            ) != result["added_units"]:
                raise ValueError("unique opportunity conservation drift")
        outcomes = {}
        for metric in stats:
            xs = stats[metric]
            q025, q50, q975 = percentiles(xs)
            obs = float(baseline_obs[metric])
            outcomes[metric] = {
                "observed": obs,
                "null_q025": q025,
                "null_median": q50,
                "null_q975": q975,
                "observed_minus_null_median": obs-q50,
                "nominal_p_ge": (1+sum(v >= obs for v in xs))/(draws+1),
                "nominal_p_le": (1+sum(v <= obs for v in xs))/(draws+1),
                "lag1_correlation": lag1(xs),
            }
        all_chains.append({
            "status": "EXPLORATORY_SPATIAL_NULL_EVALUABLE",
            "seed": seed,
            "draws": draws,
            "burn_proposals": proposals,
            "burn_accepted": accepted,
            "fraction_original_edges_changed": changed,
            "spacing_proposals_per_draw": spacing_factor*n_edges,
            "accepted_swaps_min": min(accepted_per_world),
            "accepted_swaps_median": median(accepted_per_world),
            "outcomes": outcomes,
            # Keep full primary/secondary draws for later chain mixing audit.
            "primary_null_draws": stats[MAIN_METRIC],
            "secondary_null_draws": stats[SECONDARY_METRIC],
        })
    status = (
        "EVALUABLE" if len(all_chains) == len(seeds)
        and all(x["status"]=="EXPLORATORY_SPATIAL_NULL_EVALUABLE"
                for x in all_chains)
        else "NOT_EVALUABLE"
    )
    return {
        "status": status,
        "caliper_km": caliper_km,
        "geometry": geometry,
        "chains": all_chains,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    for p in (
        "descriptors_csv", "insect_host_csv",
        "native_distribution_csv", "contemporary_distribution_csv",
        "level1_json", "unit_climate_csv", "level3_geojson", "output_json",
    ):
        ap.add_argument("--"+p.replace("_","-"), type=Path, required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[20261015, 20261016])
    ap.add_argument("--draws", type=int, default=199)
    ap.add_argument("--burn-factor", type=int, default=5)
    ap.add_argument("--spacing-factor", type=int, default=1)
    opts = ap.parse_args()
    lookup = json.loads(opts.level1_json.read_text(encoding="utf-8"))["units"]
    nrows = read_csv(opts.native_distribution_csv)
    crows = read_csv(opts.contemporary_distribution_csv)
    hrows = read_csv(opts.insect_host_csv)
    model = build_model(
        read_csv(opts.descriptors_csv), hrows, nrows, crows, lookup
    )
    focal_hosts = set(model["host_families"])
    codes = sorted({
        r["area_code_l3"] for r in nrows+crows
        if r["accepted_plant_name_id"] in focal_hosts
    })
    if len(codes) != len(model["level1"]):
        raise ValueError("resource region index mismatch")
    if any(lookup[u] != model["level1"][i] for i,u in enumerate(codes)):
        raise ValueError("hierarchy mismatch before geography")
    positions, geo_source = load_geography(
        opts.level3_geojson, codes, model["level1"]
    )
    climate, climate_source = read_unit_climate(opts.unit_climate_csv)
    if not set(codes) <= set(climate):
        raise ValueError("climate data missing any focal geographic unit")
    strata, strata_source = climate_strata(
        climate, codes, model["level1"],
        "temperature_precipitation_2x2"
    )
    if len(strata) != len(codes):
        raise ValueError("spatial climate-binning coverage incomplete")
    observed = evaluate(model, model["introduced_sets"])
    if any(observed[k] != v for k, v in EXPECTED.items()):
        raise ValueError("input reconstruction drift")
    output = {
        "schema": "chocho_spatially_local_climate_range_null_v0.1",
        "status": "EXPLORATORY_POSTHOC_GEOGRAPHIC_SPATIAL_NULL",
        "study": {
            "focal_butterflies": len(model["species"]),
            "resource_regions": len(codes),
            "introduced_plant_region_edges": sum(
                len(v) for v in model["introduced_sets"].values()
            ),
            "H": observed["introduced_host_butterfly_region_incidences"],
        },
        "source_geo": geo_source,
        "source_climate": climate_source,
        "climate_strata": strata_source,
        "observed": observed,
        "calipers": {},
        "scientific_boundary": (
            "Coarse polygon representative points proxy proximity; the "
            "null restricts plant-family×continent×climate-bin exchanges to "
            "local pairs also tethered to original introduced host geography. "
            "This is post-hoc structural redistribution, not realized plant "
            "movement or butterfly competition; chains can be non-ergodic."
        ),
    }
    for caliper in CALIPERS_KM:
        output["calipers"][str(caliper)] = run_caliper(
            model, strata, positions, caliper, opts.seeds, opts.draws,
            opts.burn_factor, opts.spacing_factor
        )
    opts.output_json.parent.mkdir(parents=True, exist_ok=True)
    opts.output_json.write_text(
        json.dumps(output, indent=2, ensure_ascii=False)+"\n", encoding="utf-8"
    )
    summary = {
        "schema":output["schema"],
        "source_geo":output["source_geo"],
        "source_climate_regions":output["source_climate"]["source_unique_valid_regions"],
        "study":output["study"],
        "calipers":{},
    }
    for k,val in output["calipers"].items():
        summary["calipers"][k] = {
            "status":val["status"],
            "geometry":val["geometry"],
            "chains":[
                {
                    "seed":c["seed"],
                    "status":c["status"],
                    "burn_accepted":c.get("burn_accepted"),
                    "fraction_original_edges_changed":c.get("fraction_original_edges_changed"),
                    "absorption":c.get("outcomes",{}).get(MAIN_METRIC),
                    "unique_added":c.get("outcomes",{}).get(SECONDARY_METRIC),
                } for c in val["chains"]
            ],
        }
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
