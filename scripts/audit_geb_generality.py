#!/usr/bin/env python3
"""Post-hoc GEB generality audit for the butterfly specialization paper.

This script does not alter any frozen ecological result. It asks whether the
already-observed host-portfolio architecture is concentrated in one butterfly
family or one native host-resource geographic context.

Geographic groups use WGSRPD Level 1 for the *native larval host-resource
footprint*. They are not butterfly occurrence realms.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean, median


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        rank = 0.5 * ((start + 1) + stop)
        for i in order[start:stop]:
            ranks[i] = rank
        start = stop
    return ranks


def pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3 or len(x) != len(y):
        return None
    mx, my = mean(x), mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    sx = sum(v * v for v in dx)
    sy = sum(v * v for v in dy)
    if sx <= 0 or sy <= 0:
        return None
    return sum(a * b for a, b in zip(dx, dy)) / math.sqrt(sx * sy)


def spearman(rows: list[dict[str, object]], xkey: str, ykey: str) -> dict[str, object]:
    pairs: list[tuple[float, float]] = []
    for row in rows:
        try:
            x = float(row[xkey])
            y = float(row[ykey])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(x) and math.isfinite(y):
            pairs.append((x, y))
    if len(pairs) < 3:
        return {"n": len(pairs), "rho": None}
    rx = average_ranks([p[0] for p in pairs])
    ry = average_ranks([p[1] for p in pairs])
    return {"n": len(pairs), "rho": pearson(rx, ry)}


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def build_family_map(rows: list[dict[str, str]]) -> dict[str, str]:
    out: dict[str, str] = {}
    for row in rows:
        family = str(row.get("Family", "")).strip()
        genus = str(row.get("Genus", "")).strip()
        species = str(row.get("Species", "")).strip()
        verbatim = str(row.get("verbatimSpecies", "")).strip()
        for name in (" ".join(x for x in (genus, species) if x), verbatim):
            if name and family:
                out.setdefault(name, family)
    return out


def family_summary(names: list[str], family_map: dict[str, str]) -> dict[str, object]:
    counts = Counter(family_map.get(name, "UNMAPPED") for name in names)
    if counts.get("UNMAPPED", 0):
        raise RuntimeError(f"unmapped butterfly names: {counts['UNMAPPED']}")
    n = len(names)
    shares = [count / n for count in counts.values()] if n else []
    eff = None if not shares else 1.0 / sum(p * p for p in shares)
    return {
        "n": n,
        "families": len(counts),
        "counts": dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))),
        "effective_family_number": eff,
    }


def collect_climate_species(obj: object, out: set[str] | None = None) -> set[str]:
    if out is None:
        out = set()
    if isinstance(obj, dict):
        if isinstance(obj.get("species"), str):
            out.add(str(obj["species"]))
        for value in obj.values():
            collect_climate_species(value, out)
    elif isinstance(obj, list):
        for value in obj:
            collect_climate_species(value, out)
    return out


def dominant_resource_continent(
    units: list[str],
    l3_to_l1: dict[str, int],
    l1_names: dict[int, str],
) -> tuple[str, int, float]:
    counts = Counter(l3_to_l1[unit] for unit in units)
    if not counts:
        raise RuntimeError("empty native resource footprint")
    best = max(counts.values())
    leaders = [code for code, value in counts.items() if value == best]
    name = "TIE" if len(leaders) != 1 else l1_names[leaders[0]]
    return name, len(counts), best / len(units)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--anthropogenic-csv", type=Path, required=True)
    ap.add_argument("--host-contribution-csv", type=Path, required=True)
    ap.add_argument("--independent-panel-json", type=Path, required=True)
    ap.add_argument("--climate-result-json", type=Path, required=True)
    ap.add_argument("--native-footprints-json", type=Path, required=True)
    ap.add_argument("--wgsrpd-level3-geojson", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    family_map = build_family_map(read_csv(args.leptraits_csv))
    anthropogenic = read_csv(args.anthropogenic_csv)
    mechanism = read_csv(args.host_contribution_csv)

    resource_eligible = anthropogenic
    adequate = [r for r in anthropogenic if truthy(r["host_taxonomy_lower_bound_adequate"])]
    expanded = [
        r for r in mechanism
        if truthy(r["host_taxonomy_lower_bound_adequate"])
        and float(r["introduced_added_units"]) > 0
    ]
    one_family = [r for r in expanded if float(r["host_family_count"]) == 1]

    panel = json.loads(args.independent_panel_json.read_text(encoding="utf-8"))
    panel_species = list(map(str, panel["panel"]["species"]))
    climate = json.loads(args.climate_result_json.read_text(encoding="utf-8"))
    climate_species = sorted(collect_climate_species(climate))

    for rows in (resource_eligible, adequate, expanded, one_family):
        for row in rows:
            row["Family"] = family_map[str(row["species"])]

    family_funnel = {
        "resource_eligible": family_summary([str(r["species"]) for r in resource_eligible], family_map),
        "host_taxonomy_adequate": family_summary([str(r["species"]) for r in adequate], family_map),
        "expanded_adequate": family_summary([str(r["species"]) for r in expanded], family_map),
        "one_family_expanded": family_summary([str(r["species"]) for r in one_family], family_map),
        "independent_panel": family_summary(panel_species, family_map),
        "climate_informative": family_summary(climate_species, family_map),
    }

    family_groups: dict[str, object] = {}
    for family in sorted({str(r["Family"]) for r in expanded}):
        subset = [r for r in expanded if r["Family"] == family]
        family_groups[family] = {
            "effective": spearman(subset, "host_family_count", "effective_contributor_number"),
            "dominance": spearman(subset, "host_family_count", "maximum_single_host_fractional_share"),
        }

    leave_one_family_out: dict[str, object] = {}
    for family in sorted({str(r["Family"]) for r in expanded}):
        subset = [r for r in expanded if r["Family"] != family]
        leave_one_family_out[family] = {
            "effective": spearman(subset, "host_family_count", "effective_contributor_number"),
            "dominance": spearman(subset, "host_family_count", "maximum_single_host_fractional_share"),
        }

    one_family_groups: dict[str, object] = {}
    for family in sorted({str(r["Family"]) for r in one_family}):
        subset = [r for r in one_family if r["Family"] == family]
        one_family_groups[family] = {
            "effective": spearman(subset, "resolved_host_species", "effective_contributor_number"),
            "dominance": spearman(subset, "resolved_host_species", "maximum_single_host_fractional_share"),
        }

    geo = json.loads(args.wgsrpd_level3_geojson.read_text(encoding="utf-8"))
    l3_to_l1 = {
        str(f["properties"]["LEVEL3_COD"]): int(f["properties"]["LEVEL1_COD"])
        for f in geo["features"]
    }
    l1_names = {
        1: "Europe",
        2: "Africa",
        3: "Asia-Temperate",
        4: "Asia-Tropical",
        5: "Australasia",
        6: "Pacific",
        7: "Northern America",
        8: "Southern America",
        9: "Antarctic",
    }

    footprint_payload = json.loads(args.native_footprints_json.read_text(encoding="utf-8"))
    footprints = footprint_payload["species"]
    presence = Counter()
    unit_totals = Counter()
    dominant: dict[str, dict[str, object]] = {}
    n_continents: list[int] = []

    for species, units in footprints.items():
        if any(unit not in l3_to_l1 for unit in units):
            missing = sorted(unit for unit in units if unit not in l3_to_l1)
            raise RuntimeError(f"{species}: unmapped WGSRPD3 units {missing}")
        counts = Counter(l3_to_l1[unit] for unit in units)
        for code, value in counts.items():
            presence[l1_names[code]] += 1
            unit_totals[l1_names[code]] += value
        name, n_l1, share = dominant_resource_continent(units, l3_to_l1, l1_names)
        dominant[species] = {"group": name, "n_level1": n_l1, "dominant_share": share}
        n_continents.append(n_l1)

    dominant_counts = Counter(v["group"] for v in dominant.values())
    for row in expanded:
        row["dominant_native_resource_continent"] = dominant[str(row["species"])]["group"]

    geography_groups: dict[str, object] = {}
    for group in sorted({str(r["dominant_native_resource_continent"]) for r in expanded}):
        subset = [r for r in expanded if r["dominant_native_resource_continent"] == group]
        geography_groups[group] = {
            "effective": spearman(subset, "host_family_count", "effective_contributor_number"),
            "dominance": spearman(subset, "host_family_count", "maximum_single_host_fractional_share"),
        }

    leave_one_geography_out: dict[str, object] = {}
    for group in sorted({str(r["dominant_native_resource_continent"]) for r in expanded}):
        subset = [r for r in expanded if r["dominant_native_resource_continent"] != group]
        leave_one_geography_out[group] = {
            "effective": spearman(subset, "host_family_count", "effective_contributor_number"),
            "dominance": spearman(subset, "host_family_count", "maximum_single_host_fractional_share"),
        }

    payload = {
        "schema": "chocho_geb_generality_audit_v0.1",
        "status": "EXPLORATORY_POSTHOC_GENERALITY_AUDIT",
        "claim_boundary": [
            "Family and native host-resource geographic replication only.",
            "WGSRPD Level-1 groups are not butterfly occurrence realms.",
            "No causal or prospective status is implied.",
        ],
        "family_funnel": family_funnel,
        "family_replication": {
            "overall": {
                "effective": spearman(expanded, "host_family_count", "effective_contributor_number"),
                "dominance": spearman(expanded, "host_family_count", "maximum_single_host_fractional_share"),
            },
            "by_family": family_groups,
            "leave_one_family_out": leave_one_family_out,
            "within_one_family": one_family_groups,
        },
        "native_resource_geography": {
            "species": len(footprints),
            "median_level1_units_per_species": median(n_continents),
            "mean_level1_units_per_species": mean(n_continents),
            "species_presence_by_level1": dict(sorted(presence.items())),
            "species_x_wgsrpd3_units_by_level1": dict(sorted(unit_totals.items())),
            "total_native_species_x_wgsrpd3_units": sum(unit_totals.values()),
            "dominant_native_resource_level1_counts": dict(
                sorted(dominant_counts.items(), key=lambda kv: (-kv[1], kv[0]))
            ),
        },
        "resource_geography_replication": {
            "definition": (
                "Species grouped by the WGSRPD Level-1 continent containing the "
                "largest share of its native larval host-resource WGSRPD3 footprint."
            ),
            "by_group": geography_groups,
            "leave_one_group_out": leave_one_geography_out,
        },
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
