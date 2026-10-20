#!/usr/bin/env python3
"""Save disjoint country/guard-region membership areas, never raw geometry.

Membership atoms let any downstream resource/climate mask compute its exact
geographic union by adding each selected atom once. Inputs are checked against
the existing source-pinned crosswalk before any derived output is written.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path

from shapely.ops import unary_union

from build_country_wgsrpd_crosswalk_20261020 import (
    NATURAL_EARTH_SHA256, dissolve_countries, equal_area_projector,
    polygonal_valid, regions_from_features, sha256,
)


RELATIVE_AREA_TOLERANCE = 1e-8
ABSOLUTE_AREA_TOLERANCE_M2 = .01


def atom_areas(country_geometry, region_geometries):
    """Return {sorted membership tuple: area} in the input projected units."""
    atoms = {}
    covered = unary_union([])
    for code, region in sorted(region_geometries.items()):
        clipped = polygonal_valid(country_geometry.intersection(region))
        if clipped.is_empty or clipped.area <= 0:
            continue
        next_atoms = {}
        for members, geometry in atoms.items():
            if not geometry.intersects(clipped):
                next_atoms[members] = geometry
                continue
            common = polygonal_valid(geometry.intersection(clipped))
            remaining = polygonal_valid(geometry.difference(clipped))
            if remaining.area > 0:
                next_atoms[members] = remaining
            if common.area > 0:
                next_atoms[tuple(sorted((*members, code)))] = common
        new_area = polygonal_valid(clipped.difference(covered))
        if new_area.area > 0:
            next_atoms[(code,)] = new_area
        atoms = next_atoms
        covered = polygonal_valid(unary_union([covered, clipped]))
    return {members: geometry.area for members, geometry in sorted(atoms.items())}


def require_same_area(actual, expected, label):
    if not math.isfinite(actual) or not math.isfinite(expected) or not math.isclose(
        actual, expected, rel_tol=RELATIVE_AREA_TOLERANCE,
        abs_tol=ABSOLUTE_AREA_TOLERANCE_M2,
    ):
        raise ValueError(f"Area identity failed for {label}: {actual} versus {expected} m2")


def build(tdwg_path, natural_earth_path, crosswalk_dir, output_dir):
    crosswalk_dir, output_dir = Path(crosswalk_dir), Path(output_dir)
    summary_path = crosswalk_dir / "country_wgsrpd_summary.json"
    links_path = crosswalk_dir / "country_wgsrpd_links.csv"
    original = json.loads(summary_path.read_text(encoding="utf-8"))
    if original.get("schema") != "country_wgsrpd_area_crosswalk_v1":
        raise ValueError("Unexpected source crosswalk schema")
    tdwg_sha, natural_sha = sha256(tdwg_path), sha256(natural_earth_path)
    if (natural_sha != NATURAL_EARTH_SHA256 or
            original["source"]["natural_earth_sha256"] != natural_sha or
            original["source"]["tdwg_sha256"] != tdwg_sha):
        raise ValueError("Input geometry hashes differ from the pinned original crosswalk")
    links_sha = sha256(links_path)
    if links_sha != original["output_links_sha256"]:
        raise ValueError("Original crosswalk link CSV hash changed")
    project = equal_area_projector()
    countries, _ = dissolve_countries(
        json.loads(Path(natural_earth_path).read_text())["features"], project)
    regions, _ = regions_from_features(
        json.loads(Path(tdwg_path).read_text())["features"], project)
    source_countries = {row["country"]: row for row in original["countries"]}
    if len(source_countries) != len(original["countries"]) or set(source_countries) != set(countries):
        raise ValueError("Country panel differs from the audited crosswalk")

    rows, audits = [], []
    for country, entry in sorted(countries.items()):
        base = source_countries[country]
        guard, material = set(base["guard_regions"]), set(base["material_regions"])
        if not material <= guard or not guard <= set(regions):
            raise ValueError(f"Invalid guard/material region sets for {country}")
        geometry = entry["geometry"]
        country_area = geometry.area
        require_same_area(country_area, base["country_area_km2"] * 1e6, country)
        region_geometries = {code: regions[code]["geometry"] for code in guard}
        atoms = atom_areas(geometry, region_geometries)
        guard_union_area = unary_union([
            geometry.intersection(g) for g in region_geometries.values()
        ]).area
        total = sum(atoms.values())
        require_same_area(total, guard_union_area, f"{country} atoms/guard union")
        require_same_area(total, base["guard_coverage"] * country_area,
                          f"{country} original guard coverage")
        material_area = sum(area for members, area in atoms.items() if material.intersection(members))
        require_same_area(material_area, base["material_coverage"] * country_area,
                          f"{country} original material coverage")
        for code, region_geometry in region_geometries.items():
            require_same_area(
                sum(area for members, area in atoms.items() if code in members),
                geometry.intersection(region_geometry).area,
                f"{country}/{code} membership",
            )
        for members, area in atoms.items():
            rows.append({"country": country, "wgsrpd3_members": "|".join(members),
                         "area_km2": area / 1e6})
        audits.append({
            "country": country, "country_area_km2": country_area / 1e6,
            "atom_count": len(atoms), "guard_union_area_km2": guard_union_area / 1e6,
            "atoms_area_km2": total / 1e6, "guard_coverage": total / country_area,
            "material_union_area_km2": material_area / 1e6,
            "material_coverage": material_area / country_area,
            "overlap_membership_area_km2": sum(
                area for members, area in atoms.items() if len(members) > 1) / 1e6,
            "overlap_excess_area_km2": sum(
                (len(members) - 1) * area for members, area in atoms.items()) / 1e6,
            "maximum_membership_size": max(map(len, atoms), default=0),
            "area_identity_checks": "PASS",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "country_area_atoms.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["country", "wgsrpd3_members", "area_km2"])
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "schema": "country_wgsrpd_disjoint_area_atoms_v1",
        "status": "PASS_SOURCE_HASHES_AND_DISJOINT_AREA_UNIONS",
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "builder_script_sha256": sha256(__file__),
        "crosswalk_summary_sha256": sha256(summary_path),
        "crosswalk_links_sha256": links_sha,
        "atoms_csv_sha256": sha256(csv_path),
        "source": {"tdwg_sha256": tdwg_sha, "natural_earth_sha256": natural_sha},
        "projection": original["projection"],
        "projection_longitude_rule": original["projection_longitude_rule"],
        "country_count": len(audits), "atom_rows": len(rows),
        "relative_area_tolerance": RELATIVE_AREA_TOLERANCE,
        "absolute_area_tolerance_square_metres": ABSOLUTE_AREA_TOLERANCE_M2,
        "partition_rule": "Split country-intersected guard regions into disjoint membership sets; combine disconnected pieces with identical membership by area. No uncovered land is added; total country land remains the denominator.",
        "selection_rule": "Count an atom once if at least one member matches the desired region mask. Restrict members to material regions where the analysis protocol requires material links.",
        "raw_geometry_serialized": False,
        "countries": audits,
    }
    (output_dir / "country_area_atoms_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdwg", required=True, type=Path)
    parser.add_argument("--natural-earth", required=True, type=Path)
    parser.add_argument("--crosswalk-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path,
                        help="Defaults to --crosswalk-dir; only derived areas and hashes are written")
    args = parser.parse_args()
    result = build(args.tdwg, args.natural_earth, args.crosswalk_dir,
                   args.output_dir or args.crosswalk_dir)
    print(json.dumps({key: result[key] for key in (
        "status", "country_count", "atom_rows", "atoms_csv_sha256",
        "crosswalk_summary_sha256", "crosswalk_links_sha256",
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
