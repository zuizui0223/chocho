#!/usr/bin/env python3
"""Derive an audited country/territory-to-WGSRPD3 area crosswalk.

Natural Earth map units are dissolved by their two-letter ISO_A2_EH code;
overseas territories therefore remain separate when they have a separate code.
No occurrence outcome is read. Exact source polygons remain outside the repo.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from pyproj import Transformer
from shapely import make_valid
from shapely.geometry import shape
from shapely.ops import transform, unary_union
from shapely.strtree import STRtree


NATURAL_EARTH_COMMIT = "117488dc884bad03366ff727eca013e434615127"
NATURAL_EARTH_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    + NATURAL_EARTH_COMMIT + "/geojson/ne_50m_admin_0_map_units.geojson"
)
NATURAL_EARTH_SHA256 = "b8d421aca6e9e08e8cdf09cc26af111cc3e0deba4fe915611d58ade71e8a4db0"
MATERIAL_MIN_SMALLER_FRACTION = .05
GUARD_MIN_SMALLER_FRACTION = .005
COUNTRY_MAPPING_FLAGS = {
    "NO": ["Natural Earth assigns Jan Mayen to NO, whereas ISO SJ includes Svalbard and Jan Mayen; no silent recoding"],
    "SJ": ["Natural Earth SJ geometry contains Svalbard but Jan Mayen is coded NO; no silent recoding"],
    "SO": ["Uncoded Somaliland map unit makes the coded Somalia footprint incomplete"],
    "CY": ["Uncoded Northern Cyprus map unit makes the coded Cyprus footprint incomplete"],
    "IN": ["Uncoded Siachen Glacier and disputed boundary interpretation; conservative country exclusion"],
    "PK": ["Uncoded Siachen Glacier and disputed boundary interpretation; conservative country exclusion"],
    "XK": ["XK is an extension rather than an officially assigned ISO 3166-1 alpha-2 code"],
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def is_two_letter(value):
    return isinstance(value, str) and re.fullmatch(r"[A-Z]{2}", value) is not None


def country_code(properties):
    """Reject conflicting real codes; use explicit EH territory mapping."""
    raw = properties.get("ISO_A2")
    extended = properties.get("ISO_A2_EH")
    if is_two_letter(raw) and is_two_letter(extended) and raw != extended:
        raise ValueError(f"Conflicting two-letter ISO fields: {raw} versus {extended}")
    if is_two_letter(extended):
        return extended
    if is_two_letter(raw):
        return raw
    return None


def polygonal_valid(geometry):
    geometry = make_valid(geometry)
    if geometry.geom_type in ("Polygon", "MultiPolygon"):
        return geometry
    if geometry.geom_type == "GeometryCollection":
        pieces = [polygonal_valid(g) for g in geometry.geoms
                  if g.geom_type in ("Polygon", "MultiPolygon", "GeometryCollection")]
        return unary_union(pieces)
    return unary_union([])


def equal_area_projector():
    # The archived MAG polygon contains -180.000026 degrees. PROJ's default
    # wrapping would send this nearly-antimeridian vertex to +179.999974,
    # creating an artificial polygon strip across the Arctic. Preserve the
    # source longitude continuity instead of wrapping individual vertices.
    transformer = Transformer.from_crs(
        "EPSG:4326", "EPSG:6933", always_xy=True, force_over=True)
    return lambda geometry: transform(transformer.transform, geometry)


def dissolve_countries(features, project):
    grouped = defaultdict(list)
    names = defaultdict(set)
    fields = defaultdict(list)
    audit = {"excluded_map_units": [], "repaired_map_units": [],
             "iso_field_fallback_map_units": []}
    for feature in features:
        p = feature["properties"]
        label = p.get("NAME", p.get("GEOUNIT", "unnamed"))
        try:
            code = country_code(p)
            reason = None if code else "No unambiguous two-letter ISO_A2 or ISO_A2_EH code"
        except ValueError as exc:
            code, reason = None, str(exc)
        if code is None:
            audit["excluded_map_units"].append({
                "name": label, "admin": p.get("ADMIN"), "reason": reason,
                "ISO_A2": p.get("ISO_A2"), "ISO_A2_EH": p.get("ISO_A2_EH"),
            })
            continue
        geometry = shape(feature["geometry"])
        if not geometry.is_valid:
            audit["repaired_map_units"].append(label)
        geometry = polygonal_valid(project(polygonal_valid(geometry)))
        if geometry.is_empty or geometry.area <= 0:
            raise ValueError(f"Empty/nonpositive mapped polygon: {label}")
        grouped[code].append(geometry)
        names[code].add(label)
        field_row = {"name": label, "ISO_A2": p.get("ISO_A2"),
                     "ISO_A2_EH": p.get("ISO_A2_EH"), "selected_country": code}
        fields[code].append(field_row)
        if p.get("ISO_A2") != code:
            audit["iso_field_fallback_map_units"].append(field_row)
    countries = {code: {
        "geometry": polygonal_valid(unary_union(parts)),
        "names": sorted(names[code]), "source_iso_fields": fields[code],
    } for code, parts in sorted(grouped.items())}
    return countries, audit


def regions_from_features(features, project):
    groups = defaultdict(list)
    levels = {}
    repaired = []
    for feature in features:
        p = feature["properties"]
        code, level = p["LEVEL3_COD"], str(p["LEVEL1_COD"])
        if code in levels and levels[code] != level:
            raise ValueError(f"Conflicting LEVEL1 assignment for {code}")
        levels[code] = level
        geometry = shape(feature["geometry"])
        if not geometry.is_valid:
            repaired.append(code)
        groups[code].append(polygonal_valid(project(polygonal_valid(geometry))))
    return {code: {"geometry": polygonal_valid(unary_union(parts)),
                   "level1": levels[code]}
            for code, parts in sorted(groups.items())}, repaired


def intersect_crosswalk(countries, regions):
    codes = sorted(regions)
    geometries = [regions[code]["geometry"] for code in codes]
    tree = STRtree(geometries)
    rows, summaries = [], []
    for country, entry in sorted(countries.items()):
        geometry = entry["geometry"]
        area = geometry.area
        if not area > 0:
            raise ValueError(f"Nonpositive country area for {country}")
        country_rows = []
        pieces = {}
        for index in tree.query(geometry, predicate="intersects"):
            code = codes[int(index)]
            region = regions[code]
            intersection = polygonal_valid(geometry.intersection(region["geometry"]))
            overlap_area = intersection.area
            if overlap_area <= 1e-8:
                continue
            region_area = region["geometry"].area
            ratio = overlap_area / min(area, region_area)
            row = {
                "country": country, "wgsrpd3": code, "level1": region["level1"],
                "intersection_km2": overlap_area / 1e6,
                "country_fraction": overlap_area / area,
                "region_fraction": overlap_area / region_area,
                "smaller_polygon_fraction": ratio,
                "material_link": int(ratio >= MATERIAL_MIN_SMALLER_FRACTION),
                "guard_link": int(ratio >= GUARD_MIN_SMALLER_FRACTION),
            }
            country_rows.append(row)
            pieces[code] = intersection
        country_rows.sort(key=lambda r: r["wgsrpd3"])
        rows.extend(country_rows)
        material = [r for r in country_rows if r["material_link"]]
        guard = [r for r in country_rows if r["guard_link"]]

        def coverage(selected):
            return unary_union([pieces[r["wgsrpd3"]] for r in selected]).area / area

        level_pieces = defaultdict(list)
        for row in material:
            level_pieces[row["level1"]].append(pieces[row["wgsrpd3"]])
        level_fractions = {
            level: unary_union(parts).area / area
            for level, parts in level_pieces.items()
        }
        dominant = max(level_fractions, key=level_fractions.get) if level_fractions else None
        flags = list(COUNTRY_MAPPING_FLAGS.get(country, []))
        material_coverage = coverage(material)
        material_sum = sum(r["country_fraction"] for r in material)
        summaries.append({
            "country": country, "names": entry["names"],
            "source_iso_fields": entry.get("source_iso_fields", []),
            "country_area_km2": area / 1e6,
            "any_coverage": coverage(country_rows),
            "material_coverage": material_coverage,
            "guard_coverage": coverage(guard),
            "material_intersection_sum_fraction": material_sum,
            "material_overlap_excess_fraction": max(0., material_sum - material_coverage),
            "material_regions": [r["wgsrpd3"] for r in material],
            "guard_regions": [r["wgsrpd3"] for r in guard],
            "material_level1_codes": sorted(level_fractions),
            "dominant_level1": dominant,
            "dominant_level1_country_fraction": level_fractions.get(dominant, 0.),
            "material_level1_country_fractions": level_fractions,
            "iso_eligible": not flags, "mapping_flags": flags,
        })
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdwg", type=Path, required=True)
    parser.add_argument("--natural-earth", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if sha256(args.natural_earth) != NATURAL_EARTH_SHA256:
        raise ValueError("Natural Earth source does not match the pinned file SHA256")
    tdwg = json.loads(args.tdwg.read_text())
    natural_earth = json.loads(args.natural_earth.read_text())
    project = equal_area_projector()
    countries, audit = dissolve_countries(natural_earth["features"], project)
    regions, repaired_regions = regions_from_features(tdwg["features"], project)
    rows, summaries = intersect_crosswalk(countries, regions)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "country_wgsrpd_links.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "country", "wgsrpd3", "level1", "intersection_km2", "country_fraction",
            "region_fraction", "smaller_polygon_fraction", "material_link", "guard_link",
        ])
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema": "country_wgsrpd_area_crosswalk_v1",
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "projection": "EPSG:6933 (WGS 84 / NSIDC EASE-Grid 2.0 Global, equal area)",
        "projection_longitude_rule": "force_over=True preserves source continuity for coordinates slightly beyond +/-180 degrees; prevents per-vertex wrapping from creating a spurious trans-Arctic MAG polygon.",
        "source": {
            "natural_earth_url": NATURAL_EARTH_URL,
            "natural_earth_commit": NATURAL_EARTH_COMMIT,
            "natural_earth_sha256": sha256(args.natural_earth),
            "tdwg_sha256": sha256(args.tdwg),
            "natural_earth_map_units": len(natural_earth["features"]),
            "tdwg_source_features": len(tdwg["features"]),
            "tdwg_unique_regions": len(regions),
        },
        "thresholds": {
            "material_min_intersection_over_smaller_polygon": MATERIAL_MIN_SMALLER_FRACTION,
            "guard_min_intersection_over_smaller_polygon": GUARD_MIN_SMALLER_FRACTION,
            "positive_area_tolerance_square_metres": 1e-8,
        },
        "country_code_rule": "Prefer two-letter ISO_A2_EH; require agreement when ISO_A2 is also two letters; dissolve all map units sharing the selected code; audit rather than impute uncoded units.",
        "coverage_rule": "Area of the union of country-intersected linked WGSRPD3 polygons divided by dissolved country map-unit area. This is coverage of the source footprint, not proof of completeness of ISO territory.",
        "country_count": len(summaries), "positive_intersection_rows": len(rows),
        "material_link_rows": sum(r["material_link"] for r in rows),
        "guard_link_rows": sum(r["guard_link"] for r in rows),
        "countries_iso_mapping_eligible": sum(r["iso_eligible"] for r in summaries),
        "repaired_tdwg_regions": repaired_regions,
        **audit,
        "countries": summaries,
        "output_links_sha256": sha256(csv_path),
        "limitations": [
            "Country-level observations cannot establish within-country co-location with a host plant.",
            "Administrative boundaries and source generalization differ; .05 material and .005 guard thresholds are fixed before ecological outcomes.",
            "WGSRPD regions can span multiple countries; country mapping inherits the regional host support without asserting presence in every local part.",
            "Uncoded or disputed map units are not silently reassigned; affected country codes carry conservative eligibility flags.",
            "Very small territories omitted by the 1:50m source remain absent from this crosswalk, not assigned zero biological occurrence.",
        ],
    }
    (out / "country_wgsrpd_summary.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in [
        "country_count", "positive_intersection_rows", "material_link_rows",
        "guard_link_rows", "countries_iso_mapping_eligible",
    ]}, indent=2))


if __name__ == "__main__":
    main()
