#!/usr/bin/env python3
"""Measured island coastline separation and present WPI location proximity.

Distances are spherical great-circle distances between densified source coast
segments, not butterfly movement costs. WPI proximity is not traffic or an
on-island port assignment. Only derived measurements are exported.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import zipfile
from pathlib import Path
from urllib.parse import urlencode

import geopandas as gpd
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import Point

from audit_gshhg_island_source_census_20261029 import dissolve_original
from audit_island_gshhg_original_source_20261029 import URL, SHA, MAX_SIZE
from audit_island_port_original_sources_20261029 import WPI, request
from audit_world_port_index_complete_20261029 import (
    EXPECTED, ORIGINAL_METADATA_SHA, ORIGINAL_IDS_SHA, validate,
)

ROOT = Path(__file__).resolve().parents[1]
R_KM = 6371.0088
# Independently reconciled with successful source Actions run 40857801018,
# artifact 12890257067, before measurement. Metadata/IDs alone do not lock coords.
PORT_PAGE_SHA = "b610760474adb1bdddd3b5c7bf154b6098864d12c958f9b3a5d57464d45a90f3"
CONTINENT_SEEDS = {
    "Eurasia": (10.0, 50.0), "Africa": (20.0, 0.0),
    "North_America": (-100.0, 45.0), "South_America": (-60.0, -15.0),
    "Australia": (135.0, -25.0),
}
ISLAND_SEEDS = {
    "Greenland": (-42.0, 72.0), "Madagascar": (47.0, -19.0),
    "Borneo": (114.0, 0.0), "Hawaii_big_island": (-155.5, 19.6),
    "Galapagos_Isabela": (-91.1, -0.7), "Izu_Oshima": (139.4, 34.74),
    "Okinawa": (127.9, 26.5), "Honshu": (138.0, 36.0),
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit_xyz(lonlat):
    """Map geographic points to a unit sphere, including across the dateline."""
    a = np.asarray(lonlat, dtype=float)
    if a.ndim != 2 or a.shape[1] != 2 or not np.isfinite(a).all():
        raise ValueError("Invalid geographic coordinates")
    if np.any(np.abs(a[:, 1]) > 90) or np.any(np.abs(a[:, 0]) > 180.000001):
        raise ValueError("Out-of-range geographic coordinates")
    lon, lat = np.radians(a).T
    return np.column_stack((np.cos(lat)*np.cos(lon),
                            np.cos(lat)*np.sin(lon), np.sin(lat)))


def chord_to_km(chord):
    return 2*R_KM*np.arcsin(np.clip(np.asarray(chord)/2, 0, 1))


def shoreline_xyz(geometry, max_segment_km):
    """Densify exterior coast arcs; remove artificial dateline cut edges.

    Both endpoints of every retained arc are included. Every point on that arc
    is within max_segment_km/2 along-arc distance of a retained sample.
    """
    if not (0 < max_segment_km <= 5):
        raise ValueError("Coast sampling step must be in (0,5] km")
    polygons = list(geometry.geoms) if geometry.geom_type == "MultiPolygon" else [geometry]
    sampled = []
    removed = 0
    for polygon in polygons:
        if polygon.geom_type != "Polygon":
            raise ValueError("Non-polygon source landmass")
        xy = np.asarray(polygon.exterior.coords, dtype=float)[:, :2]
        xyz = unit_xyz(xy)
        seam = np.isclose(np.abs(xy[:-1, 0]), 180, atol=1e-8, rtol=0) & np.isclose(
            np.abs(xy[1:, 0]), 180, atol=1e-8, rtol=0)
        removed += int(seam.sum())
        a, b = xyz[:-1][~seam], xyz[1:][~seam]
        if not len(a):
            continue
        angle = np.arctan2(np.linalg.norm(np.cross(a, b), axis=1),
                           np.einsum("ij,ij->i", a, b))
        keep = angle > 1e-12
        a, b, angle = a[keep], b[keep], angle[keep]
        if not len(a):
            continue
        if np.any(angle >= math.pi-1e-6):
            raise ValueError("Ambiguous near-antipodal source coast segment")
        divisions = np.maximum(1, np.ceil(angle*R_KM/max_segment_km)).astype(int)
        ix = np.repeat(np.arange(len(a)), divisions)
        start = np.repeat(np.cumsum(divisions)-divisions, divisions)
        t = (np.arange(len(ix))-start)/divisions[ix]
        theta = angle[ix]
        points = ((np.sin((1-t)*theta)/np.sin(theta))[:, None]*a[ix]
                  + (np.sin(t*theta)/np.sin(theta))[:, None]*b[ix])
        points /= np.linalg.norm(points, axis=1)[:, None]
        sampled.extend((points, b))
    if not sampled:
        raise ValueError("No non-cut coastline in source polygon")
    result = np.concatenate(sampled)
    return result, {"removed_dateline_cut_segments": removed,
                    "coast_sample_points": len(result)}


def coast_gap_interval(island_points, mainland_tree, step_km):
    upper = float(chord_to_km(mainland_tree.query(island_points, k=1)[0].min()))
    return max(0.0, upper-step_km), upper


def port_exposure(island_points, port_points, step_km):
    distances = chord_to_km(cKDTree(island_points).query(port_points, k=1)[0])
    upper = float(distances.min())
    result = {"nearest_wpi_coast_distance_lower_km": max(0.0, upper-step_km/2),
              "nearest_wpi_coast_distance_upper_km": upper}
    for radius in (50, 100):
        result[f"wpi_entries_within_{radius}km_certain"] = int(np.sum(distances <= radius))
        result[f"wpi_entries_within_{radius}km_possible"] = int(np.sum(distances <= radius+step_km/2))
    return result


def antarctic_cap_distance_lower(coast_points, step_km):
    """Conservative bound even when spherical source arcs bow southward."""
    southernmost_sample = float(np.arcsin(np.clip(coast_points[:, 2].min(), -1, 1)))
    sampled_cap_distance = max(0.0, southernmost_sample+math.pi/3)*R_KM
    return max(0.0, sampled_cap_distance-step_km/2)


def verify_port_pages(page_hashes, expected=PORT_PAGE_SHA):
    actual = hashlib.sha256("\n".join(sorted(page_hashes)).encode()).hexdigest()
    if actual != expected:
        raise ValueError("Original WPI source coordinates/page content changed")
    return actual


def usable_port_points(rows):
    """The numeric range check alone accepts the source's (0,0) placeholder."""
    lonlat = np.array([[float(r["attributes"]["longitude"]),
                        float(r["attributes"]["latitude"])] for r in rows])
    zerozero = np.all(lonlat == 0.0, axis=1)
    return unit_xyz(lonlat[~zerozero]), int(zerozero.sum())


def source_landmasses(cache):
    archive = cache / "gshhg-shp-2.3.7.zip"
    if not archive.exists():
        archive.write_bytes(request(URL, MAX_SIZE))
    if digest(archive) != SHA:
        raise ValueError("GSHHG source SHA256 differs from locked source")
    with zipfile.ZipFile(archive) as z:
        members = [n for n in z.namelist() if n.endswith("GSHHS_h_L1.shp")]
        if len(members) != 1:
            raise ValueError("Expected one original high-resolution L1 source")
        stem = members[0][:-4]
        for suffix in (".shp", ".shx", ".dbf", ".prj"):
            (cache / ("GSHHS_h_L1"+suffix)).write_bytes(z.read(stem+suffix))
    source = gpd.read_file(cache / "GSHHS_h_L1.shp")
    grouped = dissolve_original(source)
    grouped["area_km2"] = grouped.to_crs("EPSG:6933").area / 1e6
    if len(source) != 144749 or len(grouped) != 144744:
        raise ValueError("Original landmass census count differs")
    return grouped


def source_ports(cache):
    """Cache exact original API bytes; revalidate identities and coordinates."""
    def get(name, url, size=3_000_000):
        p = cache / name
        if not p.exists():
            p.write_bytes(request(url, size))
        return p.read_bytes()
    meta = get("wpi_metadata.json", WPI+"?f=json")
    ids_raw = get("wpi_ids.json", WPI+"/query?"+urlencode(
        {"f": "json", "where": "1=1", "returnIdsOnly": "true"}))
    if hashlib.sha256(meta).hexdigest() != ORIGINAL_METADATA_SHA or hashlib.sha256(ids_raw).hexdigest() != ORIGINAL_IDS_SHA:
        raise ValueError("WPI metadata or source roster changed")
    ids = sorted(json.loads(ids_raw)["objectIds"])
    if len(ids) != EXPECTED or len(set(ids)) != EXPECTED:
        raise ValueError("Original WPI source roster differs")
    rows, hashes = [], []
    for offset in range(0, len(ids), 200):
        chunk = ids[offset:offset+200]
        query = urlencode({"f": "json", "where": "1=1", "objectIds": ",".join(map(str, chunk)),
                           "returnGeometry": "false", "outFields": "objectid,latitude,longitude,country,index_no"})
        body = get(f"wpi_page_{offset//200:02}.json", WPI+"/query?"+query, 1_500_000)
        page = json.loads(body)
        if "error" in page or len(page.get("features", [])) != len(chunk):
            raise ValueError("Incomplete original WPI page")
        rows.extend(page["features"])
        hashes.append(hashlib.sha256(body).hexdigest())
    quality = validate(rows, ids)
    if quality["missing_or_invalid_attribute_coordinates"] != 0:
        raise ValueError("Every port must have valid original coordinates for this measurement")
    canonical_sha = verify_port_pages(hashes)
    points, excluded_zerozero = usable_port_points(rows)
    return points, {"source_uri": WPI, "metadata_sha256": ORIGINAL_METADATA_SHA,
        "ids_sha256": ORIGINAL_IDS_SHA, "page_sha256": hashes,
        "canonical_page_digest": canonical_sha, "canonical_page_digest_matches_original_actions_audit": True,
        "independent_original_audit_run": 40857801018, "independent_original_audit_artifact": 12890257067,
        "coordinate_quality": quality, "zerozero_placeholder_entries_excluded": excluded_zerozero,
        "usable_port_entries": len(points)}


def seed_indices(frame, seeds):
    result = {}
    for name, lonlat in seeds.items():
        ix = frame.index[frame.geometry.covers(Point(lonlat))].tolist()
        if len(ix) != 1:
            raise ValueError(f"Known-landform seed {name} has {len(ix)} source matches")
        result[name] = int(ix[0])
    return result


def derived_key(source_id):
    return "gshhg237_"+hashlib.sha256(f"GSHHG2.3.7:h:L1:{int(source_id)}".encode()).hexdigest()[:16]


def support_tables(rows):
    groups = {}
    for r in rows:
        area = r["area_km2"]
        area_bin = "5_to_20" if area < 20 else "20_to_100" if area < 100 else "100_to_1000" if area < 1000 else "1000_plus"
        lo, hi = r["mainland_coast_gap_lower_km"], r["mainland_coast_gap_upper_km"]
        edges = (0, 10, 50, 100, 500, float("inf"))
        bins = ("0_to_10", "10_to_50", "50_to_100", "100_to_500", "500_plus")
        bin_lo = next(i for i in range(5) if lo < edges[i+1])
        bin_hi = next(i for i in range(5) if hi < edges[i+1])
        isolation = bins[bin_hi] if bin_lo == bin_hi else "sampling_threshold_ambiguous"
        port = ("within50" if r["wpi_entries_within_50km_certain"] else
                "beyond50" if not r["wpi_entries_within_50km_possible"] else "sampling_threshold_ambiguous")
        key = (area_bin, isolation, port, r["antarctica_could_be_nearer"])
        groups[key] = groups.get(key, 0)+1
    return [{"area_class_km2": k[0], "mapped_mainland_gap_class_km": k[1],
             "wpi_coast_proximity_class_km": k[2], "antarctica_could_be_nearer": k[3],
             "islands": n} for k, n in sorted(groups.items())]


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(cache, output, step):
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    ports, port_manifest = source_ports(cache)
    print(json.dumps({"stage": "source_ports_verified", "entries": len(ports)}), flush=True)
    land = source_landmasses(cache)
    continents = seed_indices(land, CONTINENT_SEEDS)
    if len(set(continents.values())) != 5 or set(continents.values()) != set(land.index[land.area_km2 > 7_000_000]):
        raise ValueError("Five seed-defined continental landmasses disagree with prior source census")
    check_islands = seed_indices(land, ISLAND_SEEDS)
    if set(check_islands.values()) & set(continents.values()):
        raise ValueError("Known island incorrectly shares continental identity")
    islands = land[(land.area_km2 >= 5) & ~land.index.isin(continents.values())]
    if len(islands) != 8260 or int((islands.area_km2 >= 20).sum()) != 3696:
        raise ValueError("Island area census no longer matches source gate")
    mainland_points, seam_count = [], 0
    for ix in continents.values():
        points, qc = shoreline_xyz(land.loc[ix, "geometry"], step)
        mainland_points.append(points)
        seam_count += qc["removed_dateline_cut_segments"]
    mainland = np.concatenate(mainland_points)
    tree = cKDTree(mainland)
    print(json.dumps({"stage": "mainland_shoreline_ready", "points": len(mainland), "cut_segments_removed": seam_count}), flush=True)
    rows, points_total = [], 0
    for number, (ix, island) in enumerate(islands.iterrows(), 1):
        points, qc = shoreline_xyz(island.geometry, step)
        lower, upper = coast_gap_interval(points, tree, step)
        # Entire Antarctic mainland lies south of 60 S. The cap gives a
        # conservative distance lower bound, so uncertain comparisons are flagged.
        antarctic_cap_lower = antarctic_cap_distance_lower(points, step)
        row = {"island_key": derived_key(island.parent_source_ID),
               "area_km2": float(island.area_km2),
               "area_at_least20km2": int(island.area_km2 >= 20),
               "mainland_coast_gap_lower_km": lower,
               "mainland_coast_gap_upper_km": upper,
               "antarctica_could_be_nearer": int(upper > antarctic_cap_lower),
               **port_exposure(points, ports, step),
               "coast_sample_points": qc["coast_sample_points"],
               "removed_dateline_cut_segments": qc["removed_dateline_cut_segments"]}
        rows.append(row)
        points_total += len(points)
        seam_count += qc["removed_dateline_cut_segments"]
        if number % 1000 == 0:
            print(json.dumps({"islands_measured": number}), flush=True)
    bykey = {r["island_key"]: r for r in rows}
    checks = {name: bykey[derived_key(land.loc[ix, "parent_source_ID"])] for name, ix in check_islands.items()}
    tables = support_tables(rows)
    write_csv(output / "island_geography_exposures.csv", rows)
    write_csv(output / "island_geography_support.csv", tables)
    eligible = [r for r in rows if not r["antarctica_could_be_nearer"]]
    remote = [r for r in eligible if r["mainland_coast_gap_lower_km"] >= 500]
    summary = {"islands_5km2": len(rows), "islands_20km2": sum(r["area_at_least20km2"] for r in rows),
        "antarctica_comparison_flagged": len(rows)-len(eligible),
        "mainland_comparison_unambiguous_islands": len(eligible),
        "eligible_coast_gap_median_upper_km": float(np.median([r["mainland_coast_gap_upper_km"] for r in eligible])),
        "eligible_with_WPI_entry_certainly_within50km": sum(r["wpi_entries_within_50km_certain"] > 0 for r in eligible),
        "eligible_with_WPI_entry_certainly_within100km": sum(r["wpi_entries_within_100km_certain"] > 0 for r in eligible),
        "eligible_gap_certainly_at_least500km": len(remote),
        "remote500_WPI_certainly_within50km": sum(r["wpi_entries_within_50km_certain"] > 0 for r in remote),
        "remote500_WPI_certainly_beyond50km": sum(r["wpi_entries_within_50km_possible"] == 0 for r in remote),
        "remote500_WPI_50km_sampling_ambiguous": sum(r["wpi_entries_within_50km_certain"] == 0 and r["wpi_entries_within_50km_possible"] > 0 for r in remote)}
    manifest = {"schema": "chocho_island_coast_port_measurements_20261029_v1",
        "status": "GEOGRAPHIC_EXPOSURES_MEASURED_BIOLOGICAL_EFFECT_NOT_ESTIMATED",
        "source_commit": "a8416af4151e16ae7050240348765f9528d20748",
        "gshhg": {"url": URL, "zip_sha256": SHA, "resolution": "h", "level": 1,
                  "dissolved_landmasses": len(land)},
        "ports": port_manifest,
        "method": {"earth": "sphere", "radius_km": R_KM, "max_coast_segment_km": step,
                   "coast_to_coast_discretization_bound_km": step,
                   "point_to_coast_discretization_bound_km": step/2,
                   "bounds_cover": "Densification only, not GSHHG/WPI coordinate accuracy, shoreline age or spherical versus ellipsoidal Earth difference",
                   "coast_distance": "Shortest spherical arc between source coastline samples; excludes lake shores and artificial dateline cut edges",
                   "mainland_definition": "Five source landmasses independently matched to continental seeds and prior >7 million km2 census",
                   "antarctica_policy": "Absent from source L1; flag when the conservative south-of-60S cap could be nearer than a mapped continental coast. Bound uses minimum spherical sampled latitude minus half the sample spacing. Exclude flagged islands from primary support summaries.",
                   "port_interpretation": "WPI registered entries near ANY island coast, including entries on other landmasses; not an on-island-port, shipping-volume or historical-arrival measure",
                   "port_coverage": "WPI is a major-port registry, not all harbors; no nearby WPI entry does not imply no real port. The source's single (0,0) placeholder is excluded after auditing all 3689 entries.",
                   "source_ids": "Derived SHA256 keys of GSHHG version, resolution, level and sibling ID; raw source geometries and port records not redistributed"},
        "continental_seed_area_km2": {n: float(land.loc[i, "area_km2"]) for n, i in continents.items()},
        "named_island_seed_checks": checks,
        "computation_qc": {"mainland_sample_points": len(mainland), "island_sample_points": points_total,
                           "artificial_dateline_cut_segments_removed": seam_count},
        "summary": summary,
        "output_sha256": {n: digest(output/n) for n in ("island_geography_exposures.csv", "island_geography_support.csv")},
        "no_butterfly_effect_estimated": True,
    }
    (output / "island_geography_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n")
    print(json.dumps(summary, sort_keys=True), flush=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"data/exploratory/island_port_20261029")
    parser.add_argument("--max-coast-segment-km", type=float, default=1.0)
    args = parser.parse_args()
    run(args.cache_dir, args.output_dir, args.max_coast_segment_km)
