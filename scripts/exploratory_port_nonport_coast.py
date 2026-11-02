#!/usr/bin/env python3
"""Deterministic, independently reselected within-island coastal control points.

Focal WPI79 port IDs are pinned by the 2026-10-31 source audit.
GSHHG 2.3.7 sibling_id polygons and public WPI points are downloaded and
verified; this DOES NOT recover original unpublished non-port coordinates.
Generated sites are a new, pre-outcome-selection observational sample.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

import shapefile
import numpy as np
from pyproj import Geod
from shapely.geometry import Point, shape
from shapely.ops import nearest_points, unary_union

ROOT = Path(__file__).resolve().parents[1]
JOINT = ROOT / "provenance/reviewer_defenses/results/island_port_first_paper_20261031/joint_79port_WPI_plant_landcover_support.json"
PORTS_CATALOGUE = ROOT / "provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"
GSHHG_URL = "https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GSHHG_SHA256 = "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
WPI = "https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
AGENT = "chocho-coastal-matched-geometry-exploration/0.1"
GEOD = Geod(ellps="WGS84")


def http_json(url: str, params: dict) -> dict:
    request_url = url + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(request_url, headers={"User-Agent": AGENT, "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=90) as response:
        obj = json.load(response)
    if not isinstance(obj, dict) or obj.get("error"):
        raise ValueError("Source returned error: " + str(obj.get("error")))
    return obj


def source_cohort() -> tuple[list[dict], dict[int, str]]:
    raw = json.loads(JOINT.read_text(encoding="utf8"))
    rows = raw["exact_linked_source_2021_landcover_and_2001_2024_botanical_port_IDs"]
    if (len(rows), len({x["landmass_id"] for x in rows}), len({x["port_id"] for x in rows})) != (79, 36, 79):
        raise ValueError("Source-locked 79 port 36 island cohort changed")
    if raw["status"] != "ID_LINKAGE_ONLY_NOT_PLANT_SITE_CONTROL":
        raise ValueError("Old source's limitations have changed")
    with PORTS_CATALOGUE.open(newline="", encoding="utf8") as fh:
        country = {int(x["port_id"]): x["country"] for x in csv.DictReader(fh)}
    if any(int(x["port_id"]) not in country for x in rows):
        raise ValueError("WPI port ID not in source-locked port country table")
    return rows, country


def read_ports() -> dict[int, dict]:
    out = {}
    for offset in range(0, 12000, 2000):
        data = http_json(WPI, {"where": "1=1", "outFields": "OBJECTID,INDEX_NO,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
                               "returnGeometry": "false", "orderByFields": "OBJECTID",
                               "resultOffset": offset, "resultRecordCount": 2000, "f": "json"})
        features = data.get("features", [])
        for f in features:
            a = f["attributes"]
            pid = int(a["OBJECTID"])
            if pid in out:
                raise ValueError("Duplicated WPI OBJECTID")
            lat, lon = a.get("LATITUDE"), a.get("LONGITUDE")
            if lat is None or lon is None:
                continue
            lat, lon = float(lat), float(lon)
            if not (-90 < lat < 90 and -180 <= lon <= 180):
                raise ValueError("Invalid WPI coordinate")
            out[pid] = {"port_id": pid, "name": a.get("PORT_NAME"),
                        "country": a.get("COUNTRY"), "index_no": a.get("INDEX_NO"),
                        "lat": lat, "lon": lon}
        if len(features) < 2000:
            break
    if len(out) < 3000:
        raise ValueError(f"Too few WPI ports retrieved: {len(out)}")
    return out


def read_gshhg(cache: Path, island_ids: set[str]) -> tuple[dict[str, object], str]:
    cache.mkdir(parents=True, exist_ok=True)
    zip_path = cache / "gshhg-shp-2.3.7.zip"
    if not zip_path.exists():
        request = urllib.request.Request(GSHHG_URL, headers={"User-Agent": AGENT})
        with urllib.request.urlopen(request, timeout=180) as response, zip_path.open("wb") as out:
            import shutil
            shutil.copyfileobj(response, out)
    sha = hashlib.sha256(zip_path.read_bytes()).hexdigest()
    if sha != GSHHG_SHA256:
        raise ValueError("GSHHG geographic archive SHA256 mismatch")
    target = cache / "GSHHS_h_L1"
    with zipfile.ZipFile(zip_path) as z:
        for ext in (".shp", ".dbf", ".shx"):
            path = "GSHHS_shp/h/GSHHS_h_L1" + ext
            if not (target.parent / (target.name + ext)).exists():
                with z.open(path) as src, (target.parent / (target.name + ext)).open("wb") as dst:
                    import shutil
                    shutil.copyfileobj(src, dst)
    reader = shapefile.Reader(str(target), encoding="utf-8")
    field_names = [f[0] for f in reader.fields[1:]]
    if "sibling_id" not in field_names:
        raise ValueError("Wrong GSHHG level-1 coastline schema")
    sibling_idx = field_names.index("sibling_id")
    indices = defaultdict(list)
    for i, row in enumerate(reader.iterRecords()):
        k = str(row[sibling_idx])
        if k in island_ids:
            indices[k].append(i)
    missing = island_ids - set(indices)
    if missing:
        raise ValueError("GSHHG target sibling_ids missing: " + str(sorted(missing)))
    polygons = {}
    for sid, locs in indices.items():
        parts = [shape(reader.shape(i).__geo_interface__) for i in locs]
        geom = unary_union(parts) if len(parts) > 1 else parts[0]
        if not geom.is_valid:
            import shapely
            geom = shapely.make_valid(geom)
        if geom.is_empty:
            raise ValueError("Empty island polygon " + sid)
        polygons[sid] = geom
    return polygons, sha


def haversine_array(lat: float, lon: float, lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    rlat = np.radians(lat)
    lat2 = np.radians(lats)
    dl = np.radians(lons - lon)
    da = lat2 - rlat
    z = np.sin(da / 2) ** 2 + np.cos(rlat) * np.cos(lat2) * np.sin(dl / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(z, 0, 1)))


def geodesic_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    return GEOD.inv(a[0], a[1], b[0], b[1])[2] / 1000.0


def inward_point(poly, anchor: Point, *, target_m: float = 250) -> Point | None:
    """Find strictly terrestrial geodesic inward offset, no sea-side controls."""
    inside = poly.representative_point()
    azi = GEOD.inv(anchor.x, anchor.y, inside.x, inside.y)[0]
    for offset in (0, 30, -30, 60, -60, 90, -90, 120, -120, 150, -150, 180):
        lon, lat, _ = GEOD.fwd(anchor.x, anchor.y, azi + offset, target_m)
        point = Point(lon, lat)
        if poly.contains(point):
            return point
    return None


def outlines(poly):
    if poly.geom_type == "Polygon":
        return [poly.exterior]
    if poly.geom_type == "MultiPolygon":
        return [v.exterior for v in poly.geoms]
    if poly.geom_type == "GeometryCollection":
        return [v.exterior for v in poly.geoms if v.geom_type == "Polygon"]
    raise ValueError(f"Unexpected GSHHG land polygon {poly.geom_type}")


def shore_candidates(poly, *, n: int) -> list[tuple[float, float]]:
    rings = outlines(poly)
    lengths = np.array([x.length for x in rings])
    if not np.isfinite(lengths).all() or np.sum(lengths) == 0:
        return []
    pts = []
    for ring, length in zip(rings, lengths):
        count = max(8, int(round(n * length / sum(lengths))))
        for i in range(count):
            anchor = ring.interpolate((i + 0.5) / count * ring.length)
            p = inward_point(poly, anchor)
            if p is not None:
                pts.append((p.x, p.y))
    return pts


def geometry_sites(*, n_shore: int, cache: Path) -> dict:
    cohort, country = source_cohort()
    ports = read_ports()
    ids = {int(x["port_id"]) for x in cohort}
    if not ids.issubset(ports):
        raise ValueError("Missing focal WPI port IDs " + str(sorted(ids - set(ports))))
    for pid in ids:
        # The old country table holds ISO2 values, the WPI mirror may hold names.
        if not ports[pid]["country"]:
            raise ValueError(f"WPI port {pid} is missing the country label")
    islands = {str(x["landmass_id"]) for x in cohort}
    polygons, shape_sha = read_gshhg(cache, islands)
    by_island = defaultdict(list)
    for row in cohort:
        by_island[str(row["landmass_id"])].append(row)
    lats = np.array([x["lat"] for x in ports.values()])
    lons = np.array([x["lon"] for x in ports.values()])
    records = []
    for island in sorted(by_island, key=int):
        poly = polygons[island]
        candidates = shore_candidates(poly, n=n_shore)
        safe = []
        for lon, lat in candidates:
            if np.min(haversine_array(lat, lon, lats, lons)) >= 8.0:
                safe.append((lon, lat))
        chosen = []
        for row in sorted(by_island[island], key=lambda x:int(x["port_id"])):
            pid = int(row["port_id"])
            original = ports[pid]
            wpi = Point(original["lon"], original["lat"])
            port_anchor = nearest_points(wpi, poly.boundary)[1]
            port_sample = inward_point(poly, port_anchor)
            if port_sample is None:
                port_sample = wpi if poly.contains(wpi) else None
            port_center = (port_sample.x, port_sample.y) if port_sample is not None else None
            port_shift = geodesic_km((wpi.x,wpi.y), port_center) if port_center else None
            eligible = []
            for q in safe:
                if any(geodesic_km(q, old) < 2.0 for old in chosen):
                    continue
                distance = geodesic_km((wpi.x,wpi.y), q)
                if distance < 100:
                    eligible.append((distance, q))
            eligible.sort(key=lambda x:(x[0],x[1]))
            nonport = eligible[0][1] if eligible else None
            if nonport:
                chosen.append(nonport)
            quality = ("PAIRED" if nonport and port_center and port_shift is not None and port_shift <= 2.0
                       else "UNPAIRED_OR_PORT_POSITION_UNCERTAIN")
            records.append({
                "island_id": island, "port_id": pid,
                "wpi_name": original["name"], "country_original": country[pid],
                "country_mirror": original["country"], "wpi_index_no": original["index_no"],
                "wpi_longitude": wpi.x, "wpi_latitude": wpi.y,
                "port_sample_lon": port_center[0] if port_center else None,
                "port_sample_lat": port_center[1] if port_center else None,
                "port_sample_shift_from_wpi_km": port_shift,
                "nonport_lon": nonport[0] if nonport else None,
                "nonport_lat": nonport[1] if nonport else None,
                "nonport_distance_from_wpi_km": geodesic_km((wpi.x,wpi.y),nonport) if nonport else None,
                "status": quality,
                "prior_2001_2010_5km_indexed_botanical_records": int(row["original_port_2001_2010_GBIF_indexed_plants_query_count"]),
            })
        print("COAST_ISLAND", island, "ports",len(by_island[island]),
              "eligible_shore_candidates",len(safe),
              "paired",sum(r["status"] == "PAIRED" for r in records if r["island_id"] == island),
              flush=True)
    stats = {
        "focal_ports": len(records),
        "focal_islands": len(islands),
        "paired_sites": sum(x["status"]=="PAIRED" for x in records),
        "paired_islands": len({x["island_id"] for x in records if x["status"]=="PAIRED"}),
        "unpaired_cases": [dict(port_id=x["port_id"],island_id=x["island_id"],
                                reason=x["status"],shift=x["port_sample_shift_from_wpi_km"])
                           for x in records if x["status"]!="PAIRED"]
    }
    return {
        "schema":"chocho_nutrient_port_nonport_independent_coast_pilot_v1",
        "status":"INDEPENDENT_RESELECTION_OF_SITES_NOT_ORIGINAL_CONTROL_COORDINATES",
        "source_gshhg_sha256":shape_sha,
        "source_WPI_url":WPI,
        "search_radius_km_for_future_botany":1,
        "nonport_minimum_distance_from_any_indexed_WPI_port_km":8,
        "within_island_control_separation_km":2,
        "point_inside_landmass_and_offset_from_shoreline_m":250,
        "maximum_focal_port_site_WPI_shift_km":2,
        "max_focal_port_to_remote_control_km":100,
        "site_selection_uses_no_butterfly_or_plant_species_response":True,
        "historical_2019_port_proxy_not_prior_shipping":True,
        "nonport_coast_urban_vs_natural_classification_unavailable":True,
        "not_original_same_island_nonport_sites":True,
        "summary":stats,
        "rows":records,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--cache-dir",type=Path,default=Path("/tmp/chocho-gshhg"))
    p.add_argument("--shore-candidates-per-island",type=int,default=512)
    args=p.parse_args()
    receipt=geometry_sites(n_shore=args.shore_candidates_per_island,cache=args.cache_dir)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("COAST_GEOMETRY_GATE",json.dumps(receipt["summary"],ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
