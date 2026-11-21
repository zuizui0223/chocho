#!/usr/bin/env python3
"""Recheck a PRIVATELY held historical chocho source ZIP. Never republish raw ZIP.

This validator does NOT register GBIF dataset DOIs, reconstruct record-level
datasetKey attribution, or prove permission to share the supplied materials.
Only anonymized aggregate QA is written to stdout or --report-out.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

R_KM = 6371.0088
PINNED = {
 "work/barriers/ports_2019.geojson": "d580e8520dc5bd77c0ab464df3bdc9ea95653c47b00c42a510dc7c1519feaa78",
 "work/barriers/global_original_species_island_year.parquet": "8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f",
 "data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json": "58b697b86925fb6cf753e2f8e68fe646076abd28038283958accdfb63a926741",
 "outputs/port_buffer_landmass_audit.csv": "549a2414d60d921476754262bbd936857f742989a873083a099913e6b5d2a910",
 "outputs/port_alien_host_broad_all_ports.csv": "58ac8c0dd9e90f385a4f9533c8b2394450493903e3d040ca5d9c8c1847d6c1d0",
 "outputs/island_host_network_islands.csv": "bbcfaae2e10d13da0380a7de6059c0ee8f7df81e02cca7cbd255932eca8f53e8",
 "outputs/island_host_network_pairs.csv": "eb0e7623dcb88a563ff7094e485d761312d58d5838d21c6d897b1f79030f5e34",
}
EXPECTED = set(PINNED) | {"RECOVERY_MANIFEST.json", "validation/port_coordinates_receipt.json"}

def rows(z, path):
    return list(csv.DictReader(io.StringIO(z.read(path).decode("utf-8-sig"))))

def spherical_mean(coords):
    vectors = []
    for lon, lat in coords:
        lo = math.radians(float(lon))
        la = math.radians(float(lat))
        if not (-180 <= float(lon) <= 180 and -90 <= float(lat) <= 90):
            raise AssertionError("Invalid WPI longitude or latitude")
        vectors.append((math.cos(la)*math.cos(lo),
                        math.cos(la)*math.sin(lo), math.sin(la)))
    if not vectors:
        raise AssertionError("Missing eligible source port coordinates")
    v = tuple(sum(p[k] for p in vectors)/len(vectors) for k in range(3))
    norm = math.sqrt(sum(x*x for x in v))
    if norm < 1e-7:
        raise AssertionError("Degenerate spherical mean")
    return tuple(x/norm for x in v)

def check(path):
    with ZipFile(path) as z:
        if z.testzip() is not None:
            raise AssertionError("ZIP CRC failure")
        names = [a.filename for a in z.infolist()]
        if len(names) != len(set(names)) or set(names) != EXPECTED:
            raise AssertionError("Expected exactly 9 known paths with no duplicates")
        if any(n.startswith("/") or "\\" in n or ".." in n.split("/") for n in names):
            raise AssertionError("Unsafe ZIP path")
        manifest = json.loads(z.read("RECOVERY_MANIFEST.json"))
        for file in manifest["files"]:
            entry = file["archive_path"]
            if entry not in PINNED:
                raise AssertionError("Unrecognized source in manifest")
            data = z.read(entry)
            if len(data) != file["bytes"]:
                raise AssertionError("Source size differs from manifest")
            h = hashlib.sha256(data).hexdigest()
            if h != PINNED[entry] or h != file["sha256"]:
                raise AssertionError("Source bytes differ from pinned original")
        if len(manifest["files"]) != 7:
            raise AssertionError("Expected seven source files in private manifest")
        g = json.loads(z.read("work/barriers/ports_2019.geojson"))
        if g["type"] != "FeatureCollection" or len(g["features"]) != 3630:
            raise AssertionError("Expected WPI port feature collection")
        points = {}
        for f in g["features"]:
            idx = str(int(f["properties"]["OBJECTID"]))
            if f["geometry"]["type"] != "Point" or idx in points:
                raise AssertionError("Invalid or duplicate port feature")
            points[idx] = f["geometry"]["coordinates"]
        port_links = rows(z, "outputs/port_buffer_landmass_audit.csv")
        pstatus = rows(z, "outputs/port_alien_host_broad_all_ports.csv")
        island_rows = rows(z, "outputs/island_host_network_islands.csv")
        pairs = rows(z, "outputs/island_host_network_pairs.csv")
        if not (len(port_links) == len(pstatus) == len(points)):
            raise AssertionError("Port source table mismatch")
        complete = {str(int(r["port_id"])) for r in pstatus if str(r["query_complete"]).lower() == "true"}
        if len(complete) != 1596:
            raise AssertionError("Original complete-plant-query cohort changed")
        p2island = {}
        for r in port_links:
            idx = str(int(r["port_id"]))
            if (r["status"] == "SINGLE_LANDMASS" and
                int(r["islands_intersecting_5km_buffer"]) == 1 and
                idx in complete):
                if idx in p2island:
                    raise AssertionError("Duplicate port geography")
                p2island[idx] = str(r["single_landmass_id"])
        island_ids = {str(r["landmass_id"]) for r in island_rows
                      if str(r["eligible"]).lower() == "true"}
        if len(island_ids) != 51:
            raise AssertionError("Original eligible island panel drift")
        by_island = defaultdict(list)
        for port, lid in p2island.items():
            if lid in island_ids:
                by_island[lid].append(points[port])
        if set(by_island) != island_ids:
            raise AssertionError("Incomplete island geo nodes")
        nodes = {k:spherical_mean(v) for k,v in by_island.items()}
        if sum(len(v) for v in by_island.values()) != 194 or len(pairs) != 1275:
            raise AssertionError("Original 194 ports / 1275 island pair roster changed")
        errors = []
        for p in pairs:
            a, b = nodes[str(p["island_a"])], nodes[str(p["island_b"])]
            dot = sum(x*y for x,y in zip(a,b))
            km = R_KM*math.acos(max(-1,min(1,dot)))
            errors.append(abs(km-float(p["port_distance_km"])))
        maxerr = max(errors)
        if maxerr > .01:
            raise AssertionError("Original 51-island distance source does not replay")
        parquet = z.read("work/barriers/global_original_species_island_year.parquet")
        if parquet[:4] != b"PAR1" or parquet[-4:] != b"PAR1":
            raise AssertionError("Historical Parquet magic invalid")
        footlen = int.from_bytes(parquet[-8:-4],"little")
        if footlen <= 0 or footlen >= len(parquet)-8:
            raise AssertionError("Invalid Parquet Thrift footer")
        old = manifest["parquet"]
        if old["rows"] != 83573 or old["landmass_zero_rows"] != 2636:
            raise AssertionError("Private Parquet manifest population unexpectedly changed")
        new = json.loads(z.read("data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json"))
        if (new["source_islands"],new["source_butterflies"],new["cells"],len(new["shard_audits"])) != (51,225,11475,6):
            raise AssertionError("Independent public51 source receipt drift")
        # Scope guard: a source-record manifest is not underlying raw GBIF records.
        return {
          "status":"VERIFIED_PRIVATE_SOURCE_BYTES_NOT_PUBLIC_DATA_ARCHIVE",
          "known_original_source_sha256":PINNED,
          "port_points":len(points),
          "complete_plant_queries":len(complete),
          "port_coordinates_used":sum(len(v) for v in by_island.values()),
          "selected_island_landmasses":len(nodes),
          "island_pair_distances_replayed":len(pairs),
          "max_error_distance_km":maxerr,
          "historical_parquet_manifest_rows":old["rows"],
          "historical_parquet_structural_header_footer_valid":True,
          "historical_parquet_original_record_ids_not_in_metadata":True,
          "public51_index_cells":new["cells"],
          "public51_source_shards":len(new["shard_audits"]),
          "unresolved":["Original WPI PDF/file-to-GeoJSON conversion and redistribution rights",
                        "GBIF original record IDs, datasetKey/publisher attribution",
                        "Independent 155/224 analysis source data and raw replay",
                        "Private original Parquet landmass zero handling",
                        "Anonymous stable licensed reviewer deposit and author approval"],
          "warning":"No original raw files, author paths, or personal identifiers uploaded by this script."
        }

def main():
    a = argparse.ArgumentParser()
    a.add_argument("--zip",type=Path,required=True)
    a.add_argument("--report-out",type=Path)
    args = a.parse_args()
    result = check(args.zip)
    print("CHOCHO_PRIVATE_RECOVERY_REPLAY_PASS",json.dumps({
         k:v for k,v in result.items() if k != "known_original_source_sha256"
    },ensure_ascii=False))
    if args.report_out:
        args.report_out.parent.mkdir(parents=True,exist_ok=True)
        args.report_out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")

if __name__ == "__main__":
    main()
