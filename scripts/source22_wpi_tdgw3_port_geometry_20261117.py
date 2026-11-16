#!/usr/bin/env python3
"""Pin source-eligible island-port coastal endpoints to official WGSRPD3 shapes.

Only geographically unambiguous WPI points inside the TDWG polygon are
accepted. This is a strict spatial availability pilot, NOT routes or seed
transport. Offshore WPI pins may be missed; no nearest-polygon backfilling.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from shapely.geometry import Point,shape
from shapely.strtree import STRtree

def read(p):
    with p.open(newline="",encoding="utf8") as f: return list(csv.DictReader(f))

def main():
    a=argparse.ArgumentParser()
    for k in ("source-root","level3-geojson","wpi-geojson","out-dir"):
        a.add_argument("--"+k,required=True,type=Path)
    z=a.parse_args(); z.out_dir.mkdir(parents=True,exist_ok=True)
    files=list(z.source_root.rglob("within_region_eudis_native_introduced_pairs.csv"))
    if len(files)!=1: raise RuntimeError(f"Expected one original frozen EuDiS source, got {len(files)}")
    original=read(files[0])
    islands=sorted((r for r in original if r["class"]=="island"), key=lambda r:r["code"])
    if len(islands)!=22: raise RuntimeError("Original EuDiS 22 island source changed")
    eligible={r["code"]:r for r in islands}
    shapes=json.loads(z.level3_geojson.read_text())
    geometries=[]; codes=[]
    for f in shapes["features"]:
        code=str((f.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        geom=f.get("geometry")
        if code in eligible and geom:
            sh=shape(geom)
            if sh.is_empty: continue
            geometries.append(sh); codes.append(code)
    if set(codes)!=set(eligible): raise RuntimeError("Missing exact TDWG region polygons for some original 22")
    tree=STRtree(geometries)
    src=json.loads(z.wpi_geojson.read_text())
    if len(src.get("features",[]))<3000: raise RuntimeError("WPI source does not have all >3000 source ports")
    accepted=[]; multiplicity=0; unmatched=0; malformed=0
    for f in src["features"]:
        g=f.get("geometry") or {}
        coords=g.get("coordinates") or []
        if g.get("type")!="Point" or len(coords)<2: malformed+=1;continue
        lon,lat=float(coords[0]),float(coords[1])
        if not (-180<=lon<=180 and -90<=lat<=90): malformed+=1;continue
        point=Point(lon,lat)
        hits=sorted({codes[int(i)] for i in tree.query(point) if geometries[int(i)].covers(point)})
        if len(hits)!=1:
            if len(hits)>1:multiplicity+=1
            else: unmatched+=1
            continue
        p=f.get("properties") or {}
        accepted.append({
            "WGSRPD3":hits[0], "WPI_source_OBJECTID":str(p.get("OBJECTID") or ""),
            "port_name":str(p.get("PORT_NAME") or ""), "country":str(p.get("COUNTRY") or ""),
            "latitude":round(lat,8),"longitude":round(lon,8),
            "is_eligible_noaa_legacy_oscar_lat":int(abs(lat)<=80),
            "original_WPI_2019_coordinate_inside_exact_TDWG3_polygon":1
        })
    ids=[x["WPI_source_OBJECTID"] for x in accepted]
    if len(ids)!=len(set(ids)):raise RuntimeError("WPI OBJECTID repeated")
    counts=defaultdict(int)
    for x in accepted:counts[x["WGSRPD3"]]+=1
    with (z.out_dir/"source22_pinned_WPI_exact_polygon_port_points.csv").open("w",newline="",encoding="utf8") as h:
        w=csv.DictWriter(h,fieldnames=list(accepted[0]) if accepted else
           ["WGSRPD3","WPI_source_OBJECTID","port_name","country","latitude","longitude",
            "is_eligible_noaa_legacy_oscar_lat",
            "original_WPI_2019_coordinate_inside_exact_TDWG3_polygon"])
        w.writeheader();w.writerows(accepted)
    names={
        r["code"]:{
            "level1":r["level1"],
            "native_EuDiS_positive":int(r["native_trait_1"]),
            "introduced_EuDiS_positive":int(r["introduced_trait_1"]),
            "exact_WPI_points_in_polygon":counts[r["code"]]
        } for r in islands
    }
    audit={
        "status":"EXACT_TDWG3_POLYGON_WPI_SOURCE_POINTS_NO_TRANSPORT",
        "original_EuDiS_source_run":43892007427,
        "TDWG_WGSRPD3_source_SHA":"52da7828aba9d461dd133c27b3bd7a4407161f54",
        "original_eligible_regions":22,
        "WPI_GeoJSON_original_point_features":len(src["features"]),
        "accepted_exact_geometry_WPI_points":len(accepted),
        "regions_with_accepted_WPI_points":sum(n["exact_WPI_points_in_polygon"]>0 for n in names.values()),
        "regions_zero_exact_WPI_points":[c for c,n in names.items() if n["exact_WPI_points_in_polygon"]==0],
        "unmatched_or_outside_eligible_islands":unmatched,
        "multiply_covered_WPI_points":multiplicity,
        "malformed_WPI_points":malformed,
        "original_22_by_region":names,
        "critical_limits":[
          "Unmatched ports are not absent: many of the 3630 WPI ports are outside the 22 selected regions or lie offshore beyond polygon boundaries.",
          "No nearest-polygon assignment or arbitrary coastal buffer used.",
          "WPI port point is not historical ship service, island arrival, nor native plant source.",
          "OSCAR legacy current access covers 80 degrees S to 80 degrees N; no transport trajectory run."
        ]}
    (z.out_dir/"source22_exact_WPI_TDWG3_port_geometry_audit.json").write_text(
        json.dumps(audit,ensure_ascii=False,indent=2)+"\n")
    print("SOURCE22_EXACT_ISLAND_PORT_GEOGRAPHY",json.dumps({
        "wpi_source_points":len(src["features"]),
        "accepted_points":len(accepted),
        "regions_with_point":audit["regions_with_accepted_WPI_points"],
        "missing_codes":audit["regions_zero_exact_WPI_points"],
        "non_europe_with_port_point":[c for c,n in names.items() if n["level1"]!="1" and n["exact_WPI_points_in_polygon"]>0],
        "status":audit["status"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
