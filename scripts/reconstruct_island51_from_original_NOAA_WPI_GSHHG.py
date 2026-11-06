#!/usr/bin/env python3
"""Rebuild 51 island landmass port centroids from exact public WPI2019/GSHHG v2.3.7.

Preconditions:
 - original NOAA/PIRO World Port Index GeoJSON exactly 3630 point OBJECTID;
 - exact SHA256 GSHHG 2.3.7 high-resolution GSHHS_h_L1 polygon ID;
 - original fixed 1596 port GBIF-complete inventory in repository;
 - archived 51 island landmass IDs + every one of 1275 source pair distances.

No locations are inferred from a distance matrix. Reconstructed coordinates
are ACCEPTED only when all 1275 original spherical-mean inter-island distances
are reproduced to within 0.01 km. Otherwise, output diagnostic only.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import Point, box, shape
from shapely.ops import transform

from reconstruct_island51_verified_port_coordinates import spherical_mean, arc_distance

GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
GSHHG_REPAIRS=set()


def source_sha(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):d.update(chunk)
    return d.hexdigest()


def load_WPI(path):
    o=json.loads(Path(path).read_text(encoding="utf-8"))
    if o.get("type")!="FeatureCollection":raise ValueError("Not GeoJSON FeatureCollection")
    features=o.get("features",[])
    if len(features)!=3630:raise ValueError(f"Original 2019 WPI must have 3630, got {len(features)}")
    coords={}
    for feature in features:
        geom=feature["geometry"]
        if geom["type"]!="Point":raise ValueError("Original WPI point unexpectedly non-point")
        o=feature["properties"]["OBJECTID"]
        key=str(int(o))
        xy=geom["coordinates"]
        if len(xy)!=2 or not (-180<=xy[0]<=180 and -90<=xy[1]<=90):
            raise ValueError("Bad WPI latitude/longitude")
        if key in coords:raise ValueError("Duplicate WPI OBJECTID after original two page collection")
        coords[key]=(float(xy[0]),float(xy[1]))
    if len(coords)!=3630:raise ValueError("Duplicate WPI keys")
    return coords


def load_GSHHG(zip_file):
    got=source_sha(zip_file)
    if got!=GSHHG_SHA:raise ValueError("GSHHG original v2.3.7 zip checksum changed: "+got)
    with zipfile.ZipFile(zip_file) as z:
        paths=[n for n in z.namelist() if n.endswith("/h/GSHHS_h_L1.shp")]
        if len(paths)!=1:raise ValueError("GSHHS_h_L1.shp missing or ambiguous")
    # geopandas uses GDAL's vsi-zip reader; no manual ZIP unpacks necessary.
    p=f"zip://{Path(zip_file).resolve()}!{paths[0]}"
    g=gpd.read_file(p)
    if not {"sibling_id","geometry"}.issubset(g.columns):
        raise ValueError("GSHHG high-res L1 coastline changed schema "+str(list(g.columns)))
    if g.crs is None:raise ValueError("GSHHG land CRS not set")
    g=g.to_crs(4326)
    g["landmass_id"]=g.sibling_id.astype(str)
    # GSHHG sibling_id is the parent landmass identity, and can legitimately
    # occur on several polygon features. Port intersections are deduplicated by
    # sibling_id; a one-to-many GSHHG feature grouping is not a data error.
    g=g[g.geometry.notna()].copy()
    return g,got


def candidate_polygons(g,lon,lat,radius_m):
    dy=(radius_m+1000)/110500.0
    dx=dy/max(.05,math.cos(math.radians(lat)))
    if lon-dx < -180 or lon+dx > 180:
        boxes=[box(-180,max(-90,lat-dy),-180+2*dx,min(90,lat+dy)),
               box(180-2*dx,max(-90,lat-dy),180,min(90,lat+dy))]
        ids=set()
        for b in boxes:ids.update(g.sindex.query(b,predicate="intersects"))
        return sorted(ids)
    return sorted(set(g.sindex.query(box(lon-dx,max(-90,lat-dy),lon+dx,min(90,lat+dy)),
                      predicate="intersects")))


def region_overlaps(g,lon,lat,radius_m=5000):
    # Project only nearby candidate geometries into local metre units.
    # Use exact geodesic radial buffer approximation in aeqd centered on port.
    crs=CRS.from_proj4(f"+proj=aeqd +lat_0={lat:.10f} +lon_0={lon:.10f} +datum=WGS84 +units=m +no_defs")
    project=Transformer.from_crs(4326,crs,always_xy=True).transform
    circle=Point(0,0).buffer(radius_m,quad_segs=72)
    matches=[]
    nearby=candidate_polygons(g,lon,lat,radius_m)
    for i in nearby:
        geom=g.iloc[int(i)].geometry
        if geom is None:continue
        # Clipping longitude/latitude window avoids projecting a global
        # mainland polygon through a singularity on the opposite hemisphere.
        dy=(radius_m+2500)/110500.0
        dx=dy/max(.05,math.cos(math.radians(lat)))
        window=box(max(-180,lon-dx),max(-90,lat-dy),min(180,lon+dx),min(90,lat+dy))
        # GSHHG high-resolution source contains a small number of invalid
        # self-intersecting original coastline polygons. Repair their local
        # geometry without changing sibling_id landmass classification; log
        # each repaired feature for reproducibility.
        if not shapely.is_valid(geom):
            GSHHG_REPAIRS.add((str(g.iloc[int(i)].landmass_id),int(i)))
            geom=shapely.make_valid(geom)
        clipped=shapely.intersection(geom,window)
        if clipped.is_empty:continue
        projected=transform(project,clipped)
        if projected.intersects(circle):
            matches.append(str(g.iloc[int(i)].landmass_id))
    return sorted(set(matches))


def run(source_ports,land,complete,archived_islands,archived_pairs):
    ports=load_WPI(source_ports)
    polys,gsha=load_GSHHG(land)
    c=pd.read_csv(complete,dtype={"port_id":str},keep_default_na=False)
    # Archived status ledger includes ALL 3,630 indexed ports; only 1,596
    # have a fully completed GBIF botanical query. Missing entries are NOT 0.
    if len(c)!=3630 or c.port_id.nunique()!=3630:
        raise ValueError("Original port ledger must contain 3630 unique WPI ports")
    q=c[c.query_complete.astype(str).str.casefold().eq("true")]
    valid=set(str(int(x)) for x in q.port_id)
    islands=pd.read_csv(archived_islands,dtype={"landmass_id":str})
    pairs=pd.read_csv(archived_pairs,dtype={"island_a":str,"island_b":str})
    need=islands[islands.eligible.astype(str).str.lower().eq("true")].copy()
    if len(need)!=51 or len(pairs)!=1275:raise ValueError("Original 51×1275 source changed")
    wanted=set(need.landmass_id)
    available=set(polys.landmass_id)
    if not wanted.issubset(available):
        raise ValueError("The original 51 landmass sibling_id keys aren't in this GSHHG dataset: "+str(sorted(wanted-available)))
    assignment={}
    distribution=collections.Counter()
    ambiguous_source_island_candidates=[]
    for t,portid in enumerate(sorted(valid,key=int),1):
        if portid not in ports:raise ValueError("Original WPI OBJECTID not in current public layer: "+portid)
        lon,lat=ports[portid]
        overlaps=region_overlaps(polys,lon,lat,5000)
        if len(overlaps)==1:
            assignment[portid]=overlaps[0]
        elif len(overlaps)>1 and (set(overlaps)&wanted):
            ambiguous_source_island_candidates.append({
                "original_2019_WPI_port_id":portid,
                "original_coordinates_lonlat":[lon,lat],
                "all_sibling_landmass_overlaps":overlaps,
                "eligible_source51_landmass_overlaps":sorted(set(overlaps)&wanted,key=int),
                "n_intersecting_GSHHG_L1_landmasses":len(overlaps)})
        distribution[str(len(overlaps))]+=1
        if t%200==0 or t==len(valid):
            print("WPI_GSHHG_PORT_5KM_LANDMASS_PROGRESS",t,len(valid),
                  "unique",len(assignment),"source51_coverage",
                  len(set(assignment.values())&wanted),flush=True)
    assigned=collections.defaultdict(list)
    for pid,lid in assignment.items():
        if lid in wanted:assigned[lid].append(ports[pid])
    counts={lid:len(assigned[lid]) for lid in sorted(wanted,key=int)}
    old_counts={str(r.landmass_id):int(r.ports) for r in need.itertuples(index=False)}
    difference={lid:{"expected":old_counts[lid],"rebuilt":counts[lid]}
                for lid in counts if counts[lid]!=old_counts[lid]}
    present={lid for lid in wanted if counts[lid]>=1}
    nodes={lid:spherical_mean(assigned[lid]) for lid in present}
    validpairs=[]
    for row in pairs.itertuples(index=False):
        a,b=row.island_a,row.island_b
        if a in nodes and b in nodes:
            predicted=arc_distance(nodes[a]["unit_sphere_xyz"],nodes[b]["unit_sphere_xyz"])
            observed=float(row.port_distance_km)
            validpairs.append((a,b,predicted-observed))
    errors=[abs(item[2]) for item in validpairs]
    qres={
        "status":"ORIGINAL_COORDINATES_FULLY_VALIDATED" if (len(present)==51 and
                    len(validpairs)==1275 and not difference and max(errors,default=float("inf"))<=.01)
                    else "GEOSPATIAL_SOURCE_RECONSTRUCTION_MISMATCH_NOT_FOR_RASTER_MODELS",
        "expected_n_WPI":3630,"source_complete_WPI":len(valid),
        "source_GSHHG":gsha,"original_51_island_ids":sorted(wanted,key=int),
        "GSHHG_L1_polygons":len(polys),
        "source_GSHHG_invalid_geometries_repaired":len(GSHHG_REPAIRS),
        "repaired_source_landmass_ids":sorted({x[0] for x in GSHHG_REPAIRS},key=int),
        "matching_sibling_landmasses":len(wanted&available),
        "n_complete_port_assignments":len(assignment),
        "n_island_source_locations_rebuilt":len(present),
        "n_matched_original_1275_pairs":len(validpairs),
        "max_pair_distance_error_km":max(errors,default=None),
        "median_pair_distance_error_km":float(np.median(errors)) if errors else None,
        "count_mismatch_summary":difference,
        "ambiguous_source_island_port_candidates":ambiguous_source_island_candidates,
        "port_buffer_polygon_overlap_counts":dict(distribution),
        "source_files_sha256":{"WPI":source_sha(source_ports),
                               "GSHHG":gsha,"completed_index":source_sha(complete),
                               "original51":source_sha(archived_islands),
                               "original1275":source_sha(archived_pairs)},
        "diagnostics_20_largest_residuals":[
            {"island_a":a,"island_b":b,"signed_residual_km":d}
            for a,b,d in sorted(validpairs,key=lambda t:-abs(t[2]))[:20]],
        "landmass_points_only_if_validated":{
            lid:{"longitude":nodes[lid]["longitude"],"latitude":nodes[lid]["latitude"],
                 "source_complete_ports":counts[lid]}
            for lid in sorted(nodes,key=int)
        } if (len(present)==51 and len(validpairs)==1275 and not difference
              and max(errors,default=float("inf"))<=.01) else {},
        "limitations":[
            "Source method uses 5km geodesic AEQD buffer, GSHHG v2.3.7 high-res level1 land polygons and WPI2019 completed port points.",
            "A failed exact source distance reproduction prevents any AIS/wind map overlay; no invented island latitude/longitude.",
            "This is not an AIS transport model, nor actual marine corridors.",
            "Port buffer can touch islands across narrow channels; any unmatched port count is retained as a mismatch.",
            "Original output CSV port-to-landmass assignments are not themselves in the public branch, so the coast algorithm is independently reconstructed."
        ]
    }
    return qres


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("WPI-geojson","GSHHG-zip","completed-port-inventory",
              "island-source","pair-source","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    out=run(a.WPI_geojson,a.GSHHG_zip,a.completed_port_inventory,
            a.island_source,a.pair_source)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("WPI_GSHHG_51_ISLAND_COORDINATE_SOURCE_RECONSTRUCTION",
           json.dumps({k:v for k,v in out.items() if k not in (
                 "diagnostics_20_largest_residuals","landmass_points_only_if_validated",
                 "original_51_island_ids","source_files_sha256")},sort_keys=True),flush=True)
    # Nonmatching original inputs are expected to be reported, not hidden or
    # rewritten as verified island coordinates.
if __name__=="__main__":main()
