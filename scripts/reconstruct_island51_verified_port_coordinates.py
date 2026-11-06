#!/usr/bin/env python3
"""Recover 51 absolute island representative WPI port locations from ORIGINAL
source raw coordinates, refusing any geographical guesses.

Inputs:
 - original archive outputs/port_buffer_landmass_audit.csv
 - original archive work/barriers/ports_2019.geojson (WPI2019 OBJECTID)
 - archived 215/1596 port plant-inventory complete-status table
 - source-locked island_host_network_islands.csv and all 1275 original pair
   great-circle distances

Replicates original spherical-mean-of-complete-port-coordinates algorithm
from provenance/reviewer_defenses/results/island_host_network_20261030/
analyze_island_host_network.py.

No access to original port geolocations is not a license to infer coordinates
from landmass IDs or a pairwise distance matrix: reject and report missing.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import defaultdict

R_KM=6371.0088

def spherical_mean(lon_lat):
    x=y=z=0.0
    n=len(lon_lat)
    if n==0:raise ValueError("No port coordinates for a landmass")
    for lon,lat in lon_lat:
        if not (-180<=lon<=180 and -90<=lat<=90):
            raise ValueError("Impossible original WPI coordinates")
        ll=math.radians(lon)
        a=math.radians(lat)
        co=math.cos(a)
        x+=co*math.cos(ll);y+=co*math.sin(ll);z+=math.sin(a)
    x/=n;y/=n;z/=n
    norm=math.sqrt(x*x+y*y+z*z)
    if norm<1e-7:raise ValueError("Spherical average not uniquely defined")
    x/=norm;y/=norm;z/=norm
    return {"longitude":math.degrees(math.atan2(y,x)),
            "latitude":math.degrees(math.atan2(z,math.hypot(x,y))),
            "unit_sphere_xyz":[x,y,z]}

def arc_distance(a,b):
    dot=sum(x*y for x,y in zip(a,b))
    return R_KM*math.acos(max(-1,min(1,dot)))

def load_csv(p):
    with p.open(newline="",encoding="utf-8-sig") as f:
        yield from csv.DictReader(f)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("port-landmass","WPI-geojson","all-ports","islands","pairs","out"):
        ap.add_argument("--"+name,type=Path,required=True)
    args=ap.parse_args()
    files=[args.port_landmass,args.WPI_geojson,args.all_ports,args.islands,args.pairs]
    absent=[str(p) for p in files if not p.exists() or not p.is_file()]
    if absent:
        raise FileNotFoundError("Verified original source inputs must be supplied, no synthesized coordinates: "+str(absent))
    ports=list(load_csv(args.all_ports))
    complete={str(int(r["port_id"])) for r in ports if str(r["query_complete"]).lower()=="true"}
    geo={}
    for r in load_csv(args.port_landmass):
        if (r.get("status")!="SINGLE_LANDMASS" or
            int(float(r["islands_intersecting_5km_buffer"]))!=1):continue
        key=str(int(r["port_id"]))
        if key in complete:
            if key in geo:raise ValueError("Port maps to multiple landmasses")
            geo[key]=str(r["single_landmass_id"]).strip()
    source=json.loads(args.WPI_geojson.read_text(encoding="utf-8"))
    if source.get("type")!="FeatureCollection" or not isinstance(source.get("features"),list):
        raise ValueError("WPI original GeoJSON not FeatureCollection")
    coords={}
    for feature in source["features"]:
        k=str(int(feature["properties"]["OBJECTID"]))
        xy=feature["geometry"]["coordinates"]
        if len(xy)!=2:raise ValueError("WPI source contains geometry other than Point")
        if k in coords:raise ValueError("Duplicate original WPI OBJECTID")
        coords[k]=(float(xy[0]),float(xy[1]))
    missing=[k for k in geo if k not in coords]
    if missing:raise ValueError("Port IDs missing coordinates: "+str(missing[:12]))
    eligible=[r for r in load_csv(args.islands) if str(r["eligible"]).lower()=="true"]
    ids={str(r["landmass_id"]) for r in eligible}
    if len(ids)!=51:raise ValueError("Expected 51 eligible original island-landmass IDs")
    grouped=defaultdict(list)
    for port,lid in geo.items():
        if lid in ids:grouped[lid].append(coords[port])
    if set(grouped)!=ids:
        raise ValueError("Missing original verified port coordinates for selected island IDs: "+str(ids-set(grouped)))
    nodes={k:dict(spherical_mean(v),source_ports=len(v)) for k,v in grouped.items()}
    pairs=list(load_csv(args.pairs))
    if len(pairs)!=1275:raise ValueError("Expected 1275 locked island pairs")
    errors=[]
    for r in pairs:
        i=str(r["island_a"]);j=str(r["island_b"])
        if i not in nodes or j not in nodes:raise ValueError("Unknown pair endpoints")
        current=arc_distance(nodes[i]["unit_sphere_xyz"],nodes[j]["unit_sphere_xyz"])
        archived=float(r["port_distance_km"])
        errors.append(abs(current-archived))
    maxerror=max(errors)
    if not math.isfinite(maxerror) or maxerror>.01:
        raise ValueError(f"Recovered port geocentroids do not reproduce all 1275 archived distances: max km deviation {maxerror}")
    receipt={"schema":"chocho_island51_verified_2019_WPI_representative_geo_nodes",
             "status":"UNIQUE_SOURCE_COORDS_REPRODUCED_ALL_ORIGINAL_1275_GEODESICS",
             "n_source_landmasses":len(nodes),"n_source_pairs_verified":len(pairs),
             "max_original_pair_distance_error_km":maxerror,
             "sources_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
             "source_notes":[
               "Spherical mean of original WPI2019 port locations, not polygon centroids.",
               "Port query completeness and original SINGLE_LANDMASS mapping preserved.",
               "Coordinate table is derived from original verified individual port objects.",
               "Geographical positions can now be used for later AIS or MERRA-2/ERA5 raster extraction, not movement inference."],
             "island_points":{k:{kk:vv for kk,vv in t.items() if kk!="unit_sphere_xyz"}
                              for k,t in sorted(nodes.items(),key=lambda x:int(x[0]))}}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("ISLAND51_SOURCE_GEOCOORDS_VERIFIED",json.dumps({
          "n":len(nodes),"pair_distance_max_error_km":maxerror},sort_keys=True),flush=True)

if __name__=="__main__":main()
