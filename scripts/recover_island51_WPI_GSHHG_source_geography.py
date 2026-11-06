#!/usr/bin/env python3
"""Reconstruct original source GSHHG v2.3.7 landmass-ID -> NOAA WPI 2019
completed-port geometry, and independently verify all 51 island source
port-counts and 1,275 great-circle pair distances.

No synthetic island coordinates, no fitted biological response, no AIS proxy.
Original source identity: GSHHS_h_L1.sibling_id, original WPI OBJECTID.
Only allow geolocations for AIS/wind use when all source gates pass.

The original uncommitted outputs/port_buffer_landmass_audit.csv is unavailable;
this reconstructs it from original published geography and completed GBIF
port source, with several PRIOR geometric matching variants. Source Jaccards
are never used to decide which variant is correct.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from collections import Counter,defaultdict
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point,Polygon
from shapely.strtree import STRtree
from pyproj import Geod

R=6371.0088
G=Geod(ellps="WGS84")
REQUIRED_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"

def spherical_mean(locations):
    theta=np.radians(np.asarray(locations,dtype=float))
    lo=theta[:,0];la=theta[:,1]
    v=np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)]).mean(axis=1)
    norm=np.linalg.norm(v)
    if norm<1e-8:raise ValueError("Invalid original spherical average")
    return v/norm

def latlon(v):
    x,y,z=map(float,v)
    return [float(np.degrees(np.arctan2(y,x))),float(np.degrees(np.arctan2(z,np.hypot(x,y))))]

def arc(a,b):
    return float(R*np.arccos(np.clip(np.dot(a,b),-1,1)))

def is_num(v):
    try:return math.isfinite(float(v))
    except:return False

def approx_disk(lon,lat,meters=5000,n=36):
    # A true 5km geodesic circle, not 5km in degrees.
    angles=np.arange(n)*360/n
    ls,ts,_=G.fwd(np.repeat(lon,n),np.repeat(lat,n),angles,np.repeat(meters,n))
    return Polygon(zip(ls,ts))

def analyze(landfile,zipfile,geojsonfile,completefile,islandfile,pairfile):
    if hashlib.sha256(zipfile.read_bytes()).hexdigest()!=REQUIRED_SHA:
        raise ValueError("GSHHG 2.3.7 archived SHP original source sha256 changed")
    land=gpd.read_file(landfile)
    if "sibling_id" not in land.columns or "area" not in land.columns:
        raise ValueError("GSHHG polygons lack original sibling_id and area columns")
    if land.crs is None:raise ValueError("GSHHG source CRS missing")
    land=land.to_crs(4326)
    targets=pd.read_csv(islandfile,dtype={"landmass_id":str})
    targets=targets[targets.eligible.astype(str).str.casefold()=="true"].copy()
    if (len(targets),targets.landmass_id.nunique())!=(51,51):
        raise ValueError("Source 51 landmasses status changed")
    target_ids=set(targets.landmass_id)
    polygon_ids=land.sibling_id.astype(str)
    found=land[polygon_ids.isin(target_ids)].copy()
    coverage=set(found.sibling_id.astype(str))
    if len(coverage)!=51:
        raise ValueError("Original GSHHG sibling_id no longer maps to 51 archived island ids: "+str(sorted(target_ids-coverage)))
    ports=json.loads(geojsonfile.read_text(encoding="utf-8"))["features"]
    if len(ports)!=3630:raise ValueError("Original NGA WPI port count mismatch")
    loc={}
    for f in ports:
        i=str(int(f["properties"]["OBJECTID"]))
        loc[i]=tuple(map(float,f["geometry"]["coordinates"]))
    if len(loc)!=3630:raise ValueError("Nonunique original port ids")
    p=pd.read_csv(completefile,keep_default_na=False,dtype={"port_id":str})
    complete=set(p.loc[p.query_complete.astype(str).str.lower()=="true","port_id"])
    complete=complete.intersection(loc)
    if len(complete)<200:
        raise ValueError("No original completed WPI plant queries to reproduce island selection")

    # Only source 51 island shapes need refined geometry initially; this is
    # spatially efficient and the exact IDs are predefined by source.
    geom=list(found.geometry)
    code=found.sibling_id.astype(str).tolist()
    idx=STRtree(geom)
    target_hits=defaultdict(set)
    target_point_hits=defaultdict(set)
    # Reconstruct the ORIGINAL single-landmass gate by checking *all* GSHHG
    # L1 polygons, not just the subset of 51 target islands. Nearby small
    # islands and continental coastlines disqualify ambiguous 5km buffers.
    global_geoms=list(land.geometry)
    global_ids=land.sibling_id.astype(str).tolist()
    global_area=pd.to_numeric(land["area"],errors="coerce").to_numpy()
    global_tree=STRtree(global_geoms)
    globally_single_landmass=defaultdict(set)
    globally_single_island=defaultdict(set)
    target_with_point_and_unique_disk=defaultdict(set)
    target_with_unique_island_and_point=defaultdict(set)
    all_diagnostics=Counter()
    for port in sorted(complete,key=int):
        lon,lat=loc[port]
        disk=approx_disk(lon,lat)
        near=idx.query(disk)
        overlaps=[]
        for k in near:
            k=int(k)
            if shapely.intersects(geom[k],disk):
                overlaps.append(code[k])
        matches=set(overlaps)
        all_diagnostics["complete_checked"]+=1
        all_diagnostics[f"target_overlaps_{min(len(matches),3)}"]+=1
        for cid in matches:
            target_hits[cid].add(port)
        pt=Point(lon,lat)
        near_points=idx.query(pt)
        point_matches=set()
        for k in near_points:
            k=int(k)
            if shapely.intersects(geom[k],pt):
                point_matches.add(code[k])
                target_point_hits[code[k]].add(port)
        near_global=global_tree.query(disk)
        global_overlaps=set()
        global_island_overlaps=set()
        for k in near_global:
            k=int(k)
            if shapely.intersects(global_geoms[k],disk):
                global_overlaps.add(global_ids[k])
                if np.isfinite(global_area[k]) and 0<global_area[k]<7e6:
                    global_island_overlaps.add(global_ids[k])
        all_diagnostics[f"all_L1_5km_overlaps_{min(len(global_overlaps),4)}"]+=1
        all_diagnostics[f"small_islands_5km_overlaps_{min(len(global_island_overlaps),4)}"]+=1
        for cid in matches:
            if global_overlaps=={cid}:
                globally_single_landmass[cid].add(port)
                if cid in point_matches:
                    target_with_point_and_unique_disk[cid].add(port)
            if global_island_overlaps=={cid}:
                globally_single_island[cid].add(port)
                if cid in point_matches:
                    target_with_unique_island_and_point[cid].add(port)

    pairs=pd.read_csv(pairfile,dtype={"island_a":str,"island_b":str})
    if len(pairs)!=1275:raise ValueError("Expected 1275 original pair distances")
    original_counts=targets.set_index("landmass_id").ports.astype(int).to_dict()
    def variant(groups):
        matched={k:sorted(list(groups.get(k,set())),key=int) for k in sorted(target_ids,key=int)}
        counts={k:len(v) for k,v in matched.items()}
        mismatch={k:{"expected":int(original_counts[k]),"observed":counts[k],
                     "matched_port_ids":v[:60]} for k,v in matched.items()
                  if counts[k]!=original_counts[k]}
        centers={}
        for lid,portids in matched.items():
            if len(portids)>0:centers[lid]=spherical_mean([loc[k] for k in portids])
        errs=[]
        by_id=defaultdict(list)
        pairs_both_count_match=[]
        for row in pairs.itertuples(index=False):
            if row.island_a not in centers or row.island_b not in centers:
                continue
            expected=float(row.port_distance_km)
            got=arc(centers[row.island_a],centers[row.island_b])
            err=abs(expected-got)
            errs.append(err)
            by_id[row.island_a].append(err)
            by_id[row.island_b].append(err)
            if row.island_a not in mismatch and row.island_b not in mismatch:
                pairs_both_count_match.append(err)
        if not errs:errs=[float("inf")]
        top_islands=sorted((
            {"id":k,"mean_km":float(np.mean(v)),
             "max_km":max(v),"above_0p01km":sum(t>.01 for t in v)}
            for k,v in by_id.items()),key=lambda z:z["mean_km"],reverse=True)
        return {"count_matched_islands":51-len(mismatch),
             "distances_both_port_counts_matched":{
                "n_pairs":len(pairs_both_count_match),
                "max_error_km":max(pairs_both_count_match) if pairs_both_count_match else None,
                "n_pairs_error_over_0p01km":sum(x>.01 for x in pairs_both_count_match)
            },
             "top_distance_error_islands":top_islands[:15],
             "n_target_ports_overlaps":sum(counts.values()),
             "n_islands_with_coords":len(centers),
             "n_port_count_mismatch":len(mismatch),
             "port_count_mismatch":mismatch,
             "max_pair_distance_error_km":max(errs),
             "median_pair_distance_error_km":float(np.median(errs)),
             "valid_1275_source_distances":len(mismatch)==0 and len(errs)==1275 and max(errs)<=.01,
             "matched_port_id_lists":matched,
             "centroids":{k:latlon(v) for k,v in centers.items()}}

    # Algorithm candidates: point-contained vs all-target disk intersections.
    # This is purely a provenance comparison against previously frozen counts/
    # distances, NOT an exposure scan selecting on Jaccard or butterfly outcome.
    raw=variant(target_hits)
    contained=variant(target_point_hits)
    src={
      "schema":"chocho_WPI3630_GSHHG237_51island_provenance_reconstruction_v1",
      "status":"PENDING_STRICT_1275_DISTANCE_REPRODUCTION",
      "GSHHG_polygon_rows":len(land),
      "GSHHG_crs":str(land.crs),
      "GSHHG_original_sha256":REQUIRED_SHA,
      "target_51_polygon_row_count":len(found),
      "NOAA_WPI_2019_coordinate_ports":3630,
      "botanical_complete_ports_in_WPI":len(complete),
      "geometry_scan_counts":dict(all_diagnostics),
      "source_mapping_methods":{
          "geodesic_5km_buffer_intersects_target_landmass":raw,
          "WPI_point_intersects_target_landmass":contained,
          "geodesic_5km_intersects_exactly_one_global_landmass":variant(globally_single_landmass),
          "geodesic_5km_intersects_exactly_one_small_island":variant(globally_single_island),
          "WPI_point_inside_target_and_5km_exactly_one_global_landmass":variant(target_with_point_and_unique_disk),
          "WPI_point_inside_target_and_5km_exactly_one_small_island":variant(target_with_unique_island_and_point)},
      "source_51_landmass_ids":sorted(target_ids,key=int),
      "resistance_map_use_gate":"All original archived pairwise distances must match within 0.01km and all 51 completed source port counts. No imputed geography.",
      "source_code_note":"Uses GSHHG high-res L1 sibling_id and NOAA WPI OBJECTID independently of butterfly Jaccard or host response",
      "limits":["The original port-to-landmass CSV is not accessible; intersects-target alone may overinclude if 5km buffer contains multiple non-target landmasses.",
                "5km buffer circle is approximated by 36-point geodesic polygon; original exact buffer may have used another coastal projection.",
                "GSHHG dateline polygons are split; geodesic disk spanning ±180 degrees may require explicit antimeridian splitting.",
                "Source counts and distances validate reconstructed geometry; without exact agreement no ship/wind exposure must be fitted."]
    }
    passed=[name for name,v in src["source_mapping_methods"].items() if v["valid_1275_source_distances"]]
    src["passing_variants"]=passed
    if len(passed)==1:
        src["status"]="PASS_51_SOURCE_PORT_COUNTS_AND_ALL_1275_DISTANCES"
    elif len(passed)>1:
        candidates=[src["source_mapping_methods"][v]["matched_port_id_lists"] for v in passed]
        if all(c==candidates[0] for c in candidates):
            src["status"]="PASS_EQUIVALENT_GEOMETRIES_51_PORT_COUNTS_AND_DISTANCES"
        else:
            src["status"]="AMBIGUOUS_SOURCE_REPRODUCING_GEOMETRIES_DO_NOT_PROMOTE"
    return src


def main():
    a=argparse.ArgumentParser(description=__doc__)
    for field in ("gshhg-shp","gshhg-zip","wpi-geojson","plant-ports","islands","pairs","out"):
        a.add_argument("--"+field,required=True,type=Path)
    x=a.parse_args()
    try:
        result=analyze(x.gshhg_shp,x.gshhg_zip,x.wpi_geojson,
                       x.plant_ports,x.islands,x.pairs)
    except Exception as e:
        result={"status":"SOURCE_GEOMETRY_RECONSTRUCTION_ERROR",
                "exception_type":type(e).__name__,"reason":str(e)[:500]}
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,default=str)+"\n",encoding="utf-8")
    print("WPI51_GSHHG_SOURCE_GEOGRAPHY_RECONCILIATION",json.dumps({
        k:v for k,v in result.items() if k!="source_mapping_methods"},default=str)[:7000],flush=True)
    if "source_mapping_methods" in result:
        for name,v in result["source_mapping_methods"].items():
            print("WPI51_SOURCE_VARIANT",name,json.dumps({
                "matched_islands":v["count_matched_islands"],
                "n_port_count_mismatch":v["n_port_count_mismatch"],
                "max_geodesic_error_km":v["max_pair_distance_error_km"],
                "count_matched_pair_diagnostics":v["distances_both_port_counts_matched"],
                "top_distance_error_islands":v["top_distance_error_islands"][:6],
                "port_count_mismatch":v["port_count_mismatch"]
            },default=str)[:12000],flush=True)
    if result["status"].startswith("PASS_"):
        print("ISLAND51_ALL_1275_SOURCE_GREAT_CIRCLE_GEOMETRIES_VERIFIED",result["status"])
    else:
        print("ISLAND51_SOURCE_GEOGRAPHY_GATE_NOT_YET_PASSED",result["status"],flush=True)

if __name__=="__main__":main()
