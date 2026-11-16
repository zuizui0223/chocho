#!/usr/bin/env python3
"""Independent coast-identity audit of the original NOAA 2005 81 tracks.

Uses the SAME SHA-pinned GSHHG 2.3.7 intermediate L1 coastline that
produced the 0.0625deg land mask, and frozen official TDWG level-3 polygons.

Do not equate first GSHHG raster land contact to verified seed landfall.
Only match the nearest authentic polygon within 8km of the raster contact
point; preserve coarse-mask mismatch and ambiguous origin identity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import geopandas as gpd
import shapely
from pyproj import Geod
from shapely.geometry import Point,shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree

ARCHIVE_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
TDWG_SHA="52da7828aba9d461dd133c27b3bd7a4407161f54"
ORIGINAL_SOURCE_RUN=43902616033
ORIGINAL_CONTACT_RUN=43903614189
ORIGINAL_WPI_REGIONS=("AZO","CNY","CVI")
CONTACT_MAX_VECTOR_DISTANCE_KM=8.
ANCHOR_MAX_VECTOR_DISTANCE_KM=20.
WGS=Geod(ellps="WGS84")

def checksum(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(2**20),b""):
            h.update(buf)
    return h.hexdigest()

def unique(path,name):
    z=list(path.rglob(name))
    if len(z)!=1:
        raise RuntimeError(f"Expected one original artifact {name}; found {len(z)}")
    return z[0]

def read_original(root):
    original=json.loads(unique(root,"particles_diagnostic_tracks.json").read_text())
    manifest=json.loads(unique(root,"lagrangian_source_quality.json").read_text())
    if (len(original)!=81 or manifest["status"]!=
            "ACTUAL_2005_LEGACY_OSCAR_MIDPOINT_LAGRANGIAN_SOURCE_PILOT_NOT_SEED_ARRIVAL"):
        raise RuntimeError("Original 81 source trajectories were not recovered")
    expected={
       "OBSERVED_NOAA_FLOW":{"GSHHG_FINE_MASK_LAND_CONTACT":19,"COMPLETED_90_DAYS":7,
                             "SOURCE_VELOCITY_MISSING":1},
       "MATHEMATICAL_REVERSED_FLOW":{"GSHHG_FINE_MASK_LAND_CONTACT":10,"COMPLETED_90_DAYS":17},
       "STATIONARY_NO_CURRENT":{"COMPLETED_90_DAYS":27},
    }
    for name,counts in expected.items():
        got=dict(Counter(r["status"] for r in original if r.get("scenario")==name))
        if got!=counts:raise RuntimeError(f"Source {name} status drift: {got}")
    return original,manifest

def load_coast(gshhg_zip):
    sha=checksum(gshhg_zip)
    if sha!=ARCHIVE_SHA:raise RuntimeError(f"Original GSHHG zip checksum mismatch: {sha}")
    file=f"zip://{gshhg_zip.resolve()}!GSHHS_shp/i/GSHHS_i_L1.shp"
    d=gpd.read_file(file)
    if d.crs is None:raise RuntimeError("No GSHHG coordinate reference system")
    d=d.to_crs("EPSG:4326")
    if len(d)<32000 or len(d)>34000:raise RuntimeError(f"GSHHG intermediate L1 source population changed: {len(d)}")
    gg=list(d.geometry)
    if any(g is None or g.is_empty for g in gg):
        raise RuntimeError("Unexpected empty GSHHG polygon")
    # Store native feature index as an unambiguous identifier regardless of external IDs.
    ids=list(d["id"]) if "id" in d.columns else list(range(len(d)))
    if len(set(map(str,ids)))!=len(ids):
        raise RuntimeError("Duplicate GSHHG polygon id")
    sibling=list(d["sibling_id"]) if "sibling_id" in d.columns else [None]*len(d)
    return gg,ids,sibling,STRtree(gg),sha

def load_tdwg(geojson):
    d=json.loads(geojson.read_text(encoding="utf8"))
    selected={}
    for f in d["features"]:
        code=str(f.get("properties",{}).get("LEVEL3_COD") or "").strip()
        if code in ORIGINAL_WPI_REGIONS and f.get("geometry"):
            if code in selected:raise RuntimeError(f"Duplicate TDWG3 {code}")
            selected[code]=shape(f["geometry"])
    if set(selected)!=set(ORIGINAL_WPI_REGIONS):
        raise RuntimeError("Required official TDWG3 island region polygons missing")
    return selected

def dist_m(a,b):
    lo1,la1=float(a.x),float(a.y)
    lo2,la2=float(b.x),float(b.y)
    return float(WGS.inv(lo1,la1,lo2,la2)[2])

def closest_polygon(pt,gg,ids,siblings,tree):
    # Query exact intersect first to protect small polygon selection.
    idxs=[int(i) for i in tree.query(pt) if gg[int(i)].covers(pt)]
    if idxs:
        idx=min(idxs,key=lambda i:gg[i].area)
        nearest=pt
        method="VECTOR_POLYGON_CONTAINS_POINT"
    else:
        idx=int(tree.nearest(pt))
        nearest=nearest_points(pt,gg[idx])[1]
        method="NEAREST_GSHHG_VECTOR_POLYGON_NOT_CONTAINED"
    return {"poly_index":idx,"gshhg_polygon_id":str(ids[idx]),
        "gshhg_sibling_id":str(siblings[idx]) if siblings[idx] is not None else None,
        "vector_match":method,"dist_km":round(dist_m(pt,nearest)/1000,4),
        "nearest_vector_lat":float(nearest.y),"nearest_vector_lon":float(nearest.x)}

def tdwg_regions_for_point(p,tdwg):
    return sorted(c for c,g in tdwg.items() if g.covers(p))

def classify(origin,hit,origin_code,hit_region):
    if origin["dist_km"]>ANCHOR_MAX_VECTOR_DISTANCE_KM or hit["dist_km"]>CONTACT_MAX_VECTOR_DISTANCE_KM:
        return "UNRESOLVED_VECTOR_OR_PORT_COAST_MATCH"
    if origin["poly_index"]==hit["poly_index"]:
        return "SAME_ORIGIN_GSHHG_POLYGON"
    if (origin["gshhg_sibling_id"] is not None and
        origin["gshhg_sibling_id"]!="None" and
        origin["gshhg_sibling_id"]==hit["gshhg_sibling_id"]):
        return "SAME_GSHHG_SIBLING_DIFFERENT_POLYGON"
    if len(hit_region)!=1:
        return "OTHER_GSHHG_POLYGON_WITH_AMBIGUOUS_TDWG_REGION"
    if hit_region[0]==origin_code:
        return "OTHER_GSHHG_POLYGON_IN_ORIGIN_TDWG3_REGION"
    if hit_region[0] in ORIGINAL_WPI_REGIONS:
        return "GSHHG_POLYGON_IN_OTHER_FIXED_TDWG3_REGION"
    return "OTHER_GSHHG_POLYGON_OUTSIDE_THREE_FIXED_TDWG3_REGIONS"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-tracks",type=Path,required=True)
    p.add_argument("--gshhg-zip",type=Path,required=True)
    p.add_argument("--tdwg-level3-geojson",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    original,manifest=read_original(a.original_tracks)
    gg,ids,siblings,tree,archive_sha=load_coast(a.gshhg_zip)
    tdwg=load_tdwg(a.tdwg_level3_geojson)
    fixed=manifest["fixed_origins"]
    ports={}
    for code in ORIGINAL_WPI_REGIONS:
        port=fixed[code]
        pt=Point(float(port["longitude"]),float(port["latitude"]))
        g=closest_polygon(pt,gg,ids,siblings,tree)
        regions=tdwg_regions_for_point(pt,tdwg)
        if regions!=[code]:
            raise RuntimeError(f"Original anchored WPI port is not in original official TDWG: {code}, {regions}")
        ports[code]={**g,
            "port_point":[pt.x,pt.y],"name":port["name"],
            "WPI_id":port["WPI_id"],
            "origin_polygon_identity_credible":g["dist_km"]<=ANCHOR_MAX_VECTOR_DISTANCE_KM}

    collisions=[]
    for r in original:
        if r.get("status")!="GSHHG_FINE_MASK_LAND_CONTACT":continue
        row=r["trajectory"][-1]
        pt=Point(float(row["longitude"]),float(row["latitude"]))
        hit=closest_polygon(pt,gg,ids,siblings,tree)
        site=Point(hit["nearest_vector_lon"],hit["nearest_vector_lat"])
        region=tdwg_regions_for_point(site,tdwg)
        src=r["source"]
        label=classify(ports[src],hit,src,region)
        collisions.append({
            "source":src,"scenario":r["scenario"],
            "source_WPI_id":ports[src]["WPI_id"],
            "source_WPI_name":ports[src]["name"],
            "start_time":r["start_time_utc"],"launch_id":r["launch_replicate_id"],
            "NOAA_release_index":r["NOAA_source_time_index"],
            "failure_day":r["failure_day"],
            "contact_lat":float(pt.y),"contact_lon":float(pt.x),
            "nearest_vector_lat":hit["nearest_vector_lat"],
            "nearest_vector_lon":hit["nearest_vector_lon"],
            "nearest_coast_distance_km":hit["dist_km"],
            "nearest_coast_method":hit["vector_match"],
            "contact_TDWG3_regions":"|".join(region),
            "origin_GSHHG_polygon_index":ports[src]["poly_index"],
            "contact_GSHHG_polygon_index":hit["poly_index"],
            "contact_GSHHG_native_id":hit["gshhg_polygon_id"],
            "coast_identity_class":label,
            "is_external_destination_verified":False,
            "interpretation":"RASTER_LAND_MASK_CONTACT_ONLY_NOT_SEED_ARRIVAL",
        })
    if len(collisions)!=29:
        raise RuntimeError(f"Expected original 19+10 mask contacts; found {len(collisions)}")
    if Counter(r["scenario"] for r in collisions)!={
       "OBSERVED_NOAA_FLOW":19,"MATHEMATICAL_REVERSED_FLOW":10}:
        raise RuntimeError("Original land contact scenario counts no longer reproduce")

    diagnostics={}
    for scen in ("OBSERVED_NOAA_FLOW","MATHEMATICAL_REVERSED_FLOW"):
        rs=[r for r in collisions if r["scenario"]==scen]
        diagnostics[scen]={
            "original_contacts":len(rs),
            "class_counts":dict(sorted(Counter(r["coast_identity_class"] for r in rs).items())),
            "vector_coast_point_contained":sum(r["nearest_coast_method"]=="VECTOR_POLYGON_CONTAINS_POINT" for r in rs),
            "vector_nearest_but_not_contained":sum(r["nearest_coast_method"]!="VECTOR_POLYGON_CONTAINS_POINT" for r in rs),
            "within_8km_nearest_coast":sum(r["nearest_coast_distance_km"]<=CONTACT_MAX_VECTOR_DISTANCE_KM for r in rs),
            "with_other_fixed_destination_region":sum(
               r["coast_identity_class"]=="GSHHG_POLYGON_IN_OTHER_FIXED_TDWG3_REGION" for r in rs),
        }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"original_noaa2005_gshhg_land_contact_identity.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(collisions[0]));w.writeheader();w.writerows(collisions)
    audit={
       "schema":"chocho_original_NOAA2005_GSHHG_i_2p3p7_landmass_contact_identity_v1",
       "status":"ORIGINAL_GSHHG_POLYGON_SOURCE_RECONTACT_CLASSIFICATION_NOT_BIOLOGICAL_ISLAND_ARRIVAL",
       "source_particle_CI":ORIGINAL_SOURCE_RUN,
       "source_recontact_CI":ORIGINAL_CONTACT_RUN,
       "source_original_tracks":len(original),
       "GSHHG_original_archive_sha256":archive_sha,
       "GSHHG_intermediate_L1_polygons":len(gg),
       "TDWG_original_git_commit":TDWG_SHA,
       "original_WPI_coastal_anchor_polygon_gates":ports,
       "source_contamination_notice":"Ocean land contact is a raster-to-vector approximate hit, not proof a real seed landed on a particular shore",
       "coast_match_distance_threshold_km":CONTACT_MAX_VECTOR_DISTANCE_KM,
       "anchor_match_distance_threshold_km":ANCHOR_MAX_VECTOR_DISTANCE_KM,
       "diagnostics":diagnostics,
       "fixed_destinations_crossed_biological_0":False,
       "rules_and_limits":[
         "Never count GSHHG ocean land-mask crossing as actual viable seed landfall or butterfly establishment.",
         "If coarse mask first-land pixel is >8 km from vector coast, report unresolved rather than assigning landmass.",
         "Same GSHHG sibling ID can group multiple polygons; report sibling condition separately, not identical island identity.",
         "Source origin WPI coastal port identity is assigned only if <=20km from source GSHHG vector polygon.",
         "TDWG3 regional polygons contain archipelagos and are not island-level biological population boundaries.",
         "A fixed port coordinate in a destination region is not the entire coastline or historical vessel route.",
         "Classifying a hit in another fixed TDWG3 region does not establish that it was a viable seed dispersal event.",
         "Tracking scenarios differ by mathematical assumptions, not independent ecological replication.",
       ]
    }
    (a.out_dir/"original_noaa2005_gshhg_land_contact_identity_audit.json").write_text(
       json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("GSHHG_EXACT_COAST_IDENTITY_GATE",json.dumps({
      "status":audit["status"],"source_origin_anchor_credible":{
        c:x["origin_polygon_identity_credible"] for c,x in ports.items()},
      "scenarios":diagnostics,"source_polygon_count":len(gg)},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
