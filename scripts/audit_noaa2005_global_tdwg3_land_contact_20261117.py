#!/usr/bin/env python3
"""Resolve remaining original NOAA 2005 GSHHG raster contacts against ALL TDWG3 regions.

The previous source-locked source-only GSHHG vector audit classified 19 real-flow
coast contacts: 13 same source polygon, 3 a different polygon inside the original
source region and 3 on a different GSHHG polygon with unknown source/recipient region.
This script uses the entire checksum-pinned TDWG3 official map, with no outcomes.

No assertion of viable seed, port cargo transport, arrival year or butterfly colonisation.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

from pyproj import Geod
from shapely.geometry import Point,shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree

EXPECTED_GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
EXPECTED_TDWG_SHA="52da7828aba9d461dd133c27b3bd7a4407161f54"
THRESHOLD_NEAREST_TDWG_KM=5.0
THREE_FIXED=("AZO","CNY","CVI")
G=Geod(ellps="WGS84")

def source_unique(root,name):
    matches=list(root.rglob(name))
    if len(matches)!=1:raise RuntimeError(f"Expected unique archived {name}; got {len(matches)}")
    return matches[0]

def read_csv(f):
    with f.open(newline="",encoding="utf-8") as x:
        return list(csv.DictReader(x))

def read_tdwg(path):
    source=json.loads(path.read_text(encoding="utf-8"))
    polygons=[];codes=[];names=[]
    for f in source["features"]:
        p=f.get("properties") or {}
        cd=str(p.get("LEVEL3_COD") or "").strip()
        if not cd or not f.get("geometry"):continue
        shp=shape(f["geometry"])
        if shp.is_empty:continue
        polygons.append(shp);codes.append(cd)
        names.append(str(p.get("LEVEL3_NAM") or ""))
    if len(polygons)<200:
        raise RuntimeError(f"Too few source TDWG3 regions {len(polygons)}")
    if not set(THREE_FIXED).issubset(set(codes)):
        raise RuntimeError("Original three TDWG3 source islands absent")
    return polygons,codes,names,STRtree(polygons)

def km(point,other):
    return float(G.inv(point.x,point.y,other.x,other.y)[2])/1000.0

def label_point(point,polys,codes,names,tree):
    ids=sorted(set(int(i) for i in tree.query(point) if polys[int(i)].covers(point)))
    if ids:
        codes_match=sorted(set(codes[i] for i in ids))
        return {"source_status":"TDWG3_CONTAINS_ORIGINAL_GSHHG_LAND_POINT",
           "TDWG3_codes":codes_match,
           "TDWG3_names":sorted(set(names[i] for i in ids if names[i])),
           "exactly_one_region":len(codes_match)==1,
           "nearest_region_code":codes_match[0] if len(codes_match)==1 else None,
           "nearest_region_km":0.0}
    nearest=int(tree.nearest(point))
    projected=nearest_points(point,polys[nearest])[1]
    d=km(point,projected)
    return {
       "source_status":"NEAR_TDWG3_POLYGON_WITHIN_5KM_NOT_CONTAINED" if d<=THRESHOLD_NEAREST_TDWG_KM
                       else "NO_EXACT_TDWG3_REGION_AND_GEOGRAPHIC_UNCERTAINTY",
       "TDWG3_codes":[],"TDWG3_names":[],
       "exactly_one_region":False,
       "nearest_region_code":codes[nearest],
       "nearest_region_km":round(d,4)
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--original-source-dir",required=True,type=Path)
    parser.add_argument("--official-tdwg3-geojson",required=True,type=Path)
    parser.add_argument("--out-dir",required=True,type=Path)
    a=parser.parse_args()
    original=read_csv(source_unique(a.original_source_dir,"original_noaa2005_gshhg_land_contact_identity.csv"))
    audit=json.loads(source_unique(a.original_source_dir,"original_noaa2005_gshhg_land_contact_identity_audit.json").read_text())
    if (len(original)!=29 or audit["GSHHG_original_archive_sha256"]!=EXPECTED_GSHHG_SHA or
        audit["TDWG_original_git_commit"]!=EXPECTED_TDWG_SHA):
        raise RuntimeError("Original source geometry provenance mismatch")
    polys,codes,names,tree=read_tdwg(a.official_tdwg3_geojson)
    results=[]
    for r in original:
        pt=Point(float(r["nearest_vector_lon"]),float(r["nearest_vector_lat"]))
        region=label_point(pt,polys,codes,names,tree)
        original_codes=[x for x in r["contact_TDWG3_regions"].split("|") if x]
        if original_codes and sorted(original_codes)!=region["TDWG3_codes"]:
            raise RuntimeError(f"Full official TDWG3 geometry does not reproduce selected AZO/CNY/CVI overlay: {original_codes} vs {region['TDWG3_codes']}")
        orig=r["source"]
        resolved=region["exactly_one_region"]
        code=region["nearest_region_code"] if resolved else None
        classification=("WITHIN_ORIGIN_TDWG3" if resolved and code==orig else
                        "IN_OTHER_FIXED_TDWG3" if resolved and code in THREE_FIXED else
                        "IN_OTHER_GLOBAL_TDWG3_REGION" if resolved else
                        "NOT_EXACTLY_ASSIGNABLE_GLOBAL_TDWG3")
        results.append({
            "source":orig,"scenario":r["scenario"],
            "source_WPI_name":r["source_WPI_name"],
            "failure_day":int(r["failure_day"]),
            "GSHHG_polygon_native_id":r["contact_GSHHG_native_id"],
            "GSHHG_polygon_index":int(r["contact_GSHHG_polygon_index"]),
            "contact_nearest_vector_km":float(r["nearest_coast_distance_km"]),
            "GSHHG_original_identity":r["coast_identity_class"],
            "point_lon":pt.x,"point_lat":pt.y,
            "full_TDWG3_exact_codes":"|".join(region["TDWG3_codes"]),
            "full_TDWG3_names":"|".join(region["TDWG3_names"]),
            "TDWG3_match_status":region["source_status"],
            "nearest_TDWG3_code_if_not_exact":region["nearest_region_code"] if not resolved else "",
            "nearest_TDWG3_km_if_not_exact":region["nearest_region_km"] if not resolved else 0,
            "global_TDWG3_contact_class":classification,
            "is_biological_colonisation_record":False,
        })
    original_unknown=[
        r for r in results
        if r["scenario"]=="OBSERVED_NOAA_FLOW" and
           r["GSHHG_original_identity"]=="OTHER_GSHHG_POLYGON_WITH_AMBIGUOUS_TDWG_REGION"]
    if len(original_unknown)!=3:
        raise RuntimeError(f"Did not recover original exactly 3 ambiguous geography contacts: {len(original_unknown)}")
    diag={}
    for scenario in ("OBSERVED_NOAA_FLOW","MATHEMATICAL_REVERSED_FLOW"):
        pop=[r for r in results if r["scenario"]==scenario]
        diag[scenario]={
            "contacts":len(pop),
            "global_TDWG3_contact_classes":dict(Counter(r["global_TDWG3_contact_class"] for r in pop)),
            "global_TDWG3_match_statuses":dict(Counter(r["TDWG3_match_status"] for r in pop)),
            "different_TDWG3_exact_regions":[
                {"source":r["source"],"region":r["full_TDWG3_exact_codes"],
                 "day":r["failure_day"],"gshhg_poly":r["GSHHG_polygon_native_id"]}
                for r in pop if r["global_TDWG3_contact_class"]=="IN_OTHER_GLOBAL_TDWG3_REGION"
            ],
        }
    report={
        "schema":"chocho_original_noaa2005_gshhg_land_contact_full_TDWG3_region_resolution_v1",
        "status":"GLOBAL_TDWG3_SOURCE_ONLY_GEOGRAPHIC_RESOLUTION_NOT_SEED_TRANSPORT",
        "source_GSHHG_ID_CI":43906815380,
        "source_GSHHG_recheck_CI":43907339901,
        "source_particle_CI":43902616033,
        "original_GSHHG_SHA256":EXPECTED_GSHHG_SHA,
        "TDWG3_git_commit":EXPECTED_TDWG_SHA,
        "all_TDWG3_source_features":len(polys),
        "original_land_contacts":len(results),
        "original_ambiguous_actual_flow_contacts":len(original_unknown),
        "fixed_three_source_codes":THREE_FIXED,
        "source_success":True,
        "scenarios":diag,
        "original_ambiguous_resolved_or_still_unresolved":[
          {k:r[k] for k in ("source","failure_day","GSHHG_polygon_native_id",
           "contact_nearest_vector_km","full_TDWG3_exact_codes",
           "full_TDWG3_names","TDWG3_match_status",
           "nearest_TDWG3_code_if_not_exact","nearest_TDWG3_km_if_not_exact",
           "global_TDWG3_contact_class")} for r in original_unknown
        ],
        "limitations":[
          "TDWG3 boundaries represent plant geographic regions, often archipelagos and aggregates, not true source/recipient seed populations.",
          "If source point falls between GSHHG vector coastline and TDWG3 bounds, do not replace the official polygon by an invented region.",
          "Nearest TDWG3 polygon is descriptive geographic distance, not accepted as a location assignment when it does not contain the contact point.",
          "An identifiable contact on another TDWG3 region is only a modeled ocean-mask collision, never proof that viable propagules reached a different island.",
          "The historical dispersal time, original seed sources, real shipping networks, plant host usage, butterflies and colonization remain unobserved."
        ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"original_29_GSHHG_contacts_global_TDWG3_areas.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    (a.out_dir/"global_TDWG3_contact_source_gate.json").write_text(
         json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("GSHHG_GLOBAL_TDWG3_CONTACT_SOURCE_GATES",json.dumps({
       "status":report["status"],"tdwg_source_features":len(polys),
       "actual_flow":diag["OBSERVED_NOAA_FLOW"],
       "original_ambiguous":report["original_ambiguous_resolved_or_still_unresolved"]},
       ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
