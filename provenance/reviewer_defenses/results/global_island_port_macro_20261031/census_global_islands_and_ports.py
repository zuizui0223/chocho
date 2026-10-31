#!/usr/bin/env python3
"""GLOBAL island frame independent of subsequent butterfly indexing.

GSHHG high-resolution level-one polygon landmasses, >=20 km² and <1,000,000
km². Explicit area bounds omit continents. World Port Index is a source of
modern catalogued port locations only, NOT 2001-era ship traffic or imported
host plants. A catalog-absent island is not an island with no maritime access.

Before seeing any butterfly record, freeze all eligible islands and a
reproducible worldwide port/nonport size×latitude×longitude-balanced sample.
Original 17 focal islands are diagnostic only; never force-include them.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import random
import time
import zipfile
from collections import Counter,defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import URLError,HTTPError

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Point
from shapely.ops import transform, nearest_points
from pyproj import Transformer,Geod

HERE=Path(__file__).resolve().parent
GSHHG_SHA256="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
GSHHG_RELEASE="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
GEOD=Geod(ellps="WGS84")
SEED=20261031

def request(params):
    last=None
    for j in range(4):
        try:
            url=WPI+"?"+urlencode(params)
            with urlopen(Request(url,headers={"Accept":"application/json",
                "User-Agent":"chocho-global-island-port-geocensus/20261031"}),timeout=70) as f:
                data=json.load(f)
            if not isinstance(data,dict) or data.get("error"):
                raise ValueError("WPI source error "+str(data.get("error"))[:200])
            return data
        except (URLError,HTTPError,TimeoutError,ValueError) as exc:
            last=exc
            if isinstance(exc,HTTPError) and exc.code not in (408,425,429,500,502,503,504):
                break
            if j<3:time.sleep(1.5*2**j)
    raise RuntimeError("Cannot verify original WPI ports: "+repr(last))

def wpi():
    rows=[];offset=0
    while offset<12000:
        d=request({"where":"1=1","outFields":"OBJECTID,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
            "returnGeometry":"false","f":"json","orderByFields":"OBJECTID ASC",
            "resultOffset":offset,"resultRecordCount":1000})
        block=d.get("features",[])
        if not block:break
        rows.extend(x["attributes"] for x in block)
        offset+=len(block)
        if len(block)<1000 and not d.get("exceededTransferLimit",False):break
    accepted=[]
    for v in rows:
        try:
            pid=int(v["OBJECTID"]);lat=float(v["LATITUDE"]);lon=float(v["LONGITUDE"])
            if -90<=lat<=90 and -180<=lon<=180:
                accepted.append({"port_id":pid,"lon":lon,"lat":lat,
                    "name":str(v.get("PORT_NAME") or ""),
                    "country":str(v.get("COUNTRY") or "")})
        except (TypeError,KeyError,ValueError):
            continue
    if len(accepted)<3000 or len({r["port_id"] for r in accepted})!=len(accepted):
        raise ValueError("WPI source underfilled or duplicate port IDs")
    return accepted

def islands(zip_path,work,min_km2,max_km2):
    b=zip_path.read_bytes()
    digest=hashlib.sha256(b).hexdigest()
    if digest!=GSHHG_SHA256:raise ValueError("GSHHG original source SHA256 mismatch "+digest)
    work.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        matched=[p for p in z.namelist() if p.endswith((".shp",".shx",".prj",".dbf"))
                 and "/h/GSHHS_h_L1." in "/"+p]
        if len(matched)<3:raise ValueError("Missing original high-resolution L1 shapefile")
        for name in matched:z.extract(name,work)
    land=gpd.read_file(next(work.rglob("GSHHS_h_L1.shp")))
    if "sibling_id" not in land.columns:raise ValueError("Original landmass sibling_id missing")
    land["id"]=land["sibling_id"].astype("int64").astype(str)
    crs=Transformer.from_crs(4326,6933,always_xy=True)
    frame=[]
    polys=[]
    dropped=Counter()
    for iid,x in land.groupby("id",sort=False):
        geom=shapely.union_all(x.geometry.to_numpy())
        if geom.is_empty or not geom.is_valid:
            dropped["invalid"]+=1;continue
        area=transform(crs.transform,geom).area/1e6
        if area<min_km2:dropped["lt_min_area"]+=1;continue
        if area>max_km2:dropped["continent_or_gt_max_area"]+=1;continue
        center=geom.representative_point()
        if not(-180<=center.x<=180 and -90<=center.y<=90):
            dropped["invalid_geo"]+=1;continue
        frame.append({"id":iid,"area_km2":float(area),
            "centroid_lon":float(center.x),"centroid_lat":float(center.y),
            "gshhg_components":len(x),
            "port_index_count":0,"wpi_ids":""})
        polys.append(geom)
    if len(frame)<1200:
        raise ValueError("Suspiciously small GSHHG source island census "+str(len(frame)))
    return frame,polys,dropped,digest

def add_ports(frame,polys,ports,max_km):
    tree=shapely.STRtree(polys)
    counted=defaultdict(list)
    audit={"wpi_ports_total":len(ports),"ports_assigned_to_candidate_islands":0,
           "ports_ambiguous_near_coast":0,"ports_without_nearby_eligible_island":0}
    for p in ports:
        pt=Point(p["lon"],p["lat"])
        # Candidate-distance query done in degrees for bounding rectangles,
        # then actual nearest coastline assessed using WGS84 geodesic.
        dlat=max_km/110.5
        dlon=max_km/(111*max(.05,abs(math.cos(math.radians(p["lat"])))))
        possible=tree.query(shapely.box(p["lon"]-dlon,p["lat"]-dlat,
                                        p["lon"]+dlon,p["lat"]+dlat))
        matches=[]
        for k in possible:
            poly=polys[int(k)]
            if poly.covers(pt):
                distance=0.
            else:
                nearest=nearest_points(pt,poly)[1]
                _,_,meters=GEOD.inv(pt.x,pt.y,nearest.x,nearest.y)
                distance=meters/1000
            if distance<=max_km:
                matches.append((distance,int(k)))
        if not matches:
            audit["ports_without_nearby_eligible_island"]+=1;continue
        matches.sort()
        if len(matches)>1:
            audit["ports_ambiguous_near_coast"]+=1
        # Single closest candidate island. If a mainland (>1,000,000km²)
        # lies closer this can wrongly attribute a coast port to an island;
        # do not interpret counts as historical maritime access.
        counted[matches[0][1]].append(str(p["port_id"]))
        audit["ports_assigned_to_candidate_islands"]+=1
    for i,r in enumerate(frame):
        ids=sorted(set(counted[i]),key=int)
        r["port_index_count"]=len(ids)
        r["wpi_ids"]=";".join(ids)
    audit["candidate_islands_with_modern_wpi_ports"]=sum(x["port_index_count"]>0 for x in frame)
    return audit

def stratify(r):
    lat=r["centroid_lat"]
    lon=r["centroid_lon"]
    # 4 broad latitude strips × 4 longitude quadrants × 4 log-area bins.
    lat_bin=0 if lat < -30 else 1 if lat<0 else 2 if lat<30 else 3
    lon_bin=max(0,min(3,int((lon+180)/90)))
    area_bin=max(0,min(3,int(math.log10(r["area_km2"]/20)/1.1)))
    return (lat_bin,lon_bin,area_bin)

def choose(frame,n_port,n_nonport,seed=SEED):
    """Fix outcome-blind sampling within broad geography and island size.

    Random priority lists within each stratum, and use inverse-stratum-size
    weighted draws to preserve geographic dispersion instead of rich temperate
    harbor sampling. The source all-island census itself remains the true frame.
    """
    rng=random.Random(seed)
    sets={1:defaultdict(list),0:defaultdict(list)}
    for r in frame:
        key=int(r["port_index_count"]>0)
        sets[key][stratify(r)].append(r)
    selected={1:[],0:[]}
    for key,target in ((1,n_port),(0,n_nonport)):
        groups=sets[key]
        for v in groups.values():rng.shuffle(v)
        while len(selected[key])<target:
            available=sorted(s for s,g in groups.items() if g)
            if not available:break
            # Cycle once per occupied stratum in stable order, changing the
            # starting cell with deterministic randomized permutation.
            rng.shuffle(available)
            for stratum in available:
                if len(selected[key])==target:break
                if groups[stratum]:selected[key].append(groups[stratum].pop())
        if len(selected[key])<max(5,target//5):
            raise ValueError("Too few matched WPI harbor/nonharbor source islands")
    combined=selected[1]+selected[0]
    rng.shuffle(combined)
    for r in combined:r["selection_stratum"]=",".join(map(str,stratify(r)))
    return combined,{"selected_port":len(selected[1]),
         "selected_nonport":len(selected[0]),
         "not_fully_random_area_latitude_longitude_balanced":True,
         "future_butterfly_indices_not_used":True,
         "sample_seed":seed,
         "port_frame_strata":len(sets[1]),
         "nonport_frame_strata":len(sets[0])}

def write_csv(path,rows):
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--zip",type=Path,required=True)
    ap.add_argument("--work",type=Path,required=True)
    ap.add_argument("--max-area-km2",type=float,default=1_000_000)
    ap.add_argument("--min-area-km2",type=float,default=20)
    ap.add_argument("--near-wpi-km",type=float,default=3)
    ap.add_argument("--n-port",type=int,default=300)
    ap.add_argument("--n-nonport",type=int,default=300)
    args=ap.parse_args()
    frame,polys,dropped,sha=islands(args.zip,args.work,
          args.min_area_km2,args.max_area_km2)
    locations=wpi()
    audit=add_ports(frame,polys,locations,args.near_wpi_km)
    selected,sel=choose(frame,args.n_port,args.n_nonport)
    HERE.mkdir(parents=True,exist_ok=True)
    frame.sort(key=lambda x:int(x["id"]))
    selected.sort(key=lambda x:int(x["id"]))
    write_csv(HERE/"global_source_islands_20km2.csv",frame)
    write_csv(HERE/"global_outcome_blind_port_nonport_sample.csv",selected)
    receipt={
        "schema":"chocho_outcome_blind_global_GSHHG_WPI_island_port_frame_v1",
        "status":"SOURCE_VERIFIED_ISLAND_GEOGRAPHY_PORT_CENSUS_NOT_BUTTERFLY_COLONIZATION",
        "source_GSHHG":GSHHG_RELEASE,
        "source_GSHHG_sha256":sha,
        "source_WPI":WPI,
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "GSHHG_sibling_island_frame":len(frame),
        "gshhg_source_exclusions":dict(dropped),
        "minimum_km2":args.min_area_km2,
        "maximum_km2":args.max_area_km2,
        "catalogued_port_nearest_land_maximum_km":args.near_wpi_km,
        "port_source_quality":audit,
        "selection":sel,
        "sample_island_ids":[r["id"] for r in selected],
        "no_later_butterfly_record_used_in_design":True,
        "limitations":[
            "GSHHG source map of land areas is not Global Islands Database; reported island denominator may differ.",
            "Shoreline and WPI catalog locations are modern retrospective geography, not 2001–2010 vessel arrivals.",
            "Absence of WPI harbor record is not absence of historical ports, ships, human movement or planted hosts.",
            "Nearby continent/large mainland ports may be attributed to a nearby island; assignment is a location-based proxy.",
            "Island size and broad geographic spread are sample design covariates, not exact ecological propensity matching.",
            "Some very remote islands lack future butterfly records because of reporting, not actual ecological unoccupancy.",
            "No butterfly first record, HOSTS foodplant edge, historical alien status or colonization is inferred here."
        ]
    }
    (HERE/"global_island_port_frame_receipt.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"all_candidate_islands":len(frame),"sample":sel,
        "WPI":audit,"dropped":dict(dropped)},ensure_ascii=False,indent=2))

if __name__=="__main__":main()
