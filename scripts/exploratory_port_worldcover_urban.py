#!/usr/bin/env python3
"""ESA WorldCover v200 2021 geospatial audit for source-locked WPI ports.

Read the real 10m COG land-cover class pixels within 500m of both original
strict WPI78 port/non-port sampling centers, rather than an RGB WMS image.
Built-up is an OBSERVED LAND-COVER PROXY, not a shipping, historical, sewage
or nutrient exposure. Do not infer 2001-2010 urban conditions from a 2021 map.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import statistics
from contextlib import ExitStack
from collections import defaultdict
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window, from_bounds, intersection, transform as window_transform
from rasterio.errors import RasterioIOError

SOURCE="https://esa-worldcover.s3.eu-central-1.amazonaws.com"
VERSION="v200/2021/map"
CLASSES={10:"tree",20:"shrub",30:"grass",40:"cropland",50:"built",60:"bare",
         70:"snow_ice",80:"permanent_water",90:"wetland",95:"mangrove",100:"moss_lichen"}
TERRESTRIAL=set(CLASSES)-{80}
RADIUS_M=500
SOURCE_URL="https://esa-worldcover.org/en/data-access"


def tile_label(lon: float,lat: float) -> str:
    """3x3 geographic tiles indexed by the lower-left grid corner."""
    west=int(math.floor(lon/3)*3)
    south=int(math.floor(lat/3)*3)
    return f"{'N' if south>=0 else 'S'}{abs(south):02d}{'E' if west>=0 else 'W'}{abs(west):03d}"


def tile_url(tile: str) -> str:
    return f"{SOURCE}/{VERSION}/ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"


def window_tiles(lon:float,lat:float,radius_m:float) -> list[str]:
    # Local distance approximation accurate to <1% over <=1km near coasts.
    dy=radius_m/110900.0
    dx=radius_m/(111320.0*max(.10,math.cos(math.radians(lat))))
    result=set()
    for x in (lon-dx,lon+dx):
        for y in (lat-dy,lat+dy):
            result.add(tile_label(x,y))
    return sorted(result)


def disk_stat(src,lon:float,lat:float,radius_m:float) -> dict:
    if src.crs is None or src.crs.to_epsg()!=4326:
        raise ValueError("Expected native EPSG:4326 ESA WorldCover COG")
    if abs(src.transform.b)>1e-12 or abs(src.transform.d)>1e-12:
        raise ValueError("Rotated raster unsupported")
    dy=radius_m/110900.0
    dx=radius_m/(111320.0*max(.10,math.cos(math.radians(lat))))
    west,south,east,north=lon-dx,lat-dy,lon+dx,lat+dy
    if src.bounds.right <= west or src.bounds.left>=east or src.bounds.top<=south or src.bounds.bottom>=north:
        return {"counts":{},"circle_pixels":0,"tiles_overlap":False}
    bounds=(max(west,src.bounds.left),max(south,src.bounds.bottom),
            min(east,src.bounds.right),min(north,src.bounds.top))
    w=from_bounds(*bounds,transform=src.transform)
    c0=max(0,int(math.floor(w.col_off)))
    r0=max(0,int(math.floor(w.row_off)))
    c1=min(src.width,int(math.ceil(w.col_off+w.width)))
    r1=min(src.height,int(math.ceil(w.row_off+w.height)))
    if c1<=c0 or r1<=r0:
        return {"counts":{},"circle_pixels":0,"tiles_overlap":False}
    window=Window(c0,r0,c1-c0,r1-r0)
    arr=src.read(1,window=window)
    if arr.ndim!=2 or arr.size==0:raise ValueError("Empty raster COG response")
    yy,xx=np.ogrid[r0:r1,c0:c1]
    # Scale easting in metres and use 111.1 km latitude for this small area.
    pixel_lon=src.transform.c+(xx+.5)*src.transform.a
    pixel_lat=src.transform.f+(yy+.5)*src.transform.e
    circle=(
        ((pixel_lon-lon)*111320*math.cos(math.radians(lat)))**2
        +((pixel_lat-lat)*110900)**2 <= radius_m**2)
    in_circle=arr[circle]
    unique,frequency=np.unique(in_circle,return_counts=True)
    counts={int(a):int(b) for a,b in zip(unique,frequency)}
    return {"counts":counts,"circle_pixels":int(circle.sum()),"tiles_overlap":True}


class SourceTiles:
    def __init__(self):
        self.stack=ExitStack()
        self.cache={}
        self.errors={}
    def __enter__(self):return self
    def __exit__(self,*args):self.stack.close()
    def get(self,label):
        if label not in self.cache and label not in self.errors:
            try:
                src=self.stack.enter_context(rasterio.open(tile_url(label)))
                if (src.crs is None or src.crs.to_epsg()!=4326 or src.count!=1):
                    raise ValueError("Wrong raster dimensions or CRS "+label)
                self.cache[label]=src
                print("WORLDCOVER_OPEN_COG",label,src.width,src.height,flush=True)
            except (RasterioIOError,ValueError) as e:
                self.errors[label]=str(e)
        return self.cache.get(label)


def classify(site:dict, side:str,tiles:SourceTiles,radius_m:int) -> dict:
    lon=float(site["port_sample_lon"] if side=="port" else site["nonport_lon"])
    lat=float(site["port_sample_lat"] if side=="port" else site["nonport_lat"])
    expected=window_tiles(lon,lat,radius_m)
    counts=defaultdict(int)
    circle=0
    bad=[]
    for label in expected:
        src=tiles.get(label)
        if src is None:
            bad.append({"tile":label,"error":tiles.errors.get(label,"not available")[:300]})
            continue
        data=disk_stat(src,lon,lat,radius_m)
        if data["tiles_overlap"]:
            circle+=data["circle_pixels"]
            for code,n in data["counts"].items():counts[code]+=n
    valid_total=sum(n for code,n in counts.items() if code in TERRESTRIAL)
    invalid=sum(n for code,n in counts.items() if code not in TERRESTRIAL)
    result={
       "island_id":site["island_id"],"port_id":site["port_id"],
       "country":site["country_original"],"side":side,
       "latitude":lat,"longitude":lon,"radius_m":radius_m,
       "tile_names":expected,"source_tiles_missing":bad,
       "count_all_pixels_in_circle":circle,
       "count_terrestrial_nonwater_pixels":valid_total,
       "count_water_nodata_unknown_pixels":invalid,
       "nonwater_fraction_of_circle_pixels":valid_total/circle if circle else None,
       "class_pixel_counts":{str(k):counts.get(k,0) for k in sorted(CLASSES)},
       "built_fraction_among_land_pixels":counts.get(50,0)/valid_total if valid_total else None,
       "cropland_fraction_among_land_pixels":counts.get(40,0)/valid_total if valid_total else None,
       "tree_fraction_among_land_pixels":counts.get(10,0)/valid_total if valid_total else None,
       "grass_fraction_among_land_pixels":counts.get(30,0)/valid_total if valid_total else None,
       "tree_or_grass_fraction_among_land_pixels":(
           (counts.get(10,0)+counts.get(30,0))/valid_total if valid_total else None),
       "quality":"COMPLETE" if not bad and circle>=500 and valid_total>=200 else "INADEQUATE_OR_MISSING_COG",
    }
    if result["quality"]!="COMPLETE":
        for k in ("built_fraction_among_land_pixels","cropland_fraction_among_land_pixels",
                  "tree_fraction_among_land_pixels","grass_fraction_among_land_pixels",
                  "tree_or_grass_fraction_among_land_pixels"):
            result[k]=None
    return result


def focal_27_eligibility(aggregate:dict):
    if aggregate.get("site_period_facets")!=312 or aggregate.get("predeclared_source_pairs")!=78:
        raise ValueError("Focal allflora source mismatch")
    records=defaultdict(set)
    for r in aggregate["observations"]:
        if (r["plant_records"]>=30
            and r["N_scored_candidate_host_species"]>=2
            and r["N_scored_other_plant_species"]>=5):
            records[(str(r["island_id"]),int(r["port_id"]))].add((r["side"],r["period"]))
    out={k for k,vals in records.items() if len(vals)==4}
    if len(out)!=27 or len({k[0] for k in out})!=6:
        raise ValueError("Original six-island 27-pair source quality changed")
    return out


def run(geometry_path:Path,allflora_path:Path,*,radius_m=500,cohort="27"):
    geometry=json.loads(geometry_path.read_text(encoding="utf-8"))
    data=json.loads(allflora_path.read_text(encoding="utf-8"))
    rows=[x for x in geometry["rows"] if x["status"]=="PAIRED"]
    if (len(rows),len({r["island_id"] for r in rows}))!=(78,35):
        raise ValueError("Expected 78 paired source locations / 35 islands")
    eligible=focal_27_eligibility(data)
    if cohort=="27":
        rows=[r for r in rows if (str(r["island_id"]),int(r["port_id"])) in eligible]
        if len(rows)!=27:raise ValueError("Eligible source location mismatch")
    elif cohort!="78":raise ValueError("Unknown site population")
    sites=[]
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                      GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="1",
                      GDAL_HTTP_CONNECTTIMEOUT="20",GDAL_HTTP_TIMEOUT="60"):
        with SourceTiles() as tiles:
            for i,r in enumerate(rows):
                for side in ("port","nonport"):
                    v=classify(r,side,tiles,radius_m)
                    sites.append(v)
                    print("WORLD_COVER_SITE",r["island_id"],r["port_id"],side,
                          "built",v["built_fraction_among_land_pixels"],
                          "crop",v["cropland_fraction_among_land_pixels"],
                          "status",v["quality"],flush=True)
    pairs=[]
    cache={(str(r["island_id"]),int(r["port_id"]),r["side"]):r for r in sites}
    for r in rows:
        key=(str(r["island_id"]),int(r["port_id"]))
        p,n=cache[(key[0],key[1],"port")],cache[(key[0],key[1],"nonport")]
        usable=p["quality"]=="COMPLETE" and n["quality"]=="COMPLETE"
        pb=p["built_fraction_among_land_pixels"]; nb=n["built_fraction_among_land_pixels"]
        pairs.append({"island_id":key[0],"port_id":key[1],"country":r["country_original"],
             "port_built_fraction":pb,
             "nonport_built_fraction":nb,
             "port_minus_nonport_built_fraction":pb-nb if usable else None,
             "port_minus_nonport_cropland_fraction":(
                 p["cropland_fraction_among_land_pixels"]-n["cropland_fraction_among_land_pixels"] if usable else None),
             "nonport_coast_urban_atleast20pct_built":nb>=.20 if usable else None,
             "nonport_coast_rural_under5pct_built":nb<.05 if usable else None,
             "usable_both_sides":usable})
    good=[x for x in pairs if x["usable_both_sides"]]
    def summary(s):
        vals=[x["port_minus_nonport_built_fraction"] for x in s if x["usable_both_sides"]]
        return {"n":len(vals),"median":statistics.median(vals) if vals else None,
                "mean":statistics.mean(vals) if vals else None}
    return {
      "schema":"chocho_ESA_WorldCover_2021_10m_WPI_coast_landcover_v1",
      "status":"STATIC_2021_URBAN_LANDCOVER_DESCRIPTIVE_NO_HISTORICAL_PORT_CAUSAL_INFERENCE",
      "source_worldcover":SOURCE_URL,"source_data_url":SOURCE+"/"+VERSION+"/",
      "source_year":2021,"source_version":"v200","raster_native_crs":"EPSG:4326",
      "landcover_buffer_radius_m":radius_m,
      "all_values_are_pixel_fraction_among_nonwater_terrestrial_classes":True,
      "cohort":cohort,"n_target_pairs":len(rows),"n_islands_target":len({r["island_id"] for r in rows}),
      "geometry_sha256":hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
      "allflora_sha256":hashlib.sha256(allflora_path.read_bytes()).hexdigest(),
      "summary":{"pairs_with_both_sites_valid":len(good),
                 "usable_islands":len({r["island_id"] for r in good}),
                 "ports_built_mean":statistics.mean(x["port_built_fraction"] for x in good) if good else None,
                 "nonports_built_mean":statistics.mean(x["nonport_built_fraction"] for x in good) if good else None,
                 "nonport_urban_ge20percent_count":sum(x["nonport_coast_urban_atleast20pct_built"] for x in good),
                 "nonport_rural_lt5percent_count":sum(x["nonport_coast_rural_under5pct_built"] for x in good),
                 "port_more_built_count":sum(x["port_minus_nonport_built_fraction"]>0 for x in good),
                 "port_minus_nonport_built_fraction":summary(good)},
      "paired_sites":pairs,"site_buffers":sites,
      "limitations":[
        "Static ESA WorldCover 2021 cannot be interpreted as urban state in 2001-2010.",
        "ESA WorldCover classification is an estimated 10m map with site- and class-specific classification error.",
        "The 500m pixel disc is different from the 5km GBIF botanical search radius and was chosen to describe local landscape context, not exactly sampled flora habitat.",
        "Coastal land-cover denominators exclude water and nodata, which are reported separately.",
        "A 2021 urban-land-cover association cannot identify a nutrient-loading or historical shipping effect.",
        "Non-port coastal controls may not include adequately urbanized examples for separating port infrastructure from urbanization.",
        "Neither host larvae nor plant abundance was measured."
      ],
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--geometry",type=Path,required=True)
    ap.add_argument("--allflora",type=Path,required=True)
    ap.add_argument("--cohort",choices=("27","78"),default="27")
    ap.add_argument("--radius-m",type=int,choices=(250,500),default=500)
    ap.add_argument("--output",type=Path,required=True)
    q=ap.parse_args()
    result=run(q.geometry,q.allflora,radius_m=q.radius_m,cohort=q.cohort)
    q.output.parent.mkdir(parents=True,exist_ok=True)
    q.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("WORLDCOVER_2021_MATCHED_COAST_GATES",json.dumps(result["summary"],sort_keys=True),flush=True)


if __name__=="__main__":
    main()
