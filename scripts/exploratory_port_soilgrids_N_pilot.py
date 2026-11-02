#!/usr/bin/env python3
"""Feasibility of ISRIC SoilGrids 2.0 predicted topsoil total nitrogen at
source-locked butterfly port and non-port coastal sites.

ISRIC SoilGrids is a machine-learning *estimate* of stock/concentration, not
measured sewage effluent, deposition or fertilization. Missing WCS pixels are
unknown. Spatial footprint circa 250m, depth 0–5cm.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import time
import urllib.parse
import urllib.request
from pathlib import Path

from pyproj import CRS, Transformer
from rasterio.io import MemoryFile

WCS="https://maps.isric.org/mapserv"
MAPFILE="/map/nitrogen.map"
LAYER="nitrogen_0-5cm_mean"
PROJ4="+proj=igh +datum=WGS84 +units=m +no_defs"
SOURCE_DOC="https://maps.isric.org/generated/soilgrids/nitrogen.html"


def request_one(lat:float,lon:float) -> dict:
    x,y=Transformer.from_crs("EPSG:4326",CRS.from_proj4(PROJ4),
                             always_xy=True).transform(lon,lat)
    if not math.isfinite(x) or not math.isfinite(y):
        raise ValueError("Non-finite interrupted Goode Homolosine coordinates")
    # Small enough to capture the target 250m pixel and next-neighbor
    # validation, but not imply average loading over 5km GBIF buffer.
    p=[
      ("map",MAPFILE),("SERVICE","WCS"),("VERSION","2.0.1"),
      ("REQUEST","GetCoverage"),("COVERAGEID",LAYER),("FORMAT","GEOTIFF_INT16"),
      ("SUBSET",f"X({x-450:.3f},{x+450:.3f})"),
      ("SUBSET",f"Y({y-450:.3f},{y+450:.3f})"),
      ("SUBSETTINGCRS","http://www.opengis.net/def/crs/EPSG/0/152160"),
      ("OUTPUTCRS","http://www.opengis.net/def/crs/EPSG/0/152160")
    ]
    url=WCS+"?"+urllib.parse.urlencode(p)
    req=urllib.request.Request(url,headers={
        "User-Agent":"chocho-soilgrids-nitrogen-feasibility/0.1",
        "Accept":"image/tiff,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=90) as response:
        media=response.headers.get("Content-Type","")
        data=response.read(2000000)
    if not data.startswith((b"II*\x00",b"MM\x00*")):
        return {"status":"NOT_GEOTIFF","media_type":media,
                "response_prefix":data[:250].decode("utf-8",errors="replace"),
                "lon":lon,"lat":lat,"x_igh":x,"y_igh":y}
    with MemoryFile(data) as mem,mem.open() as ds:
        raw = float(next(ds.sample([(x,y)],indexes=1,masked=False))[0])
        nodata = ds.nodata
        import numpy as np
        image=ds.read(1,masked=False)
        mask=(np.isfinite(image)&(image>0))
        if nodata is not None:
            mask=mask&(image!=nodata)
        valid=(math.isfinite(raw) and (nodata is None or raw!=nodata)
               and raw>0)
        value=raw/100.0 if valid else None  # cg/kg -> g/kg
        return {"status":("VALID_POSITIVE_PIXEL" if valid else
                          "ZERO_OR_MASKED_CANNOT_INTERPRET_AS_ZERO"),
                "lon":lon,"lat":lat,"x_igh":x,"y_igh":y,
                "topsoil_total_N_g_per_kg":value,
                "raw_pixel_cg_per_kg":raw,
                "positive_pixels_in_900m_crop":int(mask.sum()),
                "zero_pixels_in_900m_crop":int(np.count_nonzero(image==0)),
                "positive_pixel_mean_g_per_kg":float(image[mask].mean()/100.0) if mask.any() else None,
                "raster_source_crs":str(ds.crs),
                "raster_crop_bounds":list(map(float,ds.bounds)),
                "pixel_nodata_value":nodata,
                "crop_width":ds.width,"crop_height":ds.height,
                "source_layer":LAYER,"source_url":SOURCE_DOC,
                "machine_predicted_not_direct_soil_sample":True}


def run(geom_path:Path):
    j=json.loads(geom_path.read_text(encoding="utf-8"))
    rows=[r for r in j["rows"] if r["status"]=="PAIRED"]
    if len(rows)!=78 or len({r["island_id"] for r in rows})!=35:
        raise ValueError("Unrecognized WPI78/35 island source")
    countries={"GB","FR","DK","SE"}
    chosen=[]
    for c in sorted(countries):
        sites=[r for r in rows if r["country_original"]==c]
        if not sites:raise ValueError("No original WPI coastal samples in "+c)
        chosen.append(sorted(sites,key=lambda r:(int(r["island_id"]),r["port_id"]))[0])
    result=[]
    for r in chosen:
        for side in ("port","nonport"):
            lat=r["port_sample_lat"] if side=="port" else r["nonport_lat"]
            lon=r["port_sample_lon"] if side=="port" else r["nonport_lon"]
            try:
                s=request_one(float(lat),float(lon))
            except Exception as ex:
                s={"status":"SOURCE_REQUEST_FAILED","error":repr(ex)}
            result.append({"island_id":r["island_id"],"port_id":r["port_id"],
                           "country":r["country_original"],"side":side,**s})
            print("SOILGRIDS_N_SITE",r["country_original"],r["port_id"],side,
                  s.get("status"),s.get("topsoil_total_N_g_per_kg"),flush=True)
            time.sleep(.3)
    positive_controls=[]
    for name,lat,lon in (("Germany_inland",50.0,8.0),
                         ("France_inland",48.0,2.5),
                         ("Switzerland_inland",47.0,9.0)):
        try:
            q=request_one(lat,lon)
        except Exception as ex:
            q={"status":"SOURCE_REQUEST_FAILED","error":repr(ex)}
        positive_controls.append({"control_label":name,**q})
        print("SOILGRIDS_N_INLAND_CONTROL",name,q.get("status"),
              q.get("topsoil_total_N_g_per_kg"),
              "positive_crop_pixels",q.get("positive_pixels_in_900m_crop"),flush=True)
    return {
       "schema":"chocho_ISRIC_soilgrids_v2_N_WCS_8point_feasibility_v1",
       "status":"FEASIBILITY_ONLY_NO_WWTP_OR_FERTILIZATION_EFFECT",
       "geometry_sha256":hashlib.sha256(geom_path.read_bytes()).hexdigest(),
       "source_wcs":WCS,"source_layer":LAYER,
       "dataset_resolution_nominal_m":250,
       "model_depth":"0-5cm",
       "total_N_concentration_units":"g N per kg model-predicted fine earth",
       "sites":result,
       "valid_pixels":sum(x["status"]=="VALID_POSITIVE_PIXEL" for x in result),
       "inland_positive_controls":positive_controls,
       "positive_inland_controls":sum(x["status"]=="VALID_POSITIVE_PIXEL" for x in positive_controls),
       "strict_source_readiness":(
           "PASS_FOR_EXPANSION" if sum(x["status"]=="VALID_POSITIVE_PIXEL" for x in result)>=6
           and sum(x["status"]=="VALID_POSITIVE_PIXEL" for x in positive_controls)>=2
           else "NOT_READY_ZERO_VALUES_REQUIRE_DIAGNOSIS"),
       "limitations":[
           "SoilGrids mapped nitrogen is machine-predicted, not soil observations at the port.",
           "The 2020-style soil-property map does not resolve 2001-2010 versus 2011-2024 temporal nutrient loading.",
           "Small coastal sites may fall in null/masked raster pixels.",
           "Soil nutrient stock is not the input flux from wastewater, atmospheric deposition, or port infrastructure.",
           "No habitat matching, plant causal response, or butterfly life-history effect is tested."
       ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--geometry",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    r=run(a.geometry)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("SOILGRIDS_TOTAL_N_WCS_FEASIBILITY",
          json.dumps({"valid":r["valid_pixels"],"total":len(r["sites"]),
                      "controls":r["positive_inland_controls"],
                      "source_readiness":r["strict_source_readiness"],
                      "errors":[x for x in r["sites"] if x["status"]!="VALID_POSITIVE_PIXEL"]},
                     ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
