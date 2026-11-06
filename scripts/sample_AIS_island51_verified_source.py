#!/usr/bin/env python3
"""Georeference actual compressed World Bank commercial AIS overview pixels to
the exactly recovered original 51 island mean NGA WPI2019 port locations.

The .tif.ovr source has no georeference of its own; derive its 2× overview
geotransform from the accompanying original .tfw world file, validate full
grid geometry and EPSG:4326 assumptions from dataset metadata.

This is a first geospatial pixel-source measurement, NOT yet path density,
circuit resistance, historical shipping route or plant-transport evidence.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import rasterio
from rasterio.windows import Window
from rasterio.transform import Affine

EXPECTED_ISLAND_STATUS="ALL_51_SPHERICAL_WPI_CENTROIDS_REPRODUCE_ORIGINAL_1275_DISTANCES"
SOURCE_ODATA="World Bank/IMF commercial AIS hourly position density Jan2015-Feb2021"


def worldfile_transform(meta,src_width,src_height):
    tfw=meta.get("original_source_world_file","")
    text=[x.strip() for x in tfw.splitlines() if x.strip()]
    if len(text)!=6:raise ValueError("Unexpected source worldfile number of lines")
    a,d,b,e,c,f=map(float,text)
    if not math.isclose(a,.005,abs_tol=.0000001) or \
       not math.isclose(e,-.005,abs_tol=.0000001) or max(abs(d),abs(b))>1e-10:
        raise ValueError("Original source 0.005 degree unrotated WGS84 geotransform not verified "+str(text))
    # The .ovr TIFF first overview of the original 0.005-degree base has
    # precisely half source width/height; the same geographic bounding box.
    # Worldfile c,f are *center* of the first native pixel, whereas GDAL's
    # affine origin is the upper-left *outer corner* of the source.
    overview=Affine(a*2,0,c-a/2,0,e*2,f-e/2)
    bottom=f-e/2+src_height*e*2
    right=c-a/2+src_width*a*2
    if abs((right-(c-a/2))-360.03)>.05 or not 168<(f-e/2-bottom)<172:
        raise ValueError("AIS overview global dimensions not compatible with original 0.005 degree worldfile")
    if abs(c+180)>1 or abs(f-85)>2:
        raise ValueError("AIS original global raster worldfile origin unexpected for worldwide WGS84")
    return overview,{"source_native_dlon_deg":a,"source_native_dlat_deg":e,
       "overview_dlon_deg":a*2,"overview_dlat_deg":e*2,
       "left":c-a/2,"right":right,
       "top":f-e/2,"bottom":bottom,
       "original_worldfile":["{:g}".format(v) for v in (a,d,b,e,c,f)]}


def sample_port_area(ds,transform,lon,lat,radius_km=5.0):
    left,top=transform.c,transform.f
    xpix=float(transform.a);ypix=-float(transform.e)
    cx=int(math.floor((lon-left)/xpix))
    cy=int(math.floor((top-lat)/ypix))
    if not (0<=cx<ds.width and 0<=cy<ds.height):
        return {"status":"OUTSIDE_AIS_SOURCE_GRID","n_valid":0}
    rx=max(1,int(math.ceil(radius_km/(111.320*math.cos(math.radians(lat))*xpix))))
    ry=max(1,int(math.ceil(radius_km/(110.900*ypix))))
    c0=max(0,cx-rx);c1=min(ds.width,cx+rx+1)
    r0=max(0,cy-ry);r1=min(ds.height,cy+ry+1)
    w=Window(c0,r0,c1-c0,r1-r0)
    arr=ds.read(1,window=w,masked=True)
    rr,cc=np.indices(arr.shape)
    cell_lon=left+(cc+c0+.5)*xpix
    cell_lat=top-(rr+r0+.5)*ypix
    disk=((cell_lon-lon)*(111.32*math.cos(math.radians(lat))))**2+\
         ((cell_lat-lat)*110.90)**2<=radius_km**2
    mask=~np.ma.getmaskarray(arr) & disk
    vals=np.asarray(arr.data,dtype=float)[mask]
    if len(vals)==0:return {"status":"SOURCE_ALL_NODATA_WITHIN_PORT_5KM","n_valid":0}
    return {
      "status":"AIS_SOURCE_VALID",
      "n_valid":int(len(vals)),"n_positive":int(np.count_nonzero(vals>0)),
      "positive_fraction":float(np.mean(vals>0)),
      "mean_AIS_broadcast_positions_per_overview_cell":float(np.mean(vals)),
      "median_AIS_broadcast_positions_per_overview_cell":float(np.median(vals)),
      "mean_log1p_AIS_positions":float(np.mean(np.log1p(vals))),
      "p90_AIS_positions":float(np.quantile(vals,.9)),
      "max_AIS_positions":float(np.max(vals)),
      "observed_overview_sample_cell_resolution_deg":xpix
    }


def run(source,verified,overview):
    meta=json.loads(source.read_text())
    if meta.get("status")!="ORIGINAL_COMMERCIAL_AIS_OVERVIEW_GDAL_PIXEL_VALUES_READ_AND_CRC_VERIFIED":
        raise ValueError("Original verified AIS raster receipt missing")
    z=json.loads(verified.read_text())
    if z.get("status")!=EXPECTED_ISLAND_STATUS or len(z.get("verified_source_51_mean_port_points",{}))!=51:
        raise ValueError("Archived 51 exact NOAA WPI island source coordinates not verified")
    if z.get("archive_original_1275_recovered_distance_max_km",1)>1e-6:
        raise ValueError("Archival distances fail exact source gate")
    with rasterio.open(overview) as ds:
        if (ds.width,ds.height,ds.count,ds.dtypes[0])!=(36003,16999,1,"int32"):
            raise ValueError("Original World Bank AIS overview dimensions/format changed")
        t,diag=worldfile_transform(meta,ds.width,ds.height)
        observations=[]
        for k,v in sorted(z["verified_source_51_mean_port_points"].items(),key=lambda x:int(x[0])):
            lon,lat=v["longitude"],v["latitude"]
            s=sample_port_area(ds,t,lon,lat)
            observations.append({"landmass_id":k,"longitude":lon,"latitude":lat,
                  "n_source_ports":v["n_complete_single_landmass_ports"],**s})
    return {
       "schema":"chocho_WPI2019_island51_WorldBankAIS2015_2021_georeferenced_port_neighborhoods_v1",
       "status":"VERIFIED_REAL_WORLDBANK_AIS_OVERVIEW_ISLAND51_LOCAL_PIXELS_NOT_ISLAND_TO_ISLAND_TRANSPORT",
       "source_WPI51_run":42173654622,"source_AIS_real_raster_run":42174287282,
       "source_WPI51_json_sha256":hashlib.sha256(verified.read_bytes()).hexdigest(),
       "source_AIS_overview_schema_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
       "WorldBank_AIS_date":"2015-01 to 2021-02",
       "WorldBank_AIS_cell_semantics":"Count of AIS positions per 0.01 degree first-level raster overview; NOT independent vessels or trips",
       "raster_georeference":diag,
       "n_target_original_islands":51,
       "n_valid_5km_port_neighbor_islands":sum(x["status"]=="AIS_SOURCE_VALID" for x in observations),
       "n_islands_with_some_AIS_positions":sum(x.get("n_positive",0)>0 for x in observations),
       "port_5km_AIS_observations":observations,
       "not_measured":"No inter-island AIS path metric, no wind vector and no Circuitscape effective ship-route resistance yet",
       "caution":"Global raster is modern shipping density, not a source of 2001–2024 introduced plants or historical colonial trade."
     }


def main():
    p=argparse.ArgumentParser()
    for name in ("source","islands","overview","out"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    output=run(a.source,a.islands,a.overview)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("WORLD_BANK_AIS_VERIFIED_51_ISLAND_PORT_DENSITY",json.dumps({
        "status":output["status"],
        "georeference":output["raster_georeference"],
        "n_islands_with_AIS_raster_coverage":output["n_valid_5km_port_neighbor_islands"],
        "n_islands_with_positive_AIS":output["n_islands_with_some_AIS_positions"],
        "sampled_origins":output["port_5km_AIS_observations"][:8]
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
