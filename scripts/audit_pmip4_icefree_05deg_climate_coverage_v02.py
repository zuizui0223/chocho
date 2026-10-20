#!/usr/bin/env python3
"""Audit if the ORIGINAL CHELSA 0/21ka rasters cover PMIP4 ICE6G_C paleo-land.

0.5 degree source-only resolution check. No plant model predictions/genotypes.
Only if 21ka climate exists on land newly emerged relative to 0ka can a grid
corridor across those shelves be interpreted without invented interpolation.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.transform import from_origin

VARS=("bio01","bio07","bio12","bio15")
MD5={0:"5ba8c0df48a12132238ffcbf675578e0",21:"9443668583800099ce9e865be90e4ad6"}

def pmip0p5(src):
    n=Path(src).name
    yr=0 if n.endswith(".0.nc") else 21 if n.endswith(".21.nc") else None
    if yr is None:raise ValueError("only exact PMIP4 0 or 21ka inputs")
    if hashlib.md5(Path(src).read_bytes()).hexdigest()!=MD5[yr]:
        raise ValueError("PMIP4 source MD5 mismatch")
    with netcdf_file(src,mmap=False) as ds:
        lat=np.asarray(ds.variables["lat"][:],float)
        lon=np.asarray(ds.variables["lon"][:],float)
        land=np.asarray(ds.variables["sftlf"][:],float)
        ice=np.asarray(ds.variables["sftgif"][:],float)
    if land.shape!=ice.shape!=(1080,2160) or lat.shape!=(1080,) or lon.shape!=(2160,):
        raise ValueError("PMIP4 1080x2160 grid changed")
    if not np.isin(land,[0,100]).all() or not np.isin(ice,[0,100]).all():
        raise ValueError("PMIP4 masks not binary percent")
    free=(land==100)&(ice==0)
    by_half_degree=free.reshape(360,3,720,3).sum(axis=(1,3))
    return np.roll(np.flipud(by_half_degree),360,axis=1)

def climate_05deg(assets):
    transform=from_origin(-180,90,.5,.5)
    collected=[]
    for asset in assets:
        with rasterio.open(asset) as src:
            if not src.crs or not src.crs.is_geographic:
                raise ValueError("CHELSA source raster geographic CRS missing")
            with WarpedVRT(src,crs="EPSG:4326",transform=transform,
                           width=720,height=360,resampling=Resampling.nearest) as vrt:
                arr=vrt.read(1,masked=True).astype(np.float32)
                result=np.asarray(arr.filled(np.nan),dtype=np.float32)
        if result.shape!=(360,720):raise ValueError("CHELSA 0.5deg output shape changed")
        collected.append(result)
    x=np.stack(collected,axis=-1)
    return np.isfinite(x).all(axis=-1)

def source_metrics(now_free,glacial_free,now_climate,glacial_climate):
    if any(x.shape!=(360,720) for x in (now_free,glacial_free,now_climate,glacial_climate)):
        raise ValueError("0.5-degree geometry does not align")
    response={"source":{},"coverage":{}}
    response["source"]["PMIP_0ka_ice_free_land_cells"]=int(now_free.sum())
    response["source"]["PMIP_21ka_ice_free_land_cells"]=int(glacial_free.sum())
    newly=glacial_free&~now_free
    retained=glacial_free&now_free
    response["source"]["newly_ice_free_land_at21ka_05deg_cells"]=int(newly.sum())
    response["coverage"]["21ka_climate_available_at_newly_icefree_land"]=int((newly&glacial_climate).sum())
    response["coverage"]["21ka_climate_missing_at_newly_icefree_land"]=int((newly&~glacial_climate).sum())
    response["coverage"]["21ka_climate_available_at_retained_icefree_land"]=int((retained&glacial_climate).sum())
    response["coverage"]["21ka_climate_missing_at_retained_icefree_land"]=int((retained&~glacial_climate).sum())
    response["coverage"]["0ka_climate_available_0ka_icefree_land"]=int((now_free&now_climate).sum())
    response["coverage"]["0ka_climate_missing_0ka_icefree_land"]=int((now_free&~now_climate).sum())
    return response

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pmip0",type=Path,required=True)
    p.add_argument("--pmip21",type=Path,required=True)
    p.add_argument("--chelsa",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    s0=pmip0p5(a.pmip0);s21=pmip0p5(a.pmip21)
    f0=climate_05deg([a.chelsa/f"CHELSA_TraCE21K_{v}_20_V1.0.tif" for v in VARS])
    f21=climate_05deg([a.chelsa/f"CHELSA_TraCE21K_{v}_-190_V1.0.tif" for v in VARS])
    cases={}
    for mode,cut in (("any_icefree_subcell",1),("majority_icefree_subcells",5),("all_nine_icefree_subcells",9)):
        m0=s0>=cut;m21=s21>=cut
        metrics={"source":{},"coverage":{}}
        vals=source_metrics(m0,m21,f0,f21)
        cases[mode]=vals
    summary={"schema":"chocho_ice6gc_CHELSA_05deg_paleoland_source_coverage_v02",
        "status":"SOURCE_FIRST_PALEOLAND_CLIMATE_COVERAGE_NOT_SPECIES_CONNECTIVITY",
        "global_regular_cells":360*720,
        "PMIP4_0_21ka_MD5":MD5,
        "cases":cases,
        "source_CHELSA_variables":list(VARS),
        "no_plant_models_or_genetics_used":True,
        "limits":["Original 10min ICE6GC fraction is coarsened to 0.5deg, so narrow corridors and small islands are unresolved",
           "New paleoland with missing CHELSA climate is NOT interpolated or classified as unsuitable",
           "Present WCVP native plant observation domains and dispersal histories are not checked by this asset compatibility audit"]}
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    np.savez_compressed(a.out/"coarse_masks_finite_05deg.npz",
      pmip0_icefree_subcells=s0.astype("u1"),pmip21_icefree_subcells=s21.astype("u1"),
      CHELSA0_fourvar_finite=f0.astype("u1"),CHELSA21_fourvar_finite=f21.astype("u1"))
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()
