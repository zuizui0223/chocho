#!/usr/bin/env python3
"""Read original 2015–2021 World Bank/IMF commercial AIS-density TIFF
at a source-explainable coarse grid, using original 0.005° TIFF overviews.

This is AIS *observation density*, not cargo volume, island-to-island voyages,
historical plant carriage, or vessel direction.
Source ZIP SHA pinned after official archive inspection.
The 1° averaging grid is a severe exploratory downscaling, and pixel mask is
never treated as observed zero ship traffic.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from affine import Affine

ORIGINAL_SHA="0cfdce41d9934982bbc8f16bd92ce7545cc088839327cbaa9b0ed1469dda31e8"
TIFF="ShipDensity_Commercial1.tif"
SOURCE_URL="https://datacatalogfiles.worldbank.org/ddh-published/0037580/5/DR0045405/shipdensity_commercial_.zip"
COARSE_SHAPE=(170,360)

def hash_file(p):
    m=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024*4),b""):m.update(b)
    return m.hexdigest()

def derive_raster(raw_zip:Path,gridout:Path,receiptout:Path):
    sha=hash_file(raw_zip)
    if sha!=ORIGINAL_SHA:raise ValueError("World Bank AIS official archive source SHA changed")
    with zipfile.ZipFile(raw_zip) as z:
        names=z.namelist()
        if TIFF not in names or TIFF+".ovr" not in names:
            raise ValueError("Original TIFF plus external pyramid must be present")
    uri=f"zip://{raw_zip.resolve()}!{TIFF}"
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      GDAL_NUM_THREADS="2",VSI_CACHE="TRUE",
                      VSI_CACHE_SIZE="50000000"):
        with rasterio.open(uri) as src:
            if str(src.crs)!="EPSG:4326" or src.count!=1 or src.width!=72006 or src.height!=33998 or src.nodata!=2147483647:
                raise ValueError("World Bank GeoTIFF original grid/NoData changed")
            overviews=src.overviews(1)
            print("AIS_SOURCE_OVERVIEW_FACTORS",overviews,flush=True)
            arr=src.read(1,out_shape=COARSE_SHAPE,masked=True,
                   resampling=Resampling.average,out_dtype="float32")
            trans=src.transform*Affine.scale(src.width/COARSE_SHAPE[1],src.height/COARSE_SHAPE[0])
            nodata=src.nodata
            source_meta={"original_width":src.width,"original_height":src.height,
                  "original_transform":tuple(src.transform),
                  "original_nodata":nodata,"original_overviews":overviews,
                  "original_dtype":src.dtypes[0],"src_description":src.descriptions}
    val=np.asarray(arr.data,dtype=np.float32)
    mask=np.asarray(np.ma.getmaskarray(arr),dtype=bool)
    invalid=~np.isfinite(val)
    mask|=invalid
    if np.any((val<0)&~mask):
        raise ValueError("Negative observed commercial AIS density")
    if len(overviews)==0:raise ValueError("Original 0.005-degree source has no efficient coarser overviews")
    count=int((~mask).sum())
    if count<40000:
        raise ValueError("Too few valid source AIS coarse cells; no interpolation of missing pixels permitted")
    vals=val[~mask]
    receipt={
      "schema":"chocho_20261106_worldbank_commercial_AIS2015_2021_1degree_derived_v1",
      "status":"SOURCE_VERIFIED_ORIGINAL_TIFF_COARSENED_NOT_ROUTE_OR_VESSEL_FLOW",
      "source_original_ZIP_url":SOURCE_URL,
      "source_original_SHA256":sha,
      "source_period":"2015-01 through 2021-02",
      "original_grid":source_meta,
      "derived_grid":{"shape":list(val.shape),"transform":list(trans)[:6],
         "crs":"EPSG:4326","cell_width_degrees":trans.a,
         "cell_height_degrees":abs(trans.e),
         "resampling":"rasterio average using original TIFF overviews, invalid NoData excluded",
         "count_valid_pixels":count,
         "count_missing_pixels":int(mask.sum()),
         "fraction_observed_zero_within_valid":float((vals==0).mean()),
         "sample_positive_quantiles_raw_AIS_index":{
           str(q):float(np.quantile(vals[vals>0],q))
           for q in (0.1,0.25,0.5,0.75,0.9,0.95,0.99)
         } if np.any(vals>0) else {},
         "count_positive_valid_pixels":int((vals>0).sum()),
         "max_positive_observed":float(np.max(vals)) if len(vals) else None,
         "metadata_not_actual_cargo_units":True
      },
      "source_nodata_excluded_not_zero":True,
      "note":"Coarse global raster for exploratory AIS-weighted marine resistance, not a direct voyage graph or measured propagule flux."
    }
    gridout.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(gridout,AIS_mean=val,source_valid_mask=(~mask),
                        raster_transform=np.array(list(trans)[:6],dtype=float))
    receiptout.parent.mkdir(parents=True,exist_ok=True)
    receiptout.write_text(json.dumps(receipt,indent=2,default=str)+"\n")
    print("AIS_SOURCE_1DEGREE_GRID_READY",json.dumps(receipt["derived_grid"],sort_keys=True)[:8500],flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser()
    for n in ("source-zip","coarse-grid","receipt"):p.add_argument("--"+n,type=Path,required=True)
    x=p.parse_args()
    derive_raster(x.source_zip,x.coarse_grid,x.receipt)
if __name__=="__main__":main()
