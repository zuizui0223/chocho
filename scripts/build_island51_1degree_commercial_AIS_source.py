#!/usr/bin/env python3
"""Average the original World Bank/IMF commercial AIS external TIFF overview
onto an exact 1° or 0.5° WGS84 grid for the independently defined ocean graph.

This is AIS *observation density*, not cargo volume, island-to-island voyages,
historical plant carriage, or vessel direction.
Source ZIP SHA pinned after official archive inspection.
The first external overview is 0.01°, already coarser than the native 0.005°
source. Its original pyramid resampling method is not supplied, so this is an
average of overview values, not a reconstructed native-resolution mean.
Missing pixels are never treated as observed zero ship traffic.
"""
from __future__ import annotations
import argparse,hashlib,json,math,shutil,tempfile,warnings,zipfile
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.errors import NotGeoreferencedWarning
from rasterio.vrt import WarpedVRT
from affine import Affine

ORIGINAL_SHA="0cfdce41d9934982bbc8f16bd92ce7545cc088839327cbaa9b0ed1469dda31e8"
TIFF="ShipDensity_Commercial1.tif"
OVERVIEW=TIFF+".ovr"
WORLDFILE="ShipDensity_Commercial1.tfw"
SOURCE_URL="https://datacatalogfiles.worldbank.org/ddh-published/0037580/5/DR0045405/shipdensity_commercial_.zip"
COARSE_SHAPE=(170,360)
COARSE_TRANSFORM=Affine(1,0,-180,0,-1,85)
MIN_VALID_FRACTION=40000/(170*360)

def hash_file(p):
    m=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024*4),b""):m.update(b)
    return m.hexdigest()


def global_grid(resolution_degrees):
    if resolution_degrees not in (.5,1.):
        raise ValueError("Supported AIS grid resolution is 0.5 or 1 degree")
    resolution=float(resolution_degrees)
    shape=(int(170/resolution),int(360/resolution))
    transform=Affine(resolution,0,-180,0,-resolution,85)
    latitude=85-resolution*(np.arange(shape[0],dtype=float)+.5)
    longitude=-180+resolution*(np.arange(shape[1],dtype=float)+.5)
    return shape,transform,latitude,longitude


def validate_source_coverage(valid):
    count=int(valid.sum())
    fraction=count/valid.size if valid.size else 0.
    if fraction<MIN_VALID_FRACTION:
        raise ValueError("Too small a valid source AIS fraction; missing pixels cannot be interpolated")
    return count,fraction


def average_overview_on_grid(overview_path, source_transform, *,
                             target_shape=None,target_transform=None,
                             resolution_degrees=1.):
    """Area-overlap average of actual overview values on an explicit WGS84 grid.

    External overview TIFFs have no CRS of their own. The caller supplies the
    transform recovered from the checksum-verified companion world file.
    """
    shape,transform,_,_=global_grid(resolution_degrees)
    if target_shape is None:target_shape=shape
    if target_transform is None:target_transform=transform
    with warnings.catch_warnings():
        # An external overview is expected to lack its own georeference. Only
        # this known warning is suppressed; the supplied transform is required.
        warnings.simplefilter("ignore",NotGeoreferencedWarning)
        source=rasterio.open(overview_path)
    with source:
        if source.count!=1 or source.nodata is None:
            raise ValueError("A one-band overview with explicit source NoData is required")
        source_width,source_height=source.width,source.height
        with WarpedVRT(source,src_crs="EPSG:4326",src_transform=source_transform,
                       crs="EPSG:4326",transform=target_transform,
                       height=target_shape[0],width=target_shape[1],
                       resampling=Resampling.average,dtype="float32",
                       src_nodata=source.nodata,nodata=np.nan,
                       warp_mem_limit=128) as target:
            values=target.read(1,masked=False)
    valid=np.isfinite(values)
    # GDAL's averaging kernel can return an edge value in a destination cell
    # that only touches the source rectangle. Require positive geographic
    # overlap; uncovered cells must not acquire apparently observed traffic.
    source_left,source_top=source_transform.c,source_transform.f
    source_right=source_left+source_width*source_transform.a
    source_bottom=source_top+source_height*source_transform.e
    target_left=target_transform.c+np.arange(target_shape[1])*target_transform.a
    target_right=target_left+target_transform.a
    target_top=target_transform.f+np.arange(target_shape[0])*target_transform.e
    target_bottom=target_top+target_transform.e
    overlap_x=np.minimum(target_right,source_right)-np.maximum(target_left,source_left)>1e-12
    overlap_y=np.minimum(target_top,source_top)-np.maximum(target_bottom,source_bottom)>1e-12
    valid&=overlap_y[:,None]&overlap_x[None,:]
    if np.any(values[valid]<0):
        raise ValueError("Negative observed commercial AIS density")
    values[~valid]=np.nan
    return values,valid


def extract_external_overview(raw_zip, destination):
    """Extract just the CRC-checked pyramid, using the original world file.

    The caller verifies the whole archive SHA before invoking this helper.
    The 9.8 GB native raster is never opened or extracted.
    """
    destination=Path(destination)
    with zipfile.ZipFile(raw_zip) as archive:
        if not {TIFF,OVERVIEW,WORLDFILE}.issubset(archive.namelist()):
            raise ValueError("Original TIFF, external overview and worldfile must be present")
        worldfile=archive.read(WORLDFILE).decode("utf-8")
        try:
            a,d,b,e,c,f=map(float,worldfile.split())
        except ValueError as exc:
            raise ValueError("Unsupported original AIS worldfile") from exc
        if (not np.isfinite([a,d,b,e,c,f]).all() or
            not math.isclose(a,.005,abs_tol=1e-12) or
            not math.isclose(e,-.005,abs_tol=1e-12) or
            abs(d)>1e-12 or abs(b)>1e-12 or abs(c+180)>1 or abs(f-85)>1):
            raise ValueError("Unsupported original AIS worldfile resolution, rotation or extent")
        native=Affine(a,b,c-a/2-b/2,d,e,f-d/2-e/2)
        transform=Affine(a*2,b*2,native.c,d*2,e*2,native.f)
        info=archive.getinfo(OVERVIEW)
        destination.parent.mkdir(parents=True,exist_ok=True)
        if shutil.disk_usage(destination.parent).free<info.file_size+128*1024*1024:
            raise RuntimeError("Insufficient scratch disk to extract original AIS external overview")
        with archive.open(info) as member,destination.open("wb") as stream:
            # Reading to EOF verifies the original ZIP member CRC32.
            shutil.copyfileobj(member,stream,length=8*1024*1024)
    if destination.stat().st_size!=info.file_size:
        raise ValueError("Extracted original AIS overview size mismatch")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",NotGeoreferencedWarning)
        source=rasterio.open(destination)
    with source:
        if source.count!=1 or source.dtypes[0]!="int32" or source.nodata!=2147483647:
            raise ValueError("Original AIS overview band, dtype or NoData changed")
        meta={
            "source_member":OVERVIEW,"source_compressed_bytes":info.compress_size,
            "source_uncompressed_bytes":info.file_size,"verified_zip_crc32":hex(info.CRC),
            "original_source_world_file":worldfile,
            "original_transform":list(native)[:6],
            "overview_transform":list(transform)[:6],
            "overview_width":source.width,"overview_height":source.height,
            "overview_nodata":source.nodata,"overview_dtype":source.dtypes[0],
            "overview_additional_pyramid_factors":source.overviews(1),
            "source_base_raster_read":False,
            "source_original_pyramid_resampling_method":"not supplied; average of overview values only",
        }
    return transform,meta

def derive_raster(raw_zip:Path,gridout:Path,receiptout:Path,*,resolution_degrees=1.):
    shape,trans,latitudes,longitudes=global_grid(resolution_degrees)
    raw_zip=Path(raw_zip)
    sha=hash_file(raw_zip)
    if sha!=ORIGINAL_SHA:raise ValueError("World Bank AIS official archive source SHA changed")
    with zipfile.ZipFile(raw_zip) as z:
        if OVERVIEW not in z.namelist() or z.getinfo(OVERVIEW).file_size!=3290970915:
            raise ValueError("Original World Bank external overview missing or changed size")
    with tempfile.TemporaryDirectory(prefix="ais-overview-",dir=raw_zip.parent) as temp:
        overview_path=Path(temp)/OVERVIEW
        source_transform,source_meta=extract_external_overview(raw_zip,overview_path)
        if (source_meta["overview_width"],source_meta["overview_height"])!=(36003,16999):
            raise ValueError("World Bank original 0.01-degree overview dimensions changed")
        print("AIS_SOURCE_EXPLICIT_EXTERNAL_OVERVIEW_READY",json.dumps(source_meta),flush=True)
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_NUM_THREADS="2",
                          GDAL_CACHEMAX=128*1024*1024):
            val,valid=average_overview_on_grid(overview_path,source_transform,
                resolution_degrees=resolution_degrees)
    count,valid_fraction=validate_source_coverage(valid)
    vals=val[valid]
    receipt={
      "schema":"chocho_20261106_worldbank_commercial_AIS2015_2021_fixed_grid_derived_v3",
      "status":"SOURCE_VERIFIED_ORIGINAL_TIFF_COARSENED_NOT_ROUTE_OR_VESSEL_FLOW",
      "source_original_ZIP_url":SOURCE_URL,
      "source_original_SHA256":sha,
      "source_period":"2015-01 through 2021-02",
      "original_grid":source_meta,
      "derived_grid":{"shape":list(val.shape),"transform":list(trans)[:6],
         "crs":"EPSG:4326","cell_width_degrees":trans.a,
         "cell_height_degrees":abs(trans.e),
         "resampling":"GDAL geographic area-overlap average of explicitly extracted 0.01-degree overview; source NoData excluded",
         "fixed_alignment":f"Exactly {resolution_degrees:g}-degree cells: west=-180, north=85, east=180, south=-85",
         "input_to_resampling":"Original extracted 0.01-degree overview, never an interpolation of a previously derived coarse grid",
         "original_native_to_overview_resampling_known":False,
         "count_valid_pixels":count,
         "fraction_valid_pixels":valid_fraction,
         "minimum_valid_fraction":MIN_VALID_FRACTION,
         "count_missing_pixels":int((~valid).sum()),
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
    np.savez_compressed(gridout,AIS_mean=val,source_valid_mask=valid,
                        raster_transform=np.array(list(trans)[:6],dtype=float),
                        lat_centers=latitudes,lon_centers=longitudes)
    receiptout.parent.mkdir(parents=True,exist_ok=True)
    receiptout.write_text(json.dumps(receipt,indent=2,default=str)+"\n")
    print("AIS_SOURCE_COARSE_GRID_READY",json.dumps(receipt["derived_grid"],sort_keys=True)[:8500],flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser()
    for n in ("source-zip","coarse-grid","receipt"):p.add_argument("--"+n,type=Path,required=True)
    p.add_argument("--resolution-degrees",type=float,choices=(.5,1.),default=1.)
    x=p.parse_args()
    derive_raster(x.source_zip,x.coarse_grid,x.receipt,resolution_degrees=x.resolution_degrees)
if __name__=="__main__":main()
