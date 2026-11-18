#!/usr/bin/env python3
"""Source-only AIS audit; does not access biological outcomes or modify inputs.

First extract the checksum-pinned source ZIP with a CRC-validating ZIP reader.
This script samples native and first-overview pixels, reconstructs coarse means
independently by exact rectangular source-pixel overlap, and records metadata.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import warnings
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import rasterio
from rasterio.errors import NotGeoreferencedWarning
from rasterio.transform import Affine
from rasterio.windows import Window

EXPECTED_SHA="0cfdce41d9934982bbc8f16bd92ce7545cc088839327cbaa9b0ed1469dda31e8"
NATIVE="ShipDensity_Commercial1.tif"


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(4*1024*1024),b""):h.update(chunk)
    return h.hexdigest()


def metadata(ds):
    return {"shape":list(ds.shape),"dtypes":list(ds.dtypes),"nodata":ds.nodata,
            "scales":list(ds.scales),"offsets":list(ds.offsets),"crs":str(ds.crs),
            "transform":list(ds.transform)[:6],"band_tags":ds.tags(1),
            "dataset_tags":ds.tags(),"image_structure":ds.tags(ns="IMAGE_STRUCTURE"),
            "overview_factors":ds.overviews(1)}


def rectangular_mean(ds,transform,left,top,step):
    """Independent exact overlap arithmetic; no GDAL resampling call."""
    right,bottom=left+step,top-step
    dx,dy=transform.a,-transform.e
    c0=max(0,int(np.floor((left-transform.c)/dx)))
    c1=min(ds.width,int(np.ceil((right-transform.c)/dx)))
    r0=max(0,int(np.floor((transform.f-top)/dy)))
    r1=min(ds.height,int(np.ceil((transform.f-bottom)/dy)))
    a=ds.read(1,window=Window(c0,r0,c1-c0,r1-r0),masked=True)
    xs=transform.c+np.arange(c0,c1)*dx
    ys=transform.f-np.arange(r0,r1)*dy
    wx=np.maximum(np.minimum(xs+dx,right)-np.maximum(xs,left),0)
    wy=np.maximum(np.minimum(ys,top)-np.maximum(ys-dy,bottom),0)
    weights=wy[:,None]*wx[None,:]
    valid=~np.ma.getmaskarray(a)
    w=np.where(valid,weights,0)
    mean=float(np.sum(np.where(valid,a.data,0)*w)/w.sum()) if w.sum()>0 else None
    return {"mean":mean,"native_pixels_in_window":int(a.size),
            "source_nodata_in_window":int((~valid).sum()),
            "valid_overlap_fraction":float(w.sum()/(step*step))}


def cases():
    locations=[("open_NPacific",-150.,30.),("open_SPacific",-120.,-20.),
        ("northAtlantic",-30.,35.),("Singapore",103.8,1.25),
        ("EnglishChannel",-2.,50.),("HongKong",114.2,22.25),
        ("IndianOcean",80.,-10.),("SouthAtlantic",-10.,-30.),
        ("TokyoBay",139.8,35.4), ("continental_zero",30.,20.)]
    seen={(x,y) for _,x,y in locations}
    for lat in (-60.,-30.,0.,30.,60.):
        for lon in (-150.,-90.,-30.,30.,90.,150.):
            if (lon,lat) not in seen:
                locations.append((f"systematic_{lon:g}_{lat:g}",lon,lat))
                seen.add((lon,lat))
    return locations


def run(archive,extracted,coarse,out):
    source_sha=digest(archive)
    if source_sha!=EXPECTED_SHA:raise ValueError("Official source ZIP SHA mismatch")
    with zipfile.ZipFile(archive) as z:
        members={x.filename:{"uncompressed_bytes":x.file_size,"CRC32":hex(x.CRC)}
                 for x in z.infolist()}
    for name in (NATIVE,NATIVE+".ovr",NATIVE+".aux.xml"):
        if (extracted/name).stat().st_size!=members[name]["uncompressed_bytes"]:
            raise ValueError("Extracted source length mismatch: "+name)
    xml=ET.parse(extracted/(NATIVE+".aux.xml"))
    native_stats={x.attrib["key"]:x.text for x in xml.findall(".//MDI")}
    histogram_total=sum(map(int,xml.findtext(".//HistCounts").split("|")))
    with np.load(coarse) as data:
        grid=np.asarray(data["AIS_mean"])
        mask=np.asarray(data["source_valid_mask"])
        coarse_transform=Affine(*data["raster_transform"])
    rows=[];deltas=[];all_native=[];all_overview=[]
    with warnings.catch_warnings(),rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                                                GDAL_CACHEMAX=64*1024*1024):
        warnings.simplefilter("ignore",NotGeoreferencedWarning)
        with rasterio.open(extracted/NATIVE) as native,rasterio.open(extracted/(NATIVE+".ovr")) as ovr:
            nm,om=metadata(native),metadata(ovr)
            if nm["overview_factors"]:raise ValueError("Native comparison must not use overview shortcuts")
            nt=native.transform
            ot=Affine(nt.a*2,0,nt.c,0,nt.e*2,nt.f)
            nearest_agreement=0;sampled_overview_pixels=0;mask_disagreements=0
            for name,lon,lat in cases():
                r,c=native.index(lon,lat)
                ro,co=r//2-32,c//2-32
                a=native.read(1,window=Window(co*2,ro*2,128,128),masked=True)
                b=ovr.read(1,window=Window(co,ro,64,64),masked=True)
                am,bm=np.ma.getmaskarray(a),np.ma.getmaskarray(b)
                nearest=a.data[::2,::2]
                match=(b.data==nearest)&(bm==am[::2,::2])
                nearest_agreement+=int(match.sum());sampled_overview_pixels+=int(b.size)
                mask_disagreements+=int(np.count_nonzero(bm!=am[::2,::2]))
                box_mean=a.reshape(64,2,64,2).mean(axis=(1,3))
                both=~bm&~np.ma.getmaskarray(box_mean)
                deltas.extend(np.abs(b.data[both]-box_mean.data[both]).tolist())
                all_native.append(a.compressed());all_overview.append(b.compressed())
                gr=int(np.floor((coarse_transform.f-lat)/(-coarse_transform.e)))
                gc=int(np.floor((lon-coarse_transform.c)/coarse_transform.a))
                left=coarse_transform.c+gc*coarse_transform.a
                top=coarse_transform.f+gr*coarse_transform.e
                nr=rectangular_mean(native,nt,left,top,coarse_transform.a)
                ov=rectangular_mean(ovr,ot,left,top,coarse_transform.a)
                observed=float(grid[gr,gc])
                error=abs(observed-ov["mean"])
                relative=error/max(abs(ov["mean"]),1.)
                rows.append({"name":name,"longitude":lon,"latitude":lat,
                    "coarse_row":gr,"coarse_col":gc,
                    "native_window_max":int(a.max()),"overview_window_max":int(b.max()),
                    "native_window_nodata":int(am.sum()),"overview_window_nodata":int(bm.sum()),
                    "overview_equals_selected_native_pixel_count":int(match.sum()),
                    "overview_pixels_compared":int(b.size),
                    "native_direct_coarse_overlap":nr,"overview_direct_coarse_overlap":ov,
                    "saved_coarse_mean":observed,"saved_coarse_valid":bool(mask[gr,gc]),
                    "saved_vs_overview_mean_absolute_error":error,
                    "saved_vs_overview_mean_relative_error":relative,
                    "saved_coarse_mean_matches_within_float32_tolerance":relative<=3e-7,
                    "overview_vs_native_coarse_mean_relative_difference":
                        (ov["mean"]-nr["mean"])/nr["mean"] if nr["mean"] else None})
    av=np.concatenate(all_native);bv=np.concatenate(all_overview);deltas=np.asarray(deltas)
    def value_profile(x):
        return {"valid_sampled_pixels":int(len(x)),"min":int(x.min()),"max":int(x.max()),
            "zero_count":int((x==0).sum()),"negative_count":int((x<0).sum()),
            "nodata_sentinel_in_valid_sample":int((x==2147483647).sum()),
            "greater_than_one_million":int((x>1_000_000).sum()),
            "quantiles":{str(q):float(np.quantile(x,q)) for q in (.5,.9,.95,.99)}}
    report={
        "schema":"chocho_published_AIS_native_overview_coarse_source_quality_v1",
        "audit_utc":datetime.now(timezone.utc).isoformat(),
        "status":"PUBLISHED_MAGNITUDES_PRESENT_IN_NATIVE_SOURCE_NO_SAMPLED_TRANSFORM_ERROR"
            if all(x["saved_coarse_mean_matches_within_float32_tolerance"] for x in rows)
            else "COARSE_TRANSFORMATION_DISCREPANCY_REQUIRES_INVESTIGATION",
        "source_ZIP_sha256":source_sha,"coarse_grid_sha256":digest(coarse),
        "official_catalogue":"https://datacatalog.worldbank.org/search/dataset/0037580/global-shipping-traffic-density",
        "official_readme":"https://datacatalogfiles.worldbank.org/ddh-published/0037580/5/DR0084213/readme.txt",
        "source_definition":"Published total hourly AIS position reports per0.005-degree cell, Jan2015-Feb2021; moving and stationary vessels.",
        "native_metadata":nm,"overview_metadata":om,"original_aux_statistics":native_stats,
        "original_aux_histogram_total":histogram_total,"native_full_pixel_count":33998*72006,
        "archive_members":members,
        "extraction":"Native and external overview sequentially extracted by Python zipfile; both member CRCs verified at EOF before audit.",
        "sampling":"10 fixed diagnostic locations plus5latitudes×6longitudes, excluding exact geographic duplicates; no biological inputs.",
        "diagnostic_locations":len(rows),
        "native_sample_profile":value_profile(av),"overview_sample_profile":value_profile(bv),
        "sampled_native_to_overview":{
            "overview_pixels_compared":sampled_overview_pixels,
            "exact_top_left_native_2x2_matches":nearest_agreement,
            "exact_match_fraction":nearest_agreement/sampled_overview_pixels,
            "nodata_mask_disagreements":mask_disagreements,
            "not_equal_to_2x2_mean_count":int((deltas>0).sum()),
            "absolute_difference_from_2x2_mean_quantiles":{
                str(q):float(np.quantile(deltas,q)) for q in (.5,.9,.95,.99,1.)}},
        "coarse_overlap_reproduction":{
            "compared_cells":len(rows),"all_within_float32_tolerance":all(x["saved_coarse_mean_matches_within_float32_tolerance"] for x in rows),
            "max_absolute_error":max(x["saved_vs_overview_mean_absolute_error"] for x in rows),
            "max_relative_error":max(x["saved_vs_overview_mean_relative_error"] for x in rows)},
        "cases":rows,
        "limits":[
            "Checks establish faithful reproduction of published numeric values, not independent physical calibration of AIS reporting counts.",
            "The native high values and scale1/offset0 provide no documented rescaling rule; no arbitrary rescaling or value clipping was applied.",
            "First overview samples are nearest-neighbor selections of native2x2 pixels, not native2x2 arithmetic means; coarse outputs are correctly labeled overview-derived indices.",
            "The source aux histogram population need not cover every native cell; its total is recorded without interpreting this as a pixel-data defect.",
            "Sampled checks do not certify every native/overview pixel or source reception coverage."
        ]}
    out.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ("status","diagnostic_locations","native_sample_profile",
        "sampled_native_to_overview","coarse_overlap_reproduction")},indent=2),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("source-zip","extracted-dir","coarse-grid","out"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args();run(a.source_zip,a.extracted_dir,a.coarse_grid,a.out)
