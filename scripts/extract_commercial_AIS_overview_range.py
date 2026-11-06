#!/usr/bin/env python3
"""Read TRUE raster pixels from official World Bank/IMF commercial AIS source.

Only fetch 124 MB compressed GeoTIFF OVERVIEW member by an exact HTTP Range
of the 480 MB ZIP, check uncompressed size and CRC32. Avoid downloading or
expanding ~9.8 GB original full-resolution TIFF. Inspect GeoTIFF overview
georeferencing, overview resampling and value semantics before use for
island pair electrical conductance. Archive only small provenance/metadata.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,os,shutil,struct,time,urllib.request,zlib,zipfile
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

from probe_commercial_AIS_ZIP_range_manifest import URL,EXPECTED_LENGTH,RemoteRangeZip,head

OVR="ShipDensity_Commercial1.tif.ovr"
TFW="ShipDensity_Commercial1.tfw"
AUX="ShipDensity_Commercial1.tif.aux.xml"


def fetch_compressed_range(out:Path,info:zipfile.ZipInfo,remote:RemoteRangeZip):
    remote.seek(info.header_offset)
    header=remote.read(30)
    values=struct.unpack("<IHHHHHIIIHH",header)
    if values[0]!=0x04034b50:raise ValueError("Official AIS ZIP local file record invalid")
    name_len,extra_len=values[-2:]
    name=remote.read(name_len).decode("utf-8")
    if name!=info.filename:raise ValueError("Archive ZIP central/local filename disagreement")
    remote.read(extra_len)
    start=remote.tell();end=start+info.compress_size-1
    req=urllib.request.Request(URL,headers={
        "Range":f"bytes={start}-{end}",
        "User-Agent":"chocho-island51-source-AIS-overview-data/1.0"})
    out.parent.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(out.parent).free < info.file_size+800_000_000:
        raise RuntimeError(f"Not enough scratch storage to decode AIS 3.3 GB overview member: {shutil.disk_usage(out.parent).free}")
    done=0
    crc=0
    total_in=0
    dec=zlib.decompressobj(-15)
    with urllib.request.urlopen(req,timeout=95) as response:
        if response.status!=206:raise ValueError(f"Server did not honor Range request: HTTP {response.status}")
        cr=response.headers.get("Content-Range","")
        if not cr.startswith(f"bytes {start}-{end}/"):
            raise ValueError("Compressed member Range returned wrong offsets "+cr)
        with out.open("wb") as stream:
            while total_in<info.compress_size:
                b=response.read(min(2*1024*1024,info.compress_size-total_in))
                if not b:raise RuntimeError("Truncated official compressed AIS overview")
                total_in+=len(b)
                data=dec.decompress(b)
                if data:
                    stream.write(data);done+=len(data);crc=zlib.crc32(data,crc)
                if total_in%(32*1024*1024)<2*1024*1024:
                    print("AIS_OFFICIAL_OVERVIEW_RANGE_DECODE_BYTES",total_in,done,flush=True)
            tail=dec.flush()
            if tail:stream.write(tail);done+=len(tail);crc=zlib.crc32(tail,crc)
    if total_in!=info.compress_size or done!=info.file_size or (crc&0xffffffff)!=info.CRC:
        raise ValueError(f"AIS ZIP original member CRC/size mismatch: {total_in}, {done}, {crc}")
    return {"official_archive_member":info.filename,"compressed_range_start":start,
        "compressed_range_end":end,"source_compressed_bytes":total_in,
        "source_uncompressed_bytes":done,"verified_zip_crc32":hex(crc & 0xffffffff),
        "source_compress_type":info.compress_type}


def raster_probe(path:Path):
    with rasterio.open(path) as ds:
        meta={
          "driver":ds.driver,"width":ds.width,"height":ds.height,
          "count":ds.count,"crs":str(ds.crs),"pixel_resolution":list(ds.res),
          "affine":list(ds.transform)[:6],
          "bounds":list(ds.bounds),
          "nodata":ds.nodata,"dtypes":list(ds.dtypes),
          "block_shapes":list(ds.block_shapes),
          "overviews_for_band1":ds.overviews(1),
          "bands_tags":[ds.tags(1)] if ds.count>=1 else [],
          "image_structure":ds.tags(ns="IMAGE_STRUCTURE"),
          "file_tags":ds.tags(),
          "sample_windows_no_geography_assumption":[]}
        # Sample a few small raw block windows to check source value types;
        # no "ship corridor" claim until georeference/aggregation verified.
        windows=[
           (ds.width//8,ds.height//4),
           (ds.width//2,ds.height//2),
           (ds.width*7//8,ds.height*3//4)]
        for c,r in windows:
            w=Window(max(0,c-3),max(0,r-3),6,6)
            x=ds.read(1,window=w,masked=True)
            mask=np.ma.getmaskarray(x)
            vals=np.asarray(x.data)[~mask]
            meta["sample_windows_no_geography_assumption"].append({
              "top_left_colrow":[int(c-3),int(r-3)],
              "n_valid":int(len(vals)),
              "n_masked":int(mask.sum()),
              "min":float(np.min(vals)) if len(vals)>0 else None,
              "max":float(np.max(vals)) if len(vals)>0 else None,
              "n_positive":int(np.sum(vals>0)) if len(vals)>0 else 0})
    return meta

def run(dest:Path,receipt:Path):
    h=head()
    if int(h.get("Content-Length","0"))!=EXPECTED_LENGTH:
        raise ValueError("Original public AIS ZIP changed size")
    remote=RemoteRangeZip(URL,EXPECTED_LENGTH)
    with zipfile.ZipFile(remote) as z:
        names={r.filename:r for r in z.infolist()}
        if not {OVR,TFW,AUX}.issubset(names):
            raise ValueError("AIS ZIP missing expected 2015-2021 overview or WGS84 raster metadata")
        wordfile=z.open(TFW).read().decode("utf-8",errors="replace")
        auxiliary=z.open(AUX).read().decode("utf-8",errors="replace")
        info=names[OVR]
        if info.file_size!=3290970915 or info.compress_size!=124257699:
            raise ValueError("Original World Bank commercial AIS overview member changed")
        extraction=fetch_compressed_range(dest,info,remote)
    grid=raster_probe(dest)
    out={
      "schema":"chocho_WorldBank_IMF_2015_2021_commercial_AIS_real_raster_overview_probe_v1",
      "status":"ORIGINAL_COMMERCIAL_AIS_OVERVIEW_GDAL_PIXEL_VALUES_READ_AND_CRC_VERIFIED",
      "official_catalogue":"https://datacatalog.worldbank.org/search/dataset/0037580/global-shipping-traffic-density",
      "source_url":URL,
      "source_zip_bytes":EXPECTED_LENGTH,
      "source_license":"CC BY 4.0",
      "source_time":"2015-01 through 2021-02 hourly AIS position broadcasts",
      "source_cell_units":"reported AIS positions per source pixel, not vessels, voyages, cargo or ecological dispersal",
      "original_source_world_file":wordfile,
      "original_raster_aux_xml":auxiliary[:6000],
      "compressed_member":extraction,
      "raster_overview":grid,
      "overview_is_not_native_0p005_degree_raster":True,
      "island_51_geographic_pixel_locations_NOT_yet_sampled":True,
      "transport_effect_NOT_yet_tested":True,
      "limits":["The ZIP member is a raster external overview; resampling method and geographic mapping need verification against base native GeoTIFF metadata before interpreting cell values.",
        "AIS satellite reception coverage, port stationarity and ship class availability influence raw counted positions.",
        "A marine conductor from averaged AIS has arbitrary transformations requiring predeclared scale and robust sensitivity.",
        "Historic plant movement predates 2015-2021 commercial AIS source, so cannot be claimed as cause.",
        "AISCargo and butterfly wind advection cannot be inferred from ship position density."]
    }
    receipt.parent.mkdir(parents=True,exist_ok=True)
    receipt.write_text(json.dumps(out,ensure_ascii=False,indent=2,default=str)+"\n",encoding="utf-8")
    print("AIS_SOURCE_REAL_OVERVIEW_PIXEL_METADATA",json.dumps({
        "source_uncompressed_bytes":extraction["source_uncompressed_bytes"],
        "driver":grid["driver"],"width":grid["width"],"height":grid["height"],
        "dtypes":grid["dtypes"],"crs":grid["crs"],"pixel_resolution":grid["pixel_resolution"],
        "bounds":grid["bounds"],"nodata":grid["nodata"],
        "overviews":grid["overviews_for_band1"],"samples":grid["sample_windows_no_geography_assumption"]
    },default=str),flush=True)
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--overview",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.overview,a.out)
if __name__=="__main__":main()
