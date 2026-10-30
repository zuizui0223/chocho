#!/usr/bin/env python3
"""Source-lock independent GSHHG 2.3.7 sea/land island topology for chocho.

Same source GSHHG version as user's separate island project. Do not assume all
Natural Earth land features or H3 land cells correspond to individual islands.
"""
import argparse,hashlib,io,json,zipfile
from pathlib import Path
from urllib.request import Request,urlopen

URL="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
SCHEMA="chocho_gshhg_237_original_island_landmass_source_gate_20261029_v01"
MAX_SIZE=160_000_000

def inspect(src):
    actual=hashlib.sha256(src).hexdigest()
    if actual!=SHA:raise ValueError("Official GSHHG release zip SHA256 mismatched")
    with zipfile.ZipFile(io.BytesIO(src)) as z:
        # The source SHA256 pins the complete official 149,157,845-byte ZIP.
        # Avoid decompressing unrelated full-resolution coastline members during source audit.
        names=z.namelist()
        l1=[n for n in names if n.replace("\\","/").endswith("GSHHS_h_L1.shp")]
        if len(l1)!=1:raise ValueError("Original high-resolution shoreline level1 not found")
        prefix=l1[0][:-4]
        for suf in (".dbf",".shx",".prj"):
            if prefix+suf not in names:raise ValueError("Incomplete original shoreline shapefile")
        if sum(x.file_size for x in z.infolist())>2_000_000_000:
            raise ValueError("GSHHG original decompression beyond bound")
        l1_shp_bytes=z.getinfo(l1[0]).file_size
        if l1_shp_bytes>150_000_000:raise ValueError("Unexpected GSHHS high L1 original shapefile size")
        readme=[n for n in names if n.lower().endswith("license.txt")]
    return {"source":"GSHHG 2.3.7 original shorelines, high resolution L1",
        "official_sha256":SHA,
        "original_archive_bytes":len(src),
        "land_polygon_source_shapefile_bytes":l1_shp_bytes,
        "gshhg_high_shoreline_level":1,
        "original_sha_verified":True,
        "original_source_ref":"https://github.com/zuizui0223/island/blob/main/config/island_source_gshhg.yml",
        "separate_island_project_thresholds_km2":[1,5,20],
        "separate_island_project_uses_5km2_primary_and_20km2_sensitivity":True,
        "continental_mainland_comparison_not_yet_identified":True}

def run(p):
    out={"schema":SCHEMA,"source_url":URL,"no_raw_island_polygon_coordinates_exported":True}
    try:
        req=Request(URL,headers={"User-Agent":"chocho-source-gshhg/1.0",
            "Accept":"application/zip,application/octet-stream"})
        with urlopen(req,timeout=130) as r:
            b=r.read(MAX_SIZE+1)
            if r.status!=200:raise ValueError("Original GSHHG source HTTP fail")
        if len(b)>MAX_SIZE:raise ValueError("Source file exceeds approved original archive size")
        out.update(inspect(b))
        out["status"]="OFFICIAL_237_GSHHG_ISLAND_SHORELINE_SOURCE_VERIFIED"
    except (OSError,ValueError,zipfile.BadZipFile) as e:
        out["status"]="OFFICIAL_ISLAND_SOURCE_INCOMPLETE"
        out["failure_type"]=type(e).__name__
        out["failure_reason"]=str(e)[:260]
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["status"],"original_archive_bytes":out.get("original_archive_bytes"),
      "failure":out.get("failure_reason")}),flush=True)
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    x=run(p.parse_args().output)
    if x["status"]!="OFFICIAL_237_GSHHG_ISLAND_SHORELINE_SOURCE_VERIFIED":raise SystemExit(2)
