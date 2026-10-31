#!/usr/bin/env python3
"""Build the global butterfly island sampling frame BEFORE using butterfly outcomes.

High-resolution GSHHG 2.3.7 level-1 landmasses with equal-area area >=20 km2,
<7,000,000 km2 (the earlier chocho continent/island separator), not coastal
points selected by later butterfly first-index events. GSHHG sibling_id is the
immutable join key shared with 17-island validation and original year sources.

Do not call every polygon 'oceanic island': Australia/mainland boundaries and
GSHHG complex/multipart landmasses must be separately reviewed. A static 2020
shoreline source is not historically annually changing island area.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

import geopandas as gpd
import shapely
from shapely.ops import transform
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
GSHHG_URL="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
SOURCE17="provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
SOURCE17_BLOB="54c4515947d2a3ee7545e4007985ef82a2af6868"
AREA_MIN=20.0
AREA_MAX=7_000_000.0
def original_17():
    raw=(ROOT/SOURCE17).read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if blob!=SOURCE17_BLOB:raise ValueError("Original 17 island source changed "+blob)
    v=json.loads(raw)
    ids={str(x["island_id"]) for x in v["islands"]}
    if len(ids)!=17:raise ValueError("Original 17 islands unexpectedly changed")
    return ids

def source_shape(archive,work):
    raw=archive.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=GSHHG_SHA:raise ValueError("GSHHG archive SHA-256 mismatch: "+digest)
    work.mkdir(exist_ok=True,parents=True)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        files=[x for x in z.namelist() if x.endswith((".shp",".shx",".dbf",".prj"))
               and "/h/GSHHS_h_L1." in "/"+x]
        if len(files)<3:raise ValueError("High-resolution level1 GSHHG shapefile missing")
        for x in files:z.extract(x,work)
    land=gpd.read_file(next(work.rglob("GSHHS_h_L1.shp")))
    if "sibling_id" not in land.columns:
        raise ValueError("GSHHG landmass source missing original sibling_id")
    land["landmass_id"]=land.sibling_id.astype(int).astype(str)
    return land,digest

def run(archive,work):
    land,sha=source_shape(archive,work)
    project=Transformer.from_crs(4326,6933,always_xy=True)
    ids=original_17()
    frame=[]
    for landmass,rows in land.groupby("landmass_id"):
        geom=shapely.union_all(rows.geometry.to_numpy())
        if geom.is_empty or not geom.is_valid:
            raise ValueError("Invalid GSHHG geometry landmass "+landmass)
        area_km2=transform(project.transform,geom).area/1e6
        if not AREA_MIN<=area_km2<AREA_MAX:continue
        point=geom.representative_point()
        llon,llat=geom.bounds[0],geom.bounds[1]
        frame.append({
            "landmass_id":landmass,
            "area_km2_equal_area":area_km2,
            "reference_lon":point.x,
            "reference_lat":point.y,
            "bbox_minlon":llon,
            "bbox_minlat":llat,
            "bbox_maxlon":geom.bounds[2],
            "bbox_maxlat":geom.bounds[3],
            "source_component_polygons":len(rows),
            "in_original_17_historical_port_botanical_cohort":str(landmass in ids),
        })
    frame.sort(key=lambda r:int(r["landmass_id"]))
    HERE.mkdir(parents=True,exist_ok=True)
    with (HERE/"global_island_frame_area20.csv").open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(frame[0]))
        writer.writeheader();writer.writerows(frame)
    present={q["landmass_id"] for q in frame}
    in17=ids&present
    areas=[q["area_km2_equal_area"] for q in frame]
    receipt={
      "schema":"chocho_global_GSHHG_island_frame_no_future_butterfly_selection_v1",
      "status":"SOURCE_VERIFIED_GLOBAL_ISLAND_SAMPLING_FRAME_NOT_BUTTERFLY_OUTCOME",
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "source_GSHHG_url":GSHHG_URL,"source_GSHHG_SHA256":sha,
      "original_17_islands_git_source":SOURCE17,
      "original_17_islands_git_blob":SOURCE17_BLOB,
      "source_gshhg_L1_polygon_rows":len(land),
      "source_unique_sibling_landmasses":int(land.landmass_id.nunique()),
      "size_threshold_km2":{"lower_inclusive":AREA_MIN,"upper_exclusive":AREA_MAX},
      "selected_islands_count":len(frame),
      "original_17_in_global_frame_count":len(in17),
      "original_17_missing_due_to_area_filter":sorted(ids-present,key=int),
      "min_area_km2":min(areas),
      "max_area_km2":max(areas),
      "candidate_butterflies_existing_GEB_island_225":225,
      "potential_species_x_island_risk_cells":225*len(frame),
      "three_period_GBIF_occurrence_searches_without_sharding":3*len(frame),
      "limits":[
        "GSHHG high-resolution Level 1 is not automatically equivalent to the 4663-island independent project geometry.",
        "Only one source-locked static coastline is used; no historical dynamic coast reconstruction.",
        "Equal-area filter does not prove oceanic isolation, continental distance or uninhabited island status.",
        "GSHHG sibling IDs are geographical locators; pre-2011 GBIF index is still needed for a valid risk denominator.",
        "Original 17-island study was a port-selected subset; global frame is not selected by future butterfly appearance.",
        "Potential species×island cells are not all climatically accessible and unindexed does NOT mean true absence.",
        "Full butterfly GBIF period records and matched historical plant resource information are separate future steps."
      ]
    }
    (HERE/"global_island_frame_receipt.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:receipt[k] for k in
        ("selected_islands_count","original_17_in_global_frame_count",
         "potential_species_x_island_risk_cells",
         "three_period_GBIF_occurrence_searches_without_sharding")},
      ensure_ascii=False,indent=2),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    args=p.parse_args()
    run(args.zip,args.work)

if __name__=="__main__":
    main()
