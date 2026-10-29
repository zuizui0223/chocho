#!/usr/bin/env python3
"""Audited GSHHG 2.3.7 global island landmass census at fixed 5/20 km².

This counts original coast polygons after original sibling/ID dissolution,
following the source rules already used by zuizui0223/island, not arbitrary
H3 centers. No individual polygon coordinates or shapefile rows are exported.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile,math
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen
from collections import Counter
import geopandas as gpd
import numpy as np
from shapely import make_valid
from shapely.ops import unary_union
from audit_island_gshhg_original_source_20261029 import URL,SHA,MAX_SIZE

SCHEMA="chocho_GSHHG_237_original_island_area_isolation_universe_20261029_v01"
ISLAND_MAINLAND_THRESHOLD_KM2=7_000_000
PRIMARY_MIN_KM2=5
SENSITIVITY_MIN_KM2=20

def dissolve_original(frame):
    """Apply GSHHG hierarchy ID grouping, no merging different source islands."""
    if frame.crs is None:raise ValueError("Original source missing CRS")
    sibling=next((c for c in ("sibling_id","SIBLING_ID") if c in frame.columns),None)
    id_col=sibling or next((c for c in ("id","ID") if c in frame.columns),None)
    if not id_col:
        raise ValueError("GSHHG original source polygon ID field missing")
    original=frame[[id_col,"geometry"]].copy()
    ids=gpd.pd.to_numeric(original[id_col],errors="coerce")
    if not (np.isfinite(ids)&(ids>=0)&(ids%1==0)).all():
        raise ValueError("Original source island identity ambiguous")
    original[id_col]=ids.astype("int64")
    members=[]
    for key,group in original.groupby(id_col,sort=True):
        geo=make_valid(unary_union(group.geometry.to_list()))
        if geo.is_empty:continue
        members.append({"parent_source_ID":int(key),"geometry":geo,"n_original_parts":len(group)})
    return gpd.GeoDataFrame(members,geometry="geometry",crs=frame.crs)

def summarize(frame):
    grouped=dissolve_original(frame)
    geo_area=grouped.to_crs("EPSG:6933").area.to_numpy()/1e6
    if not np.isfinite(geo_area).all() or np.any(geo_area<0):
        raise ValueError("Original source island areas invalid")
    mainland=geo_area>ISLAND_MAINLAND_THRESHOLD_KM2
    islands_5=(geo_area>=PRIMARY_MIN_KM2)&(~mainland)
    islands_20=(geo_area>=SENSITIVITY_MIN_KM2)&(~mainland)
    if np.sum(islands_5)<1000 or np.sum(islands_20)<500 or np.sum(mainland)<3:
        raise ValueError("GSHHG source L1 island/continent census inconsistent")
    return {
        "source_GSHHG_L1_shapefile_rows":len(frame),
        "dissolved_original_landmass_units":len(grouped),
        "original_landmass_groups_with_multiple_source_features":int(sum(x>1 for x in grouped.n_original_parts)),
        "true_landmass_classes":{
            "continental_source_landmass_area_gt7million_km2":int(mainland.sum()),
            "islands_5km2_to_mainland_threshold":int(islands_5.sum()),
            "islands_20km2_to_mainland_threshold":int(islands_20.sum()),
            "nonmainland_source_polygons_smaller_than5km2":int(((geo_area<5)&(~mainland)).sum())
        },
        "source_primary_5km2_island_area_total_km2":float(geo_area[islands_5].sum()),
        "source_sensitivity_20km2_island_area_total_km2":float(geo_area[islands_20].sum()),
        "raw_island_geometry_or_source_IDs_exported":False,
        "island_distance_to_mainland_not_yet_measured":True,
        "GSHHG_L1_omits_Antarctica_mainland_coastlines":True,
    }

def run(dest):
    with urlopen(Request(URL,headers={"User-Agent":"chocho-source-island-census/1.0"}),timeout=120) as r:
        source=r.read(MAX_SIZE+1)
    if len(source)>MAX_SIZE or hashlib.sha256(source).hexdigest()!=SHA:
        raise ValueError("Verified GSHHG 2.3.7 source ZIP changed/too large")
    with TemporaryDirectory(prefix="gshhg237-") as td:
        with zipfile.ZipFile(io.BytesIO(source)) as z:
            matches=[x for x in z.namelist() if x.replace("\\","/").endswith("GSHHS_h_L1.shp")]
            if len(matches)!=1:raise ValueError("Original GSHHG L1 shapefile missing")
            root=matches[0][:-4]
            for suf in (".shp",".dbf",".shx",".prj"):
                target=Path(td)/("GSHHS_h_L1"+suf)
                with z.open(root+suf) as stream,target.open("wb") as out:
                    for block in iter(lambda:stream.read(1024*1024),b""):
                        out.write(block)
        frame=gpd.read_file(Path(td)/"GSHHS_h_L1.shp")
        result=summarize(frame)
    output={
        "schema":SCHEMA,
        "status":"ORIGINAL_GSHHG_L1_5_AND_20_KM2_ISLAND_CENSUS_VERIFIED",
        "original_source_zip_sha256":SHA,
        "original_source_gshhg_version":"2.3.7",
        "primary_area_threshold_km2":PRIMARY_MIN_KM2,
        "sensitivity_area_threshold_km2":SENSITIVITY_MIN_KM2,
        "mainland_area_threshold_km2":ISLAND_MAINLAND_THRESHOLD_KM2,
        "derived":result,
        "no_butterfly_or_ports_arrival_effect_estimated":True,
    }
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":output["status"],"counts":result["true_landmass_classes"],"n_dissolved":result["dissolved_original_landmass_units"]},sort_keys=True),flush=True)
    return output

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    run(p.parse_args().output)
