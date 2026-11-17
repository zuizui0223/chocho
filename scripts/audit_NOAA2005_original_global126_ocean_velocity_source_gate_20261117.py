#!/usr/bin/env python3
"""NOAA old OSCAR native historical 2005 GLOBAL velocity grid source coverage gate.

Authored for source-verified original 126 WCVP islands. First check actual
old 1/3degree NOAA source at every THIRD latitude and longitude (1degree
coarsening) and THREE original 2005 release dates, not 126 trajectory arrival.
If original (2019 or GSHHG) island coastal starts do not have valid source 1deg
u/v fields nearby, record missing instead of selecting an arbitrary island.

Old NOAA OSCAR ~5-day 15m surface u/v is distinct from Earthdata daily v2.
One degree global *subsampling* is source preflight ONLY; it does not validate
the resolution of nearshore ocean current, seed buoyancy, or migration.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request,urlopen

import netCDF4
import numpy as np
from pyproj import Geod

BASE="https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.nc?"
RELEASE_SOURCE_TIME_INDICES=(907,919,931)
SOURCE_NATIVE_GRID_STRIDE=3
SOURCE_1DEG_LATS=(0,480)
SOURCE_1DEG_LONS=(0,1200)
MAX_SOURCE_BYTES=16_000_000
COAST_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
G=Geod(ellps="WGS84")
MAX_OFFSHORE_START_KM=200.

def SHA(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for v in iter(lambda:f.read(4*1024*1024),b""):h.update(v)
    return h.hexdigest()

def NOAA():
    # 2005 May/July/Sep fixed original source-grid time indices, 1deg grid
    # is literal [0:3:480][0:3:1200] of true NOAA 1/3deg velocities.
    source=",".join(
      f"{z}[907:12:931][0][0:3:480][0:3:1200]" for z in ("u","v"))
    url=BASE+quote(source,safe="[],:")
    req=Request(url,headers={"User-Agent":"chocho-WCVP126-global-NOAA-original-source/2026.11"})
    with urlopen(req,timeout=150) as rsp:
        raw=rsp.read(MAX_SOURCE_BYTES+1)
        if rsp.status!=200 or len(raw)>MAX_SOURCE_BYTES:
            raise RuntimeError("Unexpected NOAA global source data response status or size")
    if not (raw.startswith(b"CDF") or raw.startswith(b"\x89HDF\r\n\x1a\n")):
        raise RuntimeError("NOAA global 2005 query did not return a NetCDF source grid")
    return raw

def decode(raw,output):
    rawpath=output/"original_NOAA_global_source_coarse_probe.nc"
    rawpath.write_bytes(raw)
    try:
        with netCDF4.Dataset(rawpath) as f:
            must={"u","v","time","latitude","longitude"}
            if not must.issubset(f.variables):
                raise RuntimeError("NOAA NetCDF original u/v variables absent")
            def arr(k):
                a=f.variables[k][:]
                return np.array(a.filled(np.nan) if np.ma.isMaskedArray(a) else a,dtype=float)
            la=arr("latitude");lo=arr("longitude")
            u=arr("u");v=arr("v")
            ts=netCDF4.num2date(f.variables["time"][:],
                units=f.variables["time"].units,
                calendar=getattr(f.variables["time"],"calendar","standard"))
            dates=[str(t) for t in ts]
            if u.shape!=(3,1,161,401) or v.shape!=u.shape:
                raise RuntimeError(f"NOAA original 1degree subset changed: u{u.shape},v{v.shape}")
            if not (abs(la[0]-80)<1e-5 and abs(la[-1]+80)<1e-5 and
                abs(lo[0]-20)<1e-5 and abs(lo[-1]-420)<1e-5):
                raise RuntimeError("NOAA gridded coordinates are not original 1/3 data")
            if any(not x.startswith("2005-") for x in dates):
                raise RuntimeError(f"NOAA original times are not 2005: {dates}")
            if not np.all(np.diff(la)<0) or not np.all(np.diff(lo)>0):
                raise RuntimeError("NOAA original lat/lon axes not monotonic")
            ok=np.isfinite(u)&np.isfinite(v)
            return {"lat":la,"lon":lo,"u":u[:,0],"v":v[:,0],"dates":dates,
                    "finite":int(ok.sum()),"missing":int((~ok).sum())}
    finally:
        rawpath.unlink(missing_ok=True)

def fine_wet(mask,lat,lon):
    c=int(round((((lon+180)%360)-180+180)/.0625-.5))%mask.shape[1]
    r=int(round((85-lat)/.0625-.5))
    return 0<=r<mask.shape[0] and mask[r,c]>=0

def candidates(source,anchor,season_index,mask):
    lat=source["lat"];lon=source["lon"];u=source["u"];v=source["v"]
    la=float(anchor["snapped_wet_latitude"])
    lo=float(anchor["snapped_wet_longitude"])
    if not (-80<=la<=80):
        return []
    lo360=lo if lo>=20 else lo+360
    ic=int(np.argmin(abs(lon-lo360)))
    ir=int(np.argmin(abs(lat-la)))
    out=[]
    for rr in range(max(0,ir-3),min(len(lat)-1,ir+4)):
        for cc in range(max(0,ic-3),min(len(lon)-1,ic+4)):
            x=float(lat[rr]);y=float(lon[cc])
            physical_lon=((y+180)%360)-180
            if not fine_wet(mask,x,physical_lon):continue
            mat_u=u[season_index,rr:rr+2,cc:cc+2]
            mat_v=v[season_index,rr:rr+2,cc:cc+2]
            if not (np.isfinite(mat_u).all() and np.isfinite(mat_v).all()):
                continue
            _,_,d=G.inv(lo,la,physical_lon,x)
            dist_km=float(d/1000)
            if dist_km>MAX_OFFSHORE_START_KM:continue
            out.append({
              "distance_from_original_GSHHG_coastal_wet_anchor_km":round(dist_km,4),
              "latitude":x,"longitude":physical_lon,
              "NOAA_1deg_original_source_lat_index":rr,
              "NOAA_1deg_original_source_lon_index":cc,
            })
    return sorted(out,key=lambda z:(z["distance_from_original_GSHHG_coastal_wet_anchor_km"],
                                    z["NOAA_1deg_original_source_lat_index"],
                                    z["NOAA_1deg_original_source_lon_index"]))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--GSHHG-graph",type=Path,required=True)
    p.add_argument("--anchors-receipt",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    if SHA(args.GSHHG_graph)!=COAST_SHA:raise RuntimeError("Original GSHHG 0.0625degree graph SHA mismatch")
    anchor_data=json.loads(args.anchors_receipt.read_text())
    if (anchor_data["source_original_126_regions"]!=126 or
        anchor_data["source_126_regions_with_at_least_one_verified_coastal_water_anchor"]!=104):
        raise RuntimeError("Original verified 126/GSHHG shore terminal source changed")
    with np.load(args.GSHHG_graph) as g:
        mask=g["fine_cell_to_node"]
        if mask.shape!=(2720,5760):raise RuntimeError("Original GSHHG global source-mask shape changed")
    raw=NOAA()
    src=decode(raw,args.out_dir)
    src_sha=hashlib.sha256(raw).hexdigest()
    matches=[]
    for region,anchors in sorted(anchor_data["coastal_anchors_by_TDWG_region"].items()):
        if not anchors:
            matches.append({
              "region":region,"noaa_spring_summer_autumn_launch_status":"ORIGINAL_GSHHG_SHORE_ANCHOR_NOT_AVAILABLE",
              "fixed_season_validities":0,"nearest_valid_offshore_source_km_per_season":[None]*3
            })
            continue
        starts=[]
        for k,idx in enumerate(RELEASE_SOURCE_TIME_INDICES):
            opts=[]
            for original in anchors:
                opts.extend(candidates(src,original,k,mask))
            chosen=min(opts,key=lambda x:x["distance_from_original_GSHHG_coastal_wet_anchor_km"]) if opts else None
            starts.append(chosen)
        matches.append({
           "region":region,"source_original_geometry_anchor_count":len(anchors),
           "noaa_spring_summer_autumn_launch_status":
              "ALL_THREE_SOURCE_DATES_VALID" if all(starts) else
              "PARTIAL_VALID_RELEASE_DATES" if any(starts) else "NO_NOAA_VALID_COASTAL_GRID_WITHIN_200KM",
           "fixed_season_validities":sum(x is not None for x in starts),
           "nearest_valid_offshore_source_km_per_season":[
                x["distance_from_original_GSHHG_coastal_wet_anchor_km"] if x else None for x in starts],
           "original_NOAA_1deg_eligible_source_lon_lat_index_per_season":[
                [x["NOAA_1deg_original_source_lat_index"],x["NOAA_1deg_original_source_lon_index"]] if x else None
                for x in starts],
        })
    counts=Counter(r["noaa_spring_summer_autumn_launch_status"] for r in matches)
    valid_all=sum(x["fixed_season_validities"]==3 for x in matches)
    valid_any=sum(x["fixed_season_validities"]>=1 for x in matches)
    if len(matches)!=126:raise RuntimeError("Original 126 island-source statuses not reproduced")
    result={
       "schema":"chocho_WCVP126_real_NOAA2005_old_OSCAR_global_1degree_source_preflight_v1",
       "status":"GLOBAL_NOAA_OLD_OSCAR_ACTUAL_UV_SOURCE_COVERAGE_NOT_126_ISLAND_PARTICLE_TRANSPORT",
       "NOAA_dataset":"yearly_336c_0b32_9cd3",
       "NOAA_product":"OLD public NOAA OSCAR source 1/3deg originally, subsampled every third original latitude and longitude to ~1deg, 15m, ~5day",
       "NOAA_not_NASA_daily_v2":True,
       "original_source_indices":{"time_release":RELEASE_SOURCE_TIME_INDICES,
          "latitude":[0,3,480],"longitude":[0,3,1200]},
       "noaa_actual_source_release_dates":src["dates"],
       "NOAA_source_NetCDF_original_subset_sha256":src_sha,
       "NOAA_source_NetCDF_bytes":len(raw),
       "NOAA_original_global_1deg_3dates_shape":[3,161,401],
       "NOAA_original_valid_u_and_v_cell_times":src["finite"],
       "NOAA_original_missing_u_or_v_cell_times":src["missing"],
       "source_GSHHG_graph_SHA256":COAST_SHA,
       "source_GSHHG_shore_graph_valid_islands":104,
       "all_original_WCVP_islands":126,
       "fixed_valid_offshore_launch_distance_max_km":MAX_OFFSHORE_START_KM,
       "source_islands_valid_noaa_all_three_release_dates":valid_all,
       "source_islands_valid_noaa_at_least_one_release_date":valid_any,
       "original_island_source_coverage":dict(counts),
       "126_source_region_release_eligibility":matches,
       "strict_126_original_sea_current_directed_arrival_probabilities_verified":False,
       "limits":[
          "This is a source-availability probe using three actual NOAA 2005 dates; no 90-day time-dependent velocity integration or landfall was performed.",
          "Global 1degree subsampling misses original 1/3degree coastal pathways and is NOT equivalent to the earlier Atlantic native 1/3degree trajectories.",
          "An offshore starting point up to 200km from original coastline is a mathematical data-coverage diagnostic, not a real seed release.",
          "Original 22 unanchored 126-island regions are geography missingness, NOT absent ocean current or uncolonisable plants.",
          "The 126 region current direction matrix cannot be claimed unless a separate actual verified historical NOAA 90-day surface Lagrangian run is completed.",
          "A 2005 current field cannot be used as contemporary WCVP plant origination event dates or 2020 shipping voyages.",
          "No plants, viable seeds, HOSTS feeding or butterfly populations were observed by these physical grids."
       ],
    }
    (args.out_dir/"NOAA2005_global126_coastal_current_source_eligibility.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    print("NOAA126_GLOBAL_2005_UV_SOURCE_COVERAGE",json.dumps({
       "status":result["status"],"valid_all3":valid_all,
       "valid_any":valid_any,"source_dates":src["dates"],
       "source_ocean_uv_cells":src["finite"],"missing_ocean_uv_cells":src["missing"],
       "conditions":dict(counts),
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
