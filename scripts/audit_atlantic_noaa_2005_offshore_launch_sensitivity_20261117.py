#!/usr/bin/env python3
"""Post-hoc offshore-launch-distance sensitivity, actual 2005 NOAA old OSCAR.

Keep the original real 2005 current field, frozen WPI coastal source anchors
and GSHHG 0.0625 degree sea mask. Pick nearest 3 qualifying wet NOAA grid
centers within predeclared distance bands independently from outcomes.
This is a launch-assumption sensitivity, NOT observed seeds/colonization.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from run_atlantic_noaa_2005_midpoint_particles_20261117 import (
  SOURCES,TIME_RELEASE_INDICES,DAY,load_source_ports,latlon_dist,wet_pixel,
  uv_interpolate,track, ANCHORED_WPI_ID
)
# Deliberately do not recode the original <80 km nearest-three design.
BANDS=((0.,40.),(40.,80.),(80.,160.))
COAST_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
T0=898

def candidate_band(T,lat,lon,u,v,mask,release_index,anchor,low,high):
    site=anchor["latitude"],anchor["longitude"]
    cenlon=anchor["longitude"] if anchor["longitude"]>=20 else anchor["longitude"]+360
    ii=int(np.argmin(abs(lat-site[0])));jj=int(np.argmin(abs(lon-cenlon)))
    # <=160 km at all selected origins, 8 source-grid cells cover >200 km.
    sites=[]
    for i in range(max(0,ii-9),min(len(lat),ii+10)):
        for j in range(max(0,jj-9),min(len(lon),jj+10)):
            x=float(lat[i]);y=float(((lon[j]+180)%360)-180)
            d=latlon_dist(site,(x,y))
            if not (low<=d<high) or not wet_pixel(mask,x,y):continue
            if uv_interpolate(T,lat,lon,u,v,T[release_index],x,y) is None:continue
            sites.append((d,i,j,x,y))
    sites.sort(key=lambda q:(q[0],q[1],q[2]))
    chosen=[]
    for x in sites:
        if all(latlon_dist((x[3],x[4]),(z[3],z[4]))>15 for z in chosen):
            chosen.append(x)
        if len(chosen)==3:break
    return chosen

def main(a):
    r=json.loads(a.source_receipt.read_text())
    if hashlib.sha256(a.noaa_cube.read_bytes()).hexdigest()!=r["npz_sha256"]:
        raise RuntimeError("NOAA source cube does not match original run 43902616033")
    if hashlib.sha256(a.coast_npz.read_bytes()).hexdigest()!=COAST_SHA:
        raise RuntimeError("GSHHG fine-grid mask no longer original")
    with np.load(a.noaa_cube) as q:
        T=q["time_seconds_utc"].astype(float);lat=q["latitude_deg"];lon=q["longitude_deg_0_420"]
        u=q["u_m_s"].astype(float);v=q["v_m_s"].astype(float)
    with np.load(a.coast_npz) as g:
        mask=g["fine_cell_to_node"]
    if u.shape!=(53,126,201) or mask.shape!=(2720,5760):
        raise RuntimeError("Unexpected actual NOAA cube or coast-mask cohort")
    anchors=load_source_ports(a.strict_source_ports)
    records=[]
    availability=[]
    for idx in TIME_RELEASE_INDICES:
        rel=idx-T0
        assert 0<=rel<len(T)
        for src in SOURCES:
            anchor=anchors[src]
            others={k:p for k,p in anchors.items() if k!=src}
            for low,high in BANDS:
                candidates=candidate_band(T,lat,lon,u,v,mask,rel,anchor,low,high)
                band=f"{int(low)}_{int(high)}km"
                availability.append({"source":src,"NOAA_time_index":idx,
                    "distance_band":band,"valid_offshore_release_grid_cells":len(candidates)})
                for j,seed in enumerate(candidates):
                    z=track("OBSERVED_NOAA_FLOW",anchor,seed,T[rel],T,lat,lon,u,v,mask,others)
                    records.append({
                      "source":src,"NOAA_time_index":idx,"band":band,
                      "release_grid_replicate":j+1,
                      "start_dist_from_port_km":round(seed[0],3),
                      "status":z["status"],"days_computed":z["days_computed"],
                      "failure_day":z["failure_day"],
                      "endpoint_displacement_km":z["displacement_km"],
                      "endpoint_path_length_km":z["path_length_km"],
                      "within_50km_other_WPI_port_proxy":z["within_50km_of_other_WPI_port_proxy"],
                      "within_25km_other_WPI_port_proxy":z["within_25km_of_other_WPI_port_proxy"],
                      "model_control":"OFFSHORE_SOURCE_ASSUMPTION_ONLY_NOT_SEED_RELEASE",
                    })
    if not records:raise RuntimeError("Original NOAA source did not yield any offshore release")
    byband={}
    for low,high in BANDS:
        band=f"{int(low)}_{int(high)}km"
        z=[x for x in records if x["band"]==band]
        state=Counter(x["status"] for x in z)
        byband[band]={
            "potential_source_origins_times":len(SOURCES)*len(TIME_RELEASE_INDICES),
            "source_records":len(z),
            "source_grid_shortfall":sum(max(0,3-a["valid_offshore_release_grid_cells"])
                for a in availability if a["distance_band"]==band),
            "status_counts":dict(state),
            "90d_completed":sum(x["days_computed"]==90 for x in z),
            "30d_completed":sum(x["days_computed"]>=30 for x in z),
            "60d_completed":sum(x["days_computed"]>=60 for x in z),
            "within_50km_other_WPI_port_proxy":sum(x["within_50km_other_WPI_port_proxy"] for x in z),
            "median_start_distance_from_source_port_km":
                round(float(np.median([x["start_dist_from_port_km"] for x in z])),3) if z else None,
        }
    result={
      "schema":"chocho_2005_oscar_offshore_distance_band_sensitivity_v1",
      "status":"POSTHOC_OFFSHORE_LAUNCH_ASSUMPTION_SOURCE_SENSITIVITY_NO_BIOLOGICAL_ARRIVAL",
      "original_real_source_run":43902616033,"original_21-source_AZO_CNY_CVI_fixed_port_ids":ANCHORED_WPI_ID,
      "coast_mask_source_GSHHG_SHA256":COAST_SHA,
      "data_status":"No missing velocity imputed; no coastline mask altered",
      "release_distance_bands_km":[list(x) for x in BANDS],
      "within_band_selection_rule":"up to 3 nearest distinct GSHHG-wet, NOAA-bilinear-valid 1/3degree source cells for each original port and three fixed time indices; NO outcome-based selection",
      "all_source_cases":len(records),"availability":availability,
      "bands":byband,"records":records,
      "limitations":[
        "Offshore release is a hypothetical artificial-source transport assumption, NOT observed release of floating viable seed.",
        "This is an explicit post-hoc source-distance sensitivity after observing 19/27 coast-mask contacts in the nearshore pilot.",
        "Start cells are only NOAA 1/3degree resolution; distance bands may have uneven source availability.",
        "The three temporal and three spatial nearby start sites are not independent islands or propagules.",
        "Land-mask returns are not necessarily source island, and WPI distant port proximity is not full island landfall.",
        "No species-specific buoyancy duration, sea survival, actual current particle time spans beyond 90 days, verified shipping OD, or butterfly arrival."
      ],
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("NOAA_OFFSHORE_RELEASE_DISTANCE_SENSITIVITY",json.dumps({
        "status":result["status"],"source_case_count":result["all_source_cases"],
        "by_band":byband,"source_shortfalls":availability
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--noaa-cube",required=True,type=Path)
    p.add_argument("--source-receipt",required=True,type=Path)
    p.add_argument("--coast-npz",required=True,type=Path)
    p.add_argument("--strict-source-ports",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
