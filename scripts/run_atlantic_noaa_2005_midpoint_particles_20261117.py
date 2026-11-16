#!/usr/bin/env python3
"""Real NOAA OSCAR legacy 2005 ocean-only deterministic Lagrangian pilot.

Physics: 1-day midpoint RK2, interpolated actual 5-day u,v on original 1/3°
grid. Starts are explicitly offshore from source-checked WPI harbor coordinates,
NOT boat routes or botanical collection sites. GSHHG sea-mask segments terminate
when hitting resolved land. Missing u/v terminates (never zero-filled).
No currents beyond OSCAR, windage, saltwater survival, seed source/arrival,
actual coastline settlement, butterfly population or invasion inferred.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter,defaultdict

import numpy as np

SOURCES=("AZO","CNY","CVI")
ANCHORED_WPI_ID={"AZO":"3584","CNY":"428","CVI":"396"}
TIME_RELEASE_INDICES=(907,919,931)
DAY=86400.0
RUN_DAYS=90
DAY_CHECKPOINTS=(30,60,90)
PREFERRED_NEAR_PORT_KM=(25,50,100)
SOURCE_MAX_START_OFFSET_KM=80.
MAX_STEP_SUBSEGMENT_KM=2.
MASK_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
EARTH_RADIUS_KM=6371.0088

def read_csv(p):
    with p.open("r",newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def latlon_dist(a,b):
    lat1,lon1=a;lat2,lon2=b
    p,q=math.radians(lat1),math.radians(lat2)
    dp=math.radians(lat2-lat1)
    dl=math.radians((lon2-lon1+180)%360-180)
    h=math.sin(dp/2)**2+math.cos(p)*math.cos(q)*math.sin(dl/2)**2
    return 2*EARTH_RADIUS_KM*math.asin(math.sqrt(max(0,min(1,h))))

def integrate_vector(lat,lon,u,v,seconds,sign=1):
    # Coordinate step only; no latitude/longitude boundary handling here.
    dlat=sign*v*seconds/(1000*(math.pi/180)*EARTH_RADIUS_KM)
    coslat=math.cos(math.radians(lat))
    if abs(coslat)<1e-4:raise ValueError("Near-pole movement not supported")
    dlon=sign*u*seconds/(1000*(math.pi/180)*EARTH_RADIUS_KM*coslat)
    return lat+dlat,((lon+dlon+180)%360)-180

def uv_interpolate(T,lat_axis,lon_axis,u,v,second,lat,lon):
    # No interpolation across a missing value; NOAA source no longer observed
    if second<T[0] or second>T[-1]:return None
    lon360=lon if lon>=20 else lon+360
    if lat>lat_axis[0] or lat<lat_axis[-1] or lon360<lon_axis[0] or lon360>lon_axis[-1]:
        return None
    ti=max(0,min(len(T)-2,int(np.searchsorted(T,second,side="right")-1)))
    ji=int(np.searchsorted(-lat_axis,-lat,side="right")-1)
    ki=int(np.searchsorted(lon_axis,lon360,side="right")-1)
    ji=max(0,min(len(lat_axis)-2,ji))
    ki=max(0,min(len(lon_axis)-2,ki))
    td=(second-T[ti])/(T[ti+1]-T[ti])
    ld=(lat_axis[ji]-lat)/(lat_axis[ji]-lat_axis[ji+1])
    od=(lon360-lon_axis[ki])/(lon_axis[ki+1]-lon_axis[ki])
    weights=np.array([(1-ld)*(1-od),(1-ld)*od,ld*(1-od),ld*od])
    def read_component(cube):
        a=cube[ti:ti+2,ji:ji+2,ki:ki+2].reshape(2,4)
        if not bool(np.all(np.isfinite(a))):return None
        return float((1-td)*np.dot(a[0],weights)+td*np.dot(a[1],weights))
    east=read_component(u)
    if east is None:return None
    north=read_component(v)
    if north is None:return None
    if max(abs(east),abs(north))>5:raise RuntimeError("Observed source velocity is physically implausible")
    return east,north

def wet_pixel(fine_mask,lat,lon):
    if not math.isfinite(lat+lon):return False
    r=int(round((85-lat)/.0625-.5))
    c=int(round((lon+180)/.0625-.5))
    if not (0<=r<fine_mask.shape[0] and 0<=c<fine_mask.shape[1]):return False
    return bool(fine_mask[r,c]>=0)

def crosses_land(fine_mask,a,b):
    dist=latlon_dist(a,b)
    samples=max(1,int(math.ceil(dist/MAX_STEP_SUBSEGMENT_KM)))
    for i in range(1,samples+1):
        w=i/samples
        lat=a[0]*(1-w)+b[0]*w
        lon=a[1]+((b[1]-a[1]+180)%360-180)*w
        lon=((lon+180)%360)-180
        if not wet_pixel(fine_mask,lat,lon):return True,(lat,lon)
    return False,b

def valid_source_release(T,la,lo,u,v,fine_mask,timeindex,anchor):
    source_time=T[timeindex]
    candidates=[]
    near_i=int(np.argmin(abs(la-anchor["latitude"])))
    port360=anchor["longitude"] if anchor["longitude"]>=20 else anchor["longitude"]+360
    near_j=int(np.argmin(abs(lo-port360)))
    for i in range(max(0,near_i-4),min(len(la),near_i+5)):
        for j in range(max(0,near_j-4),min(len(lo),near_j+5)):
            lat=float(la[i]);lon=float(((lo[j]+180)%360)-180)
            d=latlon_dist((anchor["latitude"],anchor["longitude"]),(lat,lon))
            if d>SOURCE_MAX_START_OFFSET_KM or not wet_pixel(fine_mask,lat,lon):continue
            # Must interpolate with real surrounding water values at this time.
            if uv_interpolate(T,la,lo,u,v,source_time,lat,lon) is None:continue
            candidates.append((d,i,j,lat,lon))
    candidates.sort(key=lambda x:(x[0],x[1],x[2]))
    selected=[]
    for p in candidates:
        if all(latlon_dist((p[3],p[4]),(q[3],q[4]))>15 for q in selected):
            selected.append(p)
        if len(selected)==3:break
    return selected

def validate_fixture():
    # 2x2 finite uniform east current, and a missing corner must be terminal.
    t=np.array([0.,432000.])
    la=np.array([40.,39.])
    lo=np.array([330.,331.])
    u=np.full((2,2,2),.2)
    v=np.full((2,2,2),0.0)
    if uv_interpolate(t,la,lo,u,v,216000.,39.5,-29.5)!=(.2,0.):
        raise RuntimeError("Bilinear source interpolation control failed")
    moved=integrate_vector(39.5,-29.5,.2,0,DAY)
    if not moved[1]>-29.5 or abs(moved[0]-39.5)>1e-8:
        raise RuntimeError("Eastward 1-day advection control failed")
    u[0,0,0]=np.nan
    if uv_interpolate(t,la,lo,u,v,216000.,39.5,-29.5) is not None:
        raise RuntimeError("Masked-source interpolation incorrectly filled a missing corner")

def track(scenario,anchor,launch,start_seconds,T,la,lo,u,v,fine_mask,other_ports):
    lat,lon=launch[3],launch[4]
    start=(lat,lon)
    date0=datetime.fromtimestamp(float(start_seconds),tz=timezone.utc).isoformat()
    nearest={k:latlon_dist(start,(v["latitude"],v["longitude"]))
             for k,v in other_ports.items()}
    traj=[{"day":0,"latitude":lat,"longitude":lon,"status":"LAUNCHED_AT_VERIFIED_OFFSHORE_GRID"}]
    outcome="COMPLETED_90_DAYS"
    failure_day=None
    km=0.
    for day in range(1,RUN_DAYS+1):
        sec=float(start_seconds)+(day-1)*DAY
        if scenario=="STATIONARY_NO_CURRENT":
            end=(lat,lon)
        else:
            v0=uv_interpolate(T,la,lo,u,v,sec,lat,lon)
            if v0 is None:
                outcome="SOURCE_VELOCITY_MISSING";failure_day=day;break
            half=integrate_vector(lat,lon,*v0,DAY/2,sign=(-1 if scenario=="MATHEMATICAL_REVERSED_FLOW" else 1))
            vm=uv_interpolate(T,la,lo,u,v,sec+DAY/2,*half)
            if vm is None:
                outcome="SOURCE_MIDPOINT_MISSING";failure_day=day;break
            end=integrate_vector(lat,lon,*vm,DAY,sign=(-1 if scenario=="MATHEMATICAL_REVERSED_FLOW" else 1))
        crossed,last=crosses_land(fine_mask,(lat,lon),end)
        if crossed:
            outcome="GSHHG_FINE_MASK_LAND_CONTACT";failure_day=day
            traj.append({"day":day,"latitude":last[0],"longitude":last[1],"status":outcome})
            break
        km+=latlon_dist((lat,lon),end)
        lat,lon=end
        for name,port in other_ports.items():
            nearest[name]=min(nearest[name],latlon_dist((lat,lon),(port["latitude"],port["longitude"])))
        if day%5==0 or day in DAY_CHECKPOINTS:
            traj.append({"day":day,"latitude":lat,"longitude":lon,"status":"OFFSHORE_CONTINUED"})
    dist=latlon_dist(start,(lat,lon))
    return {
        "source":anchor["WGSRPD3"],"source_WPI_id":anchor["WPI_id"],
        "source_WPI_name":anchor["name"],
        "scenario":scenario,"start_time_utc":date0,
        "source_offshore_offset_km":round(launch[0],3),
        "launch_grid_lat_index":int(launch[1]),
        "launch_grid_lon_index":int(launch[2]),
        "start_lat":start[0],"start_lon":start[1],
        "status":outcome,"days_computed":day if outcome=="COMPLETED_90_DAYS" else (day-1),
        "failure_day":failure_day,"final_lat":lat,"final_lon":lon,
        "displacement_km":round(dist,3),"path_length_km":round(km,3),
        "minimum_distances_to_other_WPI_ports_km":{k:round(v,3) for k,v in nearest.items()},
        "within_25km_of_other_WPI_port_proxy":any(v<=25 for v in nearest.values()),
        "within_50km_of_other_WPI_port_proxy":any(v<=50 for v in nearest.values()),
        "trajectory":traj,
    }

def load_source_ports(path):
    selected={}
    for x in read_csv(path):
        k=x["WGSRPD3"]
        if k in SOURCES and x["WPI_source_OBJECTID"]==ANCHORED_WPI_ID[k] and
           x["original_WPI_2019_coordinate_inside_exact_TDWG3_polygon"]=="1":
            if k in selected:raise RuntimeError("Duplicate pinned WPI port ID")
            selected[k]={"WGSRPD3":k,"WPI_id":x["WPI_source_OBJECTID"],
                         "name":x["port_name"],"country":x["country"],
                         "latitude":float(x["latitude"]),"longitude":float(x["longitude"])}
    if set(selected)!=set(SOURCES):raise RuntimeError("Strict WPI anchors not all present")
    return selected

def main(a):
    validate_fixture()
    with np.load(a.noaa_cube) as d:
        T=d["time_seconds_utc"].astype(float)
        la=d["latitude_deg"].astype(float)
        lo=d["longitude_deg_0_420"].astype(float)
        u=d["u_m_s"].astype(float)
        v=d["v_m_s"].astype(float)
    if u.shape!=(53,126,201) or v.shape!=u.shape or len(T)!=53:
        raise RuntimeError(f"NOAA original 2005 source cube unexpectedly shaped {u.shape}")
    receipt=json.loads(a.source_receipt.read_text())
    if hashlib.sha256(a.noaa_cube.read_bytes()).hexdigest()!=receipt["npz_sha256"] or
       receipt["status"]!="AUTHENTIC_NOAA_LEGACY_2005_UV_SOURCE_CUBE_NO_TRANSPORT_INFERENCE":
       raise RuntimeError("NOAA source data checksum or stage changed")
    if hashlib.sha256(a.coast_npz.read_bytes()).hexdigest()!=MASK_SHA:
        raise RuntimeError("Prior frozen GSHHG 0.0625 degree graph identity mismatch")
    with np.load(a.coast_npz) as g:
        mask=g["fine_cell_to_node"]
        fine_lat=g["fine_latitude"];fine_lon=g["fine_longitude"]
        if mask.shape!=(2720,5760) or not np.isclose(fine_lat[0],84.96875) or
           not np.isclose(fine_lon[0],-179.96875):
            raise RuntimeError("Original GSHHG fine-grid definition mismatch")
    source_ports=load_source_ports(a.strict_source_ports)
    tracks=[]
    for i in TIME_RELEASE_INDICES:
        rel=i-898
        if rel<0 or rel>=len(T):raise RuntimeError("Release time index outside declared 2005 data")
        for src in SOURCES:
            port=source_ports[src]
            launches=valid_source_release(T,la,lo,u,v,mask,rel,port)
            if len(launches)<1:
                tracks.append({"source":src,"status":"NO_SOURCE_VALID_OFFSHORE_LAUNCH", "release_index":i})
                continue
            other={k:v for k,v in source_ports.items() if k!=src}
            for j,launch in enumerate(launches):
                for scenario in ("OBSERVED_NOAA_FLOW","MATHEMATICAL_REVERSED_FLOW","STATIONARY_NO_CURRENT"):
                    z=track(scenario,port,launch,T[rel],T,la,lo,u,v,mask,other)
                    z.update({"NOAA_source_time_index":i,"launch_replicate_id":j+1,
                              "is_biological_observation":False})
                    tracks.append(z)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"particles_diagnostic_tracks.json").open("w",encoding="utf8") as f:
        json.dump(tracks,f,ensure_ascii=False,indent=2,allow_nan=False)
    rows=[{k:z[k] for k in
         ("source","source_WPI_id","scenario","start_time_utc",
           "source_offshore_offset_km","status","days_computed",
           "failure_day","displacement_km","path_length_km",
           "within_25km_of_other_WPI_port_proxy","within_50km_of_other_WPI_port_proxy")}
          for z in tracks if "scenario" in z]
    with (a.out_dir/"lagrangian_source_only_summary.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    counts={s:dict(Counter(r["status"] for r in tracks if r.get("scenario")==s))
        for s in ("OBSERVED_NOAA_FLOW","MATHEMATICAL_REVERSED_FLOW","STATIONARY_NO_CURRENT")}
    stages={}
    for src in SOURCES:
        stages[src]={}
        for scen in counts:
            rs=[r for r in tracks if r.get("source")==src and r.get("scenario")==scen]
            stages[src][scen]={
                "tracks":len(rs),
                "90d_completed":sum(r["status"]=="COMPLETED_90_DAYS" for r in rs),
                "source_missing":sum("MISSING" in r["status"] for r in rs),
                "land_mask_contacts":sum(r["status"]=="GSHHG_FINE_MASK_LAND_CONTACT" for r in rs),
                "within_25km_other_WPI_port_proxy":sum(r["within_25km_of_other_WPI_port_proxy"] for r in rs),
                "within_50km_other_WPI_port_proxy":sum(r["within_50km_of_other_WPI_port_proxy"] for r in rs),
                "median_displacement_km":float(np.median([r["displacement_km"] for r in rs])) if rs else None,
            }
    audit={
        "status":"ACTUAL_2005_LEGACY_OSCAR_MIDPOINT_LAGRANGIAN_SOURCE_PILOT_NOT_SEED_ARRIVAL",
        "source_NetCDF":receipt["dataset"],
        "source_CI_netCDF_gate":43900899514,
        "NOAA_not_NASA_OSCAR_v2":True,
        "source_NetCDF_file_sha256":receipt["npz_sha256"],
        "source_historic_CI":43897982755,
        "original_WPI_TDWG3_ports_CI":43895414730,
        "source_GSHHG_fine_ocean_mask_SHA256":MASK_SHA,
        "source_GSHHG_fine_grid_degrees":0.0625,
        "fixed_origins":source_ports,"release_time_original_indices":TIME_RELEASE_INDICES,
        "max_float_duration_days":RUN_DAYS,"checkpoint_days":DAY_CHECKPOINTS,
        "advection_method":"1 day explicit midpoint RK2 on bilinear space/time interpolated NOAA 5-day actual u/v at 15m",
        "source_valid_launch_rule":"up to three nearest distinct wet NOAA cell centers within 80 km of original strict WPI port, bilinear source valid at release",
        "land_collision_rule":"each daily projected segment sub-sampled <=2km; stop on GSHHG fine-grid land pixel",
        "missing_current_rule":"fail-closed stop, never replace by zero",
        "positive_destination_proxy_rule":"minimum distance to other strictly polygon-matched WPI port locations <=25km/50km; NOT actual complete-island shoreline arrival",
        "scenario_counts":counts,"source_summaries":stages,
        "source_only_not_butterfly_outcomes":True,
        "limitations":[
          "WPI ports are geographic anchors, not seed collection sites or historical ship services.",
          "NOAA legacy OSCAR 15m ~5-day 1/3degree ocean currents omit windage, Stokes drift, coastal tidal eddies, land-sea subgrid effects.",
          "Reversed currents are a mathematical diagnostic, not an independent flow model or statistical placebo.",
          "Nearby 1/3 degree initial positions and season release times are not independent replicate islands.",
          "GSHHG 0.0625degree land mask approximates coastline and can miss narrow islands or create contact classification ambiguity.",
          "A <25km WPI proxy is not seed landing on an island; a missing encounter does not prove dispersal failure.",
          "A seed's floating ability, saltwater survival and habitat establishment were not observed or assumed.",
          "The original GEB butterfly +54.9%, island +118.9%, +4.12% forecasts are different units; this flow scenario does not test butterflies."
        ]
    }
    (a.out_dir/"lagrangian_source_quality.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    print("ATLANTIC_ACTUAL_OSCAR_2005_LAGRANGIAN_PILOT",json.dumps({
        "status":audit["status"],"tracks":len(tracks),
        "launch_failures":sum(z.get("status")=="NO_SOURCE_VALID_OFFSHORE_LAUNCH" for z in tracks),
        "scenario_counts":counts,"source_summaries":stages
    },ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--noaa-cube",required=True,type=Path)
    parser.add_argument("--source-receipt",required=True,type=Path)
    parser.add_argument("--coast-npz",required=True,type=Path)
    parser.add_argument("--strict-source-ports",required=True,type=Path)
    parser.add_argument("--out-dir",required=True,type=Path)
    main(parser.parse_args())
