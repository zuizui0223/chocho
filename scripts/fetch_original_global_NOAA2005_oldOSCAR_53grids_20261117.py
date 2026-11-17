#!/usr/bin/env python3
"""Retrieve source-audited NOAA old OSCAR 2005 global 1degree decimated u/v.

The original 1/3-degree NOAA mirror is sub-selected at every THIRD original
lat/lon index (1°) only for a global Lagrangian FEASIBILITY diagnostic.
Physical 15m, ~5day currents do not model windage/seed buoyancy/survival.
This source is distinct from NASA daily OSCAR v2.0 (Earthdata 401).
No missing velocity cell is silently imputed or set to zero.
"""
from __future__ import annotations
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request,urlopen

import netCDF4
import numpy as np

BASE="https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.nc?"
CHUNKS=((898,910),(911,923),(924,936),(937,950))
RELEASE_INDICES=(907,919,931)
MAX_SOURCE_CHUNK_BYTES=22_000_000
ORIGINAL_PREP_CI=44011612861

def grab(lo,hi):
    expr=",".join(f"{z}[{lo}:1:{hi}][0][0:3:480][0:3:1200]" for z in ("u","v"))
    req=Request(BASE+quote(expr,safe="[],:"),
       headers={"User-Agent":"chocho-WCVP126-native-NOAA-global-2005-1deg-Lagrangian/1.0",
                "Accept":"application/x-netcdf,*/*"})
    with urlopen(req,timeout=150) as rsp:
        payload=rsp.read(MAX_SOURCE_CHUNK_BYTES+1)
        if rsp.status!=200 or len(payload)>MAX_SOURCE_CHUNK_BYTES:
            raise RuntimeError(f"Unexpected original NOAA source size/status for [{lo},{hi}]")
    if not (payload.startswith(b"CDF") or payload.startswith(b"\x89HDF\r\n\x1a\n")):
        raise RuntimeError("NOAA source did not return authentic NetCDF")
    return payload

def decode(path):
    with netCDF4.Dataset(path) as d:
        if not {"u","v","time","latitude","longitude"}.issubset(d.variables):
            raise RuntimeError("Source griddap lacks original NOAA variables")
        def arr(k):
            v=d.variables[k][:]
            v=v.filled(np.nan) if np.ma.isMaskedArray(v) else v
            return np.asarray(v,dtype=float)
        u=arr("u");v=arr("v");lat=arr("latitude");lon=arr("longitude")
        times=netCDF4.num2date(d.variables["time"][:],
            units=d.variables["time"].units,
            calendar=getattr(d.variables["time"],"calendar","standard"))
        epoch=np.array([datetime.fromisoformat(str(t).replace(" ","T")[:19]).replace(tzinfo=timezone.utc).timestamp()
                        for t in times],dtype=float)
        if u.shape!=v.shape or u.shape[1:]!=(1,161,401):
            raise RuntimeError(f"NOAA global original selected1degree dimensions changed: u={u.shape},v={v.shape}")
        if not (np.allclose(lat[[0,-1]],[80,-80]) and
            np.allclose(lon[[0,-1]],[20,420])):
            raise RuntimeError("NOAA original 1/3degree source coordinate strides incorrect")
        if not (np.all(np.diff(lat)<0) and np.all(np.diff(lon)>0) and
                np.all(np.diff(epoch)>0)):
            raise RuntimeError("NOAA axes unexpected or timestamp repeated")
        return u[:,0].astype(np.float32),v[:,0].astype(np.float32),lat,lon,epoch

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-npz",type=Path,required=True)
    ap.add_argument("--out-receipt",type=Path,required=True)
    ap.add_argument("--previous-three-date-gate",type=Path,required=True)
    a=ap.parse_args()
    previous=json.loads(a.previous_three_date_gate.read_text())
    if (previous["status"]!="GLOBAL_NOAA_OLD_OSCAR_ACTUAL_UV_SOURCE_COVERAGE_NOT_126_ISLAND_PARTICLE_TRANSPORT" or
        previous["source_islands_valid_noaa_at_least_one_release_date"]!=94):
        raise RuntimeError("Original 2005 global NOAA release-source cohort changed")
    a.out_npz.parent.mkdir(parents=True,exist_ok=True)
    chunks=[]
    coords=None
    allu=[];allv=[];allt=[]
    for i,(start,end) in enumerate(CHUNKS):
        data=grab(start,end)
        temp=a.out_npz.with_suffix(f".part{i}.nc")
        temp.write_bytes(data)
        try:
            u,v,lat,lon,epoch=decode(temp)
        finally:
            temp.unlink(missing_ok=True)
        if len(epoch)!=end-start+1:
            raise RuntimeError("NOAA original time-axis step count differs from pinned")
        if coords is None:coords=(lat,lon)
        elif not (np.array_equal(lat,coords[0]) and np.array_equal(lon,coords[1])):
            raise RuntimeError("NOAA original 2005 chunks disagree about grid")
        allu.append(u);allv.append(v);allt.append(epoch)
        chunks.append({"time_index":[start,end],
                       "raw_NetCDF_sha256":hashlib.sha256(data).hexdigest(),
                       "raw_source_bytes":len(data),
                       "finite_both_u_and_v":int(np.sum(np.isfinite(u)&np.isfinite(v)))})
        print("NOAA2005_GLOBAL_REAL_CHUNK",json.dumps(chunks[-1]),flush=True)
    u=np.concatenate(allu);v=np.concatenate(allv);times=np.concatenate(allt)
    if u.shape!=(53,161,401) or v.shape!=u.shape or len(times)!=53 or np.any(np.diff(times)<=0):
        raise RuntimeError("Original 53 source NOAA time steps and gridded current cube not recovered")
    if (datetime.fromtimestamp(times[0],tz=timezone.utc).strftime("%Y-%m-%d")!="2005-04-12" or
        datetime.fromtimestamp(times[-1],tz=timezone.utc).strftime("%Y-%m-%d")!="2006-01-01"):
        raise RuntimeError("Original 2005/early2006 source period changed")
    dates=[datetime.fromtimestamp(times[i-898],tz=timezone.utc).strftime("%Y-%m-%d")
           for i in RELEASE_INDICES]
    prev_dates=[x[:10] for x in previous["noaa_actual_source_release_dates"]]
    if dates!=prev_dates:
        raise RuntimeError(f"Full original NOAA cube does not agree with independent three-date source preflight: {dates} vs {prev_dates}")
    np.savez_compressed(a.out_npz,time_seconds_utc=times,latitude_deg=coords[0],
       longitude_deg_0_420=coords[1],u_m_s=u,v_m_s=v)
    source={
       "schema":"chocho_global_2005_original_NOAA_legacy_OSCAR_1degree_subsample_real_u_v_source_v1",
       "status":"NOAA_OLD_OSCAR_GLOBAL_2005_53_TIME_STEPS_AUTHENTIC_NOT_FULL_1_3_NATIVE_OR_BIOLOGICAL_TRANSPORT",
       "NOAA_dataset":"yearly_336c_0b32_9cd3","NOAA_original_1over3degree_stride":3,
       "NOAA_not_NASA_daily_OSCARv2":True,"true_source_depth_m":15,
       "time_indices":[898,950],"fixed_release_indices":RELEASE_INDICES,
       "fixed_release_dates":dates,
       "source_first_time_utc":datetime.fromtimestamp(times[0],tz=timezone.utc).isoformat(),
       "source_final_time_utc":datetime.fromtimestamp(times[-1],tz=timezone.utc).isoformat(),
       "true_source_shape":[53,161,401],
       "source_chunk_checksums":chunks,
       "source_finite_both_u_v":int(np.sum(np.isfinite(u)&np.isfinite(v))),
       "source_missing_one_or_both_u_v":int(np.sum(~(np.isfinite(u)&np.isfinite(v)))),
       "source_NPZ_bytes":a.out_npz.stat().st_size,
       "source_NPZ_sha256":hashlib.sha256(a.out_npz.read_bytes()).hexdigest(),
       "source_previous_coastal_dates_preflight_CI":ORIGINAL_PREP_CI,
       "limitations":[
          "This is authentic public OLD NOAA OSCAR 1/3degree gridded data every third cell (~1degree sample), NOT native 1/3degree nearshore flow and not NASA v2.0.",
          "The 15m Eulerian velocity does not represent buoyant seeds at air-sea surface, waves, tides, windage, seed survival or nearshore habitat.",
          "A 90-day deterministic particle run using this source must be explicitly labelled a coarse physical feasibility pilot, not confirmed real seed current transport.",
          "No 126-island ecological source or receiving species distribution is inferred by this data acquisition receipt."
       ]
    }
    a.out_receipt.parent.mkdir(parents=True,exist_ok=True)
    a.out_receipt.write_text(json.dumps(source,indent=2)+"\n")
    print("NOAA2005_GLOBAL_53_SOURCE_GRID_VERIFIED",json.dumps({
       "status":source["status"],"source_shape":source["true_source_shape"],
       "fixed_dates":dates,"valid":source["source_finite_both_u_v"],
       "missing":source["source_missing_one_or_both_u_v"],
       "NPZ_bytes":source["source_NPZ_bytes"]
    }),flush=True)

if __name__=="__main__":main()
