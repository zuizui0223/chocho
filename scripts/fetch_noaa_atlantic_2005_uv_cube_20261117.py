#!/usr/bin/env python3
"""Fetch bounded historical legacy NOAA OSCAR NetCDF field for a Lagrangian pilot.

NOAA legacy OSCAR ~5-day, 1/3-degree 15-m gridded u/v. Not NASA OSCAR v2,
not wind drift, not a particle trajectory, and not measured seed arrival.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import netCDF4
import numpy as np

BASE="https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.nc?"
TIME_CHUNKS=((898,915),(916,933),(934,950))
LAT_INDICES=(105,230) # 45N to 3.333N
LON_INDICES=(865,1065) # 308.333E (-51.667) to 375E (+15), no seam
EXPECTED="yearly_336c_0b32_9cd3"
MAX_CHUNK=8_000_000

def downloaded_chunk(a,b):
    lat0,lat1=LAT_INDICES;lon0,lon1=LON_INDICES
    expr=",".join(f"{var}[{a}:1:{b}][0][{lat0}:1:{lat1}][{lon0}:1:{lon1}]" for var in ("u","v"))
    url=BASE+quote(expr,safe="[],:")
    req=Request(url,headers={"User-Agent":"chocho-atlantic-Lagrangian-source/2026.11"})
    with urlopen(req,timeout=120) as rsp:
        raw=rsp.read(MAX_CHUNK+1)
        if rsp.status!=200 or len(raw)>MAX_CHUNK:
            raise RuntimeError(f"NOAA HTTP size/status unexpected chunk {a}:{b}: http={rsp.status}, bytes={len(raw)}")
    if not (raw.startswith(b"CDF") or raw.startswith(b"\x89HDF\r\n\x1a\n")):
        raise RuntimeError(f"NOAA returned non-NetCDF for source indices {a}:{b}")
    return raw

def asfloat(v):
    z=v[:]
    if np.ma.isMaskedArray(z):z=z.filled(np.nan)
    z=np.asarray(z,dtype=np.float64)
    z[~np.isfinite(z)]=np.nan
    return z

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-npz",required=True,type=Path)
    ap.add_argument("--receipt",required=True,type=Path)
    a=ap.parse_args()
    a.out_npz.parent.mkdir(parents=True,exist_ok=True)
    tt=[];uall=[];vall=[];checks=[];lats=lons=None
    for i,(beg,end) in enumerate(TIME_CHUNKS):
        raw=downloaded_chunk(beg,end)
        tmp=a.out_npz.with_suffix(f".part{i}.nc")
        tmp.write_bytes(raw)
        try:
            with netCDF4.Dataset(tmp) as nc:
                names={"u","v","latitude","longitude","time"}
                if not names.issubset(nc.variables):raise RuntimeError(f"Source dimensions missing: {list(nc.variables)}")
                la=asfloat(nc.variables["latitude"])
                lo=asfloat(nc.variables["longitude"])
                t=nc.variables["time"]
                decoded=netCDF4.num2date(t[:],units=t.units,
                     calendar=getattr(t,"calendar","standard"))
                times=np.array([datetime.fromisoformat(str(z).replace(" ","T")[:19]).replace(tzinfo=timezone.utc).timestamp() for z in decoded],dtype=np.float64)
                u=asfloat(nc.variables["u"]);v=asfloat(nc.variables["v"])
                if u.shape != (end-beg+1,1,len(la),len(lo)) or v.shape!=u.shape:
                    raise RuntimeError(f"Wrong grid dimensions from NOAA: {u.shape}, {v.shape}")
                if len(la)!=126 or len(lo)!=201:raise RuntimeError("Grid latitude/longitude count changed")
                if np.any(np.diff(la)>=0) or np.any(np.diff(lo)<=0) or np.any(np.diff(times)<=0):
                    raise RuntimeError("Source axes are not monotonic")
                if not (np.allclose(la[0],80-LAT_INDICES[0]/3,atol=.005) and
                        np.allclose(lo[0],20+LON_INDICES[0]/3,atol=.005)):
                    raise RuntimeError("Source longitudes or latitudes shifted unexpectedly")
                if lats is None:lats,lons=la,lo
                elif not (np.array_equal(la,lats) and np.array_equal(lo,lons)):
                    raise RuntimeError("Spatial coordinates differ across exact NOAA chunks")
                tt.append(times);uall.append(u[:,0]);vall.append(v[:,0])
                checks.append({"indices":[beg,end],"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),
                         "finite_both":int(np.sum(np.isfinite(u)&np.isfinite(v)))})
        finally:
            tmp.unlink(missing_ok=True)
    ts=np.concatenate(tt);u=np.concatenate(uall).astype(np.float32);v=np.concatenate(vall).astype(np.float32)
    if len(ts)!=53 or u.shape!=(53,126,201) or v.shape!=u.shape:raise RuntimeError("Final source cube size changed")
    if np.any(np.diff(ts)<=0):raise RuntimeError("Source timestamp discontinuity or reversal")
    if not all(datetime.fromtimestamp(float(x),tz=timezone.utc).year==2005 for x in ts):
        raise RuntimeError("Source subset extended outside 2005, revise scope")
    # Old independent observed NOAA point: 2005-09-21, 30N, 150E is outside the
    # Atlantic grid, so validate 930 38.333N, 332E from the separate NetCDF IO.
    t930=930-TIME_CHUNKS[0][0];la124=124-LAT_INDICES[0];lo936=936-LON_INDICES[0]
    if datetime.fromtimestamp(ts[t930],tz=timezone.utc).strftime("%Y-%m-%d")!="2005-09-21":
        raise RuntimeError("Source time indices do not recover independent 2005-09-21 sample")
    if not (np.isfinite(u[t930,la124,lo936]) and np.isfinite(v[t930,la124,lo936])):
        raise RuntimeError("Independent test window no longer has valid source")
    np.savez_compressed(a.out_npz,time_seconds_utc=ts,latitude_deg=lats,
            longitude_deg_0_420=lons,u_m_s=u,v_m_s=v)
    d={
      "status":"AUTHENTIC_NOAA_LEGACY_2005_UV_SOURCE_CUBE_NO_TRANSPORT_INFERENCE",
      "dataset":EXPECTED,"noaa_erddap_griddap":"https://oceanwatch.pifsc.noaa.gov/erddap/griddap/"+EXPECTED,
      "product":"NOAA mirrored old OSCAR 1/3degree, ~5-day, 15m (not NASA daily 0.25degree v2)",
      "source_indices":{"time_chunks":TIME_CHUNKS,"lat":LAT_INDICES,"lon":LON_INDICES},
      "time_start_utc":datetime.fromtimestamp(ts[0],tz=timezone.utc).isoformat(),
      "time_end_utc":datetime.fromtimestamp(ts[-1],tz=timezone.utc).isoformat(),
      "grid_shape":list(u.shape),"source_NetCDF_chunks":checks,
      "finite_uv_cell_time":int(np.sum(np.isfinite(u)&np.isfinite(v))),
      "missing_uv_cell_time":int(np.sum(~(np.isfinite(u)&np.isfinite(v)))),
      "npz_sha256":hashlib.sha256(a.out_npz.read_bytes()).hexdigest(),
      "limits":["No particle trajectories have yet been computed.",
        "NOAA coarse grid can miss coastlines and nearshore conditions",
        "Gap/mask must not be imputed to zero or free advection"]
    }
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(d,indent=2)+"\n",encoding="utf8")
    print("NOAA_ATLANTIC_2005_UV_CUBE",json.dumps({
       "start":d["time_start_utc"],"end":d["time_end_utc"],"grid_shape":d["grid_shape"],
       "finite":d["finite_uv_cell_time"],"missing":d["missing_uv_cell_time"],
       "chunks":[[z["indices"],z["bytes"]] for z in checks]
    }),flush=True)

if __name__=="__main__":main()
