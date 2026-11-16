#!/usr/bin/env python3
"""Small, source-locked NOAA OLD OSCAR NetCDF grid capability gate.

The old NOAA ERDDAP 1/3-degree approximately-five-day data are not the
NASA daily OSCAR v2.0 data. No plant movement or seed transport is tested.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import netCDF4
import numpy as np

BASE = "https://oceanwatch.pifsc.noaa.gov/erddap/griddap/yearly_336c_0b32_9cd3.nc?"
BOUNDS = (930,931,124,126,935,937)
MAX_BYTES = 2_000_000

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    t0,t1,lat0,lat1,lon0,lon1=BOUNDS
    expr=",".join(f"{v}[{t0}:1:{t1}][0][{lat0}:1:{lat1}][{lon0}:1:{lon1}]" for v in ("u","v"))
    url=BASE+quote(expr,safe="[],:")
    req=Request(url,headers={"User-Agent":"chocho-oscar-source-gate-2026","Accept":"application/x-netcdf,*/*"})
    with urlopen(req,timeout=90) as response:
        raw=response.read(MAX_BYTES+1)
        content_type=response.headers.get("Content-Type","")
        if response.status!=200 or len(raw)>MAX_BYTES:
            raise RuntimeError(f"Unexpected NOAA NetCDF response: http={response.status} bytes={len(raw)}")
    if not (raw.startswith(b"CDF") or raw.startswith(b"\x89HDF\r\n\x1a\n")):
        raise RuntimeError(f"NOAA source was not actual NetCDF; content-type={content_type}")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    temp=a.out.with_suffix(".nc")
    temp.write_bytes(raw)
    try:
        with netCDF4.Dataset(temp) as ds:
            required={"u","v","time","latitude","longitude"}
            if not required.issubset(ds.variables):
                raise RuntimeError(f"Missing actual NOAA velocity dimensions: {list(ds.variables)}")
            u=np.asarray(ds.variables["u"][:].filled(np.nan),dtype=float) if hasattr(ds.variables["u"][:],"filled") else np.asarray(ds.variables["u"][:],dtype=float)
            v=np.asarray(ds.variables["v"][:].filled(np.nan),dtype=float) if hasattr(ds.variables["v"][:],"filled") else np.asarray(ds.variables["v"][:],dtype=float)
            la=np.asarray(ds.variables["latitude"][:],dtype=float)
            lo=np.asarray(ds.variables["longitude"][:],dtype=float)
            tm=np.asarray(ds.variables["time"][:])
            timevar=ds.variables["time"]
            dt=netCDF4.num2date(tm,units=timevar.units,
                calendar=getattr(timevar,"calendar","standard"))
            timestamp=[str(t) for t in dt]
            if len(la)!=3 or len(lo)!=3 or len(tm)!=2:raise RuntimeError("NetCDF source subset dimensions unexpectedly changed")
            if u.shape!=v.shape or u.shape != (2,1,3,3):raise RuntimeError(f"Expected 2x1x3x3, got {u.shape}")
            if any(int(x[:4])!=2005 for x in timestamp):raise RuntimeError("This isn't 2005 source velocity")
            finite=int(np.logical_and(np.isfinite(u),np.isfinite(v)).sum())
            if not finite:raise RuntimeError("All probed velocities missing")
            result={
                "status":"NOAA_LEGACY_OSCAR_NETCDF_GRID_READ_SUCCESS_NOT_PARTICLE_TRACKING",
                "dataset_id":"yearly_336c_0b32_9cd3",
                "year":2005,
                "original_index_window":list(BOUNDS),
                "file_bytes":len(raw),"file_sha256":hashlib.sha256(raw).hexdigest(),
                "content_type":content_type,"netcdf_dimensions":{k:len(ds.dimensions[k]) for k in ds.dimensions},
                "time":timestamp,"latitudes":la.tolist(),"longitudes":lo.tolist(),
                "matched_finite_uv_cells":finite,"u_range_m_s":[float(np.nanmin(u)),float(np.nanmax(u))],
                "v_range_m_s":[float(np.nanmin(v)),float(np.nanmax(v))],
                "source_limit":"A 3x3 Atlantic offshore point is a NetCDF IO test only; no geographic seed source, destination, coast hit or biological process inferred",
            }
    finally:
        temp.unlink(missing_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print("OSCAR_NOAA_NETCDF_IO_GATE",json.dumps(result),flush=True)

if __name__=="__main__":main()
