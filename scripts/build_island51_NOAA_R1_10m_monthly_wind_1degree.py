#!/usr/bin/env python3
"""2015-Jan to 2021-Feb publicly available NOAA NCEP/NCAR Reanalysis 1 winds.

Use monthly mean U,V wind at 10m on native 94x192 Gaussian grid.
Produces 1° U,V and mean monthly-vector speed for a scalar, UNDIRECTED
conductance proxy. NOT NASA MERRA-2, not direction-aware dispersal, nor
flight-altitude butterfly winds. Raw 10m wind U,V *vectors* are retained.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

URL_BASE="https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/"
START="2015-01-01";STOP="2021-02-28"
GRID_SHAPE=(170,360)

def hashfile(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for z in iter(lambda:f.read(1024*1024*4),b""):h.update(z)
    return h.hexdigest()

def to_grid(native,lat,lon,*,lat_centers,lon_centers):
    # NOAA Gaussian grid: descending lat, periodic 0..358.125° east lon.
    order=np.argsort(lat)
    lats=np.asarray(lat,dtype=float)[order]
    lons=np.asarray(lon,dtype=float)
    if not np.all(np.diff(lons)>0) or lons.min()<0 or lons.max()>360:
        raise ValueError("NOAA longitude coordinate unexpected")
    values=np.asarray(native,dtype=float)[order]
    if values.shape!=(len(lats),len(lons)):
        raise ValueError("NOAA Gaussian wind raster mismatched lat/lon")
    wrap=np.concatenate([values,values[:,:1]],axis=1)
    grid=RegularGridInterpolator((lats,np.r_[lons,360.]),
                    wrap,bounds_error=True,method="linear")
    yy,xx=np.meshgrid(np.asarray(lat_centers),np.mod(np.asarray(lon_centers),360),
                      indexing="ij")
    return grid(np.column_stack([yy.ravel(),xx.ravel()])).reshape(len(lat_centers),len(lon_centers))

def process(u_path,v_path,out_path,receipt_path):
    u=xr.open_dataset(u_path)
    v=xr.open_dataset(v_path)
    try:
        if "uwnd" not in u or "vwnd" not in v:
            raise ValueError("Official NOAA monthly 10m U/V wind variables absent")
        u0=u.uwnd.sel(time=slice(START,STOP))
        v0=v.vwnd.sel(time=slice(START,STOP))
        if not u0.time.equals(v0.time) or len(u0.time)!=74:
            raise ValueError("NOAA U/V months not synchronized to Jan 2015-Feb 2021")
        if tuple(u0.dims)!=("time","lat","lon") or u0.shape[1:]!=(94,192):
            raise ValueError("Unexpected 94×192 original NOAA surface_gauss wind grid")
        if tuple(v0.dims)!=tuple(u0.dims) or u0.shape!=v0.shape:
            raise ValueError("Monthly NOAA U and V dimensions different")
        if u0.attrs.get("units")!="m/s" or v0.attrs.get("units")!="m/s":
            raise ValueError("Original NOAA winds not in m/s")
        uu=u0.values.astype(float)
        vv=v0.values.astype(float)
        valid=np.isfinite(uu)&np.isfinite(vv)&(abs(uu)<150)&(abs(vv)<150)
        uclean=np.where(valid,uu,np.nan);vclean=np.where(valid,vv,np.nan)
        speed_monthly=np.sqrt(uclean**2+vclean**2)
        n_months=valid.sum(axis=0)
        if (n_months<64).any():
            raise ValueError("Missing 10+ original U/V monthly observations")
        mean_speed=np.nanmean(speed_monthly,axis=0)
        mean_u=np.nanmean(uclean,axis=0);mean_v=np.nanmean(vclean,axis=0)
        lats=np.asarray(u0.lat.values)
        lons=np.asarray(u0.lon.values)
        lat_target=np.arange(84.5,-85,-1,dtype=float)
        lon_target=np.arange(-179.5,180,1,dtype=float)
        dst={
            "speed":to_grid(mean_speed,lats,lons,lat_centers=lat_target,lon_centers=lon_target),
            "U":to_grid(mean_u,lats,lons,lat_centers=lat_target,lon_centers=lon_target),
            "V":to_grid(mean_v,lats,lons,lat_centers=lat_target,lon_centers=lon_target)
        }
        if any(a.shape!=GRID_SHAPE or not np.isfinite(a).all() for a in dst.values()):
            raise ValueError("NOAA 1 degree proxy invalid or missing")
        if not 1.<np.quantile(dst["speed"],.5)<15.:
            raise ValueError("Unphysical 10m monthly wind proxy median")
        out_path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(out_path,mean_monthly_wind_vector_speed_10m=dst["speed"].astype("float32"),
            mean_U_east_10m=dst["U"].astype("float32"),mean_V_north_10m=dst["V"].astype("float32"),
            lat_centers=lat_target,lon_centers=lon_target)
        receipt={
            "schema":"chocho_NOAA_NCEP_R1_10m_UV_monthly_2015Jan_2021Feb_1deg_v1",
            "status":"PUBLIC_NOAA_UV_REANALYSIS_SOURCE_VERIFIED_1DEG_UNDIRECTED_SPEED_PROXY_NOT_ADVECTION",
            "source_DOI_or_url":"https://www.psl.noaa.gov/data/gridded/data.ncep.reanalysis.html",
            "source_urls":{z:URL_BASE+z+".10m.mon.mean.nc" for z in ("uwnd","vwnd")},
            "source_sha256":{z:hashfile(path) for z,path in (("u",u_path),("v",v_path))},
            "source_native_resolution":"~1.875 lon × Gaussian ~1.9 lat",
            "source_original_shape":[74,94,192],
            "observed_months":74,
            "source_wind_height_m":10,
            "source_UV_units":"m/s",
            "derived_grid_degrees":1,
            "derived_grid_shape":list(GRID_SHAPE),
            "wind_speed_proxy":"Arithmetic monthly mean of magnitude of each monthly averaged U,V; NOT raw 6-hourly wind speed mean",
            "median_wind_speed_proxy_m_s":float(np.median(dst["speed"])),
            "p05_p95_wind_speed_proxy_m_s":[float(x) for x in np.quantile(dst["speed"],[.05,.95])],
            "caveats":[
              "NOAA NCEP/NCAR Reanalysis 1 (public), NOT NASA MERRA-2 or ERA5.",
              "A 10m monthly vector-mean speed cannot represent individual butterfly flight altitude, flight season or transoceanic tailwind events.",
              "Undirected scalar conductance loses true wind vector transport direction.",
              "Interpolated NOAA ~2 degree climatology onto 1 degree output does not create new observed fine-scale wind resolution.",
              "2015–2021 modern mean winds cannot explain the actual dates of historic nonnative plant introductions."
            ]}
        receipt_path.parent.mkdir(parents=True,exist_ok=True)
        receipt_path.write_text(json.dumps(receipt,indent=2)+"\n")
        print("NOAA_WIND_2015_2021_COARSE_GRID_SOURCE",json.dumps(receipt,sort_keys=True)[:13000],flush=True)
        return receipt
    finally:
        u.close();v.close()

def main():
    p=argparse.ArgumentParser()
    for x in ("u","v","out","receipt"):p.add_argument("--"+x,type=Path,required=True)
    a=p.parse_args();process(a.u,a.v,a.out,a.receipt)
if __name__=="__main__":main()
