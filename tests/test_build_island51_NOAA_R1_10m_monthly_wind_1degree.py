"""NOAA NCEP R1 10m monthly wind ingestion: geographical and vector invariants."""
import pytest
for p in ("numpy","scipy","xarray","netCDF4"): pytest.importorskip(p)
import numpy as np
import importlib.util
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"build_island51_NOAA_R1_10m_monthly_wind_1degree.py"
sp=importlib.util.spec_from_file_location("wind_source",f)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_wind_wrap_is_periodic_and_preserves_northern_southern_values():
    # native lat descending, 0/90/180/270 deg longitudes
    lats=np.array([75.,25.,-25.,-75.])
    lons=np.array([0.,90.,180.,270.])
    x=np.tile(np.array([3.,4.,5.,6.]),(4,1))+lats[:,None]*.01
    got=m.to_grid(x,lats,lons,lat_centers=[-25.,25.],
                  lon_centers=[-90.,0.,90.,179.5])
    assert got.shape==(2,4)
    assert np.isclose(got[0,0],6-.25)
    assert np.isclose(got[1,1],3+.25)
    assert np.isclose(got[1,2],4+.25)
    assert np.isclose(got[0,3],5-.25,atol=.006)

def test_wind_source_period_and_grid_are_explicit():
    assert m.START=="2015-01-01" and m.STOP=="2021-02-28"
    assert m.GRID_SHAPE==(170,360)
    assert m.URL_BASE.endswith("/surface_gauss/")
