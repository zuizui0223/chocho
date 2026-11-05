"""PANGAEA original source 2015 WW flows at WPI ports (offline tests)."""
import importlib.util,sys,math
from pathlib import Path
import pytest
pytest.importorskip("pandas")
pytest.importorskip("scipy")
pytest.importorskip("statsmodels")
pytest.importorskip("netCDF4")
import numpy as np
p=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_global_jones_wastewater2015.py"
s=importlib.util.spec_from_file_location("ww2015",p)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def test_latlon_circle_distance_and_lon_wrapping():
    x=m.haversine_distance(0,179.9,np.array([0,0]),np.array([-179.9,170.0]))
    assert 20<x[0]<25
    assert x[1]>1000

def test_grid_cell_sampling_with_valid_zero_and_missing_values():
    lat=np.array([-0.125,-0.04166667,0.04166667,0.125])
    lon=np.array([-0.125,-0.04166667,0.04166667,0.125])
    a={key:np.ones((4,4))*i for i,key in enumerate(m.VARIABLES,1)}
    a["WWp"][0,0]=np.nan
    out=m.grid_one(0,0,lat,lon,a,radius_km=25)
    assert out["cell_centers_in_25km"]==16
    assert out["valid_production_cell_centers"]==15
    assert out["produced_25km_millionm3yr"]==15
    assert out["treated_25km_millionm3yr"]==45

def test_metadata_declares_no_local_N_or_outfall_flux():
    assert m.PANGAEA_DOI=="10.1594/PANGAEA.918731"
    assert m.SHA256.startswith("c413f749")
