"""Original source WPI mean-port georeference and GSHHG geometry tests."""
import pytest
for dep in ("numpy","pandas","geopandas","pyproj","shapely"):pytest.importorskip(dep)
import importlib.util,math
from pathlib import Path
import numpy as np
from shapely.geometry import Point
path=Path(__file__).resolve().parents[1]/"scripts"/"recover_island51_WPI_GSHHG_source_geography.py"
spec=importlib.util.spec_from_file_location("wpi51",path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_exact_spherical_mean_at_date_line():
    p=m.spherical_mean([(-179.,0.),(179.,0.)])
    lon,lat=m.latlon(p)
    assert abs(abs(lon)-180.)<1e-8
    assert abs(lat)<1e-8

def test_identical_original_port_distance_formula():
    a=m.spherical_mean([(0,0)])
    b=m.spherical_mean([(90,0)])
    assert math.isclose(m.arc(a,b),m.R*math.pi/2,rel_tol=1e-9)
    c=m.spherical_mean([(30,30),(30,60)])
    lon,lat=m.latlon(c)
    assert abs(lon-30)<1e-8
    assert abs(lat-45)<1e-8

def test_geodesic_circle_is_5000_metres_not_5000_degrees():
    poly=m.approx_disk(10,42)
    assert poly.is_valid and poly.contains(Point(10,42))
    exterior=list(poly.exterior.coords)
    assert len(exterior)==37
    for x,y in exterior[:-1]:
        _,_,d=m.G.inv(10,42,x,y)
        assert abs(d-5000)<1e-3
