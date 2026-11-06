"""No invented coordinates: recover original WPI geocentres only after source lock."""
import pytest
import math,importlib.util
from pathlib import Path
root=Path(__file__).resolve().parents[1]/"scripts"
sp=importlib.util.spec_from_file_location("islandcoords",root/"reconstruct_island51_verified_port_coordinates.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_spherical_average_handles_antimeridian():
    a=m.spherical_mean([(-179,0),(179,0)])
    assert abs(abs(a["longitude"])-180)<1e-7
    assert abs(a["latitude"])<1e-9
    assert math.isclose(math.sqrt(sum(x*x for x in a["unit_sphere_xyz"])),1,abs_tol=1e-9)

def test_geodesic_reproduces_great_circle():
    a=m.spherical_mean([(0,0)])["unit_sphere_xyz"]
    b=m.spherical_mean([(90,0)])["unit_sphere_xyz"]
    assert math.isclose(m.arc_distance(a,b),m.R_KM*math.pi/2,rel_tol=1e-10)
    c=m.spherical_mean([(0,20),(0,0)])
    assert abs(c["latitude"]-10)<1e-9

def test_missing_real_coords_or_undefined_source_cannot_be_filled():
    with pytest.raises(ValueError,match="No port coordinates"):
        m.spherical_mean([])
    with pytest.raises(ValueError,match="Impossible"):
        m.spherical_mean([(500,100)])
    with pytest.raises(ValueError,match="Spherical average not uniquely defined"):
        m.spherical_mean([(0,0),(180,0)])
