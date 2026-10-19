import sys
from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import Polygon
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_ttf_host_native_observations import (
    make_matcher,classify_observation,
    NATIVE,OUTSIDE,BOUNDARY,UNMAPPED,NO_NATIVE,boundary_distance_m)

def shapes():
    geoms=[Polygon([(0,0),(2,0),(2,2),(0,2)]),
           Polygon([(2,0),(4,0),(4,2),(2,2)])]
    return geoms,["AAA","BBB"]

def classify(lat,lon,uncertainty,native):
    geoms,codes=shapes()
    return classify_observation(lat=lat,lon=lon,uncertainty_m=uncertainty,
        known_native_codes=set(native),matcher=make_matcher(geoms,codes),
        geoms=geoms,codes=codes)

def test_record_in_native_region_interior():
    assert classify(1.,1.,500.,["AAA"])==(NATIVE,"AAA")

def test_record_outside_documented_native_is_not_called_introduced():
    assert classify(1.,3.,500.,["AAA"])==(OUTSIDE,"BBB")

def test_ambiguous_touching_polygon_boundary_is_not_native():
    status,codes=classify(1.,2.,10.,["AAA"])
    assert status==BOUNDARY and "AAA" in codes and "BBB" in codes

def test_near_regional_border_position_uncertainty_excluded():
    status,_=classify(1.,1.999,300.,["AAA"])
    assert status==BOUNDARY

def test_no_native_wcvp_region_is_unknown_not_absent():
    assert classify(1.,1.,20.,[])[0]==NO_NATIVE

def test_ocean_or_no_wcvp_geometry_is_unmapped():
    assert classify(10.,10.,40.,["AAA"])==(UNMAPPED,"")

def test_quality_invalid_coordinates_fail_closed():
    for v in [np.nan,91.,-91.]:
        with pytest.raises(ValueError):classify(v,1.,0.,["AAA"])
    with pytest.raises(ValueError):classify(1.,1.,1001.,["AAA"])

def test_geodesic_boundary_distance_real_meters():
    geoms,_=shapes()
    from shapely.geometry import Point
    d=boundary_distance_m(Point(1.0,1.0),geoms[0])
    assert 100_000<d<113_000
