"""Source-free WPI harbor-to-TDWG colonial geography tests."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/join_island_ports_colonial_history.py"


def mod():
    pytest.importorskip("shapely")
    pytest.importorskip("pyproj")
    spec=importlib.util.spec_from_file_location("island_port_colonial",SCRIPT)
    a=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a)
    return a


def test_coastal_point_direct_nearest_and_unknown_not_zero():
    a=mod()
    from shapely.geometry import box
    from shapely.strtree import STRtree
    polygons=[box(-1,-1,1,1),box(2,-1,4,1)]
    tr=STRtree(polygons)
    names=["AAA","BBB"]
    code,dist,n,status=a.closest_tdgw_polygon(0,0,tr,polygons,names)
    assert (code,dist,n,status)==("AAA",0.,1,"DIRECT_COVERS")
    code,dist,n,status=a.closest_tdgw_polygon(1.05,0,tr,polygons,names)
    assert code=="AAA"
    assert 0<dist<10
    assert status=="COASTAL_NEAREST"
    code,dist,n,status=a.closest_tdgw_polygon(12,12,tr,polygons,names)
    assert code is None
    assert n==0
    assert status=="PORT_OFFSHORE_BEYOND_25KM"
    with pytest.raises(ValueError,match="invalid"):
        a.closest_tdgw_polygon(190,0,tr,polygons,names)


def test_no_unlocked_distance_cutoffs_or_missing_colonial_means():
    a=mod()
    assert a.YEARS==("years_GB","years_ESP","years_PRT","years_NL")
    assert a.LOCK.exists()
    with pytest.raises(ValueError,match="predeclared"):
        a.build(max_port_distance_km=4.0)
