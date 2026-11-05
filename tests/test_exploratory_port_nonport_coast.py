"""Sanity-check deterministic geometry helpers without network/GSHHG."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

# The submission-facing offline paper CI does not need GSHHG geospatial extras.
# This exploratory test runs fully in the dedicated coastal workflow, where
# pyshp / pyproj / shapely are installed; otherwise skip without breaking paper CI.
pytest.importorskip("shapefile")
pytest.importorskip("pyproj")
pytest.importorskip("shapely")
from shapely.geometry import Point, Polygon

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_nonport_coast.py"
spec=importlib.util.spec_from_file_location("coast_geometry",SCRIPT)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_geodesics_and_nearest_port_distance():
    assert m.geodesic_km((0,0),(0,0)) < 1e-8
    assert 110 < m.geodesic_km((0,0),(0,1)) < 112
    import numpy as np
    assert np.isclose(m.haversine_array(0,0,np.array([0.]),np.array([1.]))[0],111.195,atol=.01)


def test_shore_candidates_really_inside_original_polygon():
    square=Polygon([(-1,-1),(1,-1),(1,1),(-1,1),(-1,-1)])
    points=m.shore_candidates(square,n=24)
    assert len(points)>12
    assert all(square.contains(Point(x,y)) for x,y in points)
    near_coast=[m.geodesic_km((x,y),(-1,y)) for x,y in points if abs(x+1)<.1]
    assert len(near_coast)>0
    assert min(near_coast)<1.


def test_inland_offset_not_assumed_at_invalid_shore():
    tiny=Polygon([(0,0),(.000001,0),(.000001,.000001),(0,.000001)])
    result=m.inward_point(tiny,Point(0,0),target_m=250)
    assert result is None


def test_fixed_source_cohort_shape_and_ids():
    rows,countries=m.source_cohort()
    assert len(rows)==79
    assert len({r["landmass_id"] for r in rows})==36
    assert len({r["port_id"] for r in rows})==79
    assert all(int(r["port_id"]) in countries for r in rows)
