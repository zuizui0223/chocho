"""Pre-outcome tests of historical colonial exposure and pairwise dependence."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/"scripts/analyze_colonial_legacy_resource_beta.py"


def load():
    spec=importlib.util.spec_from_file_location("colonial_resource",CODE)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_four_empire_pair_overlap_is_not_a_route_or_calendar_overlap():
    m=load()
    a=np.array([[100,15,0,0],[50,40,20,10],[0,0,0,0]],float)
    b=np.array([[80,5,0,0],[30,15,20,10],[0,0,0,0]],float)
    result=m.shared_years(a,b,years_min=25)
    assert result.tolist()==[80,30,0]
    assert m.shared_years(a,b,years_min=50).tolist()==[80,0,0]


def test_haversine_and_fixed_rank_residualization():
    m=load()
    assert m.haversine(np.array([0.]),np.array([0.]),
                       np.array([0.]),np.array([0.]))[0]==pytest.approx(0)
    assert m.haversine(np.array([0.]),np.array([0.]),
                       np.array([0.]),np.array([1.]))[0]==pytest.approx(111.195,abs=.1)
    x=np.array([1,4,2,3,5,6],float)
    basis=m.orthogonal_controls([np.array([0,1,0,1,0,1],float)])
    residual=m.residual(x,basis)
    assert abs(np.sum(residual))<1e-8


def test_source_lock_avoids_claiming_actual_ship_transport():
    m=load()
    import json
    lock=json.loads(m.LOCK.read_text(encoding="utf8"))
    assert lock["eligibility"]["min_glonaf_area_coverage"]==.8
    assert lock["exposure"]["minimum_area_weighted_duration_years_per_region"]==25
    assert lock["randomization"]["draws"]==499
    assert "not causal" in lock["node_dependencies"].lower() or "no ordinary" in lock["node_dependencies"].lower()
