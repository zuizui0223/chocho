"""Unit tests for 2022 Lenzner naturalized-flora positive control."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/test_lenzner_alien_flora_history_control.py"


def load():
    spec=importlib.util.spec_from_file_location("alien_flora_positive",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_empire_shared_min_duration_is_not_ship_arrival_or_time_overlap():
    m=load()
    a=np.asarray([[100,10,45,0],[0,0,0,0]],float)
    b=np.asarray([[60,30,20,0],[30,60,0,0]],float)
    expected=m.years_shared(a,b,min_years=25)
    assert list(expected)==[60,0]
    with pytest.raises(ValueError):
        m.analysis([],permutations=5)


def test_pair_hash_selection_is_source_order_locked():
    m=load()
    ids=[f"site-{i:04d}" for i in range(170)]
    x,n=m.pair_sample(ids,maximum=13500,seed=31031)
    y,k=m.pair_sample(ids,maximum=13500,seed=31031)
    assert n==k==170*169//2
    assert np.array_equal(x,y)
    assert len(x)>=10000
    assert all(i<j for i,j in x)
    assert len(set(map(tuple,x.tolist())))==len(x)


def test_tied_rank_baseline():
    m=load()
    x=m.ranks(np.asarray([4,8,8,2]))
    assert list(x)==[2,3.5,3.5,1]
    assert m.correlation(x,x)==pytest.approx(1.0)
    assert m.correlation(x,-x)==pytest.approx(-1.0)
