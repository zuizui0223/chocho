"""Unit tests for 416 source-replayed native plant climatic potential connectivity.

No genotype or historical plant range is read by the test suite.
"""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from fit_native_plant_05deg_glacial_potential_graph_v02 import (
    fit_original_classifier,model_probs,grid_validity_and_analog,
    suitable_mask,count_available,GRID_N
)

def test_retrained_quadratic_native_plant_scores_are_fixed_and_deterministic():
    rng=np.random.default_rng(20)
    native=rng.normal(.2,.5,size=(60,4))
    background=rng.normal(-.3,.7,size=(75,4))
    x=rng.normal(size=(24,4))
    m,p=fit_original_classifier(native,background)
    n,q=fit_original_classifier(native,background)
    y=model_probs(m,p,x);z=model_probs(n,q,x)
    assert np.array_equal(y,z)
    assert len(y)==24 and np.all((y>0)&(y<1))
    assert p.n_output_features_==14

def test_nonanalog_grid_returns_unknown_not_occupied_for_far_outside_climate():
    climate=np.zeros((GRID_N,4))
    climate[3]=2.
    climate[5]=np.nan
    allowed=np.zeros(GRID_N,bool);allowed[:7]=True
    train=np.concatenate([np.full((35,4),-.1),np.full((40,4),.1)])
    available,analog=grid_validity_and_analog(climate,allowed,train)
    assert available.sum()==6
    assert analog.sum()==5
    assert not analog[3] and not analog[5] and not available[5]
    q=count_available(allowed,np.isfinite(climate).all(axis=1),analog)
    assert q['icefree_without_climate_n']==1
    assert q['icefree_climate_nonanalog_n']==1
    assert q['icefree_climate_analog_n']==5

def test_suitable_mask_never_mistakes_missing_or_nonanalog_for_absent():
    score=np.full(GRID_N,np.nan)
    score[0:7]=.8
    free=np.zeros(GRID_N,np.uint8);free[:7]=[0,1,4,5,8,9,9]
    analog=np.zeros(GRID_N,bool);analog[:7]=True
    finite=np.zeros(GRID_N,bool);finite[:7]=True;finite[6]=False
    sel=suitable_mask(score,.7,free,analog,finite,land_cut=5)
    assert np.flatnonzero(sel).tolist()==[3,4,5]
    permissive=suitable_mask(score,.7,free,analog,finite,land_cut=1)
    assert np.flatnonzero(permissive).tolist()==[1,2,3,4,5]
    strict=suitable_mask(score,.7,free,analog,finite,land_cut=9)
    assert np.flatnonzero(strict).tolist()==[5]
    analog[3]=False
    assert np.flatnonzero(suitable_mask(score,.7,free,analog,finite,land_cut=5)).tolist()==[4,5]

def test_immature_model_fails_closed_with_no_backfill():
    rng=np.random.default_rng(7)
    p=rng.normal(size=(29,4));b=rng.normal(size=(50,4))
    with pytest.raises(ValueError):
        fit_original_classifier(p,b)

def test_unaligned_2d_grid_and_original_climate_sources_fail_closed():
    with pytest.raises(ValueError):
        grid_validity_and_analog(np.ones((100,4)),np.ones(100),np.ones((40,4)))
    with pytest.raises(ValueError):
        suitable_mask(np.ones(100),.5,np.ones(100),np.ones(100),np.ones(100))
