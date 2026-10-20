"""Model-preserving sensitivity: record uncertainty by including modern-nonanalog climates."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from fit_native_plant_05deg_nonanalog_sensitivity_v03 import (
    suitable_mask,GRID_N,GRID_SHAPE,fit_original_classifier,model_probs
)

def test_nonanalog_inclusion_cannot_remove_climate_suitable_cells():
    probs=np.full(GRID_N,.8,dtype=float)
    land=np.ones(GRID_N,dtype=np.uint8)*5
    finite=np.ones(GRID_N,bool)
    analog=np.ones(GRID_N,bool)
    analog[10:20]=False
    without=suitable_mask(probs,.7,land,analog,finite,land_cut=5,include_nonanalog=False)
    with_extrap=suitable_mask(probs,.7,land,analog,finite,land_cut=5,include_nonanalog=True)
    assert np.all(without<=with_extrap)
    assert without.sum()==GRID_N-10 and with_extrap.sum()==GRID_N

def test_nonanalog_extrapolation_is_not_missing_climate_imputation():
    probs=np.full(GRID_N,.8,dtype=float)
    probs[7]=np.nan
    land=np.ones(GRID_N,dtype=np.uint8)*5
    finite=np.ones(GRID_N,bool);finite[7]=False
    analog=np.ones(GRID_N,bool);analog[6]=False
    x=suitable_mask(probs,.7,land,analog,finite,include_nonanalog=True)
    assert x.sum()==GRID_N-1
    assert x[6] and not x[7]

def test_original_native_sdm_fit_remains_deterministic():
    rng=np.random.default_rng(19)
    p=rng.normal(.3,.5,size=(56,4));b=rng.normal(-.5,.6,size=(75,4))
    x=rng.normal(size=(12,4))
    clf,poly=fit_original_classifier(p,b)
    ref,poly2=fit_original_classifier(p,b)
    assert np.array_equal(model_probs(clf,poly,x),model_probs(ref,poly2,x))

def test_no_biological_claim_created_by_second_model():
    s=(Path(__file__).resolve().parents[1]/"scripts"/"fit_native_plant_05deg_nonanalog_sensitivity_v03.py").read_text()
    assert 'include_nonanalog_10pct_majority_4neighbor' in s
    assert 'genetics' in s.lower()
    assert 'not real 21 ka plant occupancy' in s.lower()
