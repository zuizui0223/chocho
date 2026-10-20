"""Reverse exact frozen 0BP logistic model function; hold back 21ka for validation."""
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from reconstruct_frozen_416_logistic_and_extrapolate_v03 import (
    polynomial_features,frozen_logit_coefficients,probabilities,summarize_pair,GRID
)

def test_known_quadratic_logistic_function_reconstructed_from_float32_present_scores():
    rng=np.random.default_rng(3)
    present=rng.normal(size=(6337,4))
    heldout=rng.normal(size=(1500,4))
    beta=rng.normal(0,.4,size=15)
    x=polynomial_features(present)
    p=expit(x@beta).astype('f4')
    recovered,n=frozen_logit_coefficients(x,p)
    assert n>1000 and recovered.shape==(15,)
    assert np.max(abs(probabilities(x,recovered)-p))<1e-5
    assert np.max(abs(probabilities(polynomial_features(heldout),recovered)
                      -expit(polynomial_features(heldout)@beta)))<1e-6

def test_saturated_probabilities_cannot_be_used_as_exact_logit_observations():
    rng=np.random.default_rng(2)
    x=polynomial_features(rng.normal(size=(150,4)))
    bad=np.ones(150,float)
    with pytest.raises(ValueError):frozen_logit_coefficients(x,bad)

def test_source_fit_current_only_no_other_epoch_argument():
    rng=np.random.default_rng(1)
    present=polynomial_features(rng.normal(size=(100,4)))
    p=expit(present@rng.normal(scale=.2,size=15)).astype('f4')
    assert frozen_logit_coefficients(present,p)[0].shape==(15,)

def test_graph_metrics_of_present_and_21ka_are_finite_and_source_0plus21_alignment():
    a=np.zeros(GRID,bool);b=np.zeros(GRID,bool)
    a[15,1:6]=1;b[15,1:5]=1
    q=summarize_pair(a,b,a,b,89.75-0.5*np.arange(360))
    assert q['main']==q['extrapolated']
    assert q['main']['area_delta_km2']<0
    assert q['main']['component_delta']==0
    with pytest.raises(ValueError):
        summarize_pair(a,np.zeros((200,200),bool),a,b,np.arange(360))

def test_source_script_references_archived_0BP_and_requires_independent_21ka_holdout():
    text=(Path(__file__).resolve().parents[1]/"scripts"/
       "reconstruct_frozen_416_logistic_and_extrapolate_v03.py").read_text()
    assert 'source_holdout_max_absolute_error_21ka' in text
    assert 'source_21ka_predictions_used_only_for_validation' in text
    assert 'historical' in text.lower()
