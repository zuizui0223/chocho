"""Check exact 0/21 ka physical gating and no missing host imputation."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from pmip4_mask_host_climatic_opportunity_v01 import any_recorded_host_opportunity,delta_stats

def test_only_icefree_climatic_sites_remain_available():
    present=np.array([[.9,.8,.1,.8]])
    lgm=np.array([[.9,.8,.9,.2]])
    q=any_recorded_host_opportunity(present,lgm,np.array([.5]),
                                  np.array([1,1,0,1]),np.array([0,1,1,1]))
    assert q["climatic_suitable_fraction_0ka"]==.75
    assert q["climatic_suitable_fraction_21ka"]==.75
    assert q["icefree_climatic_fraction_0ka"]==.75
    assert q["icefree_climatic_fraction_21ka"]==.5
    assert q["icefree_climatic_delta_21ka_minus_0ka"]==pytest.approx(-.25)
    assert q["climate_suitable_but_physically_blocked_21ka_n"]==1

def test_equal_climate_scores_plus_same_terrain_means_zero_change():
    p=np.array([[.2,.8,.4],[.8,.1,.8]])
    l=np.array([1,1,0])
    q=any_recorded_host_opportunity(p,p,np.array([.5,.5]),l,l)
    assert q["climatic_delta_21ka_minus_0ka"]==0
    assert q["icefree_climatic_delta_21ka_minus_0ka"]==0

def test_plant_union_is_not_sum_of_suitability_probabilities():
    p=np.array([[.8,.1],[.9,.9]])
    q=any_recorded_host_opportunity(p,p,np.array([.6,.6]),np.ones(2),np.ones(2))
    assert q["climatic_suitable_fraction_0ka"]==1.

def test_source_alignment_fail_closed():
    with pytest.raises(ValueError):
        any_recorded_host_opportunity(np.ones((2,10)),np.ones((2,9)),
                                     np.ones(2),np.ones(10),np.ones(10))
    with pytest.raises(ValueError):
        any_recorded_host_opportunity(np.ones((1,10)),np.ones((1,10)),
                                     np.array([np.nan]),np.ones(10),np.ones(10))

def test_per_species_sign_count_including_zero():
    x=[{'x':-.2},{'x':.1},{'x':.0}]
    assert delta_stats(x,'x')=={
      'species':3,'decrease':1,'increase':1,'unchanged':1,'median_paired_delta':0.0}
