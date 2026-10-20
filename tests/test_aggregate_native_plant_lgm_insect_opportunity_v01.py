"""Tests for post-fit *potential* resource climate opportunity, not occupied ranges."""
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from aggregate_native_plant_lgm_insect_opportunity_v01 import (
    evaluate_any_host_climate,classify_insect_host_models,idstr,EXPECTED)

def test_complete_known_hosts_are_required_to_score_resource_union():
    ready={"10","20"}
    assert classify_insect_host_models({"10"},ready)=="COMPLETE_KNOWN_HOST_MODELS"
    assert classify_insect_host_models({"10","20"},ready)=="COMPLETE_KNOWN_HOST_MODELS"
    assert classify_insect_host_models({"10","30"},ready)=="PARTIAL_KNOWN_HOST_MODELS_NOT_SCORED"
    assert classify_insect_host_models({"30"},ready)=="NO_KNOWN_HOST_MODELS"
    assert classify_insect_host_models(set(),ready)=="NO_RECORDED_HOSTS"

def test_union_of_recorded_host_climate_suitability_is_not_a_probability():
    present=np.array([[.9,.3,.2,.3],[.1,.8,.2,.1]])
    glacial=np.array([[.1,.7,.2,.2],[.2,.7,.9,.1]])
    thresholds=np.array([.5,.5])
    q=evaluate_any_host_climate(present,glacial,thresholds)
    assert q["sites"]==4
    assert q["present_suitable_sites"]==2
    assert q["LGM_suitable_sites"]==2
    assert q["retained_both_sites"]==1
    assert q["lost_since_LGM_sites"]==1
    assert q["gain_in_LGM_sites"]==1
    assert q["climatic_site_Jaccard"]==pytest.approx(1/3)
    assert q["delta_LGM_minus_present"]==0.

def test_single_host_union_reproduces_underlying_plant_threshold():
    p=np.array([[.9,.1,.1,.8,.3]])
    q=np.array([[.1,.9,.2,.8,.3]])
    x=evaluate_any_host_climate(p,q,np.array([.5]))
    assert x["present_suitable_fraction"]==pytest.approx(.4)
    assert x["LGM_suitable_fraction"]==pytest.approx(.4)
    assert x["climatic_site_Jaccard"]==pytest.approx(1/3)

def test_nonfinite_plant_scores_fail_closed():
    with pytest.raises(ValueError):
        evaluate_any_host_climate(np.array([[np.nan, .2]]),np.array([[.3,.2]]),np.array([.1]))
    with pytest.raises(ValueError):
        evaluate_any_host_climate(np.ones((2,6)),np.ones((1,6)),np.array([.5,.5]))
    with pytest.raises(ValueError):
        evaluate_any_host_climate(np.ones((2,6)),np.ones((2,6)),np.array([.5]))

def test_wcvp_integer_identifiers_do_not_round():
    assert idstr("4e5")=="400000"
    for bad in ["NaN","1.2","0","-42"]:
        with pytest.raises(ValueError): idstr(bad)

def test_artifact_source_is_hash_pinned_and_retains_no_genetic_claims():
    assert len(EXPECTED)==4
    script=(Path(__file__).resolve().parents[1]/"scripts"/"aggregate_native_plant_lgm_insect_opportunity_v01.py").read_text().lower()
    assert "not occupied" in script
    assert "genetic" in script
