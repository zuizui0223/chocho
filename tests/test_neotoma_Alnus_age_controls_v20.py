"""Source-only chronology control reliability audit regression tests."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_neotoma_Alnus_age_controls_v20 import (
    chronology_object,risk_at_sample,number
)

def obj(chronologyid="14552",datasetid=24197,default=True):
    return {"data":[{"chronology":{"chronologyid":int(chronologyid),
      "agetype":"Calibrated radiocarbon years BP","default":default,
      "modelType":"linear interpolation","chronologyName":"Model 1",
      "reliableagespan":{"younger":0,"older":45548},
      "datasets":[{"datasetid":datasetid,"datasettype":"pollen"},
                  {"datasetid":24196,"datasettype":"geochronologic"}],
      "controls":[
        {"chroncontrolid":1,"age":9144,"depth":180,"controltype":"Radiocarbon, calibrated",
         "geochron":{"geochronlabnumber":"Lab-1","geochronage":{"age":8500}}},
        {"chroncontrolid":2,"age":29310,"depth":230,"controltype":"Radiocarbon, calibrated",
         "geochron":{"geochronlabnumber":"Lab-2","geochronage":{"age":25000}}}
      ]}}]}

def test_nested_neotoma_chronology_ID_dataset_link_and_control_depth():
    c=chronology_object(obj(), "14552","24197")
    assert c["isdefault"]
    assert c["source_reliable_older_BP"]==45548
    assert c["source_age_controls_n"]==2
    assert c["source_control_depth_max_cm"]==230
    assert c["control_n_with_original_14C_lab_number"]==2

def test_inside_reliable_span_and_between_two_depth_controls_are_distinct_metadata_claims():
    m=chronology_object(obj(),"14552","24197")
    r=risk_at_sample(205,19227,m)
    assert r["sample_center_inside_database_reliable_age_span"]
    assert r["relative_to_source_returned_control_depth_extent"]=="WITHIN_RETURNED_CONTROL_DEPTH_EXTENT"
    assert r["nearest_bracketing_control_age_gap_BP"]==20166
    assert not r["exact_control_at_sample_depth"]

def test_source_age_outside_trusted_span_even_if_sample_appears_calibrated():
    p=obj("13034",22683);p["data"][0]["chronology"]["reliableagespan"]={"younger":-50,"older":12540}
    p["data"][0]["chronology"]["controls"]=[{"age":1875,"depth":59,"chroncontrolid":9}]
    m=chronology_object(p,"13034","22683")
    r=risk_at_sample(392,20368,m)
    assert r["sample_center_inside_database_reliable_age_span"] is False
    assert r["relative_to_source_returned_control_depth_extent"]=="DEEPER_THAN_DEEPEST_RETURNED_CONTROL"

def test_source_age_deeper_than_last_control_is_not_automatically_false_date():
    m=chronology_object(obj(),"14552","24197")
    r=risk_at_sample(260,21000,m)
    assert r["sample_center_inside_database_reliable_age_span"] is True
    assert r["relative_to_source_returned_control_depth_extent"]=="DEEPER_THAN_DEEPEST_RETURNED_CONTROL"
    assert r["nearest_deeper_control_age_BP"] is None

def test_mismatched_original_chronology_or_site_identity_fails_closed():
    with pytest.raises(ValueError):chronology_object(obj(),"7725","24197")
    with pytest.raises(ValueError):chronology_object(obj(),"14552","14556")
    invalid=obj()
    invalid["data"][0]["chronology"]["reliableagespan"]={"younger":46000,"older":-30}
    with pytest.raises(ValueError):chronology_object(invalid,"14552","24197")

def test_reliability_bounds_must_not_replace_missing_per_sample_uncertainty():
    p=obj();m=chronology_object(p,"14552","24197")
    r=risk_at_sample(205,19227,m)
    assert "source_age_younger_bound_BP" not in r
    assert r["nearest_deeper_control_depth_cm"]==230
    assert r["nearest_shallower_control_depth_cm"]==180

def test_bad_source_numbers_are_not_trusted_age():
    for x in ["",None,True,"nan","inf","uncalibrated"]:
        assert number(x) is None
    with pytest.raises(ValueError):risk_at_sample("unverified depth","20400",chronology_object(obj(),"14552","24197"))
