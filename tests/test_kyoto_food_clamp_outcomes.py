"""Synthetic unit tests only. No butterfly outcomes from experimental observations."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import pytest

PATH=Path(__file__).resolve().parents[1]/"scripts/analyze_kyoto_food_clamp_outcomes.py"
spec=importlib.util.spec_from_file_location("anal",PATH)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def fixture(nblocks=3,with_unknown=True):
    alloc=[]; info=[]; fates=[]; visits=[]
    for b in range(nblocks):
        for arm,(competitor,clamp) in m.ARMS.items():
            plant=f"P{b}-{arm}"
            base={"plant_id":plant,"block_id":f"B{b}","treatment":arm,
                  "competitor_present":str(competitor),"food_clamp":str(clamp)}
            alloc.append({**base,"random_seed":"20261010","design_version":"kyoto_hidden_resource_bottleneck_v01"})
            info.append({**base,"native_initial_n":"2","competitor_initial_n":str(competitor*2),
                         "food_floor_cm2":"20","expected_n_resource_visits":"2"})
            for larva_i in range(2):
                # Each of four treatment arms produces a distinct outcome.
                fate={
                    "NATURAL_NO_COMPETITOR":"flight_capable_adult" if larva_i==0 else "dead_immature",
                    "NATURAL_COMPETITOR":"dead_immature",
                    "CLAMP_NO_COMPETITOR":"flight_capable_adult" if larva_i==0 else "dead_immature",
                    "CLAMP_COMPETITOR":"flight_capable_adult",
                }[arm]
                if with_unknown and b==0 and arm=="CLAMP_NO_COMPETITOR" and larva_i==1:
                    fate="lost_to_followup"
                fates.append({"plant_id":plant,"larva_id":f"{plant}:L{larva_i}",
                              "species":"Atrophaneura alcinous","fate":fate})
            for i in range(2):
                visits.append({"plant_id":plant,
                               "observation_datetime":f"2026-10-{10+i:02d}T11:00:00+09:00",
                               "native_larvae_alive_n":"2",
                               "native_late_instar_alive_n":str(i*2),
                               "accessible_leaf_area_cm2":"25" if clamp or not competitor else "12",
                               "fresh_leaf_area_added_cm2":"5" if clamp else "0",
                               "present_leaf_age_class":"mature",
                               "added_leaf_age_class":"young" if clamp else ""})
    return alloc,info,fates,visits

def test_itt_interaction_uses_all_initial_larvae_and_source_blocks():
    x=m.analyze(*fixture())
    assert x["n_independent_cages"]==12
    assert x["randomization_unit"]=="whole cage/enclosure, including all its potted plants"
    assert x["n_blocks"]==3
    assert x["ITT_missing_as_fail_interaction"]==1.0
    assert x["true_worst_case_attrition_interaction_bounds"]==[1.0-1/6,1.0]
    assert x["per_arm"]["CLAMP_NO_COMPETITOR"]["unknown_fate"]==1
    assert x["all_native_cohort_fates_known"] is False
    assert x["independent_block_bootstrap_95ci_for_missing_as_fail"] is not None
    assert x["nonresource_mechanism_causally_identified"] is False

def test_missing_original_larva_cannot_be_dropped():
    a,b,c,d=fixture()
    c.pop()
    with pytest.raises(ValueError,match="originally randomized larvae"):
        m.analyze(a,b,c,d)

def test_treatment_changes_after_randomization_rejected():
    a,b,c,d=fixture()
    b[0]["treatment"]="NATURAL_COMPETITOR"
    with pytest.raises(ValueError,match="treatment/block mismatch"):
        m.analyze(a,b,c,d)

def test_missing_resource_access_is_not_zero_and_scheduled_visits_are_required():
    a,b,c,d=fixture()
    d[0]["accessible_leaf_area_cm2"]=""
    x=m.analyze(a,b,c,d)
    assert x["per_arm"]["NATURAL_NO_COMPETITOR"]["missing_food_measurements"]==1
    assert x["per_arm"]["NATURAL_NO_COMPETITOR"]["observed_late_instar_food_visits"]==3
    assert x["per_arm"]["NATURAL_COMPETITOR"]["late_instar_visits_below_food_floor"]==3
    assert x["per_arm"]["CLAMP_COMPETITOR"]["supplemented_leaf_age_classes"]=={"young":6}
    d.pop()
    with pytest.raises(ValueError,match="scheduled resource visits"):
        m.analyze(a,b,c,d)

def test_insufficient_independent_blocks_yield_no_bootstrap_ci():
    x=m.analyze(*fixture(nblocks=2,with_unknown=False))
    assert x["independent_block_bootstrap_95ci_for_missing_as_fail"] is None
    assert x["all_native_cohort_fates_known"] is True

def test_biomass_afterward_cannot_define_a_treatment():
    a,b,c,d=fixture()
    b[0]["food_clamp"]="1"
    with pytest.raises(ValueError,match="treatment/block mismatch"):
        m.analyze(a,b,c,d)

def test_late_instar_count_cannot_exceed_total_alive_or_initial():
    a,b,c,d=fixture()
    d[1]["native_late_instar_alive_n"]="3"
    with pytest.raises(ValueError,match="late-instar count"):
        m.analyze(a,b,c,d)

def test_positive_supplement_must_have_age_class():
    a,b,c,d=fixture()
    target=next(row for row in d if row["fresh_leaf_area_added_cm2"]=="5")
    target["added_leaf_age_class"]=""
    with pytest.raises(ValueError,match="added_leaf_age_class"):
        m.analyze(a,b,c,d)

def test_late_instar_food_observation_missing_is_not_interpreted_as_zero():
    a,b,c,d=fixture()
    d[1]["accessible_leaf_area_cm2"]=""
    x=m.analyze(a,b,c,d)
    natural=x["per_arm"]["NATURAL_NO_COMPETITOR"]
    assert natural["late_instar_food_measurement_missing"]==1
    assert natural["observed_late_instar_food_visits"]==2
    assert x["n_independent_cages"]==12

def test_no_late_instar_survivors_does_not_prove_adequate_food():
    a,b,c,d=fixture()
    for row in d:
        row["native_late_instar_alive_n"]="0"
    x=m.analyze(a,b,c,d)
    assert x["per_arm"]["CLAMP_COMPETITOR"]["observed_late_instar_food_visits"]==0
    assert x["per_arm"]["CLAMP_COMPETITOR"]["late_instar_fraction_below_food_floor"] is None
    assert x["ITT_missing_as_fail_interaction"]==1.0
