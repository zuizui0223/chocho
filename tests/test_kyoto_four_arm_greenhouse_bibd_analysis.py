"""Synthetic original-cohort rows are for software validation, never butterfly observations."""
from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
P=ROOT/"scripts/analyze_kyoto_four_arm_greenhouse_bibd.py"
SPEC=importlib.util.spec_from_file_location("kyoto_bibd_fit",P)
assert SPEC and SPEC.loader
fit=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fit)
from randomize_kyoto_four_arm_greenhouse_bibd import assign, ARMS


def inventory():
    return [
        {"cage_id":f"G{g}:C{c}","greenhouse_id":f"G{g}"}
        for g in range(1,9) for c in range(1,4)
    ]


def fixture(add_greenhouse_offset=False,with_diapause=False):
    allocation=[{k:str(v) for k,v in row.items()}
                for row in assign(inventory(),seed=20261014)]
    source_probs={
        "NATURAL_NONE":.6,"NATURAL_ADD_S":.2,
        "SUPPORTED_NONE":.6,"SUPPORTED_ADD_S":.6,
    }
    info,fates,resources=[],[],[]
    for row in allocation:
        cage,gh,arm=row["cage_id"],row["greenhouse_id"],row["treatment"]
        comp,food=ARMS[arm]
        info.append({
            **{k:row[k] for k in ("cage_id","greenhouse_id","treatment",
                                   "competitor_present","food_support")},
            "focal_a_initial_n":"10","added_s_initial_n":"2" if comp else "0",
            "expected_n_food_visits":"3","food_floor_cm2":"50",
            "donor_batch_id":f"D{gh}",
        })
        n_success=round(10*(source_probs[arm]+(
            0.1 if add_greenhouse_offset and int(gh[1:])<=4 else 0.0
        )))
        for j in range(10):
            value=fit.SUCCESS if j<n_success else "dead_immature"
            if with_diapause and cage=="G1:C1" and j==9:
                value="diapause_alive"
            fates.append({
                "cage_id":cage,"larva_id":f"{cage}:A{j}",
                "species":fit.FOCAL,"cohort_role":"focal","fate":value,
            })
        for j in range(2 if comp else 0):
            fates.append({
                "cage_id":cage,"larva_id":f"{cage}:S{j}",
                "species":fit.ADDED,"cohort_role":"added","fate":"dead_immature",
            })
        for visit in range(3):
            resources.append({
                "cage_id":cage,
                "observation_datetime":f"2026-10-{10+visit:02d}T11:00:00+09:00",
                "focal_a_alive_n":"10",
                "focal_a_late_instar_alive_n":str(visit*4),
                "added_s_alive_n":"2" if comp else "0",
                "accessible_edible_leaf_area_cm2":"85" if food else "40",
                "accessible_leaf_added_area_cm2":"20" if food else "0",
                "leaf_development_stage":"young" if food else "mature",
                "added_leaf_stage":"young" if food else "",
                "leaf_source_batch_id":f"D{gh}",
                "sham_handling_done":"1",
            })
    return allocation,info,fates,resources


def test_greenhouse_adjusted_primary_reconstructs_invented_effect():
    z=fit.analyze(*fixture())
    primary=z["primary_block_adjusted"]
    assert z["design_properties"]["complete_v4_b8_k3_r6_lambda4"] is True
    assert z["primary_estimand"]=="(SUPPORTED_ADD_S-SUPPORTED_NONE)-(NATURAL_ADD_S-NATURAL_NONE)"
    assert primary["n_greenhouses"]==8
    assert primary["n_independent_whole_cages"]==24
    assert primary["n_residual_df"]==13
    assert primary["analytical_contrast_variance_multiplier"]==pytest.approx(.75)
    assert primary["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    assert primary["possible_eventual_emergence_contrast_lower"]==pytest.approx(.4)
    assert primary["possible_eventual_emergence_contrast_upper"]==pytest.approx(.4)
    assert primary["model_based_standard_error_given_homoscedastic_cage_residuals"]==pytest.approx(0,abs=1e-10)
    assert z["leave_one_greenhouse_effect_min"]==pytest.approx(.4)
    assert z["leave_one_greenhouse_effect_max"]==pytest.approx(.4)
    assert z["no_pseudo_replication_of_nested_pots_or_larvae"]
    assert z["cannot_identify_leaf_biomass_only_vs_plant_quality_vs_contact"]


def test_additive_greenhouse_effect_is_removed_from_causal_arm_contrast():
    z=fit.analyze(*fixture(add_greenhouse_offset=True))
    primary=z["primary_block_adjusted"]
    assert primary["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    assert primary["n_residual_df"]==13
    weights=primary["whole_cage_linear_contrast_weights"]
    allocation=fixture()[0]
    for gh in {row["greenhouse_id"] for row in allocation}:
        assert sum(weights[row["cage_id"]] for row in allocation
                   if row["greenhouse_id"]==gh)==pytest.approx(0)


def test_unresolved_living_pupa_does_not_become_known_mortality():
    a,b,c,d=fixture()
    x=next(row for row in c if row["cage_id"]=="G1:C1" and row["cohort_role"]=="focal"
           and row["fate"]=="dead_immature")
    x["fate"]="diapause_alive"
    z=fit.analyze(a,b,c,d)
    w=z["primary_block_adjusted"]["whole_cage_linear_contrast_weights"]["G1:C1"]
    lo=z["primary_block_adjusted"]["possible_eventual_emergence_contrast_lower"]
    hi=z["primary_block_adjusted"]["possible_eventual_emergence_contrast_upper"]
    assert z["primary_block_adjusted"]["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    assert lo<=.4<=hi
    assert hi-lo==pytest.approx(abs(w)/10)
    treatment=next(x["treatment"] for x in a if x["cage_id"]=="G1:C1")
    assert z["per_treatment_raw_arm_fates"][treatment]["living_diapause_A"]==1


def test_ambiguous_or_unidentified_fates_cannot_be_used():
    a,b,c,d=fixture()
    c[0]["fate"]="pupa_no_adult"
    with pytest.raises(ValueError,match="ambiguous pupa_no_adult"):
        fit.analyze(a,b,c,d)
    a,b,c,d=fixture()
    c.pop()
    with pytest.raises(ValueError,match="Every originally randomized"):
        fit.analyze(a,b,c,d)
    a,b,c,d=fixture()
    row=next(x for x in c if x["cohort_role"]=="added")
    row["cohort_role"]="focal"
    with pytest.raises(ValueError,match="Wrong original focal|Every originally randomized"):
        fit.analyze(a,b,c,d)


def test_donor_and_observation_conditions_cannot_be_conditioned_on_success():
    a,b,c,d=fixture()
    c0=next(x for x in a if x["treatment"]=="SUPPORTED_ADD_S")["cage_id"]
    for visit in d:
        if visit["cage_id"]==c0:
            visit["accessible_edible_leaf_area_cm2"]="0"
    z=fit.analyze(a,b,c,d)
    assert z["primary_block_adjusted"]["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    assert z["per_treatment_raw_arm_fates"]["SUPPORTED_ADD_S"]["measured_food_floor_breach_visits"]>=3
    a,b,c,d=fixture()
    d[0]["accessible_edible_leaf_area_cm2"]=""
    z=fit.analyze(a,b,c,d)
    assert sum(r["food_visits_missing_leaf_area"] for r in z["per_treatment_raw_arm_fates"].values())==1
    a,b,c,d=fixture()
    d.pop()
    with pytest.raises(ValueError,match="missing scheduled food observation"):
        fit.analyze(a,b,c,d)


def test_randomized_four_arm_incidence_and_greenhouse_ids_cannot_be_rewritten():
    a,b,c,d=fixture()
    a[0]["treatment"]=a[1]["treatment"]
    with pytest.raises(ValueError,match="assignment|incidence|Duplicate treatment|not a connected"):
        fit.analyze(a,b,c,d)
    a,b,c,d=fixture()
    b[0]["greenhouse_id"]="WRONG"
    with pytest.raises(ValueError,match="Postrandomization"):
        fit.analyze(a,b,c,d)
    a,b,c,d=fixture()
    b[0]["focal_a_initial_n"]="9"
    with pytest.raises(ValueError,match="Focal A initial cohort size"):
        fit.analyze(a,b,c,d)


def test_new_input_source_columns_and_timestamp_are_validated():
    a,b,c,d=fixture()
    d[0]["observation_datetime"]="2026-10-11"
    with pytest.raises(ValueError,match="Timezone"):
        fit.analyze(a,b,c,d)
    a,b,c,d=fixture()
    d[0]["accessible_leaf_added_area_cm2"]="12"
    with pytest.raises(ValueError,match="Accessible supplementation"):
        fit.analyze(a,b,c,d)


def test_four_arm_original_cage_status_not_six_arm_retroactive_conversion():
    a,b,c,d=fixture()
    assert len(a)==24
    assert all(x["design_version"]=="kyoto_four_arm_greenhouse_bibd_v01" for x in a)
    assert set(x["treatment"] for x in a)==set(fit.ARMS)
    assert not any("ADD_A" in x["treatment"] for x in a)
    assert sum(int(x["focal_a_initial_n"]) for x in b)==240
