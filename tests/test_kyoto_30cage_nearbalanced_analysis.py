"""Synthetic whole-cage endpoint fixtures only: NO original Kyoto adult data."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
SPEC=importlib.util.spec_from_file_location(
    "near30fit", ROOT/"scripts/analyze_kyoto_30cage_nearbalanced.py"
)
assert SPEC and SPEC.loader
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)
from randomize_kyoto_30cage_nearbalanced import assign, design_counts, ARMS


def inventory():
    return [{"cage_id":f"G{g}:C{c}","greenhouse_id":f"G{g}"}
            for g in range(1,11) for c in range(1,4)]


def fixture(additive_gh=False, source_seed=20261014, unknown=None):
    original=[{key:str(value) for key,value in r.items()}
              for r in assign(inventory(),seed=source_seed)]
    outcomes={
        "NATURAL_NONE":6,"NATURAL_ADD_S":2,
        "SUPPORTED_NONE":6,"SUPPORTED_ADD_S":6,
    }
    info,fates,food=[],[],[]
    for row in original:
        cage,gh,arm=row["cage_id"],row["greenhouse_id"],row["treatment"]
        comp,support=ARMS[arm]
        info.append({
            **{field:row[field] for field in (
                "cage_id","greenhouse_id","treatment","competitor_present","food_support"
            )},
            "focal_a_initial_n":"10",
            "added_s_initial_n":"2" if comp else "0",
            "expected_n_food_visits":"3","food_floor_cm2":"40",
            "donor_batch_id":f"D{gh}",
        })
        success_n=outcomes[arm]+(1 if additive_gh and int(gh[1:])<=5 else 0)
        for k in range(10):
            fate="flight_capable_adult" if k<success_n else "dead_immature"
            if unknown and (cage,k)==unknown:
                fate="diapause_alive"
            fates.append({
                "cage_id":cage,"larva_id":f"{cage}:FOCAL:{k}",
                "cohort_role":"focal","species":"Atrophaneura alcinous",
                "fate":fate,
            })
        for k in range(2 if comp else 0):
            fates.append({
                "cage_id":cage,"larva_id":f"{cage}:EXOTIC:{k}",
                "cohort_role":"added","species":"Sericinus montela",
                "fate":"dead_immature",
            })
        for i in range(3):
            food.append({
                "cage_id":cage,
                "observation_datetime":f"2026-10-{10+i:02d}T11:00:00+09:00",
                "focal_a_alive_n":"10",
                "focal_a_late_instar_alive_n":str(i*4),
                "added_s_alive_n":"2" if comp else "0",
                "accessible_edible_leaf_area_cm2":"60" if support else "20",
                "accessible_leaf_added_area_cm2":"20" if support else "0",
                "leaf_development_stage":"young" if support else "mature",
                "added_leaf_stage":"young" if support else "",
                "leaf_source_batch_id":f"D{gh}",
                "sham_handling_done":"1",
            })
    return original,info,fates,food


def test_30_cage_primary_it_treatment_design_and_information():
    result=m.analyze(*fixture())
    fit=result["primary_greenhouse_fixed_effect_adjusted"]
    assert result["design"]["design_valid"]
    assert result["design"]["n_whole_cages"]==30
    assert result["design"]["n_greenhouses"]==10
    assert fit["n_independent_whole_cages"]==30
    assert fit["n_residual_df"]==17
    assert fit["analytical_contrast_variance_multiplier"]==pytest.approx(.6)
    assert fit["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    assert fit["possible_eventual_emergence_contrast_lower"]==pytest.approx(.4)
    assert fit["possible_eventual_emergence_contrast_upper"]==pytest.approx(.4)
    assert result["raw_intention_to_treat_arm_mean_difference_in_differences"]==pytest.approx(.4)
    assert fit["model_based_standard_error_given_homoscedastic_cage_residuals"]==pytest.approx(0,abs=1e-9)
    assert fit["nominal_t17_95pct_interval_under_cage_error_assumptions"]==pytest.approx([.4,.4])
    assert result["cage_heteroscedasticity_HC3_sensitivity_only"]["hc3_cage_level_primary_SE"]==pytest.approx(0,abs=1e-9)
    assert result["ten_greenhouse_delete_one_jackknife_SE_descriptive_only"]==pytest.approx(0,abs=1e-9)
    assert result["jackknife_SE_is_NOT_randomization_exact_or_power_evidence"]
    assert result["no_postrandomization_food_or_survivor_selection"]
    assert result["status"]=="FUTURE_AUTHENTIC_FIELD_DATA_REQUIRED_NO_BUTTERFLY_RESULTS"


def test_greenhouse_additive_effect_cancels_fixed_effect_contrast():
    result=m.analyze(*fixture(additive_gh=True))
    fit=result["primary_greenhouse_fixed_effect_adjusted"]
    assert fit["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    assert all(
        abs(sum(weight for cage,weight in fit["whole_cage_linear_contrast_weights"].items()
                if cage.startswith(f"G{gh}:")))<1e-9
        for gh in range(1,11)
    )
    assert len(result["leave_one_greenhouse_out_description_only"])==10
    assert all(row["adjusted_primary_for_descriptive_sensitivity"]==pytest.approx(.4)
               for row in result["leave_one_greenhouse_out_description_only"])


def test_original_focal_and_added_S_cohorts_are_not_substitutable():
    a,b,c,d=fixture()
    c=[x for x in c if x["larva_id"]!="G1:C1:FOCAL:0"]
    with pytest.raises(ValueError,match="Every originally randomized larva"):
        m.analyze(a,b,c,d)
    a,b,c,d=fixture()
    c=[x for x in c if x["larva_id"]!="G1:C1:EXOTIC:0"]
    if any(x["cage_id"]=="G1:C1" and x["cohort_role"]=="added" for x in fixture()[2]):
        with pytest.raises(ValueError,match="Every originally randomized larva"):
            m.analyze(a,b,c,d)
    a,b,c,d=fixture()
    c[0]["fate"]="pupa_no_adult"
    with pytest.raises(ValueError,match="ambiguous"):
        m.analyze(a,b,c,d)


def test_living_diapause_kept_as_possible_eventual_emergence_not_failure():
    a,b,c,d=fixture()
    target=next(z for z in c if z["cage_id"]=="G1:C1" and z["fate"]=="dead_immature")
    target["fate"]="diapause_alive"
    result=m.analyze(a,b,c,d)
    fit=result["primary_greenhouse_fixed_effect_adjusted"]
    low,upper=(fit["possible_eventual_emergence_contrast_lower"],
               fit["possible_eventual_emergence_contrast_upper"])
    point=fit["greenhouse_adjusted_observed_adult_contrast"]
    assert low-1e-10<=point<=upper+1e-10
    w=fit["whole_cage_linear_contrast_weights"]["G1:C1"]
    assert upper-low==pytest.approx(abs(w)/10)
    arm=next(r["treatment"] for r in a if r["cage_id"]=="G1:C1")
    assert result["per_arm_original_cohort_outcomes"][arm]["living_diapausing_focal_A_pupae"]==1


def test_unchanged_balanced_design_cannot_hide_reassigned_cage_labels():
    a,b,c,d=fixture()
    original=[x for x in a if x["greenhouse_id"]=="G1"]
    for col in ("treatment","competitor_present","food_support"):
        original[0][col],original[1][col]=original[1][col],original[0][col]
    assert design_counts(a)["design_valid"]
    with pytest.raises(ValueError,match="seed replay mismatch"):
        m.analyze(a,b,c,d)


def test_wrong_seed_and_wrong_original_provenance_fail_closed():
    a,b,c,d=fixture()
    for row in a:
        row["random_seed"]="123456"
    with pytest.raises(ValueError,match="seed replay mismatch"):
        m.analyze(a,b,c,d)
    a,b,c,d=fixture()
    b[0]["greenhouse_id"]="G99"
    with pytest.raises(ValueError,match="greenhouse"):
        m.analyze(a,b,c,d)
    a,b,c,d=fixture()
    b[0]["donor_batch_id"]="different"
    with pytest.raises(ValueError,match="donor_batch_id"):
        m.analyze(a,b,c,d)


def test_missing_or_zero_foliage_measurements_not_used_to_select_cages():
    a,b,c,d=fixture()
    d[0]["accessible_edible_leaf_area_cm2"]=""
    r=m.analyze(a,b,c,d)
    assert sum(z["missing_food_access_measurements"] for z in r["per_arm_original_cohort_outcomes"].values())==1
    assert r["primary_greenhouse_fixed_effect_adjusted"]["greenhouse_adjusted_observed_adult_contrast"]==pytest.approx(.4)
    a,b,c,d=fixture()
    d[0]["accessible_edible_leaf_area_cm2"]="0"
    r=m.analyze(a,b,c,d)
    assert sum(z["missing_food_access_measurements"] for z in r["per_arm_original_cohort_outcomes"].values())==0
    a,b,c,d=fixture()
    d.pop()
    with pytest.raises(ValueError,match="missing scheduled food observation"):
        m.analyze(a,b,c,d)


def test_observed_count_still_30_not_thirty_larvae():
    a,b,c,d=fixture()
    r=m.analyze(a,b,c,d)
    assert r["design"]["n_whole_cages"]==30
    assert sum(x["original_assigned_focal_A_larvae"] for x in r["per_arm_original_cohort_outcomes"].values())==300
    assert sum(x["original_randomized_whole_cages"] for x in r["per_arm_original_cohort_outcomes"].values())==30
    assert sum(x["confirmed_flight_capable_focal_A_adults"] for x in r["per_arm_original_cohort_outcomes"].values())>0
    assert not r["greenhouse_specific_treatment_effects_not_identifiable_with_one_cage_per_cell"] is False


def test_heterogeneous_synthetic_cage_outcome_exposes_both_sensitivity_diagnostics():
    a,b,c,d=fixture()
    target=next(row for row in a if row["treatment"]=="NATURAL_ADD_S")
    source_cage=target["cage_id"]
    extra=next(row for row in c if row["cage_id"]==source_cage
               and row["cohort_role"]=="focal" and row["fate"]=="dead_immature")
    extra["fate"]="flight_capable_adult"
    r=m.analyze(a,b,c,d)
    hc3=r["cage_heteroscedasticity_HC3_sensitivity_only"]
    assert hc3["hc3_cage_level_primary_SE"]>0
    assert 0 < hc3["max_hat_leverage"] < 1
    assert hc3["NOT_greenhouse_cluster_robust_or_validated_biological_inference"]
    assert r["ten_greenhouse_delete_one_jackknife_SE_descriptive_only"]>0
    assert r["jackknife_SE_is_NOT_randomization_exact_or_power_evidence"]
