"""Artificial row fixtures verify arithmetic and missing-fate protections only.

Never interpret the invented adults and foliage as measured butterflies.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
p = ROOT / "scripts/analyze_kyoto_six_arm_identity_food.py"
spec = importlib.util.spec_from_file_location("kyoto6_analysis", p)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixture(nblocks=3, with_unknown=False):
    allocation, assignments, fates, resources = [], [], [], []
    arm_k = {
        "NATURAL_NONE": 3,
        "NATURAL_ADD_A": 1,
        "NATURAL_ADD_S": 0,
        "SUPPORTED_NONE": 3,
        "SUPPORTED_ADD_A": 2,
        "SUPPORTED_ADD_S": 3,
    }
    for b in range(nblocks):
        for arm, (food, add) in m.ARMS.items():
            cage = f"B{b}-{arm}"
            base = {
                "cage_id": cage, "block_id": f"B{b}", "treatment": arm,
                "food_support": str(food), "added_group": add,
            }
            allocation.append({**base, "random_seed": "20261014",
                               "design_version": m.DESIGN})
            extra_n = 0 if add == "none" else 2
            assignments.append({
                **base, "recipient_a_initial_n": "4",
                "added_larvae_initial_n": str(extra_n),
                "food_floor_cm2": "50", "expected_n_resource_visits": "3",
                "donor_source_batch_id": f"D{b}", "study_season": "PILOT_UNREAL",
            })
            for i in range(4):
                fate = m.SUCCESS if i < arm_k[arm] else "dead_immature"
                if with_unknown and b == 0 and arm == "NATURAL_ADD_A" and i == 1:
                    fate = "diapause_alive"
                fates.append({
                    "cage_id": cage, "larva_id": f"{cage}:FOCAL:{i}",
                    "cohort_role": "focal", "species": m.NAME_A, "fate": fate,
                })
            for i in range(extra_n):
                fates.append({
                    "cage_id": cage, "larva_id": f"{cage}:ADDED:{i}",
                    "cohort_role": "added",
                    "species": m.NAME_A if add == "A" else m.NAME_S,
                    "fate": "dead_immature",
                })
            for i in range(3):
                resources.append({
                    "cage_id": cage,
                    "observation_datetime": f"2026-10-{10+i:02d}T11:00:00+09:00",
                    "recipient_a_alive_n": "4", "recipient_a_late_instar_alive_n": str(i*2),
                    "added_a_alive_n": "2" if add == "A" else "0",
                    "added_s_alive_n": "2" if add == "S" else "0",
                    "accessible_leaf_area_cm2": "80" if food else "40",
                    "edible_fresh_mass_g": "10" if food else "3",
                    "accessible_leaf_area_added_cm2": "20" if food else "0",
                    "leaf_stage": "young" if food else "mature",
                    "added_leaf_stage": "young" if food else "",
                    "leaf_source_batch_id": f"D{b}",
                    "sham_handling_done": "1",
                })
    return allocation, assignments, fates, resources


def test_primary_invasive_competitor_food_rescue_and_secondary_identity_contrast():
    z = m.analyze(*fixture())
    assert z["n_randomized_cages"] == 18
    assert z["n_source_blocks"] == 3
    assert z["primary_planned_estimand"] == "S_competitor_x_food"
    assert z["primary_observed_flight_capable_adult_risk_difference"] == pytest.approx(0.75)
    assert z["other_planned_contrasts"]["identity_natural_S_minus_A"] == pytest.approx(-0.25)
    assert z["other_planned_contrasts"]["identity_supported_S_minus_A"] == pytest.approx(+0.25)
    assert z["other_planned_contrasts"]["S_competitor_natural"] == pytest.approx(-0.75)
    assert z["other_planned_contrasts"]["S_competitor_supported"] == pytest.approx(0)
    assert z["other_planned_contrasts"]["A_competitor_natural"] == pytest.approx(-0.5)
    assert z["other_planned_contrasts"]["A_competitor_supported"] == pytest.approx(-0.25)
    assert z["other_planned_contrasts"]["identity_x_food_secondary"] == pytest.approx(+0.5)
    assert z["other_planned_contrasts"]["A_competitor_x_food"] == pytest.approx(+0.25)
    assert z["worst_case_eventual_emergence_attrition_intervals"][m.PRIMARY] == pytest.approx([.75, .75])
    assert z["independent_block_bootstrap_primary_interval_if_10plus_blocks"] is None
    assert z["resource_shortage_vs_quality_vs_direct_effect_disentangled"] is False
    assert z["extra_larvae_have_same_HEADCOUNT_not_same_consumption"] is True


def test_all_randomized_focal_and_added_larvae_must_be_accounted_separately():
    allocation, assignments, fates, resources = fixture()
    missing_focal = [row for row in fates if row["larva_id"] != f"B0-NATURAL_ADD_A:FOCAL:1"]
    with pytest.raises(ValueError, match="original assigned focal"):
        m.analyze(allocation, assignments, missing_focal, resources)
    missing_extra = [row for row in fates if row["larva_id"] != f"B0-NATURAL_ADD_A:ADDED:1"]
    with pytest.raises(ValueError, match="original additional larvae"):
        m.analyze(allocation, assignments, missing_extra, resources)


def test_same_species_extra_larvae_cannot_be_mislabeled_as_focal_recipients():
    a,b,c,d = fixture()
    extra = next(x for x in c if x["cohort_role"] == "added" and x["species"] == m.NAME_A)
    extra["cohort_role"] = "focal"
    with pytest.raises(ValueError, match="original assigned focal"):
        m.analyze(a,b,c,d)


def test_assignment_and_donor_block_balance_must_not_be_rewritten_after_randomization():
    a,b,c,d = fixture()
    b[0]["added_group"] = "S"
    with pytest.raises(ValueError, match="Post-randomization"):
        m.analyze(a,b,c,d)
    a,b,c,d = fixture()
    b[0]["recipient_a_initial_n"] = "5"
    with pytest.raises(ValueError, match="original focal|must be assigned equally"):
        m.analyze(a,b,c,d)
    a,b,c,d = fixture()
    b[0]["donor_source_batch_id"] = "DIFFERENT"
    with pytest.raises(ValueError, match="donor_source_batch_id"):
        m.analyze(a,b,c,d)


def test_same_extra_headcount_across_add_A_add_S_arms_required():
    a,b,c,d = fixture()
    row = next(x for x in b if x["treatment"] == "SUPPORTED_ADD_S")
    row["added_larvae_initial_n"] = "3"
    with pytest.raises(ValueError, match="same prespecified HEADCOUNT"):
        m.analyze(a,b,c,d)


def test_missing_resources_do_not_delete_cages_or_assume_true_zero():
    a,b,c,d = fixture()
    d[0]["accessible_leaf_area_cm2"] = ""
    z = m.analyze(a,b,c,d)
    assert z["n_randomized_cages"] == 18
    assert z["per_arm"]["NATURAL_NONE"]["missing_leaf_area_visits"] == 1
    assert z["per_arm"]["NATURAL_NONE"]["late_focal_area_observed_visits"] == 6
    d.pop()
    with pytest.raises(ValueError, match="required scheduled"):
        m.analyze(a,b,c,d)


def test_postrandomization_resource_measurement_is_not_a_selection_rule():
    a,b,c,d = fixture()
    for row in d:
        if row["cage_id"].endswith("SUPPORTED_ADD_S"):
            row["accessible_leaf_area_cm2"] = "0"
    z = m.analyze(a,b,c,d)
    assert z["primary_observed_flight_capable_adult_risk_difference"] == pytest.approx(0.75)
    assert z["per_arm"]["SUPPORTED_ADD_S"]["late_focal_food_floor_breach_visits"] == 6
    assert z["primary_DOES_NOT_condition_on_achieved_food_availability"]


def test_missing_live_diapause_fate_bounds_do_not_equal_mortality():
    a,b,c,d = fixture(with_unknown=True)
    z = m.analyze(a,b,c,d)
    assert z["per_arm"]["NATURAL_ADD_A"]["living_diapausing_focal_pupae"] == 1
    assert z["per_arm"]["NATURAL_ADD_A"]["unresolved_focal_fates"] == 1
    assert z["per_arm"]["NATURAL_ADD_A"]["verified_flight_capable_adults"] == 3
    assert z["primary_observed_flight_capable_adult_risk_difference"] == pytest.approx(.75)
    lo, hi = z["worst_case_eventual_emergence_attrition_intervals"][m.PRIMARY]
    assert lo == pytest.approx(.75)
    assert hi == pytest.approx(.75)
    identity_lo, identity_hi = z["worst_case_eventual_emergence_attrition_intervals"]["identity_x_food_secondary"]
    assert identity_lo == pytest.approx(.5)
    assert identity_hi == pytest.approx(.5 + 1/12)
    assert not z["all_original_focal_fates_resolved"]
    assert z["living_diapausing_pupae_are_not_dead"]


def test_unknown_in_negative_coefficient_arm_reduces_worst_case_lower():
    a,b,c,d = fixture()
    target = next(x for x in c if x["larva_id"] == "B0-SUPPORTED_ADD_A:FOCAL:0")
    target["fate"] = "pupa_alive_pending_emergence"
    z = m.analyze(a,b,c,d)
    assert z["primary_observed_flight_capable_adult_risk_difference"] == pytest.approx(.75)
    lo, hi = z["worst_case_eventual_emergence_attrition_intervals"][m.PRIMARY]
    assert lo == pytest.approx(.75)
    assert hi == pytest.approx(.75)
    identity_lo, identity_hi = z["worst_case_eventual_emergence_attrition_intervals"]["identity_x_food_secondary"]
    assert identity_lo == pytest.approx(.5)
    assert identity_hi == pytest.approx(.5 + 1/12)


def test_invalid_extra_species_and_unresolved_pupa_code_fail():
    a,b,c,d = fixture()
    x = next(x for x in c if x["cohort_role"] == "added" and x["species"] == m.NAME_S)
    x["species"] = m.NAME_A
    with pytest.raises(ValueError, match="Added species"):
        m.analyze(a,b,c,d)
    a,b,c,d = fixture()
    c[0]["fate"] = "pupa_no_adult"
    with pytest.raises(ValueError, match="Ambiguous pupa_no_adult"):
        m.analyze(a,b,c,d)


def test_false_independent_replicates_and_nonrandomized_cages_fail():
    a,b,c,d = fixture()
    a.pop()
    with pytest.raises(ValueError, match="No original assignment|Missing randomized"):
        m.analyze(a,b,c,d)
    a,b,c,d = fixture()
    a.append(dict(a[0]))
    with pytest.raises(ValueError, match="Duplicate or blank original randomized cage"):
        m.analyze(a,b,c,d)


def test_missing_block_cannot_be_rescued_by_flattening_cages():
    a,b,c,d = fixture()
    a = [x for x in a if x["cage_id"] != "B0-NATURAL_NONE"]
    b = [x for x in b if x["cage_id"] != "B0-NATURAL_NONE"]
    c = [x for x in c if x["cage_id"] != "B0-NATURAL_NONE"]
    d = [x for x in d if x["cage_id"] != "B0-NATURAL_NONE"]
    with pytest.raises(ValueError, match="all six arms"):
        m.analyze(a,b,c,d)


def test_bootstrap_available_only_for_at_least_ten_independent_blocks():
    z = m.analyze(*fixture(nblocks=10))
    assert z["n_source_blocks"] == 10
    assert z["independent_block_bootstrap_primary_interval_if_10plus_blocks"] == pytest.approx([0.75, 0.75])
    assert z["n_randomized_cages"] == 60


def test_true_cage_zero_and_missing_area_not_confused():
    a,b,c,d = fixture()
    d[0]["accessible_leaf_area_cm2"] = "0"
    z = m.analyze(a,b,c,d)
    assert z["per_arm"]["NATURAL_NONE"]["missing_leaf_area_visits"] == 0
    assert z["per_arm"]["NATURAL_NONE"]["late_focal_food_floor_breach_visits"] == 6


def test_food_rescue_zero_with_equal_food_response_but_identity_response_still_possible():
    # Keep posthoc synthetic scenarios distinct: both additions can be equally
    # rescued by access to foliage; this cannot be described as species identity.
    risks = {
        "NATURAL_NONE": .8, "NATURAL_ADD_A": .4, "NATURAL_ADD_S": .4,
        "SUPPORTED_NONE": .8, "SUPPORTED_ADD_A": .8, "SUPPORTED_ADD_S": .8,
    }
    contrasts = m._contrasts(risks)
    assert contrasts["S_competitor_x_food"] == pytest.approx(.4)
    assert contrasts["identity_x_food_secondary"] == pytest.approx(0)
    # A nonzero total invader rescue does not imply identity-specific relief.
