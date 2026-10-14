"""Synthetic source-integrity fixtures; NO collected leaf or butterfly records."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/audit_kyoto_leaf_stage_pilot.py"
SPEC = importlib.util.spec_from_file_location("pilot_stage", PATH)
assert SPEC and SPEC.loader
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def leaves(donors=("D1", "D2", "D3")):
    return [
        {
            "donor_plant_id": donor,
            "leaf_batch_id": f"{donor}-{stage}",
            "host_species": m.HOST,
            "leaf_stage": stage,
            "sample_date": "2026-10-14",
            "shoot_node": "6" if stage == "young" else "11",
            "whole_plant_age_years": "",
            "leaf_area_cm2": "10",
            "fresh_mass_g": "1",
            "dry_mass_g": "0.2" if stage == "young" else "0.3",
        }
        for donor in donors for stage in ("young", "mature")
    ]


def feeding(families=("F1", "F2", "F3")):
    return [
        {
            "cup_id": f"{species[:1]}-{stage}-{fam}",
            "source_block_id": f"B{j}",
            "maternal_family_id": f"{species[:1]}{fam}",
            "species": species,
            "instar_class": "L1",
            "assigned_leaf_stage": stage,
            "source_donor_plant_id": ["D1", "D2", "D3"][j],
            "assay_leaf_batch_id": f"fresh-{species[:1]}-{stage}-{fam}",
            "initial_larvae_n": "2",
            "alive_end_n": "2" if stage == "young" else "1",
            "observation_minutes": "60",
            "feeding_confirmed": "yes" if stage == "young" else "unknown",
        }
        for species in sorted(m.SPECIES)
        for stage in ("young", "mature")
        for j, fam in enumerate(families)
    ]


def facilities(n=10, each=3):
    return [{"greenhouse_id": f"G{x}", "eligible_whole_cages": str(each)}
            for x in range(n)]


def test_stage_specific_plant_calibration_paired_donors_not_pseudoreplicated_leaves():
    x = m.leaf_calibration(leaves())
    assert x["independent_original_donor_plants"] == 3
    assert x["donors_paired_young_and_mature"] == 3
    assert x["structural_matched_donor_gate"]
    assert x["by_stage"]["young"]["donor_equal_mean_dry_fraction"] == pytest.approx(.2)
    assert x["by_stage"]["mature"]["donor_equal_mean_dry_fraction"] == pytest.approx(.3)
    assert x["by_stage"]["young"]["donor_equal_mean_dry_mg_per_cm2"] == pytest.approx(20)
    assert x["by_stage"]["mature"]["donor_equal_mean_dry_mg_per_cm2"] == pytest.approx(30)
    assert x["paired_within_donor_young_minus_mature_dry_fraction_mean"] == pytest.approx(-.1)
    assert x["not_a_claim_that_two_donors_support_statistical_generality"]


def test_source_matched_nochoice_assay_observations_and_unknown_denominator():
    x = m.feeding_assay(feeding(), m.leaf_calibration(leaves()))
    assert x["four_core_neonate_cells_present_with_minimum_source_repetition"]
    mature = x["by_species_instar_leaf_stage"]["Sericinus montela|L1|mature"]
    young = x["by_species_instar_leaf_stage"]["Sericinus montela|L1|young"]
    assert mature["feeding_unknown"] == 3
    assert mature["cup_feeding_confirmed_fraction_among_observed"] is None
    assert young["cup_feeding_confirmed_fraction_among_observed"] == 1
    assert young["whole_assay_cups"] == 3
    assert young["independent_maternal_family_ids"] == 3
    assert young["not_an_individual_level_or_population_level_survival_effect"]


def test_historical_three_cages_per_greenhouse_cannot_host_four_or_six_complete_arms():
    x = m.facility_audit(facilities(10, 3))
    assert x["current_eligible_total_whole_cages"] == 30
    assert x["greenhouse_local_complete_four_arm_blocks_capacity"] == 0
    assert x["greenhouse_local_complete_six_arm_blocks_capacity"] == 0
    assert x["greenhouses_supporting_at_least_four_cages"] == 0
    assert not x["four_arm_complete_greenhouse_block_structurally_available"]
    assert not x["six_arm_complete_greenhouse_block_structurally_available"]
    assert x["greenhouse_incomplete_block_design_not_analyzed_here"]
    # These numbers are a hypothetical inventory matching the historical
    # 2012 layout, NOT confirmation of present-day facility access.


def test_pilot_structural_gate_is_not_animal_experiment_or_causal_result():
    result = m.audit(leaves(), feeding(), facilities())
    assert result["operational_gate"] == "STRUCTURAL_PILOT_SUPPORTED_NOT_STUDY_APPROVED"
    assert result["true_demographic_fitness_response_estimated"] is False
    assert result["butterflies_reared_by_script"] is False
    assert result["preliminary_2023_neonate_study_reproduced"] is False
    assert result["four_arm_within_greenhouse_design_requires_redesign_or_more_capacity"]
    assert result["six_arm_within_greenhouse_design_requires_redesign_or_more_capacity"]
    assert result["GEB_main_unchanged"] is True


def test_brood_pseudoreplication_does_not_pass_source_independence_gate():
    rows = feeding(families=("F1", "F2", "F3"))
    for row in rows:
        row["maternal_family_id"] = f'{row["species"]}-ONE_FAMILY'
    result = m.audit(leaves(), rows, facilities())
    assert result["operational_gate"] == "PILOT_INPUT_INCOMPLETE_OR_INSUFFICIENT"
    assert not result["no_choice_neonate_feeding_pilot"][
        "four_core_neonate_cells_present_with_minimum_source_repetition"
    ]


def test_missing_one_foliage_class_in_pilot_prevents_stage_comparison():
    rows = [r for r in feeding() if not (
        r["species"] == "Sericinus montela" and r["assigned_leaf_stage"] == "mature"
    )]
    result = m.audit(leaves(), rows, facilities())
    assert result["operational_gate"] == "PILOT_INPUT_INCOMPLETE_OR_INSUFFICIENT"


def test_stage_source_cannot_be_fabricated_from_another_host_plant():
    rows = leaves()
    rows[0]["host_species"] = "Aristolochia kaempferi"
    with pytest.raises(ValueError, match="Aristolochia debilis"):
        m.leaf_calibration(rows)


def test_fresh_dry_water_content_sanity_and_unique_samples():
    bad = leaves()
    bad[0]["dry_mass_g"] = "2"
    with pytest.raises(ValueError, match="exceeds original fresh"):
        m.leaf_calibration(bad)
    bad = leaves()
    bad[0]["leaf_batch_id"] = bad[1]["leaf_batch_id"]
    with pytest.raises(ValueError, match="duplicate"):
        m.leaf_calibration(bad)


def test_assay_source_matches_original_donor_and_leaf_age():
    bad = feeding()
    bad[0]["source_donor_plant_id"] = "INVENTED_DONOR"
    with pytest.raises(ValueError, match="unknown plant donor"):
        m.feeding_assay(bad, m.leaf_calibration(leaves()))
    bad = feeding()
    bad[0]["assigned_leaf_stage"] = "neither"
    with pytest.raises(ValueError, match="uncalibrated"):
        m.feeding_assay(bad, m.leaf_calibration(leaves()))


def test_pilot_fate_and_feeding_observation_must_not_be_inferred_from_each_other():
    bad = feeding()
    bad[0]["alive_end_n"] = "3"
    with pytest.raises(ValueError, match="more alive"):
        m.feeding_assay(bad, m.leaf_calibration(leaves()))
    bad = feeding()
    bad[0]["feeding_confirmed"] = ""
    with pytest.raises(ValueError, match="yes/no/unknown"):
        m.feeding_assay(bad, m.leaf_calibration(leaves()))


def test_incomplete_or_duplicated_physical_greenhouse_inventory_fails_closed():
    bad = facilities()
    bad[1]["greenhouse_id"] = bad[0]["greenhouse_id"]
    with pytest.raises(ValueError, match="duplicate"):
        m.facility_audit(bad)
    bad = facilities()
    bad[0]["eligible_whole_cages"] = "-1"
    with pytest.raises(ValueError, match="invalid"):
        m.facility_audit(bad)


def test_no_repeat_cups_or_foliage_batches():
    bad = feeding()
    bad[1]["cup_id"] = bad[0]["cup_id"]
    with pytest.raises(ValueError, match="unique"):
        m.feeding_assay(bad, m.leaf_calibration(leaves()))
    bad = feeding()
    bad[1]["assay_leaf_batch_id"] = bad[0]["assay_leaf_batch_id"]
    with pytest.raises(ValueError, match="unique"):
        m.feeding_assay(bad, m.leaf_calibration(leaves()))


def test_real_capacity_rule_distinguishes_four_and_six_cage_greenhouses():
    x = m.facility_audit([
        {"greenhouse_id": "G4", "eligible_whole_cages": "4"},
        {"greenhouse_id": "G5", "eligible_whole_cages": "5"},
    ])
    assert x["four_arm_complete_greenhouse_block_structurally_available"]
    assert not x["six_arm_complete_greenhouse_block_structurally_available"]
    assert x["greenhouse_local_complete_four_arm_blocks_capacity"] == 2
    assert x["greenhouse_local_complete_six_arm_blocks_capacity"] == 0
    x = m.facility_audit([{"greenhouse_id": "G6", "eligible_whole_cages": "6"}])
    assert x["six_arm_complete_greenhouse_block_structurally_available"]
    assert x["greenhouse_local_complete_six_arm_blocks_capacity"] == 1
