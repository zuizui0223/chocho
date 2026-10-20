"""Regression tests for independent v0.8 host identity climatic patch audit."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_v08_host_identity_identifiability_v09 import (
    audit_group_row,audit_species_row,describe_species,grouped_cluster_bootstrap,
    validate_and_audit,SOURCE_SHA
)

def group(n=8,actual0=.75,actual21=.5,alt0=.5,alt21=.5):
    x={"insect_species":"Testus first","scenario":"analog","radius_km":"100",
       "other_plant_group_index":"0","original_0BP_climate_moment_gap":"0.1",
       "known_actual_host_ids":"1","other_plant_ids":"2",
       "all_original_edges_n":"12","fourway_exact_common_original_edges_n":str(n)}
    if n:
        x.update({"actual_0BP_same_fraction":str(actual0),
           "actual_21ka_same_fraction":str(actual21),
           "alternate_0BP_same_fraction":str(alt0),
           "alternate_21ka_same_fraction":str(alt21),
           "actual_paired_delta":str(actual21-actual0),
           "alternate_paired_delta":str(alt21-alt0),
           "actual_minus_alternate_paired_delta":str(actual21-actual0-alt21+alt0),
           "modern_0BP_baseline_same_patch_gap":str(actual0-alt0),
           "actual_same_to_different_n":"2","actual_different_to_same_n":"0",
           "alternate_same_to_different_n":"0","alternate_different_to_same_n":"0"})
    else:
        for field in ("actual_0BP_same_fraction","actual_21ka_same_fraction",
                      "alternate_0BP_same_fraction","alternate_21ka_same_fraction",
                      "actual_paired_delta","alternate_paired_delta",
                      "actual_minus_alternate_paired_delta",
                      "modern_0BP_baseline_same_patch_gap"):x[field]=""
    return x

def test_zero_eligible_fourway_pairs_cannot_be_assigned_effect_zero():
    assert audit_group_row(group(n=0)) is None
    q=group(n=0);q["actual_minus_alternate_paired_delta"]="0"
    with pytest.raises(ValueError):
        audit_group_row(q)

def test_actual_vs_alt_with_same_exact_edge_denominator_and_did_identity():
    p=audit_group_row(group())
    assert p["n"]==8
    assert p["diff"]==pytest.approx(-.25)
    assert p["baseline"]==pytest.approx(.25)
    q=group();q["actual_minus_alternate_paired_delta"]="0.1"
    with pytest.raises(ValueError):
        audit_group_row(q)

def test_source_one_group_species_audit_recovers_median_and_zero_missingness():
    g=group()
    row={"insect_species":"Testus first","source_frozen_other_plant_group_n":"1",
         "other_plant_groups_with_common_original_edges_n":"1",
         "other_plant_groups_with_at_least_5_common_edges_n":"1",
         "status":"COMPUTABLE_PAIRWISE_OTHER_PLANT_NEGATIVE_CONTROL",
         "median_actual_minus_otherplant_same_edge_difference":"-.25",
         "median_0BP_actual_vs_otherplant_baseline_gap":".25",
         "median_fourway_same_edge_n_per_group":"8",
         "median_difference_on_groups_with_5plus_edges":"-.25"}
    audit_species_row(row,[g])
    row["median_actual_minus_otherplant_same_edge_difference"]="0.0"
    with pytest.raises(ValueError):audit_species_row(row,[g])

def test_exactly_zero_did_is_not_missing_and_hostset_dependence_reported():
    rows=[{"insect_species":"A a","_hostset":"10",
           "median_actual_minus_otherplant_same_edge_difference":0.,
           "median_0BP_actual_vs_otherplant_baseline_gap":0.},
          {"insect_species":"A b","_hostset":"10",
           "median_actual_minus_otherplant_same_edge_difference":.1,
           "median_0BP_actual_vs_otherplant_baseline_gap":.1},
          {"insect_species":"B a","_hostset":"20",
           "median_actual_minus_otherplant_same_edge_difference":-.1,
           "median_0BP_actual_vs_otherplant_baseline_gap":-.1}]
    s=describe_species(rows)
    assert s["species"]==3 and s["genus_n"]==2
    assert s["different_known_host_combinations"]==2
    assert s["ties_no_difference_species"]==1
    assert s["actual_host_more_split_species"]==1
    assert s["other_plant_more_split_species"]==1
    assert s["median_actual_minus_otherplant_delta"]==0

def test_genus_bootstrap_deterministic_and_not_confirmatory():
    rows=[{"insect_species":"A a","difference":.1},
          {"insect_species":"A b","difference":.2},
          {"insect_species":"B a","difference":-.1}]
    x=grouped_cluster_bootstrap(rows,"difference",seed=32,draws=40)
    assert x==grouped_cluster_bootstrap(rows,"difference",seed=32,draws=40)
    assert x["genus_clusters"]==2
    assert x["confirmatory_inference"] is False
    assert x["captures_shared_plant_donors_across_genera"] is False

def test_complete_original_v08_files_292_and_2108_groups_check():
    project=Path(__file__).resolve().parents[1]
    root=project/"results/lgm_292_source_matched_otherplant_sameedge_null_v08"
    groups=root/"all_source_plant_group_fourway_contact_v08.csv"
    species=root/"all_292_insect_summary_v08.csv"
    summary=root/"summary.json"
    assert all(p.is_file() for p in (groups,species,summary))
    _,r=validate_and_audit(species,groups,summary)
    assert r["source_preselected_species"]==292
    assert r["source_frozen_donor_groups"]==2108
    assert r["source_group_scenario_radius_records"]==12648
    p=r["by_scenario_radius"]["analog_100km"]
    assert p["groups_no_fourway_common_edge"]==608
    assert p["groups_at_least_one_common_edge"]==1500
    assert p["groups_at_least_five_common_edges"]==1176
    assert p["groups_exact_zero_effect_among_contactable"]==563
    assert p["species_at_least_one_donor"]["species"]==226
    assert p["species_at_least_one_donor"]["actual_host_more_split_species"]==62
    assert p["species_at_least_one_donor"]["other_plant_more_split_species"]==79
    assert p["species_at_least_one_donor"]["ties_no_difference_species"]==85
    assert p["species_at_least_four_donors_with_five_pairs"]["species"]==176
    assert p["species_at_least_four_donors_with_five_pairs"]["median_actual_minus_otherplant_delta"]==0
