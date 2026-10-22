"""Independent frozen v08 group score and missingness audit tests."""
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_v08_otherplant_lgm_null_independent_v09 import (
    sign_counts,same_edge_recompute,number,summarize_prequalified_species
)

def fake_group():
    # Four exact common original edges:
    # actual 0BP same=3/4, actual 21ka same=2/4 (one split)
    # other 0BP same=2/4, other 21ka same=2/4 (no net change)
    return dict(all_original_edges_n="8",fourway_exact_common_original_edges_n="4",
      actual_0BP_same_fraction=".75",actual_21ka_same_fraction=".5",
      alternate_0BP_same_fraction=".5",alternate_21ka_same_fraction=".5",
      actual_paired_delta="-.25",alternate_paired_delta="0",
      actual_minus_alternate_paired_delta="-.25",
      modern_0BP_baseline_same_patch_gap=".25",
      actual_same_to_different_n="1",actual_different_to_same_n="0",
      alternate_same_to_different_n="0",alternate_different_to_same_n="0")

def test_fourway_paired_difference_independent_recalculation():
    assert same_edge_recompute(fake_group())==pytest.approx(-.25)

def test_0BP_baseline_gap_and_transition_sources_fail_closed():
    k=fake_group();k["modern_0BP_baseline_same_patch_gap"]="0.1"
    with pytest.raises(ValueError,match="baselines"):same_edge_recompute(k)
    k=fake_group();k["actual_same_to_different_n"]="0"
    with pytest.raises(ValueError,match="split/merge"):same_edge_recompute(k)
    k=fake_group();k["actual_minus_alternate_paired_delta"]="1"
    with pytest.raises(ValueError,match="difference-in-differences"):same_edge_recompute(k)

def test_incomplete_fourway_edge_remains_unevaluable_not_zero():
    k=fake_group();k["fourway_exact_common_original_edges_n"]="0"
    for f in ["actual_paired_delta","alternate_paired_delta",
              "actual_minus_alternate_paired_delta","modern_0BP_baseline_same_patch_gap"]:
        k[f]=""
    assert same_edge_recompute(k) is None
    k["actual_paired_delta"]="0"
    with pytest.raises(ValueError,match="fake genetic-climate"):
        same_edge_recompute(k)

def test_null_and_true_zero_are_distinct():
    assert number({"x":""},"x",required=False) is None
    assert number({"x":"0"},"x")==0

def test_species_signs_define_host_more_glacial_split_as_negative():
    assert sign_counts([-.2,.1,0,0])=={
        "actual_host_more_glacial_patch_loss":1,
        "other_plant_more_glacial_patch_loss":1,
        "equal_delta":2
    }

def test_original_data_source_paths_and_292_2108_12648_frozen_identifiers():
    from audit_v08_otherplant_lgm_null_independent_v09 import SOURCES
    assert len(SOURCES)==4
    assert SOURCES["result"]=="7fcd67aa367d138960aa7d0d94e6482139a5678a1c2a4be08a38dcbe4264bc55"
