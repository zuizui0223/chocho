"""Test source-locked support audit; no false host-absence or transport claim."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import math
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/audit_stage_support.py"
spec = importlib.util.spec_from_file_location("maritime_support", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_source_locked_support_and_separate_denominators():
    # The anonymous GEB reviewer-stage intentionally excludes these original
    # exploratory provenance inputs; this exception does NOT apply to paper CI.
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island inputs excluded by original anonymous GEB export policy")
    r = mod.audit(ROOT)
    all_events = r["original_coordinate_butterfly_history"]["all"]
    both = r["original_coordinate_butterfly_history"]["both_period_opportunity"]
    assert all_events["island_first_record_events"] == 12028
    assert all_events["apparent_nearest_prior_record_gap_ge500_km"] == 350
    assert both["island_first_record_events"] == 3136
    assert both["apparent_nearest_prior_record_gap_ge500_km"] == 164
    assert both["far500_nearest_prior_indexed_record_on_island"] == 104
    assert math.isclose(both["far500_nearest_prior_indexed_record_on_island_fraction"],
                        104/164, abs_tol=1e-12)
    assert both["apparent_nearest_prior_record_gap_lt100_km"] == 2529
    assert both["near100_nearest_prior_indexed_record_on_island"] == 1309

    h = r["early_port_host_link_support"]
    assert h["common_new_butterfly_species_island_pair_denominator"] == 1163
    assert h["retrospective_plant_2001_2024_matches"] == 208
    assert h["genuinely_earlier_plant_2001_2010_matches"] == 5
    assert h["earlier_host_butterfly_species"] == 5
    assert h["earlier_host_islands"] == 7
    assert h["unique_species_island_endpoints_from_pairs"] == 10
    assert h["endpoints_observed_in_separate_2011_2017_followup_cohort"] == 6
    assert h["endpoints_with_any_subsequent_record"] == 6
    assert h["endpoints_with_record_in_year_five"] == 2
    assert h["selected_endpoint_is_not_unbiased_success_comparison"]

    c = r["separate_five_year_followup_cohort"]
    assert c["original_species_island_first_record_events"] == 5665
    assert c["port_host_group_matched_but_not_dated"] == 332
    assert c["not_evaluable"] == 5333
    assert r["prior_shipping_connectivity_support"]["such_islands"] == 27
    assert r["prior_shipping_connectivity_support"]["islands_overlap_prior_51_island_host_network"] == 21
    assert r["separate_cohorts_not_joint_risk_set"]
    assert not any(r["mechanisms_identified_now"].values())
