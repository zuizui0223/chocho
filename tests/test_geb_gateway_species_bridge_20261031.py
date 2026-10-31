"""CI: species-level global GEB opportunity and island-port support joins.

No inference that a WGSRPD3 regional opportunity corresponds to a port
arrival, local larval feeding, true colonization, or observed host-driven growth.
"""
from pathlib import Path
import importlib.util,json,math
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/geb_port_island_realization_20261031"
FILE=D/"GEB_to_port_species_bridge.json"
spec=importlib.util.spec_from_file_location("gateway_source",D/"build_geb_gateway_species_bridge.py")
gate=importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)

def test_pinned_global_GEB_and_observed_port_record_support():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island-port pipeline is not part of the independent GEB anonymous review pack")
    expected=json.loads(FILE.read_text(encoding="utf8"))
    actual=gate.build(ROOT)
    assert expected.keys()==actual.keys()
    assert expected["source_git_blobs"]==actual["source_git_blobs"]
    A=expected["GEB_resource_opportunity"]
    assert (A["butterfly_species"],A["expanded_species"],A["native_butterfly_x_WGSRPD3_units"],A["introduced_added_units"])==(239,206,26530,14553)
    assert math.isclose(A["expanded_percent_relative_to_native"],100*14553/26530)
    B=expected["aligned_port_candidate_comparison"]
    assert (B["species"],B["with_global_geb_added_resource_units"],B["species_with_any_2001_2024_island_port_alien_candidate_host_record"])==(225,197,159)
    for key in ("spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts",
                "partial_rank_assoc_after_adjusting_host_count_only",
                "partial_rank_assoc_after_adjusting_log1p_host_count_native_units_and_prior_recorded_islands"):
        assert math.isclose(B[key],actual["aligned_port_candidate_comparison"][key],abs_tol=1e-10)
    assert math.isclose(B["spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"],0.6524282493600642,abs_tol=1e-10)
    assert math.isclose(B["partial_rank_assoc_after_adjusting_log1p_host_count_native_units_and_prior_recorded_islands"],0.5495755354184105,abs_tol=1e-10)
    C=expected["indexed_2011_2017_island_first_records"]
    assert (C["island_count"],C["events"],C["species"],C["events_with_positive_global_resource_gain"],C["species_with_positive_global_resource_gain"])==(17,121,64,120,63)
    assert C["any_reindexed_during_next_five_years"]==102
    assert math.isclose(expected["top10_plant_contributors_from_geb"]["fraction_of_all_added_resource_units"],0.2510105136113853,abs_tol=1e-12)
    assert expected["primary_limits"]["not_a_estimate_of_fraction_of_55_percent_realized"]
