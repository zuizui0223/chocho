"""Protect supplementary 5665-event source-locked retrospective matching.

Does not claim independently dated larval hosts or a dispersal pathway.
"""
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031"
SPEC=importlib.util.spec_from_file_location("chocho_reindex_matched",DIR/"analyze_host_reindex_matched.py")
assert SPEC is not None and SPEC.loader is not None
MOD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_original_event_support_and_matched_analysis():
    observed=MOD.analyze()
    expected=json.loads((DIR/"host_reindex_matched_results.json").read_text(encoding="utf8"))
    assert observed["butterfly_first_events_total"]==expected["butterfly_first_events_total"]==5665
    assert {k:x["events"] for k,x in observed["groups"].items()}=={
        k:x["events"] for k,x in expected["groups"].items()}
    assert sum(v["events"] for v in observed["groups"].values())==5665
    for group in observed["groups"]:
        for field in ("reindexed_within_5y","raw_reindex_rate"):
            assert math.isclose(observed["groups"][group][field],
                                expected["groups"][group][field],rel_tol=0,abs_tol=1e-12)
    for key in ("species","island","island_year"):
        a=observed["matched_comparisons"][key]
        b=expected["matched_comparisons"][key]
        assert a["paired_groups"]==b["paired_groups"]
        assert a["groups_with_positive_difference"]==b["groups_with_positive_difference"]
        assert a["groups_with_negative_difference"]==b["groups_with_negative_difference"]
        assert a["groups_with_zero_difference"]==b["groups_with_zero_difference"]
        for k in ("mean_rate_difference",):
            assert math.isclose(a[k],b[k],abs_tol=1e-12,rel_tol=0)
        for x,y in zip(a["species_or_group_bootstrap_95"],b["species_or_group_bootstrap_95"]):
            assert math.isclose(x,y,abs_tol=1e-12,rel_tol=0)
        assert a["at_least_2_events_per_group"]["paired_groups"]==b["at_least_2_events_per_group"]["paired_groups"]
        assert math.isclose(a["at_least_2_events_per_group"]["mean_rate_difference"],b["at_least_2_events_per_group"]["mean_rate_difference"],abs_tol=1e-12,rel_tol=0)
    a=observed["matched_comparisons"]["species"]["wing_size_vs_within_species_difference_spearman"]
    b=expected["matched_comparisons"]["species"]["wing_size_vs_within_species_difference_spearman"]
    assert math.isclose(a,b,abs_tol=1e-12,rel_tol=0)
    assert expected["status"]=="DESCRIPTIVE_CONFOUNDING_DIAGNOSTIC_NOT_HISTORIC_HOST_OR_MOVEMENT_EVIDENCE"


def test_reindex_diagnostic_not_recast_as_historical_experiment():
    content=(DIR/"README.md").read_text(encoding="utf8")
    assert "2001–2024" in content or "2001-2024" in content
    assert "反転" in content
