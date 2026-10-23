"""Exact rational weighting/reference regression tests for frozen 26 source cases."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
from fractions import Fraction

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("all26_weighting",
    ROOT/"scripts/analyze_gbif_all26_weighting_stability_20261023.py")
mod=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)
BASE=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"

def test_two_by_two_rational_case_and_record_tails():
    b=[{"butterfly":"A a","candidate_k":1,"gadm":{"A.1":1},"N":1,"observed_count":1},
       {"butterfly":"B b","candidate_k":1,"gadm":{"B.1":1},"N":1,"observed_count":1}]
    p=[{"accepted_id":"X","links":frozenset({"A a"}),"gadm":frozenset({"A.1"})},
       {"accepted_id":"Y","links":frozenset({"B b"}),"gadm":frozenset({"B.1"})}]
    states=mod.graph_vectors(b,p)
    assert states["states"]==2
    dist=mod.per_country(states,"case",{"A a":1,"B b":1})
    z=mod.exact_reference({"XX":dist},2)
    assert z["observed_fraction"]==1
    assert z["reference_expected_fraction"]==.5
    assert z["exact_label_upper_tail_NOT_biological_p"]==.5

def test_exact_frozen_all26_reproduction_and_weight_sensitivity():
    bot=json.loads((BASE/"gbif_all26_host_gadm_source_20261023.json").read_text(encoding="utf-8"))
    bf=json.loads((BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json").read_text(encoding="utf-8"))
    out=mod.analyse(bot,bf)
    z=out["metrics"]
    assert out["n_country_butterfly_cases"]==26
    assert out["fixed_degree_graphs"]=="13826598912"
    assert z["record"]["observed_minus_reference_pp"]==.2023547
    assert z["case"]["observed_minus_reference_pp"]==6.8580651
    assert z["case"]["exact_label_upper_tail_NOT_biological_p"]==.027714830410
    assert z["species"]["exact_label_upper_tail_NOT_biological_p"]==.099317206042
    assert z["country"]["exact_label_upper_tail_NOT_biological_p"]==.045550194376
    assert z["case_n5"]["source_unit_denominator"]==11
    assert z["case_n5"]["observed_minus_reference_pp"]==10.2756286
    omitted={x["omitted_country"]:x for x in out["posthoc_leave_one_country_out_case_metric"]}
    assert omitted["CU"]["exact_label_upper_tail_NOT_biological_p"]==.099515056794
    assert omitted["DO"]["exact_label_upper_tail_NOT_biological_p"]==.076183338991
    assert "NOT a biological p-value" in out["biology_warning"]

def test_archived_full_exact_weighting_table_is_source_exact():
    a=json.loads((BASE/"gbif_all26_weighting_stability_20261023.json").read_text(encoding="utf8"))
    bot=json.loads((BASE/"gbif_all26_host_gadm_source_20261023.json").read_text(encoding="utf8"))
    bf=json.loads((BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json").read_text(encoding="utf8"))
    out=mod.analyse(bot,bf)
    assert a["source_population"]["butterfly_country_cases"]==26
    assert a["status"]=="EXACT_POSTHOC_SOURCE_GEOMETRY_AUDIT_NOT_ECOLOGICAL_HYPOTHESIS_TEST"
    assert out["metrics"]==a["weighting_results"]
    assert out["posthoc_leave_one_country_out_case_metric"]==a["country_omission_descriptives"]
