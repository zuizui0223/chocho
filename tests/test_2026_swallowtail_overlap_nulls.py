"""Pure mathematical validation: hypothetical reference assumptions, not field data."""
from __future__ import annotations
import importlib.util,math,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(ROOT))
P=ROOT/"analyze_2026_swallowtail_overlap_nulls.py"
spec=importlib.util.spec_from_file_location("overlap",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_small_reference_exact_distribution_and_tail():
    d=m.hypergeom_pmf(2,1,1)
    assert d=={0:0.5,1:0.5}
    d=m.hypergeom_pmf(4,2,2)
    assert math.isclose(d[0],1/6)
    assert math.isclose(d[1],4/6)
    assert math.isclose(d[2],1/6)
    assert m.summary(d,0,["toy"])["lower_tail_probability_P_X_at_most_observed"]==1/6

def test_exact_convolution_preserves_probability_and_stratum_references():
    c=m.convolution([{0:.5,1:.5},{0:.5,1:.5}])
    assert c=={0:.25,1:.5,2:.25}
    assert m.quantile(c,.025)==0
    assert m.quantile(c,.975)==2
    try:
        m.hypergeom_pmf(2,3,1)
    except ValueError:pass
    else:raise AssertionError("impossible margins accepted")

def test_published_denominator_gate_and_marginal_null_expectations():
    a=m.analyze()
    assert a["print_count_gate"]=={"methods_quadrats":255,"breakdown_quadrats":215,"author_confirmation_required":True}
    assert a["ramet_simple_reference"]["observed_both"]==5
    assert math.isclose(a["ramet_simple_reference"]["expected_both_under_fixed_margins"],209*129/875,rel_tol=1e-10)
    assert a["quadrat_site_margins_reference"]["observed_both"]==11
    assert a["quadrat_month_margins_reference"]["observed_both"]==11
    assert math.isclose(a["quadrat_site_margins_reference"]["expected_both_under_fixed_margins"],48.27498452012384,rel_tol=1e-10)
    assert math.isclose(a["quadrat_month_margins_reference"]["expected_both_under_fixed_margins"],44.502574002574,rel_tol=1e-10)
    assert a["site_and_month_reference_not_independent_tests"] is True
    assert a["new_ecological_or_causal_effect_estimated"] is False

def test_environmental_sorting_nonidentifiability_constructive_example():
    x=m.two_strata_habitat_counterexample()
    assert x["total_margins"]=={"N":875,"S":209,"A":129}
    assert math.isclose(x["expected_overlap_from_environmental_sorting_alone"],200*5/400+9*124/475)
    assert 4.8<x["expected_overlap_from_environmental_sorting_alone"]<5
    assert x["kind"].startswith("HYPOTHETICAL")
    assert x["all_printed_presence_categories_exactly_reproduced"]=={
        "neither":542,"S_only":204,"A_only":124,"both":5
    }
    sunny,shady=x["strata"]
    assert sunny["contingency_cells"]=={"neither":198,"S_only":197,"A_only":2,"both":3}
    assert shady["contingency_cells"]=={"neither":344,"S_only":7,"A_only":122,"both":2}
    assert abs(sunny["both"]-sunny["within_stratum_independence_expected_both"])<=.5
    assert abs(shady["both"]-shady["within_stratum_independence_expected_both"])<=.5
    assert "not estimated Korean" in x["proof_scope"]

def test_no_collapsing_print_only_sums_into_competition_results():
    a=m.analyze()
    for k in ["ramet_simple_reference","quadrat_site_margins_reference","quadrat_month_margins_reference"]:
        assert a[k]["no_competition_inference"] is True
        assert 0<=a[k]["lower_tail_probability_P_X_at_most_observed"]<=1
        assert a[k]["observed_over_expected_ratio"]<1
    assert a["print_count_gate"]["author_confirmation_required"] is True
    assert a["GEB_PR38_untouched"] is True
