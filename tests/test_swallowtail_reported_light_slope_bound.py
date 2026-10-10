"""Mathematical unit tests, NOT inferred new Korean ramet data."""
from __future__ import annotations
import importlib.util,math
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/check_swallowtail_reported_light_slope_bound.py"
spec=importlib.util.spec_from_file_location("lightbound",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_global_expectation_bound_from_single_intercept_logit_slope():
    b=m.one_slope_expected_overlap_lower_bound()
    assert 14.6<b["lower_bound_on_EXPECTED_both"]<14.7
    assert b["lower_bound_on_EXPECTED_both"]>m.BOTH
    assert b["does_not_bound_observed_cooccurrence"] is True
    assert math.isclose(m.one_slope_expected_overlap_lower_bound(s_odds=1)["lower_bound_on_EXPECTED_both"],
                        m.S*m.A/m.N)

def test_numerical_endpoints_match_original_marginal_prevalence():
    for fraction in (.001,.25,.5,.6205,.75,.999):
        out=m.two_extreme_light_mixture(fraction)
        assert math.isclose(out["expected_S"],209,rel_tol=0,abs_tol=1e-9)
        assert math.isclose(out["expected_A"],129,rel_tol=0,abs_tol=1e-9)
        assert 0<out["p_S_at_RLI_0"]<out["p_S_at_RLI_100"]<1
        assert 0<out["p_A_at_RLI_100"]<out["p_A_at_RLI_0"]<1

def test_previous_invented_counterexample_requires_larger_light_slopes():
    x=m.original_toy_required_or()["previous_invented_two_light_strata_not_measured_Korean_RLI"]
    assert 1.040<x["implied_S_light_OR_per_one_percent"]<1.041
    assert .967<x["implied_A_light_OR_per_one_percent"]<.968
    assert x["implied_S_light_OR_per_one_percent"]>m.OR_S
    assert x["implied_A_light_OR_per_one_percent"]<m.OR_A
    assert x["reference_categories_reproduced_exactly"]=={"neither":542,"S_only":204,"A_only":124,"both":5}

def test_restrictive_0_to_100_light_response_is_not_enough_in_expectation():
    r=m.evaluate()
    out=r["two_RLI_endpoint_distribution_grid_0_to_100_numerical_minimum"]
    assert 22.8<out["expected_both_under_conditional_independence"]<23.1
    assert .60<out["fraction_at_RLI_100"]<.64
    assert r["rigorous_EXPECTED_overlap_lower_bound_under_common_intercept_assumption"]["lower_bound_on_EXPECTED_both"]>5
    assert r["observed_five_can_be_sampling_variation_even_with_higher_mean"] is True
    assert r["cannot_infer_negative_competition_or_true_light_response_as_unique_mechanism"] is True
    assert r["GEB_PR38_untouched"] is True

def test_invalid_slope_assumption_rejected():
    for data in ({"s_odds":.95},{"span":-1},{"n":0},{"a":-1}):
        try:m.one_slope_expected_overlap_lower_bound(**data)
        except ValueError:pass
        else:raise AssertionError("invalid premise accepted")
