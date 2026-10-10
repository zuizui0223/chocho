#!/usr/bin/env python3
"""Conditional mathematical sensitivity check on 2026 swallowtail light ORs.

Data are PRINTED ARTICLE MARGINS, not raw ramet light measurements.
This conditional one-intercept model is distinct from the authors' adjusted
regression. A bound on EXPECTED overlap is not a bound on observed overlap,
nor a statistical causal test of interspecific competition.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

N=875
S=209
A=129
BOTH=5
OR_S=1.009
OR_A=0.979
LIGHT_SPAN=100.0

def logistic(t):
    if t>=0:
        return 1/(1+math.exp(-t))
    exp=math.exp(t)
    return exp/(1+exp)

def logit(p):
    if not 0<p<1:raise ValueError("probability outside open unit interval")
    return math.log(p/(1-p))

def one_slope_expected_overlap_lower_bound(n=N,s=S,a=A,s_odds=OR_S,span=LIGHT_SPAN):
    """Proof: min pS >= expit(logit(mean S)-max beta*x range).

    If pS(x)=logit^-1(alpha+beta*x), beta >=0 and span is bounded,
    pSmax >= meanS. Hence pSmin >= logit^-1(logit(meanS)-beta*span).
    For ANY x histogram and ANY pA(x) >=0 with sum pA=a,
    sum pS(x)pA(x) >= a * min pS.
    This bound is conservative, but global within its explicit model class.
    """
    if not (n>0 and 0<s<n and 0<a<n and s_odds>=1 and span>=0):
        raise ValueError("invalid model assumptions")
    pmin=logistic(logit(s/n)-math.log(s_odds)*span)
    return {"lower_bound_on_EXPECTED_both":a*pmin,
            "minimum_possible_S_probability_bound":pmin,
            "derivation":"S prevalence + monotone bounded logit slope yield pS(x)>=logistic(logit(S/N)-log(OR_S)*max_RLI_span) and E_both>=A*pmin",
            "restrictive_common_intercept_and_no_unmeasured_light_site_effects":True,
            "does_not_bound_observed_cooccurrence":True}

def calibrate_intercept(mean,beta,weight_sunny):
    """Bisection matches expected prevalence for 0/100 light split."""
    if not (0<=weight_sunny<=1):raise ValueError("invalid sunny weight")
    lo,hi=-40,40
    for _ in range(100):
        mid=(lo+hi)/2
        p=(1-weight_sunny)*logistic(mid)+weight_sunny*logistic(mid+beta*LIGHT_SPAN)
        if p<mean:lo=mid
        else:hi=mid
    return (lo+hi)/2

def two_extreme_light_mixture(weight_sunny,orr_s=OR_S,orr_a=OR_A):
    if orr_s<=0 or orr_a<=0:raise ValueError("invalid odds slope")
    beta_s,beta_a=math.log(orr_s),math.log(orr_a)
    alpha_s=calibrate_intercept(S/N,beta_s,weight_sunny)
    alpha_a=calibrate_intercept(A/N,beta_a,weight_sunny)
    p_s_shade=logistic(alpha_s)
    p_s_sun=logistic(alpha_s+beta_s*LIGHT_SPAN)
    p_a_shade=logistic(alpha_a)
    p_a_sun=logistic(alpha_a+beta_a*LIGHT_SPAN)
    expected=N*((1-weight_sunny)*p_s_shade*p_a_shade
                +weight_sunny*p_s_sun*p_a_sun)
    return {"fraction_at_RLI_100":weight_sunny,
            "expected_both_under_conditional_independence":expected,
            "p_S_at_RLI_0":p_s_shade,"p_S_at_RLI_100":p_s_sun,
            "p_A_at_RLI_0":p_a_shade,"p_A_at_RLI_100":p_a_sun,
            "expected_S":N*((1-weight_sunny)*p_s_shade+weight_sunny*p_s_sun),
            "expected_A":N*((1-weight_sunny)*p_a_shade+weight_sunny*p_a_sun)}

def original_toy_required_or():
    sunny={"N":400,"S":200,"A":5,"both":3}
    shady={"N":475,"S":9,"A":124,"both":2}
    pSs=sunny["S"]/sunny["N"]
    pSh=shady["S"]/shady["N"]
    pAs=sunny["A"]/sunny["N"]
    pAh=shady["A"]/shady["N"]
    return {"previous_invented_two_light_strata_not_measured_Korean_RLI":{
                "hypothetical_RLI_for_sunny":100,"hypothetical_RLI_for_shady":0,
                "implied_S_light_OR_per_one_percent":math.exp((logit(pSs)-logit(pSh))/100),
                "implied_A_light_OR_per_one_percent":math.exp((logit(pAs)-logit(pAh))/100),
                "previous_source_light_OR":{"S":OR_S,"A":OR_A},
                "reference_categories_reproduced_exactly":{
                    "neither":542,"S_only":204,"A_only":124,"both":5}},
            "conditional_model_cannot_equate_these_implied_ORs_to_published_adjusted_effects":True}

def site_intercept_noncompetition_counterexample():
    """Show why published light slopes do NOT identify site sorting.

    The two hypothetical groups are the PREVIOUS invented sun/shade classes,
    *reinterpreted as arbitrary source groups with their own species intercepts*.
    Set a common illustrative RLI=50 to compute site intercepts; published
    within-site light ORs remain OR_S and OR_A in both groups, but unobserved
    site intercept differences can reproduce the previous exact four margins.
    These are not Korean site data and must never be treated as fitted effects.
    """
    toy=[
      {"hypothetical_source_group":"group_1_NOT_KOREA","N":400,"S":200,"A":5,"both":3},
      {"hypothetical_source_group":"group_2_NOT_KOREA","N":475,"S":9,"A":124,"both":2}
    ]
    rows=[]
    for item in toy:
        n=item["N"]
        ps=item["S"]/n
        pa=item["A"]/n
        alpha_s=logit(ps)-math.log(OR_S)*50
        alpha_a=logit(pa)-math.log(OR_A)*50
        expected=n*ps*pa
        rows.append({**item,"illustrative_shared_RLI":50,
                     "group_intercept_for_reported_OR_S":alpha_s,
                     "group_intercept_for_reported_OR_A":alpha_a,
                     "expected_both_under_conditional_independence":expected})
    sum_expected=sum(row["expected_both_under_conditional_independence"] for row in rows)
    assert 4.8<sum_expected<5
    return {
       "kind":"HYPOTHETICAL_SITE_INTERCEPT_COUNTEREXAMPLE_NOT_FITTED_FIELD_MODEL",
       "group_model":"Within each group logistic(theta_species_group + log(reported OR species)*RLI), species independent conditional on group and RLI; group intercepts unrestricted.",
       "groups":rows,
       "expected_both_in_model":sum_expected,
       "printed_2026_original_margins_exactly_reproduced_in_hypothetical_realization":{
           "N":875,"S":209,"A":129,"both":5},
       "proof":"Allowing source-group-specific intercepts even with the SAME reported light coefficients can yield expected overlap≈5. Light slope coefficients plus pooled incidence margins do not restrict unobserved site intercepts.",
       "not_a_verified_Korean_site_distribution":True,
       "not_a_causal_competition_effect":True
    }


def evaluate():
    # Grid-restricted two-endpoint simulation. A numerical minimum here is NOT
    # a theorem that all arbitrary x histograms have their optimum at endpoints.
    samples=[two_extreme_light_mixture(i/2000) for i in range(2001)]
    minimum=min(samples,key=lambda row:row["expected_both_under_conditional_independence"])
    lower=one_slope_expected_overlap_lower_bound()
    if not (lower["lower_bound_on_EXPECTED_both"]<=minimum["expected_both_under_conditional_independence"]<=S*A/N):
        raise AssertionError("mathematical scenario inconsistent")
    return {
      "schema":"chocho_swallowtail_or_constrained_rli_encounter_v01",
      "source":"Jang et al 2026, DOI 10.5141/jee.26.017",
      "source_data_type":"PUBLISHED_RAMET_MARGINS_AND_PREVIOUSLY_TRANSCRIBED_ADJUSTED_RLI_ORS_NOT_RAW_OBSERVATIONS",
      "observed_ramets":{"N":N,"S":S,"A":A,"both":BOTH},
      "reported_adjusted_RLI_OR_as_transcribed":{"S":OR_S,"A":OR_A},
      "conditional_model":"single common intercept per butterfly, logistic light slope per 1pct over 0..100, independence conditional on light; no other covariate or site/random effects",
      "baseline_global_unstratified_expected_both":S*A/N,
      "previous_toy_OR_requirement":original_toy_required_or(),
      "rigorous_EXPECTED_overlap_lower_bound_under_common_intercept_assumption":lower,
      "two_RLI_endpoint_distribution_grid_0_to_100_numerical_minimum":minimum,
      "unrestricted_source_group_intercepts_counterexample":site_intercept_noncompetition_counterexample(),
      "two_endpoint_grid_resolution_sunny_fraction":1/2000,
      "observed_five_is_not_proof_that_expected_mean_is_five":True,
      "observed_five_can_be_sampling_variation_even_with_higher_mean":True,
      "reported_OR_model_level_cannot_be_directly_applied_to_875_ramet_rows_without_raw_data":True,
      "cannot_infer_negative_competition_or_true_light_response_as_unique_mechanism":True,
      "conclusion":"Old invented sunlight-only exact four-cell match required stronger pair-specific light gradients than authors' reported ORs if interpreted at the ramet level. A simple common-intercept reported-slope model cannot have expected both=5 with these margins, but this does not infer competition: other habitat axes, nonlinearities, site effects, source scale and sampling are not identified.",
      "GEB_PR38_untouched":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    result=evaluate()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
