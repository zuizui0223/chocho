#!/usr/bin/env python3
"""Exact overlap references for PRINTED 2026 butterfly co-use tables.

No raw field data are publicly available. Reference nulls exclude sunlight,
site×month confounding, nested ramets, sampling/detection bias and behavior.
Never call a rejection of random mixing evidence for interspecific competition.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from audit_2026_swallowtail_microhabitat_published_counts import PUBLICATION,audit

def hypergeom_pmf(N,S,A):
    """Exact overlap probability at fixed row marginal counts.

    S is number of Sericinus positive units; A is Atrophaneura positive units.
    X = number of both-positive units under independent random assignment of
    fixed incidence labels to N comparable sampling units.
    """
    if any(type(v)!=int for v in (N,S,A)) or not (0<=S<=N and 0<=A<=N):
        raise ValueError("bad fixed-margin input")
    low=max(0,S+A-N)
    high=min(S,A)
    denominator=math.comb(N,A)
    vals={k:math.comb(S,k)*math.comb(N-S,A-k)/denominator for k in range(low,high+1)}
    if abs(sum(vals.values())-1)>1e-10:
        raise ValueError("hypergeometric distribution normalization failed")
    return vals

def convolution(ps):
    total={0:1.0}
    for source in ps:
        next_= {}
        for k,p in total.items():
            for j,q in source.items():
                next_[k+j]=next_.get(k+j,0.0)+p*q
        total=next_
    if abs(sum(total.values())-1)>1e-9:
        raise ValueError("convolved exact law normalization failed")
    return total

def quantile(dist,p):
    if p<0 or p>1:raise ValueError("invalid quantile")
    cumulative=0.0
    for k,v in sorted(dist.items()):
        cumulative+=v
        if cumulative>=p-1e-12:return k
    return max(dist)

def summary(pmf,observed,source_groups):
    if observed<min(pmf) or observed>max(pmf):
        raise ValueError("observed co-use infeasible with margins")
    expected=sum(k*v for k,v in pmf.items())
    var=sum((k-expected)**2*v for k,v in pmf.items())
    return {
        "stratification":source_groups,
        "observed_both":observed,
        "expected_both_under_fixed_margins":expected,
        "observed_over_expected_ratio":observed/expected if expected else None,
        "percent_below_reference_expected":100*(1-observed/expected) if expected else None,
        "lower_tail_probability_P_X_at_most_observed":sum(v for k,v in pmf.items() if k<=observed),
        "central_95pct_null_overlap_reference_interval":[quantile(pmf,.025),quantile(pmf,.975)],
        "null_sd":math.sqrt(var),
        "no_competition_inference":True
    }

def stratified_quadrat_null(groups):
    parts=[]
    components={}
    sum_observed=0
    for label,g in sorted(groups.items()):
        N=sum(g.values())
        S=g["S_only"]+g["both"]
        A=g["A_only"]+g["both"]
        x=g["both"]
        parts.append(hypergeom_pmf(N,S,A))
        components[label]={"N_printed":N,"S_positive":S,"A_positive":A,
                           "both_printed":x,"expected_given_margins":S*A/N}
        sum_observed+=x
    result=summary(convolution(parts),sum_observed,list(sorted(groups)))
    result["group_margin_breakdown"]=components
    return result

def two_strata_habitat_counterexample():
    """Hypothetical 400 sunny / 475 shady ramets.

    Within each stratum species presence is statistically independent. Marginal
    occurrences S=209, A=129, N=875 and expected overlap≈5 arise solely from
    differing habitat rates. This is mathematical non-identifiability, NOT a
    reconstruction of the unshared Korean field observations.
    """
    rows=[
        {"hypothetical_microhabitat":"sunny","N":400,"S":200,"A":5},
        {"hypothetical_microhabitat":"shady","N":475,"S":9,"A":124},
    ]
    total={"N":sum(r["N"] for r in rows),
           "S":sum(r["S"] for r in rows),
           "A":sum(r["A"] for r in rows)}
    expected=sum(r["S"]*r["A"]/r["N"] for r in rows)
    return {
        "kind":"HYPOTHETICAL_CONSTRUCTIVE_COUNTEREXAMPLE_NOT_OBSERVED_DATA",
        "strata":[{**r,"expected_both_under_within_stratum_independence":r["S"]*r["A"]/r["N"]} for r in rows],
        "total_margins":total,
        "expected_overlap_from_environmental_sorting_alone":expected,
        "observed_paper_overlap":5,
        "proof_scope":"The printed global margins + overlap cannot alone distinguish competitive exclusion from environmental sorting. This hypothetical example is not estimated Korean light strata."
    }

def analyze(source=PUBLICATION):
    original_audit=audit(source)
    sites=source["table1_quadrat_by_site"]
    months=source["table1_quadrat_by_month"]
    ramets=source["table1_ramet_status"]
    ramet_n=sum(ramets.values())
    ramet_s=ramets["S_only"]+ramets["both"]
    ramet_a=ramets["A_only"]+ramets["both"]
    rs=summary(hypergeom_pmf(ramet_n,ramet_s,ramet_a),
               ramets["both"],["unstratified_ramets"])
    observed_site=sum(x["both"] for x in sites.values())
    observed_month=sum(x["both"] for x in months.values())
    if observed_site!=observed_month or original_audit["source_based_plausible_quadrat_denominator"]!=215:
        raise ValueError("printed quadrat source consistency invariant failed")
    scenario=two_strata_habitat_counterexample()
    if scenario["total_margins"]!={"N":ramet_n,"S":ramet_s,"A":ramet_a}:
        raise ValueError("hypothetical margin-alignment error")
    return {
        "schema":"chocho_swallowtail_2026_printed_overlap_exact_null_v01",
        "source_publication":"Jang et al 2026 J Ecol Environ DOI 10.5141/jee.26.017",
        "source_kind":"printed_published_aggregates_only_no_raw_ramet_or_quadrat_rows",
        "print_count_gate":{"methods_quadrats":source["declared_total_quadrats_in_methods"],
                            "breakdown_quadrats":original_audit["source_based_plausible_quadrat_denominator"],
                            "author_confirmation_required":True},
        "ramet_simple_reference":rs,
        "quadrat_site_margins_reference":stratified_quadrat_null(sites),
        "quadrat_month_margins_reference":stratified_quadrat_null(months),
        "sun_shade_counterexample":scenario,
        "site_and_month_reference_not_independent_tests":True,
        "within_site_month_light_and_host_size_not_adjustable_without_original_data":True,
        "new_ecological_or_causal_effect_estimated":False,
        "GEB_PR38_untouched":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    result=analyze()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
