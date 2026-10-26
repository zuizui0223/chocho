#!/usr/bin/env python3
"""Audit Fisher 2026 published survival-proportion test against raw source.

Original R code applies chisq.test(c(48.8,60.8)); this is a goodness-of-fit
on two percentages, NOT a survival comparison retaining 240 larvae per arm.
Two contradictory 20241 source-label scenarios are held open. Original authors'
2026 Insects experiment, not a novel new biological result or misconduct claim.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import chisquare, chi2_contingency, fisher_exact

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"provenance/reviewer_defenses/results/monarch_leaf_age_pupation_stage_conversion_20261026.json"
SOURCE_TREE="acbb93c75e1f4810355c6ce03e6076e2303b2a98"
PAPER_DOI="10.3390/insects17020215"
EXPECTED_BLOBS={
    "121125_Fisher2025_MonarchMilkweedNvO_AnalysisCode.R":"de410164236730994c559bcf115a1e64e187be18",
    "SurviveToPupation.csv":"33f269e78f235cd1039394162a6871370db5e020",
    "2324_GD_PupalDays.csv":"7066610f45f0baace4cfb679c5b7ee3355ce881f",
    "2324_GD_PupalMass.csv":"bf2978fb1d7070b8916378813898f771d39e79d6",
}
PUBLISHED={
    "day6_top_percent":75.4,
    "day6_bottom_percent":74.2,
    "day6_chisq":0.0096257,
    "day6_p":0.9218,
    "pupation_top_percent":60.8,
    "pupation_bottom_percent":48.8,
    "pupation_chisq":1.3139,
    "pupation_p":0.2517,
}


def gitsha(raw:bytes)->str:
    return hashlib.sha1(
        b"blob "+str(len(raw)).encode("ascii")+b"\0"+raw
    ).hexdigest()


def original_sources(path):
    hashes={}
    for name,sha in EXPECTED_BLOBS.items():
        raw=(path/name).read_bytes()
        computed=gitsha(raw)
        if computed!=sha:
            raise ValueError(f"Original source differs: {name}: {computed}")
        hashes[name]=sha
    code=(path/"121125_Fisher2025_MonarchMilkweedNvO_AnalysisCode.R").read_text(
        encoding="utf-8-sig")
    literal_six="chisq.test(c(74.2,75.4))"
    literal_pupa="chisq.test(c(48.8,60.8))"
    if code.count(literal_six)!=1 or code.count(literal_pupa)!=1:
        raise ValueError("Published-author survival goodness-of-fit on percentages source lines changed")
    with (path/"SurviveToPupation.csv").open(newline="",encoding="utf-8-sig") as f:
        data=list(csv.DictReader(f))
    if len(data)!=8:
        raise ValueError("Original 8 summary source rows missing")
    summary={}
    for d in data:
        key=(d["YearTrial"],d["Treatment"])
        if key in summary:
            raise ValueError("Duplicate original survival cell")
        summary[key]=(int(d["#Survive"]),int(d["Total"]))
    if any(n!=60 for a,n in summary.values()):
        raise ValueError("Original 60 assigned larvae per treatment/trial changed")
    return hashes,summary,{"day6":literal_six,"pupation":literal_pupa}


def check_original_receipt():
    z=json.loads(RECEIPT.read_text(encoding="utf-8"))
    if (
        z["status"]!="SOURCE_TREATMENT_LABEL_INDETERMINATE_TWO_SCENARIOS_RETAINED"
        or z["source_assigned_larvae"]!=480
        or z["reported_pupation_successes_across_blocks"]!=262
        or len(z["original_trial_label_disagreements"])!=2
    ):
        raise ValueError("Original four-block treatment contradiction qualification missing")
    return z


def decimal_percent(n,d):
    return 100.0*n/d


def compare_case(obj,name):
    a=int(obj["top_survived"])
    n1=int(obj["top_assigned"])
    c=int(obj["bottom_survived"])
    n0=int(obj["bottom_assigned"])
    if n1!=240 or n0!=240 or a+c!=262:
        raise ValueError("An original source trust sensitivity changed assigned animals")
    table=np.array([[a,n1-a],[c,n0-c]],dtype=int)
    chi2,chi_p,dof,expect=chi2_contingency(table,correction=False)
    yf,yp,dof_y,expect_y=chi2_contingency(table,correction=True)
    fe=fisher_exact(table,alternative="two-sided")
    estimated=float(obj["MH_common_OR_top_vs_bottom"])
    exact=float(obj["conditional_fixed_margin_exact"]["two_sided_probability_order_p"])
    ci=obj["trial_cluster_bootstrap"]["common_odds_ratio_bootstrap_95"]
    return {
        "scenario":name,
        "source_observations":{"Top_pupated":a,"Top_assigned":n1,
                               "Bottom_pupated":c,"Bottom_assigned":n0},
        "Top_rate_percent":decimal_percent(a,n1),
        "Bottom_rate_percent":decimal_percent(c,n0),
        "difference_percentage_points":decimal_percent(a,n1)-decimal_percent(c,n0),
        "unstratified_fisher_exact_two_sided_p":float(fe.pvalue),
        "unstratified_fisher_exact_OR":float(fe.statistic),
        "unstratified_2x2_chi2_no_Yates_p":float(chi_p),
        "unstratified_2x2_chi2_no_Yates_statistic":float(chi2),
        "unstratified_2x2_chi2_with_Yates_p":float(yp),
        "trial_stratified_Mantel_Haenszel_OR":estimated,
        "trial_stratified_four_block_fixed_margin_exact_p":exact,
        "trial_cluster_four_block_bootstrap_95pct_OR":ci,
        "source_integrity_unresolved":True,
        "pooled_tests_ignore_trial_cluster_correlation":True,
    }


def run(original):
    hashes,source_counts,literals=original_sources(original)
    z=check_original_receipt()
    data=z["experimental_survival_scenarios"]
    s0=compare_case(data["S0_original_survival_summary"],"S0_original_survival_summary")
    s1=compare_case(data["S1_pupal_individual_count_hypothesis"],"S1_hypothetical_stage_records")
    expected_20241={
        ("20241","Top"):(30,60),
        ("20241","Bottom"):(36,60),
    }
    for k,x in expected_20241.items():
        if source_counts.get(k)!=x:
            raise ValueError("Original 20241 treatment discrepancy source changed")
    goodness_6=chisquare([
        PUBLISHED["day6_bottom_percent"],PUBLISHED["day6_top_percent"]
    ])
    goodness_p=chisquare([
        PUBLISHED["pupation_bottom_percent"],PUBLISHED["pupation_top_percent"]
    ])
    if not (math.isclose(float(goodness_6.statistic),PUBLISHED["day6_chisq"],abs_tol=1e-5)
            and math.isclose(float(goodness_6.pvalue),PUBLISHED["day6_p"],abs_tol=2e-4)
            and math.isclose(float(goodness_p.statistic),PUBLISHED["pupation_chisq"],abs_tol=1e-4)
            and math.isclose(float(goodness_p.pvalue),PUBLISHED["pupation_p"],abs_tol=2e-4)):
        raise ValueError("Literal paper R percentage-as-count statistics do not reproduce published results")
    if not (math.isclose(s1["Top_rate_percent"],60.8,abs_tol=.05)
            and math.isclose(s1["Bottom_rate_percent"],48.3333333,abs_tol=.01)
            and math.isclose(s0["Top_rate_percent"],58.33333,abs_tol=.01)):
        raise ValueError("Original independently reconstructed group proportions changed")

    report={
        "schema":"chocho_fisher2026_published_percentage_chisq_and_source_discrepancy_audit_v01",
        "status":"ORIGINAL_PUBLISHED_TEST_DENOMINATOR_MISMATCH_AND_SOURCE_LABEL_DISCORDANCE_VERIFIED",
        "original_article_doi":PAPER_DOI,
        "article_section":"2.6 proportion-test Methods and 3.2 larval survival Results",
        "article_reported":PUBLISHED,
        "original_author_code_git_commit":SOURCE_TREE,
        "source_git_blob_sha1":hashes,
        "source_literal_R_test_calls":literals,
        "recomputed_on_two_published_percentages":{
            "six_day_chisquare_goodness_of_fit":float(goodness_6.statistic),
            "six_day_chisquare_goodness_of_fit_p":float(goodness_6.pvalue),
            "pupation_chisquare_goodness_of_fit":float(goodness_p.statistic),
            "pupation_chisquare_goodness_of_fit_p":float(goodness_p.pvalue),
            "statistical_warning":"chisq.test(c(48.8,60.8)) uses percentages as two observed frequency values; it omits 240 assigned larvae per treatment and does not test the 2×2 survival-outcome table.",
        },
        "published_percent_and_raw_count_denominator_check":{
            "assigned_per_leaf_treatment":240,
            "published_top_60_8_most_plausible_integer_count":146,
            "published_bottom_48_8_rounded_integer_count_for_240":117,
            "stage_data_bottom_positive_actual":116,
            "stage_data_bottom_percentage":116/240*100,
            "summary_data_top_positive_actual":140,
            "summary_data_bottom_positive_actual":122,
            "published_bottom_percentage_is_neither_original_summary_nor_individual_record_exact_percentage":True,
        },
        "unresolved_20241_original_source_treatment_counts":{
            "summary_top":30,"summary_bottom":36,
            "pupal_days_top":36,"pupal_days_bottom":30,
            "pupal_mass_top":36,"pupal_mass_bottom":30,
        },
        "two_scenarios_no_trust_selection":[s0,s1],
        "scientific_claim":"No unique hypothesis-testing conclusion about pupation success is justified until original YearTrial20241 assignment and published percentage discrepancy are resolved by study authors.",
        "not_implied":[
            "Authorship wrongdoing, misconduct or intentional misreporting",
            "A new causal milkweed phenology mechanism discovered by this source audit",
            "Adult eclosion, field population growth, host plant introduction benefit, competition or global resource conversion rate",
        ],
        "author_clarification_questions":[
            "Was YearTrial20241 Top=30 and Bottom=36 in SurviveToPupation.csv, or Top=36 and Bottom=30 in both 262-row pupal files, the correct treatment-specific pupation record?",
            "What underlying numerator and denominator yield the published 48.8% Bottom survival to pupation (versus 116/240=48.3% recorded in the two pupal files)?",
            "Should the chi-square analysis on hard-coded percentages be replaced by one using 4×2 original pupation successes/failures, with explicit trial structure?",
        ],
        "source_original_animal_ID_rows_not_republished":True,
    }
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    x=p.parse_args()
    r=run(x.source_repo)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":r["status"],
        "published_p":r["article_reported"]["pupation_p"],
        "reproduced_p":r["recomputed_on_two_published_percentages"]["pupation_chisquare_goodness_of_fit_p"],
        "scenarios":[{
            "n":z["scenario"],"survival_top":z["Top_rate_percent"],
            "survival_bottom":z["Bottom_rate_percent"],
            "fisher_p":z["unstratified_fisher_exact_two_sided_p"],
            "trial_conditional_exact_p":z["trial_stratified_four_block_fixed_margin_exact_p"]
        } for z in r["two_scenarios_no_trust_selection"]],
        "published_count_disagreement":r["published_percent_and_raw_count_denominator_check"],
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
