#!/usr/bin/env python3
"""Fisher 2026 monarch: real Day-6 mass ID linked to observed pupation.

Observational individual prognosis *within recorded trial×leaf group*;
not proof of growth causality or a resolved 20241 treatment label.
20,000 conditional label permutations preserve each source stratum's n/events.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_COMMIT="acbb93c75e1f4810355c6ce03e6076e2303b2a98"
PREVIOUS=ROOT/"provenance/reviewer_defenses/results/monarch_leaf_age_pupation_stage_conversion_20261026.json"
BLOBS={
    "2324_GD_AllMass.csv":"3905f64bfe2667ce566a09d291e7c848378a11bf",
    "2324_GD_6DayMass.csv":"4f390c96c27d1635e01d6f0877bd81db84a05a0b",
    "2324_GD_PupalMass.csv":"bf2978fb1d7070b8916378813898f771d39e79d6",
    "SurviveToPupation.csv":"33f269e78f235cd1039394162a6871370db5e020",
}
BLOCKS=("20232","20241","20242")
LEAVES=("Top","Bottom")
PERMUTATIONS=20000
SEED=20261026


def source_file(repo,name):
    file=repo/name
    raw=file.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=BLOBS[name]:
        raise ValueError(f"Original experiment file source drift: {name}")
    with file.open(newline="",encoding="utf-8-sig") as fd:
        rows=list(csv.DictReader(fd))
    if not rows:
        raise ValueError("Empty original source file")
    return rows


def parse_masses(repo):
    both=source_file(repo,"2324_GD_AllMass.csv")
    stage6=source_file(repo,"2324_GD_6DayMass.csv")
    pupa=source_file(repo,"2324_GD_PupalMass.csv")
    source_file(repo,"SurviveToPupation.csv")
    prev=json.loads(PREVIOUS.read_text(encoding="utf-8"))
    if (prev["status"]!="SOURCE_TREATMENT_LABEL_INDETERMINATE_TWO_SCENARIOS_RETAINED"
        or len(prev["original_trial_label_disagreements"])!=2):
        raise ValueError("Original published treatment-label discrepancy not documented")

    stages={}
    for r in both:
        st=r["Stage"].strip()
        if st not in {"Day 6","Pupa"}:
            raise ValueError("Unrecognised original measured juvenile stage")
        key=(r["YearTrial"].strip(),r["Treatment"].strip(),r["ID"].strip())
        if key[0] not in set(BLOCKS)|{"20231"} or key[1] not in LEAVES or not key[2]:
            raise ValueError("Original trial, leaf source or larva ID missing")
        stagekey=(*key,st)
        if stagekey in stages:
            raise ValueError("Duplicate larva×measured stage ID in source")
        val=float(r["Mass"])
        if not math.isfinite(val) or val<=0:
            raise ValueError("Invalid positive originally measured larval/pupal mass")
        stages[stagekey]=val
    if len(stages)!=532:
        raise ValueError("Original 532 measured stage records changed")

    for name,expected_stage,measure_col,expected_n in [
        ("2324_GD_6DayMass.csv","Day 6","SixDayMass",270),
        ("2324_GD_PupalMass.csv","Pupa","PupalMass",262),
    ]:
        rows=stage6 if expected_stage=="Day 6" else pupa
        if len(rows)!=expected_n:
            raise ValueError("Original source individual measurement count changed")
        other={}
        for r in rows:
            k=(r["YearTrial"].strip(),r["Treatment"].strip(),r["ID"].strip(),expected_stage)
            if k in other:
                raise ValueError("Duplicate original individual/stage in "+name)
            other[k]=float(r[measure_col])
        indexed={k:v for k,v in stages.items() if k[3]==expected_stage}
        if set(indexed)!=set(other) or any(
            not math.isclose(v,other[k],rel_tol=0,abs_tol=1e-9) for k,v in indexed.items()
        ):
            raise ValueError("Independent source stage extracts do not reproduce original AllMass.csv")
    data6={k[:3]:v for k,v in stages.items() if k[3]=="Day 6"}
    data_pupa={k[:3]:v for k,v in stages.items() if k[3]=="Pupa"}
    if len(data6)!=270 or len(data_pupa)!=262:
        raise ValueError("Original stage count source drift")
    if {k[0] for k in data6}!=set(BLOCKS):
        raise ValueError("Only three named Day6 source trials are eligible")
    first_trial={k for k in data_pupa if k[0]=="20231"}
    matched={k for k in data6 if k in data_pupa}
    if (len(first_trial)!=72 or len(matched)!=190 or
        any(k not in data6 for k in data_pupa if k[0] in BLOCKS)):
        raise ValueError("Original pupal observations are not wholly linked to Day6 IDs in eligible source trials")
    return data6,data_pupa,{
        "AllMass_rows":len(stages),
        "Day6_individuals_measured":len(data6),
        "Pupa_individuals_measured":len(data_pupa),
        "Day6_and_later_pupation_ID_matched":len(matched),
        "Day6_measured_not_in_pupal_record":len(data6)-len(matched),
        "pupa_only_unmeasured_Day6_source_first_trial_20231":len(first_trial),
        "original_fate_treatment_count_discrepancy_20241_unresolved":True,
        "all_day6_eligible_trial_pupae_match_original_individual_IDs":True,
    }


def in_group_auc(values,success):
    n=len(values)
    pos=int(np.sum(success))
    neg=n-pos
    if pos<=0 or neg<=0:
        raise ValueError("No within-stratum positive-negative pupation comparison")
    ranks=rankdata(values,method="average")
    summed=float(np.sum(ranks[success]))
    auc=(summed-pos*(pos+1)/2)/(pos*neg)
    if not -1e-10<=auc<=1+1e-10:
        raise ValueError("Monarch within-stratum rank concordance invalid")
    return float(np.clip(auc,0,1)),ranks,pos,neg


def summarize(data6,pupae,swapped=False):
    strata=[]
    for block in BLOCKS:
        for label in LEAVES:
            original_label=("Bottom" if label=="Top" else "Top") if (
                swapped and block=="20241") else label
            vals=[(k[2],m,k in pupae) for k,m in data6.items()
                  if k[0]==block and k[1]==original_label]
            if not vals:
                raise ValueError("Missing original Day6 mass stratum")
            mass=np.array([x[1] for x in vals],dtype=float)
            success=np.array([x[2] for x in vals],dtype=bool)
            auc,ranks,npup,nfail=in_group_auc(mass,success)
            strata.append({
                "YearTrial":block,
                "recorded_label":label,
                "original_source_label_under_hypothesis":original_label,
                "n_Day6_measured":len(mass),
                "n_later_pupated":npup,
                "n_not_observed_pupating":nfail,
                "later_pupation_fraction_given_day6_measured":npup/len(mass),
                "Day6_mass_mean_later_pupated":float(np.mean(mass[success])),
                "Day6_mass_mean_without_later_pupal_record":float(np.mean(mass[~success])),
                "within_stratum_prognostic_mass_AUC":auc,
                "ranks":ranks,
                "n_positive":npup,
            })
    if sum(z["n_Day6_measured"] for z in strata)!=270 or sum(z["n_later_pupated"] for z in strata)!=190:
        raise ValueError("Original source cohort changed under stratification")
    return strata


def perm_test(strata):
    r=np.random.default_rng(SEED)
    values=np.array([x["within_stratum_prognostic_mass_AUC"] for x in strata])
    observed=float(np.mean(values))
    null=np.empty(PERMUTATIONS)
    for j in range(PERMUTATIONS):
        ranks=[]
        for z in strata:
            n=z["n_Day6_measured"]
            pos=z["n_positive"]
            indices=r.choice(n,size=pos,replace=False)
            a=(float(np.sum(z["ranks"][indices]))-pos*(pos+1)/2)/(pos*(n-pos))
            ranks.append(a)
        null[j]=np.mean(ranks)
    if not np.isfinite(null).all():
        raise ValueError("Nonfinite null ranks")
    q=np.quantile(null,[.025,.5,.975])
    obs_dev=abs(observed-.5)
    p=float((1+np.sum(np.abs(null-.5)>=obs_dev-1e-12))/(PERMUTATIONS+1))
    trial_values={
        block:float(np.mean([x["within_stratum_prognostic_mass_AUC"] for x in strata
                             if x["YearTrial"]==block])) for block in BLOCKS
    }
    loo={block:float(np.mean([x["within_stratum_prognostic_mass_AUC"] for x in strata
                             if x["YearTrial"]!=block])) for block in BLOCKS}
    return {
        "observed_equal_6_strata_Mann_Whitney_AUC":observed,
        "interpretation":"Pr(a randomly chosen later-pupating Day6-measured larva has higher Day6 mass than an otherwise same-trial same-recorded-leaf measured larva not recorded pupating); ties 0.5",
        "null":"randomize outcome identity within each YearTrial×recorded leaf keeping n and #pupae fixed",
        "permutation_draws":PERMUTATIONS,
        "seed":SEED,
        "null_q025_median_q975":q.tolist(),
        "two_sided_deviation_from_half_empirical_fraction":p,
        "trial_mean_AUC":trial_values,
        "leave_one_trial_out_AUC":loo,
        "source_selection":"Only individuals surviving and weighed on Day6; no causal inference for assigned 480 larvae",
    }


def run(source):
    day6,pupae,gate=parse_masses(source)
    strata=summarize(day6,pupae)
    sensitivity=summarize(day6,pupae,swapped=True)
    original_value=float(np.mean([x["within_stratum_prognostic_mass_AUC"] for x in strata]))
    alternative=float(np.mean([x["within_stratum_prognostic_mass_AUC"] for x in sensitivity]))
    if not math.isclose(original_value,alternative,rel_tol=0,abs_tol=1e-13):
        raise ValueError("Within-block leaf group rename must not change observed rank AUC")
    inference=perm_test(strata)
    export=[]
    for z in strata:
        export.append({k:v for k,v in z.items() if k not in {"ranks","n_positive"}})
    return {
        "schema":"chocho_individual_monarch_day6_mass_pupation_rank_AUC_v01_20261026",
        "status":"ORIGINAL_INDIVIDUAL_STAGE_LINK_CONFIRMED_SOURCE_CONDITIONAL_PROGNOSIS_ONLY",
        "author_paper_DOI":"10.3390/insects17020215",
        "author_original_Git_commit":ORIGINAL_COMMIT,
        "exact_input_source_Git_blobs":BLOBS,
        "individual_stage_link_quality":gate,
        "observed_six_original_strata":export,
        "fixed_prognostic_analysis":inference,
        "yearTrial20241_entire_recorded_TopBottom_group_rename_AUC_invariant":True,
        "source_label_swap_AUC_change":alternative-original_value,
        "no_larval_individual_ID_or_source_raw_value_persisted":True,
        "limitations":[
            "Posthoc original experiment reanalysis, not new butterfly experiment or across-species ecological law",
            "Rank AUC is prognostic observational relationship conditional on surviving to measured Day6",
            "One of four original survival trials lacks Day6 measurements; sample is 270/480 assigned",
            "The 80 Day6-measured animals without Pupa records are not assigned causes of loss or adult fate",
            "YearTrial20241 original Top/Bottom labels conflict in treatment success table; biological leaf-treatment effect remains unresolved",
            "A high AUC does not establish a milkweed age mediator or causal mass-survival intervention",
            "No local introduced host arrival, native/alien comparison, host competition or field butterfly population increase",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=run(args.source_repo)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    a=result["fixed_prognostic_analysis"]
    print(json.dumps({
        "status":result["status"],
        "fate_link":result["individual_stage_link_quality"],
        "auc":a["observed_equal_6_strata_Mann_Whitney_AUC"],
        "null":a["null_q025_median_q975"],
        "two_sided_tail":a["two_sided_deviation_from_half_empirical_fraction"],
        "individual_strata_AUC":[{
            "block":z["YearTrial"],"leaf_source_label":z["recorded_label"],
            "n":z["n_Day6_measured"],"pupae":z["n_later_pupated"],
            "AUC":z["within_stratum_prognostic_mass_AUC"]
        } for z in result["observed_six_original_strata"]],
        "source_label_swap_unchanged":result["yearTrial20241_entire_recorded_TopBottom_group_rename_AUC_invariant"],
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
