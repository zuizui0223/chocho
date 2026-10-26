#!/usr/bin/env python3
"""Original experimentally measured monarch leaf-age stage conversion audit.

Source Fisher et al. (2026), DOI 10.3390/insects17020215; GitHub source pinned
at acbb93c75e1f4810355c6ce03e6076e2303b2a98. Existing paper's individual
source evidence, not a new host relocation, wild competition or demographic
mechanism claim.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

TREE="acbb93c75e1f4810355c6ce03e6076e2303b2a98"
DOI="10.3390/insects17020215"
BLOBS={
    "SurviveToPupation.csv":"33f269e78f235cd1039394162a6871370db5e020",
    "2324_GD_6DayMass.csv":"4f390c96c27d1635e01d6f0877bd81db84a05a0b",
    "2324_GD_PupalDays.csv":"7066610f45f0baace4cfb679c5b7ee3355ce881f",
    "2324_GD_PupalMass.csv":"bf2978fb1d7070b8916378813898f771d39e79d6",
    "LeafCharacteristics.csv":"745bed22c3f3ef38736ed33effbbb7453ec154e6",
    "PercentNitrogen.csv":"8207f3d8e2acc8845c359b2f999ec697f02171b2",
}
BLOCKS=("20231","20232","20241","20242")
RNG_SEED=20261026
DRAWS=10000


def read_original(root,name):
    path=root/name
    raw=path.read_bytes()
    git_sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()
    if git_sha!=BLOBS[name]:
        raise ValueError(f"Source SHA drifted for {name}: {git_sha}")
    with path.open(newline="",encoding="utf-8-sig") as f:
        records=list(csv.DictReader(f))
    if not records:
        raise ValueError("Original source file missing rows "+name)
    return records


def source_survival(rows):
    table={}
    for r in rows:
        block=str(r["YearTrial"]).strip()
        tr=str(r["Treatment"]).strip()
        if block not in BLOCKS or tr not in {"Top","Bottom"} or (block,tr) in table:
            raise ValueError("Unexpected/duplicated trial by leaf-age assignment")
        a=int(r["#Survive"])
        n=int(r["Total"])
        if n!=60 or a<0 or a>n:
            raise ValueError("Original sample assignment n=60 changed")
        expected=round(100*a/n,7)
        if abs(float(r["%Survive"])-expected)>1e-5:
            raise ValueError("Original reported survival fraction inconsistent")
        year=int(r["Year"])
        trial=int(r["Trial"])
        if block!=str(year)+str(trial):
            raise ValueError("Original YearTrial labels inconsistent")
        table[(block,tr)]={"survive":a,"total":n}
    if len(rows)!=8 or len(table)!=8:
        raise ValueError("Expected 4 matched blocks × 2 treatments")
    if sum(x["total"] for x in table.values())!=480:
        raise ValueError("Not original 480 assigned larvae")
    return table


def mh_or(blocks):
    numerator=denom=0.0
    for a,b,c,d in blocks:
        n=a+b+c+d
        numerator+=a*d/n
        denom+=b*c/n
    return numerator/denom if denom>0 else None


def two_sided_fixed_margin_hypergeometric(blocks):
    # A_i = top pupation successes conditional on each experimental
    # YearTrial total pupated and the original 60 top / 60 bottom assignment.
    pmf=np.array([1.0],dtype=float)
    minsum=0
    for a,b,c,d in blocks:
        M=a+b+c+d
        K=a+c
        top_n=a+b
        low=max(0,top_n-(M-K))
        high=min(top_n,K)
        xs=np.arange(low,high+1)
        probs=hypergeom.pmf(xs,M,K,top_n)
        if not np.isclose(probs.sum(),1.0,atol=1e-10):
            raise ValueError("Hypergeometric trial null distribution does not sum to 1")
        pmf=np.convolve(pmf,probs)
        minsum+=low
    observed=int(sum(a for a,b,c,d in blocks))
    point=observed-minsum
    if point<0 or point>=len(pmf):
        raise ValueError("Observed top survivorship outside fixed-margin support")
    obsprob=pmf[point]
    two=float(np.sum(pmf[pmf<=obsprob+1e-12]))
    high=float(np.sum(pmf[point:]))
    low=float(np.sum(pmf[:point+1]))
    return {
        "test":"Four-block fixed-margin exact hypergeometric convolution",
        "conditional_observed_top_pupae":observed,
        "conditional_tail_top_high":high,
        "conditional_tail_top_low":low,
        "two_sided_probability_order_p":min(1.0,two),
        "max_supported_top_counts":minsum+len(pmf)-1,
        "min_supported_top_counts":minsum,
        "limitation":"Conditional treatment randomization does NOT estimate between-trial population variance",
    }


def robust_trial_bootstrap(block_data,draws=DRAWS):
    rng=np.random.default_rng(RNG_SEED)
    log_or=[]
    risk_deltas=[]
    for _ in range(draws):
        k=rng.integers(0,len(block_data),size=len(block_data))
        chosen=[block_data[i] for i in k]
        pooled=mh_or(chosen)
        if pooled is None or not math.isfinite(pooled) or pooled<=0:
            continue
        log_or.append(math.log(pooled))
        risk_deltas.append(float(np.mean([(a/(a+b))-(c/(c+d)) for a,b,c,d in chosen])))
    if len(log_or)!=draws:
        raise ValueError("Unexpected bootstrap degeneracy in four trial blocks")
    logci=np.quantile(log_or,[.025,.5,.975],method="linear")
    diffci=np.quantile(risk_deltas,[.025,.5,.975],method="linear")
    return {
        "trial_clusters":len(block_data),
        "draws":draws,"successful":len(log_or),"seed":RNG_SEED,
        "common_odds_ratio_bootstrap_95":[float(math.exp(logci[0])),float(math.exp(logci[2]))],
        "median_OR_bootstrap":float(math.exp(logci[1])),
        "paired_risk_difference_trial_bootstrap_95":[float(diffci[0]),float(diffci[2])],
        "paired_risk_difference_bootstrap_median":float(diffci[1]),
        "bootstrap_top_better_fraction":sum(v>0 for v in risk_deltas)/len(risk_deltas),
        "qualification":"Only 4 real experimental blocks, not general field-year/site replicates",
    }


def numeric_stats(records,name,column,labels=("Top","Bottom")):
    required={"YearTrial","Treatment",column}
    if not required.issubset(records[0]):
        raise ValueError(f"{name} missing {sorted(required-set(records[0]))}")
    out={}
    individs=set()
    for row in records:
        block=str(row["YearTrial"]).strip()
        leaf=str(row["Treatment"]).strip()
        if block not in BLOCKS or leaf not in labels:
            raise ValueError(f"Unexpected trial or leaf treatment in {name}")
        if "ID" in row:
            token=(block,leaf,str(row["ID"]).strip())
            if token in individs:
                raise ValueError("Duplicate individual ID in "+name)
            individs.add(token)
        try:
            v=float(row[column])
        except (ValueError,TypeError):
            continue
        if not math.isfinite(v) or v<0:
            raise ValueError("Nonpositive/outlier source should be audited; no imputation in "+name)
        out.setdefault((block,leaf),[]).append(v)
    stats={}
    for block in BLOCKS:
        stats[block]={}
        for leaf in labels:
            arr=out.get((block,leaf),[])
            stats[block][leaf]={
                "n_measured":len(arr),
                "mean":float(np.mean(arr)) if arr else None,
                "median":float(np.median(arr)) if arr else None,
            }
    pairmeans=[]
    for block in BLOCKS:
        top=stats[block]["Top"]["mean"]
        bottom=stats[block]["Bottom"]["mean"]
        if top is not None and bottom is not None:
            pairmeans.append(top-bottom)
    return {
        "file":name,
        "original_rows":len(records),
        "value":column,
        "per_YearTrial":stats,
        "available_blocks_for_pair_comparison":len(pairmeans),
        "equal_block_mean_Top_minus_Bottom_conditional_on_measured":float(np.mean(pairmeans)) if pairmeans else None,
        "not_recruitment_or_survival_estimate":True,
    }


def leaf_source_stats(root):
    f=read_original(root,"LeafCharacteristics.csv")
    g=read_original(root,"PercentNitrogen.csv")
    out={}
    for name,rows,label,fields in [
        ("LeafCharacteristics.csv",f,"Leaf",("Force","Thickness","Length","Width")),
        ("PercentNitrogen.csv",g,"Leaf",("NitrogenPerc",)),
    ]:
        groups=defaultdict(list)
        for r in rows:
            leaf=str(r[label]).strip()
            if leaf not in {"Top","Bottom"}:
                raise ValueError("Unexpected plant leaf position/age label in "+name)
            groups[leaf].append(r)
        summary={}
        for leaf,recs in groups.items():
            stats={}
            for key in fields:
                vals=[float(r[key]) for r in recs if r.get(key) and r[key]!="NA"]
                if len(vals)>0 and not all(math.isfinite(v) for v in vals):
                    raise ValueError("Invalid original plant trait value")
                stats[key]={"n":len(vals),"mean":float(np.mean(vals)) if vals else None}
            summary[leaf]={"rows":len(recs),"traits":stats}
        out[name]={"original_rows":len(rows),"groups":summary,
                   "no_individual_caterpillar_linkage":True}
    return out


def survival_scenario(table, source_label):
    blocks=[]
    pblock={}
    for block in BLOCKS:
        top=table[(block,"Top")]
        bottom=table[(block,"Bottom")]
        a=top["survive"]
        bb=top["total"]-a
        c=bottom["survive"]
        d=bottom["total"]-c
        blocks.append((a,bb,c,d))
        pblock[block]={
            "top_pupation":a,"top_denominator":a+bb,
            "bottom_pupation":c,"bottom_denominator":c+d,
            "Top_minus_Bottom_probability":a/(a+bb)-c/(c+d),
            "within_block_OR_top_bottom":(a*d/(bb*c)) if bb*c else None,
        }
    common=mh_or(blocks)
    if common is None:
        raise ValueError("Four source trial margins fail common odds computation")
    loo=[mh_or(blocks[:i]+blocks[i+1:]) for i in range(4)]
    a=sum(s[0] for s in blocks)
    b=sum(s[1] for s in blocks)
    c=sum(s[2] for s in blocks)
    d=sum(s[3] for s in blocks)
    return {
        "scenario":source_label,
        "top_survived":a,"top_assigned":a+b,
        "bottom_survived":c,"bottom_assigned":c+d,
        "top_pupation_rate":a/(a+b),
        "bottom_pupation_rate":c/(c+d),
        "pooled_risk_difference":a/(a+b)-c/(c+d),
        "MH_common_OR_top_vs_bottom":common,
        "within_4_trial_blocks":pblock,
        "leave_one_trial_out_MH_OR_range":[float(min(loo)),float(max(loo))],
        "conditional_fixed_margin_exact":two_sided_fixed_margin_hypergeometric(blocks),
        "trial_cluster_bootstrap":robust_trial_bootstrap(blocks),
        "source_status":"ORIGINAL_SUMMARY_AS_RECORDED" if source_label=="S0_survival_table_first"
                        else "HYPOTHETICAL_PUPAL_RECORD_COUNT_ALTERNATIVE_NOT_AUTHOR_CORRECTION",
    }


def run(root):
    raw=read_original(root,"SurviveToPupation.csv")
    original_table=source_survival(raw)
    days=numeric_stats(read_original(root,"2324_GD_PupalDays.csv"),
                       "2324_GD_PupalDays.csv","PupaDays")
    mass=numeric_stats(read_original(root,"2324_GD_PupalMass.csv"),
                       "2324_GD_PupalMass.csv","PupalMass")
    early=numeric_stats(read_original(root,"2324_GD_6DayMass.csv"),
                        "2324_GD_6DayMass.csv","SixDayMass")
    source_success_total=sum(r["survive"] for r in original_table.values())
    if days["original_rows"]!=source_success_total or mass["original_rows"]!=source_success_total:
        raise ValueError("Overall pupal stage data mismatch original summary total")
    disagreements=[]
    alt={k:dict(v) for k,v in original_table.items()}
    for block in BLOCKS:
        for leaf in ("Top","Bottom"):
            summary=original_table[(block,leaf)]["survive"]
            day_n=days["per_YearTrial"][block][leaf]["n_measured"]
            mass_n=mass["per_YearTrial"][block][leaf]["n_measured"]
            if day_n!=mass_n:
                raise ValueError("PupalDays and PupalMass source records disagree with one another")
            if summary!=day_n:
                disagreements.append({
                    "YearTrial":block,"original_leaf_label":leaf,
                    "summary_survived":summary,
                    "pupal_days_individuals":day_n,
                    "pupal_mass_individuals":mass_n,
                })
                alt[(block,leaf)]["survive"]=day_n
    if (len(disagreements)!=2
        or {(d["YearTrial"],d["original_leaf_label"]) for d in disagreements} !=
        {("20241","Top"),("20241","Bottom")}):
        raise ValueError("Source label contradictions differ from strictly documented 20241 swap")
    if sum(v["survive"] for v in alt.values())!=source_success_total:
        raise ValueError("Alternative incorrectly changes total number of pupae")
    s0=survival_scenario(original_table,"S0_survival_table_first")
    s1=survival_scenario(alt,"S1_pupal_day_and_mass_record_first_hypothesis")
    return {
        "schema":"chocho_monarch_leaf_age_stage_conversion_discrepant_treatments_v02_20261026",
        "source_doi":DOI,
        "original_author_git_commit":TREE,
        "exact_source_git_blobs":BLOBS,
        "status":"SOURCE_TREATMENT_LABEL_INDETERMINATE_TWO_SCENARIOS_RETAINED",
        "source_assigned_larvae":480,
        "source_experimental_blocks":len(BLOCKS),
        "reported_pupation_successes_across_blocks":source_success_total,
        "survival_outcome":"successful_pupation_not_adult_emergence",
        "original_trial_label_disagreements":disagreements,
        "contradiction_provenance":"Only 2024 trial1 Top/Bottom source labels inverted between summary and BOTH PupalDays/PupalMass; no source edited",
        "experimental_survival_scenarios":{
            "S0_original_survival_summary":s0,
            "S1_pupal_individual_count_hypothesis":s1,
        },
        "secondary_measured_stage_outcomes":{
            "Day6_mass":early,
            "pupal_mass_among_successes":mass,
            "pupal_days_among_successes":days,
        },
        "secondary_stage_20241_treatment_label_ambiguous":True,
        "plant_leaf_phenotypes":leaf_source_stats(root),
        "scientific_status":"POSTHOC_REANALYSIS_SOURCE_PROVENANCE_ERROR_NOT_NEW_RECRUITMENT_MECHANISM",
        "limitations":[
            "Published source label conflict remains unresolved pending author correction",
            "Original experimental work belongs to Fisher et al. 2026; do not claim novel treatment biology",
            "Four blocks, conditional exact test and bootstrap identify different uncertainties",
            "Asclepias syriaca leaf age is not host plant relocation or alien-native status",
            "Pupation is not adult eclosion; cannot infer new butterfly population recruitment",
            "Growth or development values may be selection-conditioned stage outcomes",
            "No matched larvae/habitat to global 239 butterfly resource map",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    r=run(a.source_repo)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":r["status"],
        "disagreements":r["original_trial_label_disagreements"],
        "assigned":r["source_assigned_larvae"],
        "scenario_results":{k:{
            "top_rate":v["top_pupation_rate"],
            "bottom_rate":v["bottom_pupation_rate"],
            "MH_OR":v["MH_common_OR_top_vs_bottom"],
            "bootstrap":v["trial_cluster_bootstrap"],
            "exact_p":v["conditional_fixed_margin_exact"]["two_sided_probability_order_p"],
        } for k,v in r["experimental_survival_scenarios"].items()},
        "stages":{k:{"n":v["original_rows"],
            "measured_Top_minus_Bottom":v["equal_block_mean_Top_minus_Bottom_conditional_on_measured"]}
            for k,v in r["secondary_measured_stage_outcomes"].items()},
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
