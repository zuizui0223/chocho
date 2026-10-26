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


def run(root):
    raw=read_original(root,"SurviveToPupation.csv")
    table=source_survival(raw)
    blocks=[]
    pblock={}
    for block in BLOCKS:
        t=table[(block,"Top")]
        b=table[(block,"Bottom")]
        a=t["survive"]; bb=t["total"]-a; c=b["survive"]; d=b["total"]-c
        blocks.append((a,bb,c,d))
        pblock[block]={
            "top_pupation":a,"top_denominator":a+bb,
            "bottom_pupation":c,"bottom_denominator":c+d,
            "Top_minus_Bottom_probability":a/(a+bb)-c/(c+d),
            "within_block_OR_top_bottom":(a*d/(bb*c)) if bb*c>0 else None,
        }
    odds=mh_or(blocks)
    if odds is None:
        raise ValueError("Four source blocks have no estimable common odds")
    loo=[mh_or(blocks[:i]+blocks[i+1:]) for i in range(4)]
    a=sum(s[0] for s in blocks)
    b=sum(s[1] for s in blocks)
    c=sum(s[2] for s in blocks)
    d=sum(s[3] for s in blocks)
    days=numeric_stats(read_original(root,"2324_GD_PupalDays.csv"),
                       "2324_GD_PupalDays.csv","PupaDays")
    mass=numeric_stats(read_original(root,"2324_GD_PupalMass.csv"),
                       "2324_GD_PupalMass.csv","PupalMass")
    first_mass=numeric_stats(read_original(root,"2324_GD_6DayMass.csv"),
                             "2324_GD_6DayMass.csv","SixDayMass")
    if days["original_rows"] != a+c or mass["original_rows"]!=a+c:
        raise ValueError("Source success stages counts do not equal original pupation counts")
    for block in BLOCKS:
        for leaf in ("Top","Bottom"):
            survived=table[(block,leaf)]["survive"]
            if days["per_YearTrial"][block][leaf]["n_measured"]!=survived:
                raise ValueError("PupalDays block treatment n differs from source pupation survivors")
            if mass["per_YearTrial"][block][leaf]["n_measured"]!=survived:
                raise ValueError("PupalMass block treatment n differs from source pupation survivors")
    return {
        "schema":"chocho_monarch_leaf_age_trial_block_pupation_stage_conversion_v01_20261026",
        "source_doi":DOI,
        "original_author_git_commit":TREE,
        "exact_source_git_blobs":BLOBS,
        "scientific_status":"POSTHOC_ORIGINAL_EXPERIMENT_CALIBRATION_NOT_NOVEL_FIELD_BUTTERFLY_RECRUITMENT",
        "source_assigned_larvae":a+b+c+d,
        "source_experimental_blocks":len(BLOCKS),
        "survival_outcome":"successful_pupation_not_adult_emergence",
        "top_total_survival":a,"top_total_assigned":a+b,
        "bottom_total_survival":c,"bottom_total_assigned":c+d,
        "top_pupation_rate":a/(a+b),
        "bottom_pupation_rate":c/(c+d),
        "pooled_top_minus_bottom_risk_difference":a/(a+b)-c/(c+d),
        "primary_Mantel_Haenszel_trial_stratified_pupation_OR_top_vs_bottom":odds,
        "within_4_real_trial_blocks":pblock,
        "leave_one_YearTrial_out_OR_range":[float(min(loo)),float(max(loo))],
        "exact_trial_conditional_test":two_sided_fixed_margin_hypergeometric(blocks),
        "YearTrial_cluster_uncertainty":robust_trial_bootstrap(blocks),
        "secondary_measured_stage_outcomes":{
            "Day6_mass":first_mass,
            "pupal_mass_among_successes":mass,
            "pupal_days_among_successes":days,
        },
        "plant_leaf_phenotypes":leaf_source_stats(root),
        "limitations":[
            "Existing authors' controlled original 2026 leaf-age experiment, not independent new ecological causal effect",
            "Asclepias syriaca leaf position top/bottom is not alien/native plant establishment",
            "Pupation does not imply adult eclosion or offspring recruitment",
            "Measured growth/development conditional on source measured/surviving cohort",
            "Four trial blocks only: individual caterpillars cannot substitute for replicated independent environments",
            "No plant leaf area source map or climate-season region link to chocho",
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
        "status":r["scientific_status"],
        "assigned":r["source_assigned_larvae"],
        "blocks":r["within_4_real_trial_blocks"],
        "top_rate":r["top_pupation_rate"],
        "bottom_rate":r["bottom_pupation_rate"],
        "stratified_OR":r["primary_Mantel_Haenszel_trial_stratified_pupation_OR_top_vs_bottom"],
        "trial_bootstrap":r["YearTrial_cluster_uncertainty"],
        "exact_fixed_block_p":r["exact_trial_conditional_test"]["two_sided_probability_order_p"],
        "stages":{k: {"n":v["original_rows"],
                        "conditional_Top_minus_Bottom":v["equal_block_mean_Top_minus_Bottom_conditional_on_measured"]}
                     for k,v in r["secondary_measured_stage_outcomes"].items()}
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
