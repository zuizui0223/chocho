#!/usr/bin/env python3
"""Source-checked descriptive regional resampling of independent GBIF butterfly endpoint.

The original four models, leave-one-WGSRPD3-region-out fits, and three effort-
eligibility cohorts are FROZEN. This script does NOT recover dated effort cells,
refit the original model, perform occupancy detection correction or model ecology.
It audits source outcomes and describes uncertainty of the already saved paired
regional cross-validation scores. Uses Python stdlib only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
CASES = {
    "primary": ("gbif_global224_observed_beta_region_holdout_v1.csv",
                "gbif_global224_observed_beta_v1.json", 169),
    "relaxed": ("gbif_global224_observed_beta_relaxed_region_holdout_v1.csv",
                "gbif_global224_observed_beta_relaxed_v1.json", 194),
    "stringent": ("gbif_global224_observed_beta_stringent_region_holdout_v1.csv",
                  "gbif_global224_observed_beta_stringent_v1.json", 133),
}
MODELS = ("M0_all_flora_effort_prior", "M1_known_host",
          "M2_other_aliens", "M3_both")
DRAW_COUNT = 4999
SEED = 20261119
NUMERIC_TOLERANCE = 1e-8

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def mean(values):
    vals = list(values)
    if not vals:
        raise RuntimeError("An empty set cannot be scored")
    return sum(vals) / len(vals)

def quantile(values, q):
    a = sorted(values)
    x = (len(a)-1)*q
    i = math.floor(x)
    j = math.ceil(x)
    return a[i]+(x-i)*(a[j]-a[i])

def score(rows):
    m = {key: mean(r[key] for r in rows) for key in MODELS}
    if m[MODELS[2]] <= 0 or m[MODELS[0]] <= 0:
        raise RuntimeError("Invalid MSE denominator")
    return {
        "H_given_other_percent": 100*(m[MODELS[2]]-m[MODELS[3]])/m[MODELS[2]],
        "H_given_base_percent": 100*(m[MODELS[0]]-m[MODELS[1]])/m[MODELS[0]],
        "other_given_base_percent": 100*(m[MODELS[0]]-m[MODELS[2]])/m[MODELS[0]],
        "means": m,
    }

def world_equal_score(rows):
    buckets = defaultdict(list)
    for row in rows:
        buckets[row["world_division"]].append(row)
    if not buckets:
        raise RuntimeError("No world divisions")
    m2 = mean(mean(x[MODELS[2]] for x in grp) for grp in buckets.values())
    m3 = mean(mean(x[MODELS[3]] for x in grp) for grp in buckets.values())
    return 100*(m2-m3)/m2

def read_original_case(root, case, filename, receiptname, expected_n):
    source_file = root / filename
    receipt_file = root / receiptname
    receipt = json.loads(receipt_file.read_text())
    data = defaultdict(dict)
    with source_file.open(newline="", encoding="utf8") as f:
        reader = csv.DictReader(f)
        required = {"model","wgsrpd3_code","world_division",
                    "test_pairs","training_pairs","mse","mae"}
        if set(reader.fieldnames or ()) != required:
            raise RuntimeError("Source heldout region CSV columns drift")
        count=0
        for r in reader:
            count += 1
            model,code = r["model"],r["wgsrpd3_code"]
            if model not in MODELS or not code:
                raise RuntimeError("Model/code not from 224 source")
            reg=data[code]
            if model in reg:
                raise RuntimeError("Duplicate source heldout model×region")
            if "world_division" in reg and reg["world_division"]!=r["world_division"]:
                raise RuntimeError("Model-specific regional world-division drift")
            reg["world_division"]=r["world_division"]
            if not (len(r["world_division"])==1 and r["world_division"].isdigit()):
                raise RuntimeError("Bad world-division source")
            # A regional holdout includes n-1 test pairs and all other pair training.
            if (int(r["test_pairs"])!=expected_n-1 or
                int(r["training_pairs"])!=(expected_n-1)*(expected_n-2)//2):
                raise RuntimeError("Original held-out geographic pair count drift")
            mse=float(r["mse"])
            mae=float(r["mae"])
            if not (math.isfinite(mse) and mse>0 and math.isfinite(mae) and mae>=0):
                raise RuntimeError("Nonfinite original heldout scoring")
            reg[model]=mse
    if count!=expected_n*len(MODELS) or len(data)!=expected_n:
        raise RuntimeError(f"Original 224 GBIF {case} cohort size drift")
    codes=sorted(data)
    rows=[data[c] for c in codes]
    if any(len(r)!=5 for r in rows):
        raise RuntimeError("Missing held-out region model")
    if int(receipt["cohort"]["region_count"])!=expected_n:
        raise RuntimeError("Independent GBIF source receipt region count changed")
    if int(receipt["cohort"]["global_source_mapped_candidate_species"])!=224:
        raise RuntimeError("Independent GBIF candidate species count changed")
    if int(receipt["cohort"]["evaluated_region_pairs"])!=expected_n*(expected_n-1)//2:
        raise RuntimeError("Independent GBIF paired universe changed")
    values=score(rows)
    expected=receipt["heldout_region_model_scores"]
    for model in MODELS:
        if abs(values["means"][model]-expected[model]["region_equal_MSE"])>NUMERIC_TOLERANCE:
            raise RuntimeError(f"Original empirical 224 region MSE mismatch: {case} {model}")
    if abs(values["H_given_other_percent"]-receipt["heldout_region_model_comparisons"]["known_host_after_other_aliens"]["MSE_reduction_percent"])>NUMERIC_TOLERANCE:
        raise RuntimeError(f"Original observed H|other estimate drift: {case}")
    return codes,rows,receipt,source_file,receipt_file

def boot(rows, seed):
    groups=defaultdict(list)
    for r in rows:
        groups[r["world_division"]].append(r)
    worlds=sorted(groups)
    n=len(rows)
    rng=random.Random(seed)
    all_draw=[]
    world_strat=[]
    world_equal=[]
    for _ in range(DRAW_COUNT):
        sampled=[rows[rng.randrange(n)] for _ in range(n)]
        strat=[]
        for world in worlds:
            g=groups[world]
            strat.extend(g[rng.randrange(len(g))] for _ in range(len(g)))
        all_draw.append(score(sampled)["H_given_other_percent"])
        world_strat.append(score(strat)["H_given_other_percent"])
        world_equal.append(world_equal_score(strat))
    return {
        "draw_count":DRAW_COUNT,
        "seed":seed,
        "region_resample_unstratified_descriptive_95": [quantile(all_draw,.025),quantile(all_draw,.975)],
        "within_world_stratified_descriptive_95": [quantile(world_strat,.025),quantile(world_strat,.975)],
        "equal_world_weight_descriptive_95": [quantile(world_equal,.025),quantile(world_equal,.975)],
        "resampling_unit_note":"These regional CV errors share training data and dependent test pairs. Percentile intervals are descriptive, not source-independent 95% confidence intervals.",
    }

def analyse(root):
    out={}
    for offset,(case,(filename,receiptname,n)) in enumerate(CASES.items()):
        codes,rows,receipt,f,rf=read_original_case(root,case,filename,receiptname,n)
        worlds=sorted({r["world_division"] for r in rows})
        levels=[]
        for w in worlds:
            sub=[r for r in rows if r["world_division"]==w]
            other=[r for r in rows if r["world_division"]!=w]
            levels.append({
                "world_division":w,
                "n_source_regions":len(sub),
                "same_original_folds_within_world_H_given_other_pct":score(sub)["H_given_other_percent"],
                "remove_world_from_existing_region_errors_H_given_other_pct":score(other)["H_given_other_percent"],
                "not_retrained":True,
            })
        v=score(rows)
        out[case]={
            "cohort":{"n_regions":n,
                      "n_unique_region_pairs":n*(n-1)//2,
                      "frozen_GBIF_candidate_butterflies":224,
                      "region_labels":codes,
                      "world_division_counts":{w:sum(r["world_division"]==w for r in rows) for w in worlds}},
            "source":{"region_score_csv_sha256":sha(f),"original_GBIF_readout_sha256":sha(rf),
                      "original_GBIF_doi":receipt["sources"]["GBIF_doi"],
                      "original_annual_species_region_year_sha256":receipt["sources"]["annual_GBIF_species_region_year_sha256"],
                      "original_annual_effort_sha256":receipt["sources"]["annual_GBIF_region_year_effort_sha256"]},
            "effects":{"H_given_other_percent":v["H_given_other_percent"],
                       "H_given_base_percent":v["H_given_base_percent"],
                       "other_given_base_percent":v["other_given_base_percent"],
                       "regions_H_given_other_improve":sum(r[MODELS[3]]<r[MODELS[2]] for r in rows),
                       "regions_H_given_other_worse":sum(r[MODELS[3]]>r[MODELS[2]] for r in rows),
                       "equal_world_weight_H_given_other_percent":world_equal_score(rows),
                       "source_original_mean_MSE":v["means"]},
            "regional_resampling":boot(rows,SEED+offset),
            "original_region_errors_world_influence":levels,
            "original_observational_gate":receipt["cohort"]["region_recording_threshold"],
        }
    return {
        "schema":"chocho_GBIF224_observed_region_paired_resampling_20261119_v1",
        "status":"SOURCE_CHECKED_POST_HOC_DESCRIPTIVE_REANALYSIS_NOT_DETECTION_CORRECTION",
        "audit_note":"No refitted original GloNAF/HOSTS model, new historical plant timestamps or GBIF occurrence rows; source-original held-out region CSV scores checked exactly against original report.",
        "cases":out,
        "limits":[
            "All three effort-eligibility cohorts overlap and are not independent studies.",
            "Region bootstrap samples reused original dependent pairwise leave-one-region-out predictions and training sets; descriptive sensitivity only.",
            "World influence analysis removes worlds from stored region errors and DOES NOT retrain a model on held-out world geography.",
            "Dataset×calendar-day counts were used only in source-original region thresholds. No matching or detection correction to underlying daily GBIF records was performed.",
            "Original independent GBIF assemblage was recorded rather than complete occupancy. Missing indexed records are not absence.",
            "WCVP/HOSTS/GloNAF botanical distribution inputs are contemporary and cannot identify historical arrivals or causal plant/ship/consumer pathways.",
            "Source original 169/194/133 region eligibility changes population and recording effort simultaneously; threshold sensitivity is not a causal effort adjustment.",
        ],
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,default=ROOT)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    data=analyse(args.source_root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    for k,v in data["cases"].items():
        print(k, json.dumps({"n":v["cohort"]["n_regions"],
                  "H_given_other_pct":v["effects"]["H_given_other_percent"],
                  "n_better":v["effects"]["regions_H_given_other_improve"],
                  "regional_bootstrap":v["regional_resampling"]["region_resample_unstratified_descriptive_95"],
                  "within_world_bootstrap":v["regional_resampling"]["within_world_stratified_descriptive_95"],
                  "equal_world_weight":v["effects"]["equal_world_weight_H_given_other_percent"]}))

if __name__=="__main__":
    main()
