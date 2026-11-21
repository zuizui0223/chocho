#!/usr/bin/env python3
"""Compare source-identical original GBIF224 recorded-assemblage ecological model refit.

This validates original model coefficients and ALL 169 leave-region-out folds,
not merely arithmetic on a published MSE aggregate. No original GBIF
coordinates, personal records, or third-party flora are copied to this receipt.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

DOI="10.15468/dl.umqmsf"
MODEL_SCHEMA="chocho_global224_GBIF_observed_WGSRPD3_2018_2024_beta_v1"
EXPECTED_REGIONS=169
EXPECTED_PAIRS=14196
EXPECTED_SPECIES=224
EXPECTED_FOLD_ROWS=676
SRC_COMMIT="dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0"


def num_tree(a,b,path="root",out=None):
    if out is None:out=[]
    if type(a)!=type(b):
        out.append({"path":path,"reason":"type_mismatch"});return out
    if isinstance(a,dict):
        if a.keys()!=b.keys():
            out.append({"path":path,"reason":"keys_mismatch"});return out
        for k in a:num_tree(a[k],b[k],path+"."+str(k),out)
    elif isinstance(a,list):
        if len(a)!=len(b):
            out.append({"path":path,"reason":"array_length"});return out
        for i,(x,y) in enumerate(zip(a,b)):num_tree(x,y,path+"."+str(i),out)
    elif isinstance(a,(float,int)) and not isinstance(a,bool):
        if not math.isclose(float(a),float(b),rel_tol=1e-6,abs_tol=1e-9):
            out.append({"path":path,"reason":"numeric_mismatch","original":a,"recomputed":b})
    elif a!=b:
        out.append({"path":path,"reason":"value_mismatch"})
    return out


def readfold(path):
    with path.open("r",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))


def validate(original:Path,refit:Path,oldfold:Path,newfold:Path)->dict:
    a=json.loads(original.read_text(encoding="utf-8"))
    b=json.loads(refit.read_text(encoding="utf-8"))
    for d in (a,b):
        assert d["schema"]==MODEL_SCHEMA
        assert d["sources"]["GBIF_doi"]==DOI
        assert d["cohort"]["region_count"]==EXPECTED_REGIONS
        assert d["cohort"]["evaluated_region_pairs"]==EXPECTED_PAIRS
        assert d["cohort"]["global_source_mapped_candidate_species"]==EXPECTED_SPECIES
    # Full aggregate statistical results and adverse held-world outcomes
    # should agree. Version tags and paths may contain nonnumeric differences.
    names=("cohort","outcome_Jaccard_summary","partial_rank",
           "heldout_region_model_scores","heldout_region_model_comparisons",
           "heldout_world_division_results")
    diffs=[]
    for key in names:num_tree(a[key],b[key],key,diffs)
    oa,ob=readfold(oldfold),readfold(newfold)
    if len(oa)!=EXPECTED_FOLD_ROWS or len(ob)!=EXPECTED_FOLD_ROWS:
        raise AssertionError(f"Expected {EXPECTED_FOLD_ROWS} held-region model-fold records")
    oa.sort(key=lambda row:(row["model"],row["wgsrpd3_code"]))
    ob.sort(key=lambda row:(row["model"],row["wgsrpd3_code"]))
    for i,(rowa,rowb) in enumerate(zip(oa,ob)):
        for k in ("model","wgsrpd3_code","world_division","test_pairs","training_pairs"):
            if rowa[k]!=rowb[k]:
                diffs.append({"path":f"fold[{i}].{k}","reason":"fold_design_mismatch"})
        for k in ("mse","mae"):
            if not math.isclose(float(rowa[k]),float(rowb[k]),rel_tol=1e-6,abs_tol=1e-9):
                diffs.append({"path":f"fold[{i}].{k}","reason":"fold_score_mismatch"})
    result={
        "schema":"chocho_GEB_v05_original_GBIF224_source_model_refit_QA_v1",
        "status":("SOURCE224_OBSERVED_MODEL_PRIMARY_REFIT_MATCH"
                  if not diffs else "SOURCE224_OBSERVED_MODEL_REFIT_DIFFERENCE"),
        "source_model_commit":SRC_COMMIT,
        "source_original_gbif_doi":DOI,
        "source_cohort":{"butterflies":EXPECTED_SPECIES,"regions":EXPECTED_REGIONS,
                         "region_pairs":EXPECTED_PAIRS,"held_region_model_rows":EXPECTED_FOLD_ROWS},
        "compared_result_sections":list(names),
        "relative_tolerance":1e-6,
        "absolute_tolerance":1e-9,
        "difference_count":len(diffs),
        "first_differences":diffs[:20],
        "original_host_after_nonlisted_alien_MSE_reduction_pct":
            a["heldout_region_model_comparisons"]["known_host_after_other_aliens"]["MSE_reduction_percent"],
        "refit_host_after_nonlisted_alien_MSE_reduction_pct":
            b["heldout_region_model_comparisons"]["known_host_after_other_aliens"]["MSE_reduction_percent"],
        "scope_limitations":[
            "Original 224 species observed-assemblage primary model from pinned regional/year aggregate inputs, not an independent GBIF occurrence archive.",
            "Raw occurrence ZIP-to-annual aggregation is audited separately; this model does not prove butterfly establishment or host causality.",
            "The later 29x2 effort-matched simulations and separately sourced 155-species RF model are NOT validated by this original-model primary refit.",
            "No source-level HOSTS/GloNAF or individual GBIF data are redistributed here.",
            "No anonymous peer-review deposit or archival DOI is created by this receipt."
        ],
    }
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original",type=Path)
    ap.add_argument("--refit",type=Path)
    ap.add_argument("--original-folds",type=Path)
    ap.add_argument("--refit-folds",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--self-test",action="store_true")
    z=ap.parse_args()
    if z.self_test:
        assert not num_tree({"a":[1.,2]},{"a":[1.+1e-9,2.]})
        assert num_tree({"a":1},{"a":2})
        assert EXPECTED_SPECIES==224 and EXPECTED_PAIRS==14196
        print("GBIF224_REFIT_NUMERIC_COMPARISON_SELFTEST_PASS")
        return
    if any(getattr(z,k) is None for k in ("original","refit","original_folds","refit_folds","out")):
        ap.error("All model input and output arguments required")
    r=validate(z.original,z.refit,z.original_folds,z.refit_folds)
    z.out.parent.mkdir(parents=True,exist_ok=True)
    z.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":r["status"],"difference_count":r["difference_count"],
                      "cohort":r["source_cohort"]},ensure_ascii=False))
    if r["difference_count"]:raise SystemExit("GBIF224_ORIGINAL_MODEL_REFIT_DIFFERENCE")


if __name__=="__main__":
    main()
