#!/usr/bin/env python3
"""Re-fit the exact original GBIF224 M2/M3 predictors on intensity-thinned GBIF response.

Source construction and original LOIO model come from immutable archived
chocho commit dfd49cfc. No butterfly nonrecord is a biological absence.
Synthetic effort-thinning is not estimation of detection probability.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np

FROZEN_SOURCE_SHA = "dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0"
SOURCE_DIR = "provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031"
N_REGIONS = 169
N_PAIRS = 14196
N_DRAWS = 29
SEEDS = {"dataset_key_calendar_day": 20261119, "indexed_candidate_record_count": 20261120}
ORIGINAL_GAIN = 0.09129078138747806
ORIGINAL_EARLY_MEAN = 0.11421295035628568
ORIGINAL_LATE_MEAN = 0.16788018941797136
MODEL_M2 = "M2_other_aliens"
MODEL_M3 = "M3_both"
MODEL_NAMES = ("M0_all_flora_effort_prior", "M1_known_host", MODEL_M2, MODEL_M3)
TOL = 1e-8

def mean(xs):
    return float(sum(xs)/len(xs))

def quantile(xs, p):
    a=sorted(xs)
    if not a:
        raise RuntimeError("Empty quantile sample")
    x=(len(a)-1)*p
    i=math.floor(x)
    j=math.ceil(x)
    return float(a[i]+(x-i)*(a[j]-a[i]))

def jaccard_on_original_pairs(masks, ii, jj):
    values=np.zeros(len(ii),dtype=float)
    n_undefined=0
    for row,(i,j) in enumerate(zip(ii,jj)):
        a,b=masks[int(i)],masks[int(j)]
        union=(a|b).bit_count()
        if union:
            values[row]=(a&b).bit_count()/union
        else:
            n_undefined+=1
    return values,n_undefined

def mean_mse(y, model_arrays, ia, ib, fit_fn, lam):
    # Exactly as original region-held-out scoring: every held region has
    # 168 test dyads, all other 14,028 pairs in training.
    values={model:[] for model in (MODEL_M2,MODEL_M3)}
    for held in range(N_REGIONS):
        test=(ia==held)|(ib==held)
        train=~test
        if int(test.sum())!=168 or int(train.sum())!=14028:
            raise RuntimeError("Original region-held-out 169-region fold drift")
        for model in (MODEL_M2,MODEL_M3):
            pred=fit_fn(train,test,y,model_arrays[model],lam)
            if not np.isfinite(pred).all():
                raise RuntimeError("Nonfinite source model predictions")
            err=pred-y[test]
            values[model].append(float(np.mean(err*err)))
    m2=mean(values[MODEL_M2])
    m3=mean(values[MODEL_M3])
    if m2<=0:
        raise RuntimeError("Invalid source M2 heldout MSE")
    return {
        "mean_MSE_M2":m2,
        "mean_MSE_M3":m3,
        "H_after_O_MSE_gain_percent":100*(m2-m3)/m2,
        "held_regions_M3_better":sum(a>b for a,b in zip(values[MODEL_M2],values[MODEL_M3])),
        "held_regions_M3_worse":sum(a<b for a,b in zip(values[MODEL_M2],values[MODEL_M3])),
    }

def analyse(args):
    frozen_root=args.frozen_root.resolve()
    current_scripts=args.current_scripts.resolve()
    if not (frozen_root / "scripts/analyze_gbif_global224_observed_beta.py").is_file():
        raise RuntimeError("Frozen original source scripts absent")
    # Keep old original imports resolving to the original archived scripts.
    sys.path.insert(0,str(current_scripts))
    sys.path.insert(0,str(frozen_root/"scripts"))
    import audit_gbif224_annual_effort_thinning_20261119 as thinning
    import analyze_gbif_global224_observed_beta as original
    from analyze_gbif_31_observed_host_beta import fit_ridge_predict

    if thinning.PIN!=FROZEN_SOURCE_SHA or thinning.N_REGIONS!=N_REGIONS:
        raise RuntimeError("Prior source audit provenance changed")
    if original.LOCK.read_bytes() != (frozen_root / SOURCE_DIR / "gbif_global224_observed_beta_lock_v1.json").read_bytes():
        raise RuntimeError("Frozen original model lock changed")

    original_folder=frozen_root / SOURCE_DIR
    original_csv=original_folder / "gbif_global224_observed_beta_region_holdout_v1.csv"
    original_receipt=json.loads((original_folder/"gbif_global224_observed_beta_v1.json").read_text())
    capture={}
    reference_scoring=original.score_held_out
    def capture_scoring(y,models,ia,ib,region_codes,world,lam=10.0):
        if capture:
            raise RuntimeError("Original model construction called more than once")
        capture.update({
            "raw_late_y":np.asarray(y,float).copy(),
            "models":{name:np.column_stack(values).astype(float,copy=True) for name,values in models.items()},
            "ia":np.asarray(ia,dtype=int).copy(),
            "ib":np.asarray(ib,dtype=int).copy(),
            "regions":list(region_codes),
            "world":list(world),
            "lambda":lam
        })
        return reference_scoring(y,models,ia,ib,region_codes,world,lam)
    original.score_held_out=capture_scoring
    try:
        original_result=original.analyse(eligibility="primary")
    finally:
        original.score_held_out=reference_scoring

    if set(capture["models"])!=set(MODEL_NAMES) or len(capture["regions"])!=169:
        raise RuntimeError("Original four model/source region set changed")
    if len(capture["ia"])!=N_PAIRS or len(capture["raw_late_y"])!=N_PAIRS:
        raise RuntimeError("Original 169-region pair model changed")
    if len(set(capture["regions"]))!=N_REGIONS:
        raise RuntimeError("Original region identities not unique")
    if capture["lambda"]!=10.:
        raise RuntimeError("Original regularization lambda changed")
    if original_result["cohort"]["region_count"]!=N_REGIONS:
        raise RuntimeError("Original independent observed cohort changed")
    if abs(original_result["heldout_region_model_comparisons"]["known_host_after_other_aliens"]["MSE_reduction_percent"] - ORIGINAL_GAIN)>TOL:
        raise RuntimeError("Recomputed original H-after-O MSE does not match source")
    if abs(original_receipt["heldout_region_model_comparisons"]["known_host_after_other_aliens"]["MSE_reduction_percent"] - ORIGINAL_GAIN)>TOL:
        raise RuntimeError("Original immutable GBIF readout mismatch")

    fold_by_name={}
    with original_csv.open(encoding="utf8",newline="") as f:
        for row in csv.DictReader(f):
            key=(row["model"],row["wgsrpd3_code"])
            if key in fold_by_name:
                raise RuntimeError("Original source region score duplicate")
            fold_by_name[key]=float(row["mse"])
    if len(fold_by_name)!=4*169:
        raise RuntimeError("Original full four-model score table missing")
    for model in MODEL_NAMES:
        for row in original_result["folds"][model]:
            code=row["wgsrpd3_code"]
            if abs(row["mse"]-fold_by_name[(model,code)])>TOL:
                raise RuntimeError(f"Original model did not reproduce source regional MSE: {model} {code}")

    source_args=argparse.Namespace(
        cells=original_folder/"gbif_global224_annual_v1.species_region_year.csv.gz",
        effort=original_folder/"gbif_global224_annual_v1.region_year_effort.csv",
        audit=original_folder/"gbif_global224_annual_v1.audit.json",
        region_scores=original_csv
    )
    codes,species,budgets,sp,source_audit=thinning.load_sources(source_args)
    if len(species)!=224 or set(codes)!=set(capture["regions"]):
        raise RuntimeError("Original annual species or region identities not identical")
    sorted_regions,per_region,_=thinning.construct(codes,species,budgets,sp)
    original_masks=thinning.original_masks(sorted_regions,per_region)
    lookup={code:i for i,code in enumerate(sorted_regions)}
    reorder=[lookup[code] for code in capture["regions"]]
    raw_prior=[original_masks["prior"][i] for i in reorder]
    raw_late=[original_masks["later"][i] for i in reorder]
    ia,ib=capture["ia"],capture["ib"]
    prior_y,undef_prior=jaccard_on_original_pairs(raw_prior,ia,ib)
    later_y,undef_late=jaccard_on_original_pairs(raw_late,ia,ib)
    if undef_prior or undef_late:
        raise RuntimeError("Original butterfly pair Jaccard unexpectedly undefined")
    if not np.allclose(later_y,capture["raw_late_y"],atol=1e-12,rtol=0):
        raise RuntimeError("Source pairwise current butterfly Jaccard != original frozen model")
    if abs(np.mean(prior_y)-ORIGINAL_EARLY_MEAN)>1e-12 or abs(np.mean(later_y)-ORIGINAL_LATE_MEAN)>1e-12:
        raise RuntimeError("Source original prior/late Jaccard mean mismatch")
    base_model=capture["models"][MODEL_M2]
    if base_model.shape[1]!=23 or capture["models"][MODEL_M3].shape[1]!=24:
        raise RuntimeError(f"Original M2/M3 feature dimensions changed: {[x.shape for x in capture['models'].values()]}")
    # Original fixed source controls: 9 geography, one history, whole flora,
    # two whole-flora richness, THEN prior butterfly similarity at index 13.
    prior_col=13
    if not np.allclose(base_model[:,prior_col],prior_y,atol=1e-12,rtol=0):
        raise RuntimeError("Cannot locate original prior butterfly covariate exactly")
    if not np.allclose(capture["models"][MODEL_M3][:,prior_col],prior_y,atol=1e-12,rtol=0):
        raise RuntimeError("M3 prior butterfly control not source-identical")

    valid_fold_keys=set(zip(np.minimum(ia,ib),np.maximum(ia,ib)))
    if len(valid_fold_keys)!=N_PAIRS:
        raise RuntimeError("Original pair source has duplicates")

    base_m2=base_model.copy()
    base_m3=capture["models"][MODEL_M3].copy()
    results={}
    for proxy_idx,(proxy,seed) in enumerate(SEEDS.items()):
        source_key=("days" if proxy=="dataset_key_calendar_day" else "records")
        rates,vectors=thinning.make_sampling_vectors(sorted_regions,per_region,budgets,source_key)
        rng=random.Random(seed)
        scores=[]
        for draw in range(N_DRAWS):
            masks=thinning.simulate(sorted_regions,vectors,rng)
            earlier=[masks["prior"][i] for i in reorder]
            later=[masks["later"][i] for i in reorder]
            yy,undef_l=jaccard_on_original_pairs(later,ia,ib)
            xx,undef_p=jaccard_on_original_pairs(earlier,ia,ib)
            X2=base_m2.copy()
            X3=base_m3.copy()
            X2[:,prior_col]=xx
            X3[:,prior_col]=xx
            value=mean_mse(yy,{MODEL_M2:X2,MODEL_M3:X3},ia,ib,fit_ridge_predict,capture["lambda"])
            value.update({
                "draw":draw+1,
                "prior_Jaccard_mean":float(xx.mean()),
                "later_Jaccard_mean":float(yy.mean()),
                "late_minus_prior":float(yy.mean()-xx.mean()),
                "prior_both_empty_dyads":int(undef_p),
                "later_both_empty_dyads":int(undef_l),
                "prior_empty_regions":int(sum(k==0 for k in earlier)),
                "later_empty_regions":int(sum(k==0 for k in later)),
            })
            scores.append(value)
        gain=[r["H_after_O_MSE_gain_percent"] for r in scores]
        rawdelta=[r["late_minus_prior"] for r in scores]
        results[proxy]={
            "effort_proxy":source_key,
            "seed":seed,
            "draw_count":len(scores),
            "mean_period_effort_region_retention": {
                "prior":mean([rates[(c,"prior")] for c in sorted_regions]),
                "later":mean([rates[(c,"later")] for c in sorted_regions])
            },
            "H_after_O_held_region_MSE_reduction_pct":{
                "median":quantile(gain,.5),
                "mean":mean(gain),
                "min":min(gain),
                "max":max(gain),
                "MC_central_95":[quantile(gain,.025),quantile(gain,.975)],
                "fraction_positive":sum(g>0 for g in gain)/len(gain),
            },
            "thinned_late_minus_prior_Jaccard":{
                "median":quantile(rawdelta,.5),
                "MC_central_95":[quantile(rawdelta,.025),quantile(rawdelta,.975)]
            },
            "num_draws_with_undefined_outcome_pairs":sum(r["later_both_empty_dyads"]>0 for r in scores),
            "num_draws_with_undefined_prior_pairs":sum(r["prior_both_empty_dyads"]>0 for r in scores),
            "first_five_draws":scores[:5],
            "per_draw":scores,
        }
        print("THINNED_HOST_GAIN",json.dumps({
            "proxy":proxy,
            "n_draws":len(scores),
            "original_gain":ORIGINAL_GAIN,
            "gain_summary":results[proxy]["H_after_O_held_region_MSE_reduction_pct"],
            "Jaccard_summary":results[proxy]["thinned_late_minus_prior_Jaccard"]
        },ensure_ascii=False),flush=True)

    source={
        "GBIF_DOI":"10.15468/dl.umqmsf",
        "pinned_original_repo_commit":FROZEN_SOURCE_SHA,
        "source_audit":source_audit,
        "all_4_original_models_exactly_reproduced":True,
        "same_169_original_regions":True,
        "same_14196_original_region_pairs":True,
        "exact_raw_prior_and_late_Jaccard_for_each_pair":True,
        "original_source_H_after_O_pct":ORIGINAL_GAIN,
        "original_raw_prior_Jaccard_mean":float(prior_y.mean()),
        "original_raw_later_Jaccard_mean":float(later_y.mean()),
        "n_original_source_predictor_columns_M2":int(base_model.shape[1]),
        "n_original_source_predictor_columns_M3":int(base_m3.shape[1]),
        "original_prior_Jaccard_source_column":prior_col,
    }
    result={
        "schema":"chocho_GBIF224_observed_HO_on_effort_thinned_communities_20261119_v1",
        "status":"SOURCE_REPRODUCED_EXPLORATORY_RETROSPECTIVE_OBSERVATION_INTENSITY_SIMULATION_NOT_OCCUPANCY",
        "sample":{"butterflies":224,"regions":169,"original_pairs":14196,
                  "leave_one_region_out_folds":169},
        "methods":{"draws_per_proxy":N_DRAWS,
                   "predictors":"frozen original M2 and M3 GloNAF/HOSTS + regional native/climate/colonial/effort covariates; only prior butterfly Jaccard re-simulated",
                   "outcome":"GBIF later observed butterfly Jaccard after source-conditional indexed-record thinning",
                   "standardisation":"original training-only ridge lambda=10, clipped [0,1]"},
        "source":source,
        "results":results,
        "limits":[
            "Indexed GBIF records are not standardized repeat visits, and no nonrecord means a proven ecological absence.",
            "Thinning is an assumption-dependent, positive-record simulation, not fitted detection/occupancy correction.",
            "Underlying historical GloNAF/HOSTS data are contemporary and may postdate GBIF observations.",
            "The original late-period effort/record controls are held fixed while synthetic late species presence is thinned; the resulting model is retrospective and can be misspecified.",
            "Model M2/M3 is evaluated on 29 correlated simulations per proxy; MC intervals describe conditional variability rather than biological confidence.",
            "Original 239 butterfly potential-resource increase and independent larval records are distinct outcomes; host-future colonisation, cargo shipping and reproductive performance are not identified."
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--frozen-root",type=Path,required=True)
    ap.add_argument("--current-scripts",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    # Original model paths and GloNAF/TDWG sources are source-root-relative.
    original_cwd=Path.cwd()
    try:
        os.chdir(args.frozen_root.resolve())
        result=analyse(args)
    finally:
        os.chdir(original_cwd)
    print("SOURCE_REPLAY_PASS",json.dumps({
        "n":result["sample"],
        "H_after_O_original":result["source"]["original_source_H_after_O_pct"],
        "proxies":{p:r["H_after_O_held_region_MSE_reduction_pct"] for p,r in result["results"].items()}
    },ensure_ascii=False))

if __name__=="__main__":
    main()
