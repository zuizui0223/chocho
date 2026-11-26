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


def transfer_by_entire_world(y,model_arrays,ia,ib,region_world,fit_fn,lam):
    """Full unseen TDWG Level1: no held-world region endpoint occurs in training."""
    worlds=[str(x) for x in region_world]
    if sorted(set(worlds))!=list("12345678") or len(worlds)!=N_REGIONS:
        raise RuntimeError("Original eight TDWG Level1 groups changed")
    pairs={}
    for w in sorted(set(worlds)):
        held=np.asarray([x==w for x in worlds],dtype=bool)
        test=held[ia]|held[ib]
        train=~test
        k=int(held.sum())
        expected_test=N_PAIRS-((N_REGIONS-k)*(N_REGIONS-k-1)//2)
        if int(test.sum())!=expected_test or int(train.sum())!=N_PAIRS-expected_test:
            raise RuntimeError("Original held-world complete dyad population changed")
        if int(test.sum())<30 or int(train.sum())<100:
            raise RuntimeError("Original held-world fold too small")
        scores={}
        preds={}
        for model in (MODEL_M2,MODEL_M3):
            pred=fit_fn(train,test,y,model_arrays[model],lam)
            if not np.isfinite(pred).all():
                raise RuntimeError("Nonfinite held-world model prediction")
            scores[model]=float(np.mean((pred-y[test])**2))
            preds[model]=pred
        within=(held[ia]&held[ib])[test]
        cross=~within
        partition={
            "within_held_world":score_pair_partition(y[test],preds[MODEL_M2],preds[MODEL_M3],within),
            "held_world_to_other":score_pair_partition(y[test],preds[MODEL_M2],preds[MODEL_M3],cross),
        }
        if sum(x["n_pairs"] for x in partition.values())!=int(test.sum()):
            raise RuntimeError("Held-world pair classes not exhaustive")
        for model,col in ((MODEL_M2,"MSE_M2"),(MODEL_M3,"MSE_M3")):
            reconstructed=sum(z["n_pairs"]*z[col] for z in partition.values())/int(test.sum())
            if abs(reconstructed-scores[model])>1e-12:
                raise RuntimeError("Held-world prediction MSE does not reconstruct from within/cross pair classes")
        mse2,mse3=scores[MODEL_M2],scores[MODEL_M3]
        if mse2<=0:
            raise RuntimeError("Invalid held-world source mean squared error")
        pairs[w]={
            "world_division":w,
            "regions_held_out":k,
            "source_train_pairs":int(train.sum()),
            "source_test_pairs":int(test.sum()),
            "MSE_M2":mse2,
            "MSE_M3":mse3,
            "H_after_O_MSE_improvement_pct":100*(mse2-mse3)/mse2,
            "dyad_partitions":partition,
        }
    units=list(pairs.values())
    e2=mean([r["MSE_M2"] for r in units])
    e3=mean([r["MSE_M3"] for r in units])
    pair_m2=sum(r["MSE_M2"]*r["source_test_pairs"] for r in units)/sum(r["source_test_pairs"] for r in units)
    pair_m3=sum(r["MSE_M3"]*r["source_test_pairs"] for r in units)/sum(r["source_test_pairs"] for r in units)
    return {
        "world_equal_mean_MSE_gain_pct":100*(e2-e3)/e2,
        "test_pair_weighted_mean_MSE_gain_pct":100*(pair_m2-pair_m3)/pair_m2,
        "mean_of_eight_signed_world_pct":mean([r["H_after_O_MSE_improvement_pct"] for r in units]),
        "n_worlds_improve":sum(r["H_after_O_MSE_improvement_pct"]>0 for r in units),
        "per_world":pairs,
        "not_iid_and_not_pooled_independent_pairs":True,
    }

def summarize_strict_world_scores(scores,raw):
    aggregates={
        "world_equal_mean_MSE_gain_pct",
        "test_pair_weighted_mean_MSE_gain_pct",
        "mean_of_eight_signed_world_pct",
    }
    record={"raw_original":raw, "n_simulations":len(scores)}
    for name in sorted(aggregates):
        arr=[s["world_transfer"][name] for s in scores]
        record[name]={
            "median":quantile(arr,.5),
            "mean":mean(arr),
            "min":min(arr),"max":max(arr),
            "MC_central_95":[quantile(arr,.025),quantile(arr,.975)],
            "n_positive":sum(x>0 for x in arr),
            "n_draws":len(arr),
        }
    n=[s["world_transfer"]["n_worlds_improve"] for s in scores]
    record["number_of_eight_worlds_improving"]={
        "median":quantile(n,.5),
        "min":min(n),"max":max(n),
        "n_simulations":len(n),
    }
    world_out={}
    for w in sorted(raw["per_world"]):
        vals=[s["world_transfer"]["per_world"][w]["H_after_O_MSE_improvement_pct"] for s in scores]
        world_out[w]={
            "original_raw":raw["per_world"][w],
            "median_H_after_O_MSE_gain_pct":quantile(vals,.5),
            "MC_central_95":[quantile(vals,.025),quantile(vals,.975)],
            "min":min(vals),"max":max(vals),
            "n_positive":sum(x>0 for x in vals),
            "n_draws":len(vals),
        }
    record["per_world"]=world_out
    return record


def diagnostics_distribution(x):
    a=np.asarray(x,dtype=float)
    if a.ndim!=1 or len(a)==0 or not np.isfinite(a).all():
        raise RuntimeError("Missing, nonfinite or multidimensional diagnostics array")
    return {
        "n":int(len(a)),
        "mean":float(a.mean()),
        "sd":float(a.std()),
        "median":float(np.median(a)),
        "q05":float(np.quantile(a,.05)),
        "q25":float(np.quantile(a,.25)),
        "q75":float(np.quantile(a,.75)),
        "q95":float(np.quantile(a,.95)),
        "fraction_exact_zero":float(np.mean(a==0)),
    }

def score_pair_partition(ytest,p2,p3,loc):
    """Disaggregate the ORIGINAL world held-out MSE by source dyad category."""
    z=np.asarray(loc,dtype=bool)
    if not z.any():
        raise RuntimeError("No dyads in held-world partition")
    y,p,q=ytest[z],p2[z],p3[z]
    e2=float(np.mean((p-y)**2))
    e3=float(np.mean((q-y)**2))
    if e2<=0:
        raise RuntimeError("No held-world baseline MSE")
    delta=q-p
    algebraic=float(2*np.mean((y-p)*delta)-np.mean(delta*delta))
    if abs((e2-e3)-algebraic)>1e-12:
        raise RuntimeError("Predictive change cannot be decomposed into calibration/error change")
    return {
        "n_pairs":int(len(y)),
        "MSE_M2":e2,
        "MSE_M3":e3,
        "H_after_O_MSE_improvement_pct":100*(e2-e3)/e2,
        "mean_recorded_butterfly_Jaccard":float(y.mean()),
        "mean_predicted_butterfly_Jaccard_M2":float(p.mean()),
        "mean_predicted_butterfly_Jaccard_M3":float(q.mean()),
        "mean_prediction_change_M3_minus_M2":float(delta.mean()),
        "fraction_M3_predicts_higher_Jaccard":float(np.mean(delta>0)),
        "squared_error_difference_M2_minus_M3":e2-e3,
        "exact_algebraic_error_difference":algebraic,
    }

def source_plant_and_effort_support(capture,calls,budgets):
    """Read original source GloNAF H/O masks, not a new host-list selection."""
    if len(calls)!=7:
        raise RuntimeError("Source original Jaccard construction changed from 4 temporal+3 plant calls")
    for name,idx,model,col in (
        ("whole",4,MODEL_M2,10),
        ("H",5,MODEL_M3,22),
        ("O",6,MODEL_M2,22),
        ("O_M3",6,MODEL_M3,23),
    ):
        observed=calls[idx]["jaccard"]
        if (len(observed)!=N_PAIRS or
            not np.allclose(observed,capture["models"][model][:,col],rtol=0,atol=1e-12)):
            raise RuntimeError("Original source plant "+name+" Jaccard cannot be tied to frozen predictor column")
    if not np.allclose(calls[0]["jaccard"],capture["models"][MODEL_M2][:,13],rtol=0,atol=1e-12):
        raise RuntimeError("Frozen original prior butterfly Jaccard covariate mismatch")
    if not np.allclose(calls[2]["jaccard"],capture["raw_late_y"],rtol=0,atol=1e-12):
        raise RuntimeError("Frozen later GBIF butterfly similarity source mismatch")
    masks={name:[int(x) for x in calls[idx]["masks"]]
           for name,idx in (("whole",4),("H",5),("O",6))}
    if any(len(x)!=N_REGIONS for x in masks.values()):
        raise RuntimeError("Frozen 169-region GloNAF masks have changed")
    if any((h|o)!=whole or (h&o)!=0
           for whole,h,o in zip(masks["whole"],masks["H"],masks["O"])):
        raise RuntimeError("H/O source masks are not a partition of all GloNAF aliens")
    codes=capture["regions"]
    worlds=[str(x) for x in capture["world"]]
    original_regions=[]
    for k,(code,w) in enumerate(zip(codes,worlds)):
        H=masks["H"][k].bit_count()
        O=masks["O"][k].bit_count()
        whole=masks["whole"][k].bit_count()
        days_pr=budgets[code]["prior"]["days"]
        days_la=budgets[code]["later"]["days"]
        rec_pr=budgets[code]["prior"]["records"]
        rec_la=budgets[code]["later"]["records"]
        if min(days_pr,days_la,rec_pr,rec_la,whole)<=0:
            raise RuntimeError("Missing original source effort or GloNAF alien species")
        original_regions.append({
            "world":w,
            "source_H_alien_taxa":H,
            "source_O_alien_taxa":O,
            "source_H_fraction_of_all_alien_taxa":H/whole,
            "source_all_alien_taxa":whole,
            "prior_dataset_days":days_pr,
            "later_dataset_days":days_la,
            "later_vs_prior_dataset_days":days_la/days_pr,
            "prior_indexed_butterfly_records":rec_pr,
            "later_indexed_butterfly_records":rec_la,
            "later_vs_prior_indexed_records":rec_la/rec_pr,
        })
    source_region_summary={}
    groups={
        "South_America_Level1_8": [r for r in original_regions if r["world"]=="8"],
        "outside_South_America_Level1_8": [r for r in original_regions if r["world"]!="8"],
        "North_America_Level1_7": [r for r in original_regions if r["world"]=="7"],
        "all_169":original_regions,
    }
    for g,rs in groups.items():
        source_region_summary[g]={
            "n_original_regions":len(rs),
            "source_registered_host_alien_plant_richness":diagnostics_distribution([x["source_H_alien_taxa"] for x in rs]),
            "source_unlisted_alien_plant_richness":diagnostics_distribution([x["source_O_alien_taxa"] for x in rs]),
            "source_H_fraction_of_all_alien_plants":diagnostics_distribution([x["source_H_fraction_of_all_alien_taxa"] for x in rs]),
            "source_total_alien_plant_richness":diagnostics_distribution([x["source_all_alien_taxa"] for x in rs]),
            "prior_dataset_day_effort":diagnostics_distribution([x["prior_dataset_days"] for x in rs]),
            "later_dataset_day_effort":diagnostics_distribution([x["later_dataset_days"] for x in rs]),
            "late_over_prior_dataset_day_ratio":diagnostics_distribution([x["later_vs_prior_dataset_days"] for x in rs]),
            "late_over_prior_indexed_butterfly_record_ratio":diagnostics_distribution([x["later_vs_prior_indexed_records"] for x in rs]),
        }
    ia,ib=capture["ia"],capture["ib"]
    held=np.asarray([w=="8" for w in worlds],dtype=bool)
    training=~(held[ia]|held[ib])
    within=held[ia]&held[ib]
    cross=(held[ia]|held[ib])&(~within)
    test=within|cross
    if (int(training.sum())!=11026 or int(within.sum())!=190 or
        int(cross.sum())!=2980 or int(test.sum())!=3170):
        raise RuntimeError("Original South American 20-region dyads changed")
    descriptors={
        "original_H_alien_Jaccard":calls[5]["jaccard"],
        "original_O_alien_Jaccard":calls[6]["jaccard"],
        "original_whole_alien_Jaccard":calls[4]["jaccard"],
        "original_prior_recorded_butterfly_Jaccard":calls[0]["jaccard"],
    }
    descriptor_summary={}
    for name,raw in descriptors.items():
        values=np.asarray(raw,dtype=float)
        ref=diagnostics_distribution(values[training])
        std=ref["sd"]
        descriptor_summary[name]={
            "training_source_excludes_all_South_America":ref,
            "South_America_test_all_3170_pairs":diagnostics_distribution(values[test]),
            "South_America_test_within_190_pairs":diagnostics_distribution(values[within]),
            "South_America_test_cross_2980_pairs":diagnostics_distribution(values[cross]),
            "all_test_mean_minus_training_mean_scaled_by_training_sd":
                float((values[test].mean()-ref["mean"])/std) if std>1e-9 else None,
            "all_test_fraction_outside_training_q05_q95":
                float(np.mean((values[test]<ref["q05"])|(values[test]>ref["q95"]))),
        }
    return {
        "source_plant_ID_masks_exactly_match_original_frozen_model":True,
        "source_H_means_catalogued_HOSTS_larval_alien_plant_not_all_true_hosts":True,
        "source_169_original_region_gloNAF_H_O_and_effort_summary":source_region_summary,
        "South_America_heldout_dyad_source_support":descriptor_summary,
        "South_America_train_n":int(training.sum()),
        "South_America_within_n":int(within.sum()),
        "South_America_cross_n":int(cross.sum()),
        "South_America_all_test_n":int(test.sum()),
    }

def summarize_partition(simulations, original, world):
    groups=("within_held_world","held_world_to_other")
    out={"original_unthinned":original["per_world"][world]["dyad_partitions"],
         "n_simulations":len(simulations),"by_pair_class":{}}
    for part in groups:
        v=[item["world_transfer"]["per_world"][world]["dyad_partitions"][part] for item in simulations]
        dist={}
        for k in ("H_after_O_MSE_improvement_pct","mean_prediction_change_M3_minus_M2",
                  "mean_recorded_butterfly_Jaccard","squared_error_difference_M2_minus_M3"):
            vals=[x[k] for x in v]
            dist[k]={
                "median":quantile(vals,.5),
                "MC_central_95":[quantile(vals,.025),quantile(vals,.975)],
                "min":min(vals),"max":max(vals),
                "n_positive":sum(x>0 for x in vals),
            }
        dist["source_pairs"]=v[0]["n_pairs"]
        out["by_pair_class"][part]=dist
    return out

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
    jaccard_source_calls=[]
    reference_jaccard=original.jaccard
    def capture_jaccard(masks,ia,ib):
        values=reference_jaccard(masks,ia,ib)
        jaccard_source_calls.append({
            "masks":list(masks),
            "jaccard":np.asarray(values[0],dtype=float).copy(),
            "undefined":int(values[1])
        })
        return values
    original.jaccard=capture_jaccard
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
        original.jaccard=reference_jaccard

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

    raw_world=transfer_by_entire_world(
        later_y,
        {MODEL_M2:base_model,MODEL_M3:capture["models"][MODEL_M3]},
        ia,ib,capture["world"],fit_ridge_predict,capture["lambda"])
    archived_world=original_result["heldout_world_division_results"]["by_world_division"]
    if set(archived_world)!=set(raw_world["per_world"]):
        raise RuntimeError("Original world divisions differ from archive")
    for w,row in raw_world["per_world"].items():
        earlier=archived_world[w]
        if (row["regions_held_out"]!=earlier["regions_heldout"]
            or row["source_test_pairs"]!=earlier["test_pairs"]
            or row["source_train_pairs"]!=earlier["training_pairs"]):
            raise RuntimeError("Original strict geographic pair fold mismatch for "+w)
        for model,score_key in ((MODEL_M2,"MSE_M2"),(MODEL_M3,"MSE_M3")):
            if abs(row[score_key]-earlier[model])>TOL:
                raise RuntimeError("Cannot reproduce original world-holdout MSE: "+w+"/"+model)
        if abs(row["H_after_O_MSE_improvement_pct"]-earlier["known_host_given_other_pct"])>TOL:
            raise RuntimeError("Original held-world signed host increment drift: "+w)
    print("SOURCE_WORLD_TRANSFER_REPRODUCED",json.dumps({
        "raw_original_equal_world_gain_pct":raw_world["world_equal_mean_MSE_gain_pct"],
        "n_worlds_improve":raw_world["n_worlds_improve"],
        "original_world_8_pct":raw_world["per_world"]["8"]["H_after_O_MSE_improvement_pct"]
    }),flush=True)

    source_diagnostics=source_plant_and_effort_support(capture,jaccard_source_calls,budgets)
    print("ORIGINAL_GLO_NAF_SUPPORT_AUDITED",json.dumps({
        "South_America_known_H_fraction":source_diagnostics["source_169_original_region_gloNAF_H_O_and_effort_summary"]["South_America_Level1_8"]["source_H_fraction_of_all_alien_plants"]["median"],
        "outside_South_America_known_H_fraction":source_diagnostics["source_169_original_region_gloNAF_H_O_and_effort_summary"]["outside_South_America_Level1_8"]["source_H_fraction_of_all_alien_plants"]["median"],
        "held_test_H_jaccard":source_diagnostics["South_America_heldout_dyad_source_support"]["original_H_alien_Jaccard"]["South_America_test_all_3170_pairs"]["median"]
    }),flush=True)

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
            value["world_transfer"]=transfer_by_entire_world(
                yy,{MODEL_M2:X2,MODEL_M3:X3},
                ia,ib,capture["world"],fit_ridge_predict,capture["lambda"])
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
            "South_America_within_vs_cross_dyad_calibration":summarize_partition(scores,raw_world,"8"),
            "North_America_within_vs_cross_dyad_calibration":summarize_partition(scores,raw_world,"7"),
            "strict_entire_TDWG_Level1_world_holdout":summarize_strict_world_scores(scores,raw_world),
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
        "all_eight_original_entire_world_holdout_scores_reproduced":True,
        "original_raw_entire_world_transfer":raw_world,
        "source_H_O_plant_and_observation_effort_support_audit":source_diagnostics,
        "original_raw_prior_Jaccard_mean":float(prior_y.mean()),
        "original_raw_later_Jaccard_mean":float(later_y.mean()),
        "n_original_source_predictor_columns_M2":int(base_model.shape[1]),
        "n_original_source_predictor_columns_M3":int(base_m3.shape[1]),
        "original_prior_Jaccard_source_column":prior_col,
    }
    result={
        "schema":"chocho_GBIF224_South_America_host_transfer_failure_diagnostics_20261119_v1",
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
            "Original 239 butterfly potential-resource increase and independent larval records are distinct outcomes; host-future colonisation, cargo shipping and reproductive performance are not identified.",
            "Entire-world heldout dyads are dependent: cross-world dyads occur in two world-fold tests, and the eight world groups are highly uneven (Level1 6 has only two regions).",
            "World-balanced aggregate can mask negative transfer in individual held-out biogeographic regions; no universality claim is licensed."
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
    # Resolve paths BEFORE switching to frozen original working directory.
    args.frozen_root=args.frozen_root.resolve()
    args.current_scripts=args.current_scripts.resolve()
    args.output=args.output.resolve()
    try:
        os.chdir(args.frozen_root)
        result=analyse(args)
    finally:
        os.chdir(original_cwd)
    print("SOUTH_AMERICA_FAILURE_DECOMPOSITION",json.dumps({
        p:{
            "within_MSE_gain":v["South_America_within_vs_cross_dyad_calibration"]["by_pair_class"]["within_held_world"]["H_after_O_MSE_improvement_pct"],
            "cross_MSE_gain":v["South_America_within_vs_cross_dyad_calibration"]["by_pair_class"]["held_world_to_other"]["H_after_O_MSE_improvement_pct"],
            "cross_M3_minus_M2_pred":v["South_America_within_vs_cross_dyad_calibration"]["by_pair_class"]["held_world_to_other"]["mean_prediction_change_M3_minus_M2"]
        } for p,v in result["results"].items()},ensure_ascii=False),flush=True)
    print("ENTIRE_WORLD_TRANSFER_SUMMARY",json.dumps({
        p:{
          "world_equal":v["strict_entire_TDWG_Level1_world_holdout"]["world_equal_mean_MSE_gain_pct"],
          "n_worlds_improving":v["strict_entire_TDWG_Level1_world_holdout"]["number_of_eight_worlds_improving"],
          "world_8":v["strict_entire_TDWG_Level1_world_holdout"]["per_world"]["8"]
        } for p,v in result["results"].items()},ensure_ascii=False),flush=True)
    print("SOURCE_REPLAY_PASS",json.dumps({
        "n":result["sample"],
        "H_after_O_original":result["source"]["original_source_H_after_O_pct"],
        "proxies":{p:r["H_after_O_held_region_MSE_reduction_pct"] for p,r in result["results"].items()}
    },ensure_ascii=False))

if __name__=="__main__":
    main()
