#!/usr/bin/env python3
"""Three-chain exact gain-margin null for alien flora / butterfly resource overlap.

Source-native bits, per-region additions, per-butterfly additions and the
Level1×island/mainland/mixed butterfly gain footprint stay fixed. Fixed WCVP
alien flora predictor and known-source 5fold held-both-islands transfer.
Exploratory source-circular calibration, NOT independent butterfly ecology.
"""
from __future__ import annotations
import argparse,csv,json,math,random,sys
from collections import Counter
from pathlib import Path
import numpy as np

import audit_island104_alien_flora_to_butterfly_resource_20261118 as bridge

SEEDS=(20261119,20261120,20261121)
DRAWS_PER_CHAIN=75
ATTEMPTS_PER_DRAW=72765
BURN_MULTIPLIER=5
E_OR_S=("E","S")
ORIGINAL={"E":26.36323764850138,"S":33.76782049915366}

def csvrows(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
def pair_sim(x,y):
    denom=(x|y).bit_count()
    if denom==0:raise RuntimeError("Source-current resource pair suddenly undefined after exact margins")
    return (x&y).bit_count()/denom

def setup(original_root,entries,area,region_class_file):
    lib=original_root/"scripts"
    sys.path.insert(0,str(lib))
    import analyze_global_resource_gain_fixed_margins as source
    root=original_root/"provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031"
    matrix=root/"beta_baseline_v3.matrix.csv.gz"
    receipt=root/"beta_baseline_v3.json"
    if bridge.sha(matrix)!=bridge.SOURCE_SHA or bridge.sha(receipt)!=bridge.REPORT_SHA:
        raise RuntimeError("Original 239×369 resource source SHA no longer pinned")
    codes,taxa,native,added,current,counts=source.read_matrix(matrix,receipt)
    if len(codes)!=369 or len(taxa)!=239 or counts!={"native":26530,"added":14553,"current":41083}:
        raise RuntimeError("Source 239×369 cell margins changed")
    lv=json.loads((original_root/"data/frozen/wgsrpd3_level1_lookup_v0.1.json").read_text())["units"]
    if len(lv)!=369 or set(lv)!=set(codes):
        raise RuntimeError("Source source Level1 mask invalid")
    cls={r["area_code_l3"]:r["class"] for r in csvrows(region_class_file)}
    if len(cls)!=369 or set(cls)!=set(codes) or Counter(cls.values())!={"island":126,"mainland":237,"mixed":6}:
        raise RuntimeError("Original island/mainland/mixed geography classification changed")
    strata=[str(lv[c])+"|"+cls[c] for c in codes]
    ix={c:i for i,c in enumerate(codes)}
    for e in entries:
        if e["island1"] not in ix or e["island2"] not in ix:
            raise RuntimeError("Unmatched 104-island resource-source code")
    source_ij=[(ix[e["island1"]],ix[e["island2"]]) for e in entries]
    expected_byregion,expected_byspecies=source.marginal_counts(native,added,len(taxa))
    if sum(expected_byregion)!=14553 or sum(expected_byspecies)!=14553:
        raise RuntimeError("Original resource gain margin source invalid")
    return source,native,added,codes,strata,source_ij,expected_byregion,expected_byspecies

def score_y(y,X,folds):
    results={}
    for name,(b,p) in X.items():
        b_err=[];p_err=[]
        for fold,tr,te in folds:
            pred_b=bridge.src.ridge(b[tr],y[tr],b[te])
            pred_p=bridge.src.ridge(p[tr],y[tr],p[te])
            b_err.append(float(np.mean((y[te]-pred_b)**2)))
            p_err.append(float(np.mean((y[te]-pred_p)**2)))
        b0=float(np.mean(b_err));p0=float(np.mean(p_err))
        if b0<=0:raise RuntimeError("Null baseline has zero MSE")
        results[name]={"gain":float(100*(b0-p0)/b0),
                       "fold_gains":list(map(float,[100*(x-y)/x for x,y in zip(b_err,p_err)]))}
    return results

def chain_diagnostic(x):
    arr=np.asarray(x,dtype=float);n=len(arr)
    m=arr.mean()
    err=arr-m
    var=float(np.dot(err,err)/n)
    if var==0:
        return {"lag1":None,"positive_pair_window_lags":0,"approx_ess":None}
    acf=[]
    for lag in range(1,min(n//2,32)):
        acf.append(float((np.dot(err[:-lag],err[lag:])/(n-lag))/var))
    sm=0.0;lags=0
    for k in range(0,len(acf)-1,2):
        pair=acf[k]+acf[k+1]
        if pair<=0:break
        sm+=pair;lags+=2
    tau=max(1.,1.+2*sm)
    return {"lag1":acf[0],"positive_pair_window_lags":lags,"approx_ess":float(n/tau),
            "first_half_mean":float(arr[:n//2].mean()),
            "second_half_mean":float(arr[n//2:].mean())}

def main():
    ap=argparse.ArgumentParser()
    for k in ("archived-root","plants","geometry","island-area","region-class","out-dir"):
        ap.add_argument("--"+k,type=Path,required=True)
    a=ap.parse_args()
    e,area,original=bridge.setup(a.archived_root,a.plants,a.geometry,a.island_area)
    if original["source_pairs"]!=5355 or len(e)!=5353:
        raise RuntimeError("Original 104 island 5355 botanical/5353 resource-defined pair source drift")
    if bridge.src.RIDGE_ALPHA!=20:
        raise RuntimeError("Original ridge alpha changed")
    folds=bridge.fold_cohort(e,native_only=False)
    if len(folds)!=5 or sum(len(te) for _,_,te in folds)!=1088:
        raise RuntimeError("Original double-unseen 1088 test pair source changed")
    X={
        name:(bridge.features(e,area,include_current_margins=(name=="S"),include_alien_similarity=False),
              bridge.features(e,area,include_current_margins=(name=="S"),include_alien_similarity=True))
        for name in E_OR_S
    }
    y_obs=np.array([r["response"] for r in e])
    observed=score_y(y_obs,X,folds)
    for name in E_OR_S:
        if abs(observed[name]["gain"]-ORIGINAL[name])>1e-8:
            raise RuntimeError("Source bridge gain not exactly reproduced")
    source,native,added,regions,strata,indices,regional,species=setup(
        a.archived_root,e,area,a.region_class)
    # For each source dyad, fixed plant composition and botanical geography;
    # only potential 239-consumer resource *placements* change.
    x_alien=np.array([r["alien_flora_Jaccard"] for r in e])
    summary={"schema":"chocho_original239_flora_to_resource_104islands_three_chain_null_v1",
             "status":"SOURCE_LOCKED_EXPLORATORY_WCVP_SHARED_SOURCE_MARGIN_NULL",
             "protocol":"results/island104_alien_flora_resource_margin_null_20261118/PROTOCOL.md",
             "source_pr118_commit":"bb4d8f00fbdb97fa2d90f606a8c9e7e0620646f2",
             "source_original_matrix_sha256":bridge.SOURCE_SHA,
             "source_original104_plant_pairs":5355,
             "source_resource_defined_pairs":len(e),"source_double_unseen_pairs":1088,
             "source_n_butterflies":239,"source_resource_counts":original["source_units"],
             "fixed_source_baseline_ES":ORIGINAL,
             "original_es_model_scores":observed,
             "seeds":list(SEEDS),"draws_per_chain":DRAWS_PER_CHAIN,
             "burnin_proposals_per_chain":BURN_MULTIPLIER*ATTEMPTS_PER_DRAW,
             "proposals_per_retained_draw":ATTEMPTS_PER_DRAW,
             "null_holds_native_and_added_region_species_and_Level1_island_class_margins":True,
             "chains":[]}
    combined=[]
    for chain_no,seed in enumerate(SEEDS):
        sampled=added.copy()
        edges=source.additions_list(sampled)
        rng=random.Random(seed)
        burn_accept=source.swap_new_resources(native,sampled,edges,rng,
                       BURN_MULTIPLIER*ATTEMPTS_PER_DRAW,region_strata=strata)
        if burn_accept<1000:
            raise RuntimeError("Cannot mix original source eligible 2×2 margin swaps")
        rows=[]
        for draw in range(DRAWS_PER_CHAIN):
            accepted=source.swap_new_resources(native,sampled,edges,rng,
                         ATTEMPTS_PER_DRAW,region_strata=strata)
            if accepted<200:raise RuntimeError("Too few original constrained resource swaps")
            now_region,now_species=source.marginal_counts(native,sampled,len(species))
            if now_region!=regional or now_species!=species:
                raise RuntimeError("Original resource gain count preservation failed")
            if any(n&a for n,a in zip(native,sampled)):
                raise RuntimeError("A null gain overwrote a native known-host resource")
            cur=[n|x for n,x in zip(native,sampled)]
            target=np.array([pair_sim(cur[i],cur[j]) for i,j in indices],float)
            fits=score_y(target,X,folds)
            rr={"chain":chain_no+1,"seed":seed,"draw":draw+1,
                "accepted":accepted,
                "gain_E":fits["E"]["gain"],
                "gain_S":fits["S"]["gain"],
                "raw_spearman":bridge.corr_spearman(x_alien,target)}
            rows.append(rr)
            combined.append(rr)
        ss={
            "chain":chain_no+1,"seed":seed,"accepted_burnin_swaps":burn_accept,
            "retained_draws":len(rows),"accepted_swaps_min":min(x["accepted"] for x in rows),
            "accepted_swaps_median":float(np.median([x["accepted"] for x in rows])),
            "gain_E_mean":float(np.mean([x["gain_E"] for x in rows])),
            "gain_S_mean":float(np.mean([x["gain_S"] for x in rows])),
            "gain_S_diagnostic":chain_diagnostic([x["gain_S"] for x in rows])
        }
        summary["chains"].append(ss)
        print("FIXED_SOURCE_MARGIN_NULL_CHAIN_FINISHED",json.dumps(ss),flush=True)
    if len(combined)!=len(SEEDS)*DRAWS_PER_CHAIN:
        raise RuntimeError("Expected source resource marginal null draw number changed")
    results={}
    for name in E_OR_S:
        observed_gain=observed[name]["gain"]
        vals=np.array([z["gain_"+name] for z in combined])
        results[name]={
            "observed_gain_pct":observed_gain,
            "null_median_gain_pct":float(np.median(vals)),
            "null_mean_gain_pct":float(np.mean(vals)),
            "null_q025_gain_pct":float(np.quantile(vals,.025)),
            "null_q975_gain_pct":float(np.quantile(vals,.975)),
            "null_min_gain_pct":float(vals.min()),
            "null_max_gain_pct":float(vals.max()),
            "observed_minus_null_median_gain_points":float(observed_gain-np.median(vals)),
            "descriptive_high_tail":float((1+int(np.sum(vals>=observed_gain)))/(1+len(vals))),
            "n_null_scores":len(vals)
        }
    summary["observed_vs_null"]=results
    summary["limits"]=[
        "All 3 source-randomized chains are simulated under a post-observed-results design and share one WCVP/HOSTS panel; retained states are serially correlated.",
        "2×2 switches preserve butterfly/region counts and island/Level1 strata but not plant phylogeny, family matching, climatic similarity or actual species dispersal.",
        "The hypothetical resource set is computed from the SAME contemporary WCVP plant distributions and fixed HOSTS; a positive excess is neither independent trophic evidence nor a direct mechanism.",
        "The S predictor uses test-island contemporary butterfly resource counts, which are known only from source-derived response margins; S is retrospective.",
        "GBIF occurrence and endemics' true larval use are not measured by the null."
    ]
    a.out_dir.mkdir(parents=True,exist_ok=True)
    summary_path=a.out_dir/"margin_null_flora_resource_source_summary.json"
    summary_path.write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    with (a.out_dir/"margin_null_225_source_scores.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(combined[0]));w.writeheader();w.writerows(combined)
    print("WCVP104_ALIEN_FLORA_RESOURCE_EXACT_GAIN_MARGIN_NULL",json.dumps({
        "observed_vs_null":results,"source_pairs":len(e),"strict_pairs":1088,
        "seeds":list(SEEDS),"chains":len(summary["chains"])
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
