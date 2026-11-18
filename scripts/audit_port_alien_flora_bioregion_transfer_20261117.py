#!/usr/bin/env python3
"""Hold out every TDWG Level1 region to test port×island-context predictive value.

A post-hoc OUT-OF-REGION falsification of the single-paper WCVP botanical
port/moderation claims. This is cross-sectional botanical geography (no
seed dispersal, vessel OD, larval host performance, butterflies or colonisation).

Evaluate separate models with and without TWO port×standardized island context
interaction terms. Because a totally unseen botanical Level1 has no estimable
fixed-effect intercept, BOTH trained models exclude Level1 dummies by design.
They receive the SAME train-only normalization and covariates.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import audit_island_port_alien_flora_native_area_moderators_20261117 as source

SEED=20261117+3001
N_BOOT=2499

def get(rows,key):
    return np.array([float(r[key]) for r in rows])

def arrays(train,test,*,flexible,effort_key=None):
    # Avoid leakage of training normalization from test-geography responses or
    # test-geography covariate distributions.
    def normalize(col):
        lo=get(train,col)
        mu=float(np.mean(lo));sd=float(np.std(lo,ddof=1))
        if sd<1e-10:raise RuntimeError(f"Cannot normalize train-only {col} in TDWG holdout")
        return (lo-mu)/sd,(get(test,col)-mu)/sd
    pa,pe=normalize("ln2ports")
    aa,ae=normalize("lnarea")
    na,ne=normalize("lnnative")
    a_train=[np.ones(len(train)),pa,aa,na]
    a_test=[np.ones(len(test)),pe,ae,ne]
    names=["intercept","port","area","native"]
    if flexible:
        a_train += [pa*pa,aa*aa,na*na,aa*na]
        a_test += [pe*pe,ae*ae,ne*ne,ae*ne]
        names += ["port2","area2","native2","area_native"]
    if effort_key:
        eff_a,eff_e=normalize(effort_key)
        a_train.append(eff_a);a_test.append(eff_e);names.append(effort_key)
    x_train=np.column_stack(a_train);x_test=np.column_stack(a_test)
    # Only extra columns are predeclared 2 moderators:
    # log(ports)×log(area), log(ports)×log(native richness).
    x_train_mod=np.column_stack([x_train,pa*aa,pa*na])
    x_test_mod=np.column_stack([x_test,pe*ae,pe*ne])
    if np.linalg.matrix_rank(x_train)<x_train.shape[1] or np.linalg.matrix_rank(x_train_mod)<x_train_mod.shape[1]:
        raise RuntimeError("Leave-region-out train feature matrix rank deficient")
    return x_train,x_test,x_train_mod,x_test_mod,names

def beta(X,Y):
    return np.linalg.lstsq(X,Y,rcond=None)[0]

def bootstrap(level_scores,rng):
    levs=list(level_scores)
    a=np.array([level_scores[l]["mse_base"] for l in levs])
    b=np.array([level_scores[l]["mse_with_moderators"] for l in levs])
    out=np.empty(N_BOOT,dtype=float)
    for i in range(N_BOOT):
        sel=rng.integers(0,len(levs),size=len(levs))
        m0=float(np.mean(a[sel]));m1=float(np.mean(b[sel]))
        out[i]=(m0-m1)/m0*100.
    return [float(np.quantile(out,.025)),float(np.quantile(out,.975))]

def one_cohort(rows,name,*,flexible,effort_key=None):
    levels=sorted({r["level1"] for r in rows},key=int)
    if len(levels)!=9:raise RuntimeError(f"Original 9 botanical biogeographical groups not retained in {name}: {levels}")
    scored=[]
    group={}
    for lev in levels:
        tr=[r for r in rows if r["level1"]!=lev]
        te=[r for r in rows if r["level1"]==lev]
        if len(tr)<65 or not te:raise RuntimeError("Too few out-of-region training or held examples")
        y_train=get(tr,"lnalien")
        y_test=get(te,"lnalien")
        x0,xe,x1,xex,names=arrays(tr,te,flexible=flexible,effort_key=effort_key)
        pb=xe@beta(x0,y_train)
        pm=xex@beta(x1,y_train)
        if not (np.isfinite(pb).all() and np.isfinite(pm).all()):
            raise RuntimeError("Unstable OOS predictions")
        e0=(pb-y_test)**2
        e1=(pm-y_test)**2
        scored.extend([{"region":lev,"code":r["code"],"squared_base":float(q0),
            "squared_moderated":float(q1)} for r,q0,q1 in zip(te,e0,e1)])
        group[lev]={
          "n_islands":len(te),
          "mse_base":float(np.mean(e0)),
          "mse_with_moderators":float(np.mean(e1)),
          "pct_improvement_if_positive":float((np.mean(e0)-np.mean(e1))/np.mean(e0)*100),
          "is_better_moderated":bool(np.mean(e1)<np.mean(e0))
        }
    scores=np.array([[r["squared_base"],r["squared_moderated"]] for r in scored])
    pooled=np.mean(scores,axis=0)
    macro=np.mean([[v["mse_base"],v["mse_with_moderators"]] for v in group.values()],axis=0)
    rng=np.random.default_rng(SEED+len(rows)+len(name))
    return {
      "cohort":name,"n_islands":len(rows),"n_held_out_TDWG_Level1":len(levels),
      "main_effects":"flexible_quadratic_port_area_native" if flexible else "linear_port_area_native",
      "effort_nuisance":effort_key,
      "held_out_region_fixed_effects_intentionally_not_in_model":True,
      "training_norm_source":"only training Level1 islands in each fold, never the held out region",
      "response":"ln(1+WCVP accepted introduced flora)",
      "predictors_base":"WPI port count, island area, accepted native plant richness"+(", nonlinear squares+area×native" if flexible else "")+(", "+effort_key if effort_key else ""),
      "predictors_added":"two port×standardized area/native interactions, fit only on train source",
      "pooled_test_mse_no_interaction":float(pooled[0]),
      "pooled_test_mse_both_interactions":float(pooled[1]),
      "pooled_test_relative_MSE_improvement_pct":float((pooled[0]-pooled[1])/pooled[0]*100),
      "macro_equal_level1_MSE_no_interaction":float(macro[0]),
      "macro_equal_level1_MSE_with_both_interactions":float(macro[1]),
      "macro_equal_level1_relative_MSE_improvement_pct":float((macro[0]-macro[1])/macro[0]*100),
      "descriptive_level1_resample_95pct_relative_MSE_improvement_pct":bootstrap(group,rng),
      "n_Level1_regions_improved":sum(v["is_better_moderated"] for v in group.values()),
      "level1_results":group,
      "all_test_islands_scored_once":len(scored),
      "limits":[
        "Nonrandom botanical geography and WPI historical association; test is out-of-region prediction, not new biological field validation.",
        "Nine Level1 held areas share all remaining Level1s in their training data, hence folds are not independent replication.",
        "Original main regressions included Level1 fixed effects, which cannot be used unchanged to extrapolate to an unseen Level1: both OOS models intentionally omit that FE.",
        "A positive OOS improvement supports associational regional transfer of a moderation term only; it does not identify seed transport, host quality or butterfly establishment.",
        "The bootstrap resamples just 9 held-out Level1 strata and is descriptive, not a prospective inference."
      ]
    }

def main():
    a=argparse.ArgumentParser()
    for n in ("flora","hosts","effort","out-dir"):a.add_argument("--"+n,required=True,type=Path)
    v=a.parse_args()
    for n,p in (("flora",v.flora),("hosts",v.hosts),("effort",v.effort)):
        if source.gitblob(p)!=source.SOURCE_SHA[n]:
            raise RuntimeError("Original Git-blob source changed "+n)
    full=source.vec_data(source.rows(v.flora),source.rows(v.hosts),source.rows(v.effort))
    hosts=[r for r in full if r["has_host"]]
    efforts=[r for r in hosts if r["has_effort"]]
    candidates=[
      ("all126_linear",full,False,None),
      ("all126_flexible",full,True,None),
      ("matched123_linear",hosts,False,None),
      ("matched123_flexible",hosts,True,None),
      ("botanical98_flexible_current_effort",efforts,True,"lneffort"),
      ("botanical98_flexible_pre1900_effort",efforts,True,"lnpre1900effort"),
    ]
    res=[]
    for name,items,flex,effort in candidates:
        q=one_cohort(items,name,flexible=flex,effort_key=effort)
        res.append(q)
        print("PORT_ALIEN_FLORA_OOS_LEVEL1",json.dumps({
            "name":name,"n":len(items),"region_macro_pct":q["macro_equal_level1_relative_MSE_improvement_pct"],
            "pooled_pct":q["pooled_test_relative_MSE_improvement_pct"],
            "regions_improved":q["n_Level1_regions_improved"],
            "descriptive95":q["descriptive_level1_resample_95pct_relative_MSE_improvement_pct"]
        }),flush=True)
    report={
      "schema":"chocho_123_126_98_WCVP_port_area_native_held_TDWG_Level1_predictive_check_v1",
      "status":"SOURCE_VERIFIED_OUT_OF_BIOREGION_MODERATOR_PREDICTION_DIAGNOSTIC_NOT_CAUSAL",
      "source_pr173_CI":44003165119,
      "source_nonlinear_historical_CI":44005033127,
      "source_Git_blob_shas":source.SOURCE_SHA,
      "no_fold_uses_heldout_TDWG_Level1_as_training":True,
      "sampling_unit":"WCVP native/introduced flora × WGSRPD3 island region, not butterfly or voyage",
      "source_samples":[126,123,98],
      "n_test_variants":len(res),
      "results":res,
      "limits":[
       "Cross-regional split checks model association transfer and does not test whether the historic WPI port caused introduction.",
       "All model classes are post-hoc decided after PR173 discovery, no tuning or selective hold-out choice based on test sign.",
       "The same present-day WCVP flora and current WPI ports underlie all six comparisons; do not cite them as six biological replications.",
       "A 9-biogeographic-group bootstrap is a descriptive uncertainty display, not an independent confirmatory experiment.",
       "Neither larval host feeding, the date of plant establishment nor butterfly range expansion is evaluated in this response."
      ]
    }
    v.out_dir.mkdir(parents=True,exist_ok=True)
    (v.out_dir/"held_biogeographic_level1_moderation_prediction.json").write_text(
      json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("PORT_ALIEN_FLORA_OOS_SOURCE_LOCKED_GATE",json.dumps({
       "n_variants":len(res),"n_samples":report["source_samples"],
       "macro_improvement_pct":{q["cohort"]:q["macro_equal_level1_relative_MSE_improvement_pct"] for q in res}
    }),flush=True)

if __name__=="__main__":main()
