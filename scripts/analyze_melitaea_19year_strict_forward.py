#!/usr/bin/env python3
"""Strict year-forward, predeclared host-loss prediction (not causal inference).

Training for test year t includes ONLY earlier index years k<t with outcomes
observed by the start of test year t. No future data train the test-year model.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline

import audit_melitaea_19year_fallback_support as source
import analyze_melitaea_19year_host_fallback_prediction as old

FIRST_YEAR=2010
LAST_YEAR=2017

def temporal_predictions(panel,increment=False,start=FIRST_YEAR,end=LAST_YEAR):
    arr=np.array([[np.nan if v is None else v for v in r["base"]]+
                  ([r["interaction"]] if increment else []) for r in panel],dtype=float)
    y=np.array([r["y"] for r in panel],dtype=int)
    year=np.array([r["year"] for r in panel],dtype=int)
    pred=np.full(len(panel),np.nan)
    train_counts={}
    for test_year in range(start,end+1):
        test=year==test_year
        train=year<test_year
        if not test.any() or len(set(y[train]))!=2:
            raise ValueError(f"invalid prequential calendar year {test_year}")
        if np.any(year[train]>=test_year):
            raise ValueError("future information leaked into training")
        mod=make_pipeline(SimpleImputer(strategy="median",add_indicator=True,keep_empty_features=True),
                          StandardScaler(),
                          LogisticRegression(C=1,solver="lbfgs",max_iter=1000))
        mod.fit(arr[train],y[train])
        pred[test]=mod.predict_proba(arr[test])[:,1]
        train_counts[str(test_year)]=int(train.sum())
    return np.clip(pred,1e-9,1-1e-9),train_counts

def evaluate(panel):
    primary=[r for r in panel if FIRST_YEAR<=r["year"]<=LAST_YEAR]
    if not primary:
        raise ValueError("no forward risk set")
    ap,train_count=temporal_predictions(panel,False)
    bp,_=temporal_predictions(panel,True)
    years=np.array([r["year"] for r in panel])
    select=(years>=FIRST_YEAR)&(years<=LAST_YEAR)
    p=ap[select];q=bp[select]
    y=np.array([r["y"] for r in primary])
    if not (np.isfinite(p).all() and np.isfinite(q).all()):
        raise ValueError("missing prospective prediction")
    losses=lambda pred:-(y*np.log(pred)+(1-y)*np.log(1-pred))
    diff=losses(p)-losses(q)
    yr=[r["year"] for r in primary]
    per_year={str(year):{"n":int(sum(t==year for t in yr)),
                        "outcome_zero_count":int(sum(r["y"] for r in primary if r["year"]==year)),
                        "paired_logloss_gain":float(diff[np.array(yr)==year].mean())}
              for year in range(FIRST_YEAR,LAST_YEAR+1)}
    by_patch=old.cluster_interval(diff,[r["patch"] for r in primary],seed=20261010,iterations=2000)
    by_year=old.cluster_interval(diff,yr,seed=20261011,iterations=2000)
    good_years=sum(a["paired_logloss_gain"]>0 for a in per_year.values())
    gain=float(np.mean(diff))
    passed=(gain>0 and by_patch[0]>0 and by_year[0]>0 and good_years>=6)
    cells={}
    for decline in (0,1):
        for high in (0,1):
            rows=[r for r in primary if r["decline"]==decline and r["backup_high"]==high]
            cells[f"decline_{decline}_backup_{high}"]={
                "n":len(rows),"next_year_no_nest_count":sum(r["y"] for r in rows)}
    return {
      "test_years":[FIRST_YEAR,LAST_YEAR],"forward_test_n":len(primary),
      "forward_test_patches":len({r["patch"] for r in primary}),
      "training_rows_by_year":train_count,
      "primary_outcome":"next autumn no larval nests (not proven local extinction)",
      "baseline_logloss":float(np.mean(losses(p))),
      "augmented_logloss":float(np.mean(losses(q))),
      "forward_logloss_gain":gain,
      "forward_brier_gain":float(np.mean((p-y)**2-(q-y)**2)),
      "patch_cluster_95ci":by_patch,
      "year_cluster_95ci":by_year,
      "positive_test_years":good_years,"per_year":per_year,"cell_support":cells,
      "frozen_forward_gate_pass":bool(passed),
      "decision":("FORWARD_PREDICTION_ASSOCIATION_ONLY" if passed
                  else "NO_ROBUST_FORWARD_HOST_BACKUP_GAIN"),
      "retrospective_host_memory_or_causality_identified":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    tsv,net,sha=old.get_source()
    data,metadata=source.ingest(tsv)
    gate=source.support(data)
    if not gate["data_support_pass"]:
        raise ValueError("original source support gate failed")
    panel=old.build_panel(data,old.read_static(net))
    if len(panel)!=gate["three_year_consecutive_occupied_risk_set"]:
        raise ValueError("panel not equal to predeclared risk set")
    effect=evaluate(panel)
    result={
      "schema":"chocho_melitaea_19year_strict_forward_v01",
      "original_source":"Schulz et al 2020 Zenodo 4987060",
      "source_sha256":sha,
      "protocol":"docs/exploratory/MELITAEA_19YEAR_STRICT_FORWARD_VALIDATION_V01.json",
      "source_gate":gate,
      "strict_forward":effect,
      "host_removal_not_randomized":True,
      "botanical_introduction_not_tested":True,
      "GEB_PR38_unchanged":True
    }
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
