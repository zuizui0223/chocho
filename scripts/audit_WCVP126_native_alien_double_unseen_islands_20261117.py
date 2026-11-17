#!/usr/bin/env python3
"""STRICT double-unseen island transfer diagnostic on source-matched WCVP pair flora.

Read the exact prior 5,355-pair matrix and original independent GSHHG marine
geometry + WPI source tables. Deterministically hash the 104 eligible island
codes into five bins, never based on floristic outcome, source WPI or sea
connectivity. Within each bin predict island pairs with BOTH endpoints held out,
using ONLY pairs with both endpoints in the OTHER bins as training.

The entire held fold has zero occurrence of any test-island as a training
endpoint. Cross-boundary island pairs are excluded from BOTH train/test in the
particular fold. Test model features use CURRENT marginal species richness of
held-out islands (conditional association, not prospective plant immigration).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 as src

SEED=20261117+4242
N_BIN=5
MODELS=("B","G","P","GP")
NBOOT=1999

def label(code):
    return int(hashlib.sha256(("chocho-WCVP126-double-unseen-v1|"+code).encode()).hexdigest()[:10],16)%N_BIN

def evaluate(d,flora,area_by):
    all_ids=sorted({v["island1"] for v in d}|{v["island2"] for v in d})
    grouped=defaultdict(set)
    for code in all_ids:grouped[label(code)].add(code)
    if len(grouped)!=N_BIN or min(map(len,grouped.values()))<12:
        raise RuntimeError(f"Deterministic hash fold size untenable: {[(k,len(v)) for k,v in grouped.items()]}")
    target=np.array([r[flora] for r in d],dtype=float)
    features=src.make_X(d,flora,area_by)
    rows=[]
    used_test=set()
    for fold in range(N_BIN):
        held=grouped[fold]
        trainidx=np.array([j for j,x in enumerate(d)
            if x["island1"] not in held and x["island2"] not in held],dtype=int)
        testidx=np.array([j for j,x in enumerate(d)
            if x["island1"] in held and x["island2"] in held],dtype=int)
        if len(trainidx)<2500 or len(testidx)<70:
            raise RuntimeError(f"Insufficient double-unseen source support fold={fold}, train={len(trainidx)} test={len(testidx)}")
        assert not any(d[j]["island1"] in held or d[j]["island2"] in held for j in trainidx)
        assert all(d[j]["island1"] in held and d[j]["island2"] in held for j in testidx)
        if set(map(int,testidx)) & used_test:
            raise RuntimeError("Test pair reused across held two-island folds")
        used_test.update(map(int,testidx))
        row={
            "fold":fold,"n_unseen_islands":len(held),
            "n_train_pairs":len(trainidx),"n_test_double_unseen_pairs":len(testidx),
            "held_island_codes":sorted(held)}
        for m in MODELS:
            pred=src.ridge(features[m][trainidx],target[trainidx],features[m][testidx])
            row["MSE_"+m]=float(np.mean((pred-target[testidx])**2))
        rows.append(row)
    n_test=sum(r["n_test_double_unseen_pairs"] for r in rows)
    b=np.array([[r["MSE_"+m] for m in MODELS] for r in rows])
    macro=b.mean(axis=0)
    weighted=np.average(b,weights=np.array([r["n_test_double_unseen_pairs"] for r in rows]),axis=0)
    def compare(v):
        B,G,P,GP=v
        return {
          "G_vs_B_pct":float((B-G)/B*100),
          "P_vs_B_pct":float((B-P)/B*100),
          "GP_vs_B_pct":float((B-GP)/B*100),
          "G_added_after_P_pct":float((P-GP)/P*100),
          "P_added_after_G_pct":float((G-GP)/G*100)}
    improvements=compare(macro)
    rng=np.random.default_rng(SEED+(1 if flora=="native" else 3))
    boot={k:[] for k in improvements}
    for _ in range(NBOOT):
        idx=rng.integers(0,N_BIN,N_BIN)
        q=compare(b[idx].mean(axis=0))
        for k,v in q.items():boot[k].append(v)
    return {
      "outcome":flora+"_Jaccard",
      "n_source_islands":len(all_ids),"n_source_pairs_with_both_Jaccard_defined":len(d),
      "n_double_unseen_test_pairs":n_test,
      "fraction_source_pairs_strictly_tested":n_test/len(d),
      "folds":rows,
      "macro_equal_fold_mse":dict(zip(MODELS,map(float,macro))),
      "pooled_pair_weighted_mse":dict(zip(MODELS,map(float,weighted))),
      "macro_equal_fold_relative_mse_improvement_percent":improvements,
      "pooled_pair_weighted_relative_mse_improvement_percent":compare(weighted),
      "descriptive_5fold_bootstrap95_not_independent_validation":{
        k:[float(np.quantile(v,.025)),float(np.quantile(v,.975))]
        for k,v in boot.items()},
      "n_folds_G_improved":sum(r["MSE_G"]<r["MSE_B"] for r in rows),
      "n_folds_P_improved":sum(r["MSE_P"]<r["MSE_B"] for r in rows),
      "test_islands_absent_from_train_endpoints":True,
      "no_crossfold_pairs_in_test":True,
      "no_folds_optimized_on_WCVP_response":True,
      "model_uses_current_richness_margins_of_unseen_islands":True
    }

def main():
    p=argparse.ArgumentParser()
    for k in ("plants","geometry","island-area","sea-audit","shipping-audit","out-dir"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    sea=json.loads(a.sea_audit.read_text())
    ais=json.loads(a.shipping_audit.read_text())
    if sea["source_valid_original_island_pair_geometric_shortest_sea"]!=5356:
        raise RuntimeError("Original 5356 geographic sea paths source changed")
    if ais["strict126_unordered_plant_pairs_with_one_or_more_actual_source_2020_port_OD"]!=7:
        raise RuntimeError("Observed 2020 AIS source pair support changed")
    data,unknown=src.input_pair(a)
    source_ids={v["island1"] for v in data}|{v["island2"] for v in data}
    if len(source_ids)!=104 or len(data)!=5355 or len(unknown)!=1:
        raise RuntimeError("Original WCVP 5355 common native/alien pairs, one unknown Jaccard, 104 islands source changed")
    area_data=src.rows(a.island_area)
    amap={r["code"]:float(r["area_km2"]) for r in area_data}
    if len(amap)!=126:raise RuntimeError("Original 126 source region area metadata changed")
    native=evaluate(data,"native",amap)
    alien=evaluate(data,"alien",amap)
    if not native["n_double_unseen_test_pairs"]==alien["n_double_unseen_test_pairs"]:
        raise RuntimeError("Strict cross native/alien source pair denominators diverged")
    z={
      "schema":"chocho_native_alien_wcvp_104_island_double_unseen_5hashfold_20261117_v1",
      "status":"STRICT_DOUBLE_UNSEEN_ISLAND_GEOMETRY_PORT_PAIR_CV_NOT_CURRENTS_OR_ACTUAL_SEED_IMPORT",
      "original_accepted_WCVP_native_alien_source_CI":44010011288,
      "original_GSHHG_only_sea_geometry_source_CI":44010673580,
      "original_source_verified_AIS2020_pair_support_CI":44010307551,
      "parent_original_leave_one_island_source_CI":44011154673,
      "source_common_islands":104,"source_common_pairs":5355,
      "excluded_original_one_Jaccard_undefined_pair":unknown,
      "fold_allocation":"SHA256 literal chocho-WCVP126-double-unseen-v1|TDWG3 code modulo5; no outcomes used",
      "test_pair_must_have_two_unseen_island_endpoints":True,
      "training_pairs_exclude_every_test_island_endpoint":True,
      "native_flora":native,"alien_flora":alien,
      "source_limits":[
       "G is same original GSHHG coastline-constrained UNDIRECTED weighted shortest water path, not actual directional OSCAR sea current or Circuitscape resistance.",
       "P is contemporaneous WPI port-count product, NOT measured vessel OD; only seven original 2020 shipping pair links are verified within the full 126 island set.",
       "Source similarity is contemporary native/introduced WCVP accepted flora; species richness from held island is treated as known at prediction time.",
       "A geographic CV using existing botanical records is not temporally independent ecological replication or causal evidence of sea currents introducing plants.",
       "The five hashed island folds are dependent via overlapping training data and sample only within-fold pairs; fold bootstrap confidence intervals are descriptive."
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"original_WCVP126_native_alien_double_unseen_5fold.json").write_text(
       json.dumps(z,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("WCVP126_DOUBLE_UNSEEN_NATIVE_ALIEN_CROSS",
       json.dumps({
         "status":z["status"],
         "n_source_pairs":len(data),
         "n_strict_test_pairs":native["n_double_unseen_test_pairs"],
         "native":native["macro_equal_fold_relative_mse_improvement_percent"],
         "alien":alien["macro_equal_fold_relative_mse_improvement_percent"],
         "native_descriptive_CI":native["descriptive_5fold_bootstrap95_not_independent_validation"],
         "alien_descriptive_CI":alien["descriptive_5fold_bootstrap95_not_independent_validation"],
       },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
