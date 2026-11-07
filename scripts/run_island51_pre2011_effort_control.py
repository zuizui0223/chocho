#!/usr/bin/env python3
"""Post-hoc butterfly host/circuit sensitivity to PRE-2011 all-plant *recording effort*."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from exploratory_island51_circuit_network_holdout import source, train_predict
from run_island51_host_circuit_joint import (
    INPUT_PATHS, TARGETS, load_original_port_counts, load_resistance_matrix,
    make_joint_features, validate_input_lock,
)
from run_island51_observed_ocean_circuits import sha256

PARENT = ROOT / "results/island51_host_circuit_joint_20261106"
ATLAS_LOCK = ROOT / "results/island51_host_circuit_atlas_20261106/input_lock.json"
ROSTER = ROOT / INPUT_PATHS["islands"]
PAIRS = ROOT / INPUT_PATHS["pairs"]
PARENT_FOLDS = PARENT / "model_fold_losses.csv"
MODELS = ("CR","HR","HCR","CRE","HRE","HCRE")
COMPARISONS = (
    ("CR_to_HCR","CR","HCR","prior_host_identity_increment"),
    ("HR_to_HCR","HR","HCR","prior_marine_circuit_increment"),
    ("CR_to_CRE","CR","CRE","early_botanical_effort_beyond_circuits"),
    ("CRE_to_HCRE","CRE","HCRE","primary_identity_beyond_historical_effort"),
    ("HR_to_HRE","HR","HRE","early_botanical_effort_beyond_hosts"),
    ("HRE_to_HCRE","HRE","HCRE","secondary_circuits_beyond_historical_effort"),
    ("HCR_to_HCRE","HCR","HCRE","early_botanical_effort_beyond_joint"),
)

HIST_DIR = ROOT / "provenance/reviewer_defenses/results/island_dated_host_match_20261030"
HIST_GIT_BLOB_SHA = {"island_port_early_inventory_status.csv":"30ffb98ac91aa3cc56451f5202000bb10ccdb023","island_port_early_plant_availability.csv":"0a01831056c2cd59a4bd41bd1dbd4c84d77cf922","island_port_early_inventory_audit.json":"95a8ea7dc2096611837ebbd66aa27025b7fc9c96"}
PARENT_RICHNESS_MSE = {
    "new_butterfly_jaccard":{"CR":0.005849988979,"HR":0.005954970620,"HCR":0.005716559115},
    "earlier_butterfly_jaccard":{"CR":0.010607946157,"HR":0.010007384512,"HCR":0.010104237739},
}
PARENT_RICHNESS_MAE = {
    "new_butterfly_jaccard":{"CR":0.053338303691,"HR":0.053883752394,"HCR":0.052531444456},
    "earlier_butterfly_jaccard":{"CR":0.068154599964,"HR":0.066458025189,"HCR":0.066369520658},
}

def git_blob_sha(raw):
    import hashlib
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def early_botanical_features(q, ids, parent_ports):
    file_receipts = {}
    for name,expected in HIST_GIT_BLOB_SHA.items():
        path=HIST_DIR/name
        raw=path.read_bytes()
        if git_blob_sha(raw)!=expected:
            raise ValueError("Historical source file changed: "+name)
        file_receipts[name]={"git_blob_sha":expected,"sha256":sha256(path)}
    state=pd.read_csv(HIST_DIR/"island_port_early_inventory_status.csv",dtype={"landmass_id":str})
    inventory=pd.read_csv(HIST_DIR/"island_port_early_plant_availability.csv",dtype={"landmass_id":str})
    audit=json.loads((HIST_DIR/"island_port_early_inventory_audit.json").read_text())
    expected={"status":"COMPLETE_DATED_INVENTORY","planned_ports":194,"complete_ports":194,
              "live_query_completed_ports":141,"cached_zero_ports":53,
              "completely_resolved_islands":51,"changed_count_ports":0}
    if any(audit.get(k)!=v for k,v in expected.items()):
        raise ValueError("Early botanical audit scope changed")
    if len(state)!=194 or state.port_id.nunique()!=194 or len(inventory)!=194 or inventory.port_id.nunique()!=194:
        raise ValueError("Historical port cohort incomplete")
    merged=state.merge(inventory,on=["port_id","landmass_id"],validate="one_to_one",
                       suffixes=("","_plan"),indicator=True)
    if len(merged)!=194 or not merged["_merge"].eq("both").all():
        raise ValueError("Historical original port join failed")
    if not merged.complete.astype(str).str.lower().eq("true").all():
        raise ValueError("Historical original port query incomplete")
    if merged.status.value_counts().to_dict()!={"COMPLETE_EARLY_FACETS":141,
                                                "CACHED_ZERO_EARLY_RECORDS":53}:
        raise ValueError("Historical 141/53 source status changed")
    if not ((merged.count_change==0).all() and
            (merged.early_plant_records==merged.early_plant_records_plan).all() and
            (merged.original_cached_count==merged.early_plant_records).all()):
        raise ValueError("Historical count cache mismatch")
    if ((merged.early_plant_records<0).any() or
        (merged.loc[merged.status.eq("CACHED_ZERO_EARLY_RECORDS"),"early_plant_records"]!=0).any()):
        raise ValueError("Invalid historical zero-record evidence")
    grouped=merged.groupby("landmass_id").agg(
        port_count=("port_id","count"),early_record_sum=("early_plant_records","sum"),
        ports_with_record=("early_plant_records",lambda s:int((s>0).sum())))
    if set(grouped.index)!=set(ids) or {k:int(grouped.loc[k,"port_count"]) for k in ids}!=parent_ports:
        raise ValueError("Historical port grouping differs from parent")
    counts={k:int(grouped.loc[k,"early_record_sum"]) for k in ids}
    values=np.log1p(np.array([counts[k] for k in ids],dtype=float))
    left,right=values[q.i.to_numpy(dtype=int)],values[q.j.to_numpy(dtype=int)]
    features={"sum_log1p_2001_2010_allplant_buffer_records":left+right,
              "abs_difference_log1p_2001_2010_allplant_buffer_records":np.abs(left-right)}
    coverage=[{"island_id":k,"original_port_count":int(grouped.loc[k,"port_count"]),
               "early_plant_buffer_record_count_sum":counts[k],
               "ports_with_early_plant_records":int(grouped.loc[k,"ports_with_record"])} for k in ids]
    receipt={"islands_with_early_plant_records":sum(v>0 for v in counts.values()),
             "islands_without_early_plant_records":sum(v==0 for v in counts.values()),
             "sum_port_buffer_early_plant_records":sum(counts.values()),
             "files":file_receipts,"zero_record_not_absence":True,
             "overlapping_port_buffers_can_duplicate_observations":True,
             "historical_host_species_identities_not_included":True}
    return features,coverage,receipt




def host_counts(path, ids):
    frame = pd.read_csv(path, dtype={"landmass_id": str})
    eligible = frame.loc[frame.eligible.astype(str).str.lower().eq("true")]
    if len(eligible) != 51 or eligible.landmass_id.nunique() != 51 or set(eligible.landmass_id) != set(ids):
        raise ValueError("Host counts do not cover the fixed 51 islands")
    numbers = {}
    for row in eligible.itertuples(index=False):
        x = float(row.alien_host_species)
        if not np.isfinite(x) or x < 0 or x != int(x):
            raise ValueError("Nonnegative integer host counts required")
        numbers[str(row.landmass_id)] = int(x)
    return numbers


def richness_features(q, ids, counts):
    if set(counts) != set(ids):
        raise ValueError("Host richness island identity mismatch")
    log = np.log1p(np.array([counts[k] for k in ids], dtype=float))
    l, r = log[q.i.to_numpy(dtype=int)], log[q.j.to_numpy(dtype=int)]
    values = {
        "sum_log1p_candidate_host_taxa": l + r,
        "abs_difference_log1p_candidate_host_taxa": np.abs(l - r),
    }
    if any(v.shape != (len(q),) or not np.isfinite(v).all() for v in values.values()):
        raise ValueError("Unaligned richness predictors")
    return values


def compare_folds(folds):
    if len(folds) != len(TARGETS) * len(MODELS) * 51 or folds.duplicated(["target", "model", "heldout_island_id"]).any():
        raise ValueError("Incomplete 51-island model-fold support")
    rows = []
    for target in TARGETS:
        f = folds.loc[folds.target.eq(target)]
        for name, ref, aug, role in COMPARISONS:
            a = f.loc[f.model.eq(ref)].set_index("heldout_island_id")
            b = f.loc[f.model.eq(aug)].set_index("heldout_island_id")
            if len(a) != 51 or len(b) != 51 or set(a.index) != set(b.index):
                raise ValueError("Unmatched fold identities")
            b = b.reindex(a.index)
            g = a.MSE.to_numpy() - b.MSE.to_numpy()
            g_mae = a.MAE.to_numpy() - b.MAE.to_numpy()
            x, y = float(a.MSE.mean()), float(b.MSE.mean())
            xa, ya = float(a.MAE.mean()), float(b.MAE.mean())
            rows.append({
                "target": target, "comparison": name, "role": role,
                "reference_model": ref, "augmented_model": aug,
                "islands": 51, "heldout_endpoint_predictions": 2550,
                "MSE_reference": x, "MSE_augmented": y,
                "MSE_reduction_percent": 100 * (x-y) / x,
                "islands_mse_improved": int(np.sum(g > 0)),
                "islands_mse_worse": int(np.sum(g < 0)),
                "islands_mse_tied": int(np.sum(g == 0)),
                "MAE_reference": xa, "MAE_augmented": ya,
                "MAE_reduction_percent": 100 * (xa-ya) / xa,
                "islands_mae_improved": int(np.sum(g_mae > 0)),
                "islands_mae_worse": int(np.sum(g_mae < 0)),
            })
    return pd.DataFrame(rows)


def run(out_dir):
    out = Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT / p).resolve()) for p in
           ("data", "provenance", "results/island51_host_circuit_joint_20261106",
            "results/island51_host_circuit_atlas_20261106")):
        raise ValueError("Cannot overwrite source, parent, or atlas artifacts")
    validated = validate_input_lock(PARENT / "input_lock.json")
    atlas = json.loads(ATLAS_LOCK.read_text())
    pin = {item["path"]: item["sha256"] for item in atlas["files"]}
    relative_fold = str(PARENT_FOLDS.relative_to(ROOT))
    if pin.get(relative_fold) != sha256(PARENT_FOLDS):
        raise ValueError("Parent fold checksum changed from atlas-pinned version")
    q, ids, _ = source(PAIRS, ROSTER)
    if len(ids) != 51 or len(q) != 1275:
        raise ValueError("Original island/pair support changed")
    prep = json.loads((ROOT / INPUT_PATHS["preparation"]).read_text())
    ports = load_original_port_counts(ROSTER, ids, prep)
    counts = host_counts(ROSTER, ids)
    richness = richness_features(q, ids, counts)
    early,early_coverage,early_receipt=early_botanical_features(q,ids,ports)
    matrices = {key: load_resistance_matrix(ROOT / INPUT_PATHS[key], ids) for key in ("uniform", "ais", "wind")}
    parent = pd.read_csv(PARENT_FOLDS, dtype={"heldout_island_id": str}).set_index(
        ["target", "model", "heldout_island_id"], verify_integrity=True)
    folds, predictions, features_by_target = [], [], {}
    max_mse = max_mae = 0.0
    ii, jj = q.i.to_numpy(dtype=int), q.j.to_numpy(dtype=int)
    for target in TARGETS:
        old = make_joint_features(q, ids, ports, matrices, target)
        cr=dict(old["C"],**richness)
        hr=dict(old["H"],**richness)
        hcr=dict(old["HC"],**richness)
        features={"CR":cr,"HR":hr,"HCR":hcr,
                  "CRE":dict(cr,**early),"HRE":dict(hr,**early),
                  "HCRE":dict(hcr,**early)}
        features_by_target[target] = {m: list(features[m]) for m in MODELS}
        design = {m: np.column_stack(list(features[m].values())) for m in MODELS}
        y = q[target].to_numpy(dtype=float)
        for i, island in enumerate(ids):
            test = (ii == i) | (jj == i)
            train = ~test
            if (int(test.sum()), int(train.sum())) != (50, 1225):
                raise ValueError("Invalid leave-one-island-out fold")
            for oldmodel in ("C","H","HC"):
                original_design=np.column_stack(list(old[oldmodel].values()))
                original_pred=train_predict(original_design,y,train,test)
                p=parent.loc[(target,oldmodel,island)]
                dx=abs(float(np.mean((y[test]-original_pred)**2))-float(p.MSE))
                dy=abs(float(np.mean(np.abs(y[test]-original_pred)))-float(p.MAE))
                if dx>1e-12 or dy>1e-12:
                    raise ValueError("Parent C/H/HC fold no longer reproduces")
                max_mse,max_mae=max(max_mse,dx),max(max_mae,dy)
            for m in MODELS:
                estimate = train_predict(design[m], y, train, test)
                obs = y[test]
                if len(estimate) != 50 or not np.isfinite(estimate).all():
                    raise ValueError("Invalid prediction")
                mse = float(np.mean((obs-estimate)**2))
                mae = float(np.mean(np.abs(obs-estimate)))
                folds.append({"target": target, "model": m, "heldout_island_id": island,
                              "train_pairs": 1225, "test_pairs": 50, "MSE": mse, "MAE": mae})
                for pidx, value, actual in zip(np.flatnonzero(test), estimate, obs):
                    predictions.append({"target": target, "model": m, "heldout_island_id": island,
                                        "source_pair_index": int(pidx), "observed": float(actual),
                                        "predicted": float(value),
                                        "squared_error": float((actual-value)**2),
                                        "absolute_error": float(abs(actual-value))})
    if len(predictions) != 30600 or len(folds) != 612:
        raise ValueError("Missing per-endpoint predictions or model folds")
    fold_frame = pd.DataFrame(folds)
    for target in TARGETS:
        for model,expected in PARENT_RICHNESS_MSE[target].items():
            original=float(fold_frame[(fold_frame.target==target)&(fold_frame.model==model)].MSE.mean())
            if abs(original-expected)>1e-9:
                raise ValueError("PR130 original richness-control MSE not reproduced: "+target+"/"+model)
            mae_expected=PARENT_RICHNESS_MAE[target][model]
            mae_observed=float(fold_frame[(fold_frame.target==target)&(fold_frame.model==model)].MAE.mean())
            if abs(mae_observed-mae_expected)>1e-9:
                raise ValueError("PR130 original richness-control MAE not reproduced: "+target+"/"+model)
    summary = compare_folds(fold_frame)
    if len(summary) != 14:
        raise ValueError("Missing fixed model comparisons")
    out.mkdir(parents=True, exist_ok=True)
    fold_frame.to_csv(out/"model_fold_losses.csv", index=False, float_format="%.17g")
    summary.to_csv(out/"comparison_summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(predictions).to_csv(out/"heldout_predictions.csv", index=False, float_format="%.17g")
    pd.DataFrame(early_coverage).to_csv(out/"historical_botanical_coverage_by_island.csv",index=False)
    receipt = {
        "schema": "chocho_island51_pre2011_allplant_effort_sensitivity_v1",
        "status": "ALL_51_ISLANDS_HISTORICAL_BOTANICAL_EFFORT_CONTROL_COMPLETE",
        "status_role": "exploratory_posthoc_sensitivity_not_independent_validation",
        "source_lock_status": validated["status"], "parent_folds_sha256": sha256(PARENT_FOLDS),
        "max_parent_mse_difference": max_mse, "max_parent_mae_difference": max_mae,
        "island_count": 51, "unique_island_pairs": 1275, "original_port_count": sum(ports.values()),
        "candidate_host_counts": {k: counts[k] for k in ids},
        "new_marginal_terms": list(richness), "early_botanical_effort_terms":list(early), "early_plant_record_scope":early_receipt, "feature_names": features_by_target,
        "models": list(MODELS), "fixed_comparisons": [list(x) for x in COMPARISONS],
        "fold_rows": 612, "endpoint_prediction_rows": 30600,
        "estimator": "train-only standardized OLS with intercept, predictions clipped to [0,1]",
        "folds_independent": False,
        "plant_inventory_pooled_2001_2024": True,
        "historical_plant_observation_proxy_not_historical_host_jaccard":True,
        "causal_colonization_or_transport_identified": False,
        "output_sha256": {p.name: sha256(p) for p in sorted(out.glob("*.csv"))},
    }
    (out/"result_receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2)+"\n")
    print("PRE2011_EFFORT_RESULT_JSON="+summary.to_json(orient="records", double_precision=12), flush=True)
    print("PRE2011_EFFORT_CHECK_JSON="+json.dumps({
        "parent_max_mse_error": max_mse, "parent_max_mae_error": max_mae,
        "model_fold_rows": len(folds), "predictions": len(predictions),
    }, sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, required=True)
    run(ap.parse_args().out_dir)
