#!/usr/bin/env python3
"""Source-ranked port island training-influence audit on fixed double-unseen folds."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import math
import numpy as np
import pandas as pd

from exploratory_island51_circuit_network_holdout import source, train_predict
from run_island51_host_circuit_joint import (
    INPUT_PATHS, TARGETS, load_original_port_counts, load_resistance_matrix,
    validate_input_lock,
)
from run_island51_latitude_controls import (
    ROOT, ROSTER, PAIRS, GEOGRAPHY, BASE_MODELS, ORIGINAL_COMPARISONS,
    make_features, latitude_pair_controls, partition_geographic_blocks,
    holdout_plan, CATEGORIES,
)
from run_island51_pre2011_effort_control import git_blob_sha
from run_island51_observed_ocean_circuits import sha256

MODEL_NAMES=tuple(model+"_LAT" for model in BASE_MODELS)
COMPARE=tuple((name+"_LAT",ref+"_LAT",aug+"_LAT",role)
              for name,ref,aug,role in ORIGINAL_COMPARISONS)
SCENARIO_NAMES=("FULL","DROP_15","DROP_9","DROP_15_9")
SCENARIO_DROPS={"FULL":(),"DROP_15":("15",),"DROP_9":("9",),
                "DROP_15_9":("15","9")}
PARENT=ROOT/"results/island51_latitude_controls_20261107/contrast_scope_summary.csv"
PARENT_BLOB="0d7189c9a7f538116c9ba4eaabe52644126f77c1"

def verify_port_roster(ids,port_counts):
    if len(ids)!=51 or set(port_counts)!=set(ids) or sum(port_counts.values())!=194:
        raise ValueError("The exact 51-island, 194-port roster is required")
    ranked=sorted(ids,key=lambda key:(-port_counts[key],int(key)))
    top=[(key,int(port_counts[key])) for key in ranked[:2]]
    if top!=[("15",39),("9",21)]:
        raise ValueError("Largest two original port groups do not match source-only selection")
    return pd.DataFrame([{"original_port_rank":i+1,"island_id":key,
       "original_port_count":int(port_counts[key]),
       "excluded_by_fixed_port_training_scenarios":key in {"15","9"}}
       for i,key in enumerate(ranked)])

def scenario_plan(q,ids,folds,source_plan):
    i,j=q.i.to_numpy(dtype=int),q.j.to_numpy(dtype=int)
    rows=[];masks=[]
    if len(folds)!=15 or len(source_plan)!=15:
        raise ValueError("Must preserve 15 original double-unseen block folds")
    spec=source_plan.set_index("fold_id",verify_integrity=True)
    coverage={scenario:np.zeros(len(q),dtype=np.uint8) for scenario in SCENARIO_NAMES}
    for scenario in SCENARIO_NAMES:
        drop_ids=SCENARIO_DROPS[scenario]
        drop=set(drop_ids)
        for fid,base_train,test in folds:
            original=spec.loc[fid]
            held=set(json.loads(original.heldout_island_ids))
            if not held or not held.issubset(set(ids)):
                raise ValueError("Original geographic block heldout identities invalid")
            excluded=held|drop
            excluded_indices=np.asarray([k for k,island in enumerate(ids)
                                         if island in excluded],dtype=int)
            expected_train=(~np.isin(i,excluded_indices))&(~np.isin(j,excluded_indices))
            train=base_train&(~np.isin(i,excluded_indices))&(~np.isin(j,excluded_indices))
            if (not np.array_equal(train,expected_train) or np.any(test&train) or
                int(test.sum())!=int(original.test_pairs) or
                int(train.sum())!=math.comb(51-len(excluded),2) or
                np.any(np.isin(i[train],excluded_indices)) or
                np.any(np.isin(j[train],excluded_indices))):
                raise ValueError("Training island exclusion, double-unseen holdout or support invalid")
            coverage[scenario]+=test.astype(np.uint8)
            masks.append((scenario,fid,train,test))
            rows.append({
                "scenario":scenario, "fold_id":fid,
                "block_a":int(original.block_a),"block_b":int(original.block_b),
                "pair_class":str(original.pair_class),
                "source_heldout_island_ids":json.dumps(sorted(held,key=int)),
                "additional_training_excluded_ids":json.dumps(list(drop_ids)),
                "total_training_excluded_island_ids":json.dumps(sorted(excluded,key=int)),
                "original_fold_train_pairs":int(original.train_pairs),
                "train_pairs":int(train.sum()),"test_pairs":int(test.sum()),
                "training_islands":51-len(excluded),
            })
    if len(masks)!=60 or len(rows)!=60 or any(
        not np.array_equal(cover,np.ones(len(q),dtype=np.uint8))
        for cover in coverage.values()):
        raise ValueError("Every original pair must be tested once in every scenario")
    return masks,pd.DataFrame(rows)

def compare_tables(pred,fold):
    scope=[]
    contrast_folds=[]
    for scenario in SCENARIO_NAMES:
        for target in TARGETS:
            p=pred[pred.scenario.eq(scenario)&pred.target.eq(target)]
            f=fold[fold.scenario.eq(scenario)&fold.target.eq(target)]
            for name,reference,augmented,role in COMPARE:
                x=p[p.model.eq(reference)].set_index("source_pair_index",verify_integrity=True)
                y=p[p.model.eq(augmented)].set_index("source_pair_index",verify_integrity=True).reindex(x.index)
                if len(x)!=1275 or len(y)!=1275 or y.predicted.isna().any():
                    raise ValueError("Source pair prediction support invalid for "+scenario)
                joined=pd.DataFrame({"pair_class":x.pair_class.to_numpy(),
                    "r_mse":x.squared_error.to_numpy(dtype=float),
                    "a_mse":y.squared_error.to_numpy(dtype=float),
                    "r_mae":x.absolute_error.to_numpy(dtype=float),
                    "a_mae":y.absolute_error.to_numpy(dtype=float)})
                if not np.array_equal(x.pair_class.to_numpy(),y.pair_class.to_numpy()):
                    raise ValueError("Contrasted models have different geographic test pairs")
                for category,expected in [("all_pairs",1275),
                                           ("within_longitude_block",235),
                                           ("between_longitude_blocks",1040)]:
                    k=joined if category=="all_pairs" else joined[joined.pair_class.eq(category)]
                    if len(k)!=expected:
                        raise ValueError("Unexpected fixed test pair scope count")
                    r_mse,a_mse=float(k.r_mse.mean()),float(k.a_mse.mean())
                    r_mae,a_mae=float(k.r_mae.mean()),float(k.a_mae.mean())
                    if r_mse<=0 or r_mae<=0:
                        raise ValueError("Reference loss not strictly positive")
                    scope.append({"scenario":scenario,"target":target,"contrast":name,
                        "role":role,"pair_class":category,"unique_evaluated_pairs":len(k),
                        "MSE_reference":r_mse,"MSE_augmented":a_mse,
                        "MSE_reduction_percent":100*(r_mse-a_mse)/r_mse,
                        "MAE_reference":r_mae,"MAE_augmented":a_mae,
                        "MAE_reduction_percent":100*(r_mae-a_mae)/r_mae,
                        "pairs_with_lower_squared_error":int((k.r_mse>k.a_mse).sum()),
                        "pairs_with_higher_squared_error":int((k.r_mse<k.a_mse).sum()),
                        "pairs_with_lower_absolute_error":int((k.r_mae>k.a_mae).sum())})
                a=f[f.model.eq(reference)].set_index("fold_id",verify_integrity=True)
                b=f[f.model.eq(augmented)].set_index("fold_id",verify_integrity=True).reindex(a.index)
                if len(a)!=15 or len(b)!=15 or b.MSE.isna().any():
                    raise ValueError("Missing source block-pair folds")
                for fid,row in a.iterrows():
                    other=b.loc[fid]
                    if (row.test_pairs!=other.test_pairs or row.train_pairs!=other.train_pairs or
                        row.pair_class!=other.pair_class):
                        raise ValueError("Compared models use different original fold support")
                    contrast_folds.append({
                        "scenario":scenario,"target":target,"comparison":name,"fold_id":fid,
                        "pair_class":row.pair_class,"train_pairs":int(row.train_pairs),
                        "test_pairs":int(row.test_pairs),
                        "MSE_reference":float(row.MSE),"MSE_augmented":float(other.MSE),
                        "MSE_gain":float(row.MSE-other.MSE),
                        "MAE_reference":float(row.MAE),"MAE_augmented":float(other.MAE),
                        "MAE_gain":float(row.MAE-other.MAE),
                        "MSE_improved":bool(row.MSE>other.MSE)})
    if len(scope)!=120 or len(contrast_folds)!=600:
        raise ValueError("Require exactly 120 scope summaries and 600 comparison folds")
    return pd.DataFrame(scope),pd.DataFrame(contrast_folds)

def reproduce_parent(scope):
    if git_blob_sha(PARENT.read_bytes())!=PARENT_BLOB:
        raise ValueError("PR134 parent comparison scores changed unexpectedly")
    old=pd.read_csv(PARENT)
    expected=old[old.contrast.isin([c[0] for c in COMPARE])]
    current=scope[scope.scenario.eq("FULL")]
    if len(old)!=72 or len(expected)!=30 or len(current)!=30:
        raise ValueError("Missing original latitude-adjusted comparisons from PR134")
    max_error=0.0
    for row in current.itertuples(index=False):
        p=expected[(expected.target==row.target)&
                   (expected.contrast==row.contrast)&
                   (expected.pair_class==row.pair_class)]
        if len(p)!=1 or int(p.iloc[0].unique_evaluated_pairs)!=row.unique_evaluated_pairs:
            raise ValueError("Changed original target/contrast/geographic scope")
        for field in ("MSE_reference","MSE_augmented","MAE_reference","MAE_augmented"):
            d=abs(float(p.iloc[0][field])-float(getattr(row,field)))
            max_error=max(max_error,d)
            if d>1e-9:
                raise ValueError("Original PR134 latitude model no longer reproduces: "+field)
    return max_error

def run(out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/p).resolve()) for p in ("data","provenance","results")):
        raise ValueError("New sensitivity output must not overwrite sources/results")
    lock=validate_input_lock(ROOT/"results/island51_host_circuit_joint_20261106/input_lock.json")
    q,ids,_=source(PAIRS,ROSTER)
    if len(ids)!=51 or len(q)!=1275:
        raise ValueError("Original 51-island/1275-pair source not retained")
    preparation=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    original_ports=load_original_port_counts(ROSTER,ids,preparation)
    ranked=verify_port_roster(ids,original_ports)
    matrices={name:load_resistance_matrix(ROOT/INPUT_PATHS[name],ids)
              for name in ("uniform","ais","wind")}
    blocks,locations=partition_geographic_blocks(GEOGRAPHY,ids)
    original_folds,original_spec=holdout_plan(q,ids,blocks)
    masks,manifest=scenario_plan(q,ids,original_folds,original_spec)
    features,names,historical=make_features(q,ids,original_ports,matrices)
    latitude=latitude_pair_controls(q,ids)
    designs={}
    for target in TARGETS:
        designs[target]={}
        for model in BASE_MODELS:
            design=np.column_stack([features[target][model]]+
                                   [latitude[col] for col in latitude])
            if design.shape[0]!=1275 or not np.isfinite(design).all():
                raise ValueError("Latitude-adjusted model design is not finite")
            designs[target][model+"_LAT"]=design
    predictions=[];losses=[]
    ii=q.i.to_numpy(dtype=int); jj=q.j.to_numpy(dtype=int)
    source_a=q.island_a.to_numpy(dtype=str);source_b=q.island_b.to_numpy(dtype=str)
    for scenario,fid,train,test in masks:
        meta=manifest[manifest.scenario.eq(scenario)&manifest.fold_id.eq(fid)]
        if len(meta)!=1:
            raise ValueError("Scenario/fold identity missing")
        kind=str(meta.iloc[0].pair_class)
        source_idx=np.flatnonzero(test)
        for target in TARGETS:
            response=q[target].to_numpy(dtype=float)
            observed=response[test]
            for model in MODEL_NAMES:
                predicted=train_predict(designs[target][model],response,train,test)
                if predicted.shape!=(len(source_idx),) or not np.isfinite(predicted).all():
                    raise ValueError("Model prediction shape or finiteness invalid")
                sq=(observed-predicted)**2
                ae=np.abs(observed-predicted)
                losses.append({
                    "scenario":scenario,"target":target,"model":model,"fold_id":fid,
                    "pair_class":kind,"train_pairs":int(train.sum()),
                    "test_pairs":int(test.sum()),
                    "MSE":float(np.mean(sq)),"MAE":float(np.mean(ae))})
                for idx,y,p,sqerr,abserr in zip(source_idx,observed,predicted,sq,ae):
                    predictions.append({
                        "scenario":scenario,"target":target,"model":model,
                        "fold_id":fid,"pair_class":kind,"source_pair_index":int(idx),
                        "island_a":source_a[idx],"island_b":source_b[idx],
                        "observed":float(y),"predicted":float(p),
                        "squared_error":float(sqerr),"absolute_error":float(abserr)})
    p=pd.DataFrame(predictions);fold=pd.DataFrame(losses)
    if (len(p)!=40800 or len(fold)!=480 or
        p.duplicated(["scenario","target","model","source_pair_index"]).any() or
        fold.duplicated(["scenario","target","model","fold_id"]).any()):
        raise ValueError("Original full source pair/model/scenario support incomplete")
    scope,contrast_folds=compare_tables(p,fold)
    max_parent_diff=reproduce_parent(scope)
    out.mkdir(parents=True,exist_ok=True)
    frames={
       "original_port_roster_rank.csv":ranked,
       "source_longitude_island_blocks.csv":locations,
       "scenario_fold_specification.csv":manifest,
       "heldout_pair_predictions.csv":p,
       "model_fold_losses.csv":fold,
       "contrast_fold_losses.csv":contrast_folds,
       "contrast_scope_summary.csv":scope,
    }
    for name,table in frames.items():
        table.to_csv(out/name,index=False,float_format="%.17g")
    receipt={
       "schema":"chocho_original51_top_port_training_influence_v1",
       "status":"ALL_1275_DOUBLE_UNSEEN_PAIRS_EVALUATED_IN_FOUR_TRAINING_POPULATIONS",
       "parent_results_known_before_sensitivity":True,
       "exact_original_source_islands":51,"original_biological_ports":194,
       "original_pair_count":1275,"original_source_holdout_folds":15,
       "scenarios":SCENARIO_DROPS,"top_two_port_groups":[["15",39],["9",21]],
       "models":MODEL_NAMES,"comparisons":[x[0] for x in COMPARE],
       "saved_predictions":len(p),"model_fold_rows":len(fold),
       "contrast_fold_rows":len(contrast_folds),"contrast_scope_rows":len(scope),
       "training_population_changed_but_all_test_pairs_unchanged":True,
       "parent_30_latitude_comparisons_reproduced_max_difference":max_parent_diff,
       "source_input_lock_status":lock["status"],
       "archived_historical_plant_effort":historical,
       "original_source_matrices_unchanged":True,
       "host_jaccard_pooled_2001_2024_may_postdate_butterfly_target":True,
       "source_port_groups_not_ship_routes_or_butterfly_colonization":True,
       "not_an_independent_validation_population":True,
       "no_naive_independent_fold_p_values":True,
       "output_sha256":{key:sha256(out/key) for key in frames},
    }
    (out/"execution_receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("TOP_PORT_INFLUENCE_RESULT_JSON="+scope.to_json(orient="records",double_precision=12),flush=True)
    print("TOP_PORT_INFLUENCE_RECEIPT_JSON="+json.dumps({
        "islands":51,"pairs_per_scenario":1275,"scenarios":len(SCENARIO_NAMES),
        "predictions":len(p),"model_fold_rows":len(fold),
        "contrast_fold_rows":len(contrast_folds),"scope_rows":len(scope),
        "parent_max_absolute_score_difference":max_parent_diff},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)
