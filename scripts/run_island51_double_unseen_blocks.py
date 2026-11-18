#!/usr/bin/env python3
"""Fixed 5-longitude-block, both-endpoints-unseen butterfly similarity prediction."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from exploratory_island51_circuit_network_holdout import source,train_predict
from run_island51_host_circuit_joint import (
    INPUT_PATHS,TARGETS,load_original_port_counts,load_resistance_matrix,
    make_joint_features,validate_input_lock,
)
from run_island51_richness_control import host_counts,richness_features
from run_island51_pre2011_effort_control import early_botanical_features,git_blob_sha
from run_island51_observed_ocean_circuits import sha256

GEOGRAPHY=ROOT/"results/island51_host_circuit_atlas_20261106/island_geography.csv"
GEOGRAPHY_GIT_BLOB="f0c5b1139e54b4872d28682f0991f89f2368ab43"
ROSTER=ROOT/INPUT_PATHS["islands"]
PAIRS=ROOT/INPUT_PATHS["pairs"]
PREVIOUS=ROOT/"results/island51_pre2011_effort_control_20261107/comparison_summary.csv"
MODELS=("CRE","HRE","HCRE")
COMPARISONS=(("CRE_to_HCRE","CRE","HCRE","host_composition_beyond_circuits"),
             ("HRE_to_HCRE","HRE","HCRE","circuits_beyond_host_composition"))
BLOCK_SIZES=(11,10,10,10,10)
CATEGORIES=("all_pairs","within_longitude_block","between_longitude_blocks")
TARGET_COUNT=1275

def partition_geographic_blocks(path,ids):
    if git_blob_sha(Path(path).read_bytes())!=GEOGRAPHY_GIT_BLOB:
        raise ValueError("Audited representative island geography is not the pinned file")
    frame=pd.read_csv(path,dtype={"island_id":str,"landmass_id":str})
    if len(frame)!=51 or frame.island_id.nunique()!=51 or set(frame.island_id)!=set(ids):
        raise ValueError("Audited island geography must cover all original51 IDs once")
    if (frame.island_id.astype(str)!=frame.landmass_id.astype(str)).any():
        raise ValueError("Geographic island_id and landmass_id differ")
    x=frame.representative_longitude.to_numpy(dtype=float)
    y=frame.representative_latitude.to_numpy(dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(x < -180) or np.any(x > 180) or np.any(y < -90) or np.any(y > 90):
        raise ValueError("Invalid original geographic representative longitude/latitude")
    ordered=frame.assign(numeric_id=frame.island_id.astype(int)).sort_values(
        ["representative_longitude","numeric_id"],kind="mergesort").reset_index(drop=True)
    labels={}
    rows=[]
    start=0
    for block,n in enumerate(BLOCK_SIZES):
        for item in ordered.iloc[start:start+n].itertuples(index=False):
            sid=str(item.island_id)
            if sid in labels:
                raise ValueError("Duplicate island in longitude partition")
            labels[sid]=block
            rows.append({"island_id":sid,"longitude_block":block,
                         "representative_longitude":float(item.representative_longitude),
                         "representative_latitude":float(item.representative_latitude),
                         "group_example_port":str(item.example_port_name)})
        start+=n
    if start!=51 or set(labels)!=set(ids) or sorted(list(labels.values())).count(0)!=11:
        raise ValueError("Longitude partition counts or identity wrong")
    return labels,pd.DataFrame(rows)

def holdout_plan(q,ids,blocks):
    a=q.i.to_numpy(dtype=int);b=q.j.to_numpy(dtype=int)
    gid=np.asarray([blocks[id] for id in ids],dtype=int)
    ba=gid[a];bb=gid[b]
    masks=[]
    rows=[]
    used=np.zeros(len(q),dtype=int)
    for k in range(5):
        for h in range(k,5):
            held=np.where((gid==k)|(gid==h))[0]
            train=(~np.isin(a,held)) & (~np.isin(b,held))
            test= ((ba==k)&(bb==h))|((ba==h)&(bb==k))
            n=int(test.sum());ntrain=int(train.sum())
            expected= (BLOCK_SIZES[k]*(BLOCK_SIZES[k]-1)//2 if k==h else BLOCK_SIZES[k]*BLOCK_SIZES[h])
            remain=51-len(held)
            if (n!=expected or ntrain!=remain*(remain-1)//2 or
                np.any(train&test) or np.any(np.isin(a[train],held)) or np.any(np.isin(b[train],held))):
                raise ValueError("A longitude fold has invalid heldout endpoint support")
            used+=test.astype(int)
            fid=f"B{k}_B{h}"
            rows.append({"fold_id":fid,"block_a":k,"block_b":h,
                         "pair_class":"within_longitude_block" if k==h else "between_longitude_blocks",
                         "heldout_islands":len(held),"training_islands":remain,
                         "train_pairs":ntrain,"test_pairs":n,
                         "heldout_island_ids":json.dumps([ids[z] for z in held.tolist()])})
            masks.append((fid,train,test))
    if len(masks)!=15 or not np.array_equal(used,np.ones(len(q),dtype=int)) or len(q)!=TARGET_COUNT:
        raise ValueError("All original island pairs must appear in exactly one double-unseen fold")
    return masks,pd.DataFrame(rows)

def make_features(q,ids,port_counts,matrices):
    history,coverage,receipt=early_botanical_features(q,ids,port_counts)
    richness=richness_features(q,ids,host_counts(ROSTER,ids))
    designs={}
    names={}
    for target in TARGETS:
        prior=make_joint_features(q,ids,port_counts,matrices,target)
        cre=dict(prior["C"],**richness,**history)
        hre=dict(prior["H"],**richness,**history)
        hcre=dict(prior["HC"],**richness,**history)
        m={"CRE":cre,"HRE":hre,"HCRE":hcre}
        designs[target]={model:np.column_stack(list(v.values())) for model,v in m.items()}
        names[target]={model:list(v.keys()) for model,v in m.items()}
    return designs,names,receipt

def summarize(predictions,folds):
    result=[]
    comparison_folds=[]
    for target in TARGETS:
        p=predictions[predictions.target.eq(target)]
        f=folds[folds.target.eq(target)]
        for label,ref,aug,role in COMPARISONS:
            x=p[p.model.eq(ref)].set_index("source_pair_index")
            y=p[p.model.eq(aug)].set_index("source_pair_index").reindex(x.index)
            if len(x)!=TARGET_COUNT or len(y)!=TARGET_COUNT or y.predicted.isna().any():
                raise ValueError("Incomplete pair predictions for comparison")
            joined=pd.DataFrame({"source_pair_index":x.index.to_numpy(dtype=int),
                "pair_class":x.pair_class.to_numpy(),
                "reference_squared_error":x.squared_error.to_numpy(dtype=float),
                "augmented_squared_error":y.squared_error.to_numpy(dtype=float),
                "reference_absolute_error":x.absolute_error.to_numpy(dtype=float),
                "augmented_absolute_error":y.absolute_error.to_numpy(dtype=float)})
            for cat in CATEGORIES:
                z=joined if cat=="all_pairs" else joined[joined.pair_class.eq(cat)]
                n=len(z)
                expected={"all_pairs":1275,"within_longitude_block":235,"between_longitude_blocks":1040}[cat]
                if n!=expected:
                    raise ValueError("Wrong geographic comparison sample")
                r_mse=float(z.reference_squared_error.mean())
                a_mse=float(z.augmented_squared_error.mean())
                r_mae=float(z.reference_absolute_error.mean())
                a_mae=float(z.augmented_absolute_error.mean())
                if r_mse<=0 or r_mae<=0:
                    raise ValueError("Reference error must be positive for ratio")
                d=z.reference_squared_error-z.augmented_squared_error
                dam=z.reference_absolute_error-z.augmented_absolute_error
                result.append({"target":target,"contrast":label,"role":role,"pair_class":cat,
                    "unique_evaluated_pairs":n,
                    "MSE_reference":r_mse,"MSE_augmented":a_mse,
                    "MSE_reduction_percent":100*(r_mse-a_mse)/r_mse,
                    "MAE_reference":r_mae,"MAE_augmented":a_mae,
                    "MAE_reduction_percent":100*(r_mae-a_mae)/r_mae,
                    "pairs_with_lower_squared_error":int((d>0).sum()),
                    "pairs_with_higher_squared_error":int((d<0).sum()),
                    "pairs_with_lower_absolute_error":int((dam>0).sum())})
            left=f[f.model.eq(ref)].set_index("fold_id")
            right=f[f.model.eq(aug)].set_index("fold_id").reindex(left.index)
            if len(left)!=15 or len(right)!=15 or right.MSE.isna().any():
                raise ValueError("Missing comparison fold")
            for fid,row in left.iterrows():
                other=right.loc[fid]
                if (row.test_pairs!=other.test_pairs or row.train_pairs!=other.train_pairs or
                    row.pair_class!=other.pair_class):
                    raise ValueError("Mismatched geographic fold")
                comparison_folds.append({"target":target,"comparison":label,"fold_id":fid,
                    "pair_class":row.pair_class,"test_pairs":int(row.test_pairs),
                    "MSE_reference":float(row.MSE),"MSE_augmented":float(other.MSE),
                    "MSE_gain":float(row.MSE-other.MSE),
                    "MAE_reference":float(row.MAE),"MAE_augmented":float(other.MAE),
                    "MAE_gain":float(row.MAE-other.MAE),
                    "MSE_improved":bool(row.MSE>other.MSE)})
    if len(result)!=12 or len(comparison_folds)!=60:
        raise ValueError("Expected all 12 scope summaries and 60 contrast folds")
    return pd.DataFrame(result),pd.DataFrame(comparison_folds)

def run(out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/p).resolve()) for p in ("data","provenance","results")):
        raise ValueError("Do not overwrite sources or previous results")
    lock=validate_input_lock(ROOT/"results/island51_host_circuit_joint_20261106/input_lock.json")
    q,ids,_=source(ROOT/INPUT_PATHS["pairs"],ROOT/INPUT_PATHS["islands"])
    if len(q)!=1275 or len(ids)!=51:
        raise ValueError("Original 51-island cohort changed")
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    ports=load_original_port_counts(ROOT/INPUT_PATHS["islands"],ids,prep)
    matrices={key:load_resistance_matrix(ROOT/INPUT_PATHS[key],ids) for key in ("uniform","ais","wind")}
    blocks,locations=partition_geographic_blocks(GEOGRAPHY,ids)
    folds,fold_spec=holdout_plan(q,ids,blocks)
    designs,feature_names,historical=make_features(q,ids,ports,matrices)
    records=[];losses=[]
    for target in TARGETS:
        response=q[target].to_numpy(dtype=float)
        for fid,train,test in folds:
            item=fold_spec.set_index("fold_id").loc[fid]
            cls=item.pair_class
            for model in MODELS:
                predicted=train_predict(designs[target][model],response,train,test)
                observed=response[test]
                if predicted.shape!=(int(test.sum()),) or not np.isfinite(predicted).all():
                    raise ValueError("Prediction shape/values invalid")
                errs=(observed-predicted)**2
                abses=np.abs(observed-predicted)
                losses.append({"target":target,"model":model,"fold_id":fid,
                  "pair_class":cls,"train_pairs":int(train.sum()),"test_pairs":int(test.sum()),
                  "MSE":float(np.mean(errs)),"MAE":float(np.mean(abses))})
                for idx,pred,obs,square,ab in zip(np.flatnonzero(test),predicted,observed,errs,abses):
                    records.append({"target":target,"model":model,"fold_id":fid,
                        "pair_class":cls,"source_pair_index":int(idx),
                        "island_a":str(q.iloc[idx].island_a),"island_b":str(q.iloc[idx].island_b),
                        "observed":float(obs),"predicted":float(pred),
                        "squared_error":float(square),"absolute_error":float(ab)})
    pred=pd.DataFrame(records);model_folds=pd.DataFrame(losses)
    if len(pred)!=7650 or len(model_folds)!=90 or pred.duplicated(["target","model","source_pair_index"]).any():
        raise ValueError("Incomplete both-unseen prediction count")
    summary,contrast_folds=summarize(pred,model_folds)
    out.mkdir(parents=True,exist_ok=True)
    tables={"longitude_island_blocks.csv":locations,
            "fold_specification.csv":fold_spec,
            "heldout_pair_predictions.csv":pred,
            "model_fold_losses.csv":model_folds,
            "contrast_fold_losses.csv":contrast_folds,
            "contrast_scope_summary.csv":summary}
    for path,frame in tables.items():
        frame.to_csv(out/path,index=False,float_format="%.17g")
    receipt={"schema":"chocho_original51_double_unseen_longitude_blocks_v1",
       "status":"ALL_ORIGINAL_1275_PAIRS_EVALUATED_WITH_BOTH_ENDPOINTS_ABSENT",
       "prior_LOIO_results_already_seen":True,"new_double_unseen_results_known_when_frozen":False,
       "target_independent_longitude_blocks":list(BLOCK_SIZES),
       "source_graph_and_host_jaccard_unmodified":True,
       "source_lock_status":lock["status"],"geography_sha":sha256(GEOGRAPHY),
       "islands":51,"unique_source_pairs":1275,"folds":15,
       "within_block_pairs":235,"between_block_pairs":1040,
       "models":list(MODELS),"comparisons":[c[0] for c in COMPARISONS],
       "saved_predictions":len(pred),"model_fold_rows":len(model_folds),
       "contrast_fold_rows":len(contrast_folds),"contrast_scope_rows":len(summary),
       "historical_botanical_recording_proxy":historical,
       "feature_names":feature_names,
       "pooled_host_jaccard_2001_2024_may_postdate_target":True,
       "historical_host_identity_not_reconstructed":True,
       "longitude_blocks_not_bioregions_or_independent_validation_population":True,
       "no_independent_fold_p_values_or_causal_claims":True,
       "outputs_sha256":{key:sha256(out/key) for key in tables}}
    (out/"execution_receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("DOUBLE_UNSEEN_RESULT_JSON="+summary.to_json(orient="records",double_precision=12),flush=True)
    print("DOUBLE_UNSEEN_RECEIPT_JSON="+json.dumps({
       "islands":51,"pairs":1275,"folds":15,"within":235,"between":1040,
       "predictions":len(pred),"fold_rows":len(model_folds),
       "contrast_fold_rows":len(contrast_folds)},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)
