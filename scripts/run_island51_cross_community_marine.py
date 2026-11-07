#!/usr/bin/env python3
"""Cross-community marine prediction on the exact 990 historical-plant-indexed island pairs.

Uses a single fully archived 194-port source, the original 51-island butterfly pairs,
and existing 15 geography-only both-endpoints-unseen splits. No GBIF network requests.
"""
from __future__ import annotations
import argparse
import json
import gzip
from pathlib import Path
import numpy as np
import pandas as pd

from exploratory_island51_circuit_network_holdout import source,train_predict
from run_island51_host_circuit_joint import (
    INPUT_PATHS,load_original_port_counts,load_resistance_matrix,
    make_joint_features,validate_input_lock,
)
from run_island51_richness_control import host_counts,richness_features
from run_island51_pre2011_effort_control import early_botanical_features
from run_island51_latitude_controls import (
    ROOT,ROSTER,PAIRS,GEOGRAPHY,latitude_pair_controls,
    partition_geographic_blocks,holdout_plan
)
from run_island51_pre2011_general_flora_prediction import (
    early_general_flora_source,general_flora_features
)
from run_island51_observed_ocean_circuits import sha256

SOURCE=ROOT/"results/island51_pre2011_public_general_flora_20261107/source"
MODEL_NAMES=("BASE","AIS","WIND","BOTH")
PANELS=(
    ("early_general_flora","early_general_flora_jaccard"),
    ("earlier_butterflies","earlier_butterfly_jaccard"),
    ("later_butterflies","new_butterfly_jaccard"),
    ("later_butterflies_retro_hosts","new_butterfly_jaccard"),
)
COMPARISONS=(
    ("BASE_to_AIS","BASE","AIS","marine_AIS_increment"),
    ("BASE_to_WIND","BASE","WIND","mean_windspeed_increment"),
    ("AIS_to_BOTH","AIS","BOTH","wind_beyond_AIS"),
    ("WIND_to_BOTH","WIND","BOTH","AIS_beyond_wind"),
    ("BASE_to_BOTH","BASE","BOTH","joint_marine_increment"),
)
SCOPES=("all_both_plants_indexed","within_longitude_block","between_longitude_blocks","both_plants_at_least_10_specieskeys")

def validate_source(q,ids,root_source=SOURCE):
    plant_sets,source_port_coverage,source_island_coverage,source_receipt=early_general_flora_source(root_source,ids)
    controls,plant_jaccard,masks=general_flora_features(q,ids,plant_sets)
    counts=np.asarray([len(plant_sets[k]) for k in ids],dtype=int)
    good=(counts[q.i.to_numpy(dtype=int)]>0)&(counts[q.j.to_numpy(dtype=int)]>0)
    high=(counts[q.i.to_numpy(dtype=int)]>=10)&(counts[q.j.to_numpy(dtype=int)]>=10)
    if (len(q)!=1275 or len(ids)!=51 or sum(counts>0)!=45 or
        int(np.sum(good))!=990 or int(np.sum(counts>=10))!=39 or
        int(np.sum(high))!=741 or
        not np.array_equal(good,masks["both_early_species_indexed"])):
        raise ValueError("The original 45 indexed-plant islands / 990 pairs / 741 high-index pairs differ")
    if not np.isfinite(plant_jaccard).all():
        raise ValueError("Non-finite source early general flora Jaccard")
    return plant_sets,source_port_coverage,source_island_coverage,source_receipt,controls,plant_jaccard,good,high

def build_common(q,ids,original_port_counts,matrices,source_island_coverage,plant_controls):
    # Only source geometry and source-plant observation opportunity; no butterfly outcome
    # and no historical general-flora Jaccard in the plant outcome baseline.
    g=make_joint_features(q,ids,original_port_counts,matrices,"earlier_butterfly_jaccard")["G"]
    if "earlier_butterfly_Jaccard" in g or "host_jaccard" in g:
        raise ValueError("Butterfly target or contemporary host Jaccard leaked into plant baseline")
    archive_effort,coverage,historical=early_botanical_features(q,ids,original_port_counts)
    latitude=latitude_pair_controls(q,ids)
    rows=source_island_coverage.set_index("island_id")
    record=np.array([float(rows.loc[k,"sum_live_original_port_5km_records_NOT_unique_occurrences"])
                     for k in ids],dtype=float)
    if np.any(record<0) or not np.isfinite(record).all():
        raise ValueError("Live source botanical record index counts unavailable")
    transformed=np.log1p(record)
    left,right=transformed[q.i.to_numpy(dtype=int)],transformed[q.j.to_numpy(dtype=int)]
    indexed_species_richness={
        key:value for key,value in plant_controls.items()
        if key in (
          "sum_log1p_2001_2010_general_plant_specieskey_richness",
          "abs_difference_log1p_2001_2010_general_plant_specieskey_richness"
        )
    }
    if len(indexed_species_richness)!=2:
        raise ValueError("Original indexed plant richness controls not present")
    live_effort={
      "sum_log1p_live_2001_2010_port_buffer_indexed_plant_records":left+right,
      "abs_difference_log1p_live_2001_2010_port_buffer_indexed_plant_records":np.abs(left-right),
    }
    common=dict(g,**latitude,**archive_effort,**indexed_species_richness,**live_effort)
    if len(common)!=15:
        raise ValueError("Fixed common geography + plant opportunity covariates differ: "+str(list(common)))
    if any(np.asarray(v).shape!=(1275,) or not np.isfinite(v).all() for v in common.values()):
        raise ValueError("Incomplete source common marine-control data")
    return common,{"old_record_effort_receipt":historical,
                  "new_botanical_effort_controls":list(live_effort)}

def model_designs(q,ids,port_counts,matrices,source_island_coverage,plant_controls,plant_jaccard):
    baseline,source_receipt=build_common(q,ids,port_counts,matrices,source_island_coverage,plant_controls)
    idx_a=q.i.to_numpy(dtype=int);idx_b=q.j.to_numpy(dtype=int)
    from island51_ocean_circuit import resistance_feature
    marine={name:resistance_feature(matrices[name],idx_a,idx_b) for name in ("ais","wind")}
    original_counts=host_counts(ROSTER,ids)
    original_richness=richness_features(q,ids,original_counts)
    designs={};names={}
    for panel,outcome in PANELS:
        features=dict(baseline)
        if panel.startswith("later"):
            features["previous_2001_2010_butterfly_jaccard"]=q.earlier_butterfly_jaccard.to_numpy(dtype=float)
        if panel=="later_butterflies_retro_hosts":
            features["retrospective_host_jaccard_pooled_2001_2024"]=q.host_jaccard.to_numpy(dtype=float)
            features.update(original_richness)
        if "early_general_flora_jaccard" in features:
            raise ValueError("Plant target composition used as its own predictor")
        families={
           "BASE":dict(features),
           "AIS":dict(features,native_AIS_circuit_resistance=marine["ais"]),
           "WIND":dict(features,mean_windspeed_circuit_resistance=marine["wind"]),
           "BOTH":dict(features,native_AIS_circuit_resistance=marine["ais"],
                       mean_windspeed_circuit_resistance=marine["wind"]),
        }
        designs[panel]={model:np.column_stack(list(ff.values())) for model,ff in families.items()}
        names[panel]={model:list(ff) for model,ff in families.items()}
        if any(not np.isfinite(mat).all() or mat.shape[0]!=1275 for mat in designs[panel].values()):
            raise ValueError("Invalid indexed source panel covariates")
    targets={
      "early_general_flora":plant_jaccard,
      "earlier_butterflies":q.earlier_butterfly_jaccard.to_numpy(dtype=float),
      "later_butterflies":q.new_butterfly_jaccard.to_numpy(dtype=float),
      "later_butterflies_retro_hosts":q.new_butterfly_jaccard.to_numpy(dtype=float),
    }
    return designs,names,targets,source_receipt

def summarize(predictions,folds,scopes):
    contrasts=[]
    per_fold=[]
    for panel,_ in PANELS:
        p=predictions[predictions.panel.eq(panel)]
        f=folds[folds.panel.eq(panel)]
        for label,ref,aug,role in COMPARISONS:
            refrows=p[p.model.eq(ref)].set_index("source_pair_index",verify_integrity=True).sort_index()
            augrows=p[p.model.eq(aug)].set_index("source_pair_index",verify_integrity=True).sort_index()
            if len(refrows)!=990 or len(augrows)!=990 or not np.array_equal(refrows.index,augrows.index):
                raise ValueError("Matched target source pairs or comparison populations differ")
            for scope in SCOPES:
                indices=np.flatnonzero(scopes[scope])
                keep=np.intersect1d(refrows.index.to_numpy(dtype=int),indices)
                if scope=="all_both_plants_indexed" and len(keep)!=990:
                    raise ValueError("Matched 990 original indexed-plant pairs must be preserved")
                if scope=="both_plants_at_least_10_specieskeys" and len(keep)!=741:
                    raise ValueError("Matched 741 original high-botanical-richness pairs must be retained")
                if len(keep)==0:
                    raise ValueError("Source-only comparison scope empty")
                x=refrows.loc[keep];y=augrows.loc[keep]
                r_mse=float(x.squared_error.mean());a_mse=float(y.squared_error.mean())
                r_mae=float(x.absolute_error.mean());a_mae=float(y.absolute_error.mean())
                if r_mse<=0 or r_mae<=0:
                    raise ValueError("Reference source panel loss is zero")
                contrasts.append({
                    "panel":panel,"contrast":label,"role":role,
                    "pair_class":scope,"unique_source_pairs":len(keep),
                    "MSE_reference":r_mse,"MSE_augmented":a_mse,
                    "MSE_gain_percent":100*(r_mse-a_mse)/r_mse,
                    "MAE_reference":r_mae,"MAE_augmented":a_mae,
                    "MAE_gain_percent":100*(r_mae-a_mae)/r_mae,
                    "pairs_with_lower_squared_error":int(np.sum(x.squared_error>y.squared_error)),
                    "pairs_with_higher_squared_error":int(np.sum(x.squared_error<y.squared_error)),
                    "pairs_with_lower_absolute_error":int(np.sum(x.absolute_error>y.absolute_error)),
                })
            x=f[f.model.eq(ref)].set_index("fold_id",verify_integrity=True)
            y=f[f.model.eq(aug)].set_index("fold_id",verify_integrity=True).reindex(x.index)
            if len(x)!=15 or len(y)!=15 or y.MSE.isna().any():
                raise ValueError("Matched source target model-fold comparisons missing")
            for fid,row in x.iterrows():
                other=y.loc[fid]
                if int(row.train_pairs)!=int(other.train_pairs) or int(row.test_pairs)!=int(other.test_pairs):
                    raise ValueError("Contrasted models use different source train/test folds")
                per_fold.append({
                    "panel":panel,"contrast":label,"fold_id":fid,
                    "pair_class":row.pair_class,"train_pairs":int(row.train_pairs),
                    "test_pairs":int(row.test_pairs),
                    "MSE_reference":float(row.MSE),"MSE_augmented":float(other.MSE),
                    "MSE_gain":float(row.MSE-other.MSE),
                    "MAE_reference":float(row.MAE),"MAE_augmented":float(other.MAE),
                    "MAE_gain":float(row.MAE-other.MAE),
                })
    if len(contrasts)!=80 or len(per_fold)!=300:
        raise ValueError("Fixed four-panel, five-contrast, four-scope source support differs")
    return pd.DataFrame(contrasts),pd.DataFrame(per_fold)

def run(out_dir,source_dir=SOURCE):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/k).resolve()) for k in ("data","provenance","results")):
        raise ValueError("New source-locked comparison must not overwrite original inputs or results")
    original_lock=validate_input_lock(ROOT/"results/island51_host_circuit_joint_20261106/input_lock.json")
    q,ids,_=source(PAIRS,ROSTER)
    if len(ids)!=51 or len(q)!=1275:
        raise ValueError("Original source butterfly 51-island panel altered")
    plants,port_coverage,island_cov,plant_receipt,controls,plant_jaccard,eligible,high=validate_source(q,ids,source_dir)
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    counts=load_original_port_counts(ROSTER,ids,prep)
    matrices={key:load_resistance_matrix(ROOT/INPUT_PATHS[key],ids) for key in ("uniform","ais","wind")}
    geo,locations=partition_geographic_blocks(GEOGRAPHY,ids)
    original_folds,original_spec=holdout_plan(q,ids,geo)
    source_pair_class=np.array(["within_longitude_block" if geo[str(a)]==geo[str(b)]
        else "between_longitude_blocks" for a,b in zip(q.island_a,q.island_b)],dtype=object)
    scopes={
       "all_both_plants_indexed":eligible,
       "within_longitude_block":eligible&(source_pair_class=="within_longitude_block"),
       "between_longitude_blocks":eligible&(source_pair_class=="between_longitude_blocks"),
       "both_plants_at_least_10_specieskeys":high,
    }
    if sum(int(scopes[x].sum()) for x in ("within_longitude_block","between_longitude_blocks"))!=990:
        raise ValueError("Botanical-indexed pair geography categories do not partition 990")
    designs,names,targets,receipt= model_designs(q,ids,counts,matrices,island_cov,controls,plant_jaccard)
    source_pair_ledger=pd.DataFrame({
       "source_pair_index":np.arange(len(q),dtype=int),
       "island_a":q.island_a.to_numpy(dtype=str),
       "island_b":q.island_b.to_numpy(dtype=str),
       "early_general_flora_jaccard_imputed":plant_jaccard,
       "early_general_flora_both_islands_indexed":eligible,
       "early_general_flora_both_islands_at_least_10_specieskeys":high,
       "earlier_butterfly_jaccard":q.earlier_butterfly_jaccard.to_numpy(dtype=float),
       "later_butterfly_jaccard":q.new_butterfly_jaccard.to_numpy(dtype=float),
       "longitude_pair_class":source_pair_class,
       "source_port_distance_km":q.port_distance_km.to_numpy(dtype=float),
    })
    fold_manifest=[]; masks=[]
    tested=np.zeros(len(q),dtype=np.uint8)
    for fid,original_train,original_test in original_folds:
        train=original_train&eligible
        test=original_test&eligible
        spec=original_spec.loc[original_spec.fold_id.eq(fid)]
        if len(spec)!=1 or int(test.sum())<1 or int(train.sum())<20:
            raise ValueError("Source-only eligible original block-pair holdout lacks support")
        held=set(json.loads(spec.iloc[0].heldout_island_ids))
        if any(a in held or b in held for a,b in q.loc[train,["island_a","island_b"]].itertuples(index=False,name=None)):
            raise ValueError("Training contains an island in a heldout geographic block")
        tested+=test.astype(np.uint8)
        masks.append((fid,train,test,str(spec.iloc[0].pair_class)))
        fold_manifest.append({
            "fold_id":fid,"pair_class":str(spec.iloc[0].pair_class),
            "source_heldout_blocks":str(spec.iloc[0].block_a)+":"+str(spec.iloc[0].block_b),
            "source_heldout_island_ids":json.dumps(sorted(held,key=int)),
            "train_source_pairs":int(train.sum()),"test_source_pairs":int(test.sum()),
            "original_unfiltered_train_pairs":int(original_train.sum()),
            "original_unfiltered_test_pairs":int(original_test.sum()),
        })
    if not np.array_equal(tested,eligible.astype(np.uint8)) or len(masks)!=15:
        raise ValueError("Missing/duplicated source botanical-eligible heldout test pair")
    pred=[];losses=[]
    for panel,_ in PANELS:
        observed=targets[panel]
        for fid,train,test,kind in masks:
            orig_indices=np.flatnonzero(test)
            for model in MODEL_NAMES:
                output=train_predict(designs[panel][model],observed,train,test)
                if not np.isfinite(output).all() or output.shape!=(len(orig_indices),):
                    raise ValueError("Incomplete model source pair predictions")
                source_values=observed[test]
                squared=(source_values-output)**2
                absolute=np.abs(source_values-output)
                losses.append({
                  "panel":panel,"model":model,"fold_id":fid,"pair_class":kind,
                  "train_pairs":int(train.sum()),"test_pairs":int(test.sum()),
                  "MSE":float(np.mean(squared)),"MAE":float(np.mean(absolute)),
                })
                for index,prediction,value,se,ae in zip(orig_indices,output,source_values,squared,absolute):
                    pred.append({
                       "panel":panel,"model":model,"fold_id":fid,
                       "pair_class":kind,"source_pair_index":int(index),
                       "island_a":str(q.iloc[index].island_a),
                       "island_b":str(q.iloc[index].island_b),
                       "observed":float(value),"predicted":float(prediction),
                       "squared_error":float(se),"absolute_error":float(ae),
                    })
    preds=pd.DataFrame(pred);fold=pd.DataFrame(losses)
    if (len(preds)!=15840 or len(fold)!=240 or
        preds.duplicated(["panel","model","source_pair_index"]).any() or
        fold.duplicated(["panel","model","fold_id"]).any()):
        raise ValueError("Saved 15,840 source pair model predictions/folds incomplete")
    summaries,comparison_folds=summarize(preds,fold,scopes)
    out.mkdir(parents=True,exist_ok=True)
    frames={
        "source_pair_community_panel.csv":source_pair_ledger,
        "source_eligible_geo_folds.csv":pd.DataFrame(fold_manifest),
        "model_fold_losses.csv":fold,
        "comparison_fold_losses.csv":comparison_folds,
        "comparison_scope_summary.csv":summaries,
    }
    for name,frame in frames.items():
        frame.to_csv(out/name,index=False,float_format="%.17g")
    with gzip.open(out/"heldout_predictions.csv.gz","wt",newline="",encoding="utf-8",compresslevel=6) as stream:
        preds.to_csv(stream,index=False,float_format="%.17g")
    result={
      "schema":"chocho_original51_2001_2010_plants_butterflies_marine_common_pairs_v1",
      "status":"FOUR_COMMON_PAIR_PANELS_15_FOLDS_ALL_PREDICTIONS_COMPLETE",
      "original_source_islands":51,
      "historical_plant_indexed_islands":int(sum(len(plants[k])>0 for k in ids)),
      "matched_common_source_pairs":int(eligible.sum()),
      "quality_min10_specieskeys_pairs":int(high.sum()),
      "source_folds":15,
      "original_WPI_ports":194,"marine_qualified_ports":182,
      "input_lock_status":original_lock["status"],
      "original_plant_snapshot_sha256":sha256(Path(source_dir)/"GBIF_original194_port_2001_2010_speciesKey.csv.gz"),
      "original_gbif_plant_receipt_sha256":sha256(Path(source_dir)/"source_receipt.json"),
      "source_botanical_effort_provenance":receipt,
      "panels":[p[0] for p in PANELS],"models":list(MODEL_NAMES),
      "source_panel_covariates":names,
      "saved_model_predictions":len(preds),"model_fold_rows":len(fold),
      "contrast_fold_rows":len(comparison_folds),"comparison_scope_rows":len(summaries),
      "source_scope_sizes":{k:int(v.sum()) for k,v in scopes.items()},
      "historical_plant_speciesKey_not_butterfly_host_or_verified_introduced":True,
      "original_AIS_2015_2021_postdates_2001_2010_plants_and_earlier_butterflies":True,
      "retrospective_host_panel_2001_2024_not_chronologically_safe":True,
      "no_causal_routes_or_colonization_claim":True,
      "out_sha256":{p.name:sha256(p) for p in out.iterdir() if p.is_file() and p.name!="execution_receipt.json"}
    }
    (out/"execution_receipt.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("CROSS_COMMUNITY_RESULT_JSON="+summaries.to_json(orient="records",double_precision=12),flush=True)
    print("CROSS_COMMUNITY_COUNTS_JSON="+json.dumps({
       "plants":45,"matched_pairs":990,"high_quality_pairs":741,"panels":4,
       "predictions":len(preds),"model_folds":len(fold),"contrast_folds":len(comparison_folds),
       "scope_rows":len(summaries),"scope_sizes":result["source_scope_sizes"]},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,required=True)
    ap.add_argument("--source-dir",type=Path,default=SOURCE)
    opt=ap.parse_args()
    run(opt.out_dir,opt.source_dir)
