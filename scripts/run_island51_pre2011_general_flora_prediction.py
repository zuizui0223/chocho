#!/usr/bin/env python3
"""Source-pinned all-51-island historical GENERAL Plantae composition holdout prediction.

Runs on COMPLETE original194 WPI/GBIF speciesKey acquisition, not restricted HOSTS edges.
Do not interpret early all-plant indexed species as larval-host presence/absence.
"""
from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from exploratory_island51_circuit_network_holdout import source, train_predict
from run_island51_host_circuit_joint import (
    INPUT_PATHS, TARGETS, load_original_port_counts,
    load_resistance_matrix, validate_input_lock
)
from run_island51_latitude_controls import (
    ROOT, ROSTER, PAIRS, GEOGRAPHY, make_features,
    latitude_pair_controls, partition_geographic_blocks, holdout_plan
)
from run_island51_pre2011_effort_control import git_blob_sha
from run_island51_observed_ocean_circuits import sha256
from acquire_island51_pre2011_public_flora import (
    original_source_rows, WPI_BLOB, ARCHIVED_BLOB,
)

MODELS=("HRE_LAT","HRE_AIS_LAT","HB","HBA","HBF","HBFA")
CONTRASTS=(
    ("HRE_to_HRE_AIS_LAT","HRE_LAT","HRE_AIS_LAT","original_AIS_reproduction"),
    ("HRE_LAT_to_HB","HRE_LAT","HB","recorded_early_flora_richness_and_coverage"),
    ("HB_to_HBF","HB","HBF","early_general_flora_composition_without_AIS"),
    ("HBA_to_HBFA","HBA","HBFA","early_general_flora_composition_conditional_AIS"),
    ("HB_to_HBA","HB","HBA","AIS_conditional_early_flora_richness"),
    ("HBF_to_HBFA","HBF","HBFA","AIS_conditional_early_general_flora_composition"),
)
SCOPES=("all_pairs","within_longitude_block","between_longitude_blocks","both_early_species_indexed")
PARENT=ROOT/"results/island51_latitude_controls_20261107/contrast_scope_summary.csv"
PARENT_BLOB="0d7189c9a7f538116c9ba4eaabe52644126f77c1"

def early_general_flora_source(source_dir,original_ids):
    path=Path(source_dir)
    old_geo,old_counts=original_source_rows()
    raw=path/"GBIF_original194_port_2001_2010_speciesKey.csv.gz"
    site=path/"GBIF_original194_port_2001_2010_coverage.csv"
    islands=path/"GBIF_original51_early_general_flora_islands.csv"
    receipt_file=path/"source_receipt.json"
    receipt=json.loads(receipt_file.read_text(encoding="utf-8"))
    if (receipt.get("schema")!="chocho_original194_direct_public_GBIF_Plant_speciesKey_2001_2010_v1" or
        receipt.get("status")!="COMPLETE_ALL_194_SOURCE_PORTS_INDEPENDENT_GBIF_LIVE_ACQUISITION" or
        receipt.get("source_WPI_Git_blob")!=WPI_BLOB or
        receipt.get("original_archive_Git_blob")!=ARCHIVED_BLOB or
        receipt.get("source_port_original_geo_count")!=194 or
        receipt.get("source_eligible_islands")!=51):
        raise ValueError("Full original source-pinned public GBIF plant history unavailable")
    for p in (raw,site,islands):
        if sha256(p)!=receipt["source_sha256"].get(p.name):
            raise ValueError("Saved live GBIF source file content checksum mismatch "+p.name)
    coverage=pd.read_csv(site,dtype={"island_id":str,"port_id":int})
    indexed=pd.read_csv(raw,dtype={"island_id":str,"port_id":int,"GBIF_speciesKey":int,"indexed_records_2001_2010":int})
    island_table=pd.read_csv(islands,dtype={"island_id":str})
    if (len(coverage)!=194 or coverage.port_id.nunique()!=194 or
        {int(k) for k in coverage.port_id}!={int(k) for k in old_geo} or
        len(island_table)!=51 or island_table.island_id.nunique()!=51 or
        set(island_table.island_id)!=set(original_ids)):
        raise ValueError("Original complete GBIF source port/island cohort not retained")
    for r in coverage.itertuples(index=False):
        p=int(r.port_id)
        if (r.island_id!=old_geo[p]["island_id"] or
            int(r.original_early_port_buffer_records)!=old_counts[p] or
            int(r.live_GBIF_2001_2010_plant_records)<0):
            raise ValueError("GBIF plant count panel changed source port identity or recorded count")
    if (len(indexed)!=int(receipt["total_port_speciesKey_rows"]) or
        indexed.duplicated(["port_id","GBIF_speciesKey"]).any() or
        (indexed.GBIF_speciesKey<=0).any() or
        (indexed.indexed_records_2001_2010<=0).any()):
        raise ValueError("Malformed / incomplete original194 GBIF speciesKey facets")
    checked=indexed.merge(coverage[["port_id","island_id","live_GBIF_2001_2010_plant_records"]],
        on=["port_id","island_id"],how="left",validate="many_to_one",indicator=True)
    if checked["_merge"].ne("both").any():
        raise ValueError("Source historic plant taxon outside original WPI port/island")
    species_counts=checked.groupby("port_id")["indexed_records_2001_2010"].sum().to_dict()
    for r in coverage.itertuples(index=False):
        if species_counts.get(r.port_id,0)>int(r.live_GBIF_2001_2010_plant_records):
            raise ValueError("GBIF species-facet count exceeds total occurrence index for a source port")
    plants={i:set() for i in original_ids}
    for row in indexed.itertuples(index=False):
        plants[row.island_id].add(int(row.GBIF_speciesKey))
    richness={sid:len(plants[sid]) for sid in original_ids}
    by_port=coverage.groupby("island_id").port_id.count().to_dict()
    for r in island_table.itertuples(index=False):
        if (r.island_id not in plants or
            int(r.original_port_count)!=int(by_port[r.island_id]) or
            int(r.unique_early_indexed_general_plant_speciesKeys)!=richness[r.island_id]):
            raise ValueError("Full source island early flora species union/richness mismatch")
    if sum(len(v)==0 for v in plants.values())!=int(receipt["islands_with_no_speciesKey_not_true_absence"]):
        raise ValueError("General historical plant no-index support changed")
    return plants,coverage,island_table,receipt

def general_flora_features(q,ids,sets):
    if len(ids)!=51 or set(ids)!=set(sets):
        raise ValueError("General early flora coverage must preserve all original51 island identities")
    species_count=np.array([len(sets[k]) for k in ids],dtype=float)
    la=species_count[q.i.to_numpy(dtype=int)]
    lb=species_count[q.j.to_numpy(dtype=int)]
    log_a=np.log1p(la);log_b=np.log1p(lb)
    empty_a=la==0;empty_b=lb==0
    jaccard=[]
    for a,b in zip(q.i.to_numpy(dtype=int),q.j.to_numpy(dtype=int)):
        one=sets[ids[a]];two=sets[ids[b]]
        den=len(one|two)
        jaccard.append(float(len(one&two)/den) if den else 0.0)
    jaccard=np.asarray(jaccard,dtype=float)
    controls={
       "sum_log1p_2001_2010_general_plant_specieskey_richness":log_a+log_b,
       "abs_difference_log1p_2001_2010_general_plant_specieskey_richness":np.abs(log_a-log_b),
       "only_one_endpoint_no_2001_2010_plant_specieskey":np.logical_xor(empty_a,empty_b).astype(float),
       "both_endpoints_no_2001_2010_plant_specieskey":(empty_a&empty_b).astype(float),
    }
    masks={"both_early_species_indexed":(~empty_a)&(~empty_b)}
    if (any(not np.isfinite(x).all() or x.shape!=(1275,) for x in [jaccard,*controls.values()]) or
        np.any(jaccard<0) or np.any(jaccard>1)):
        raise ValueError("General flora historical pair covariates missing or out-of-range")
    return controls,jaccard,masks

def build_models(q,ids,port_counts,matrices,plant_sets):
    original,names,historical=make_features(q,ids,port_counts,matrices)
    latitude=latitude_pair_controls(q,ids)
    new,floral_jaccard,masks=general_flora_features(q,ids,plant_sets)
    designs={}
    feature_names={}
    for target in TARGETS:
        base=np.column_stack([original[target]["HRE"]]+[latitude[k] for k in latitude])
        ais=np.column_stack([original[target]["HRE_AIS"]]+[latitude[k] for k in latitude])
        new_cols=np.column_stack(list(new.values()))
        names0=names[target]["HRE"]+list(latitude)
        namesA=names[target]["HRE_AIS"]+list(latitude)
        if (len(names0)!=base.shape[1] or len(namesA)!=ais.shape[1] or
            namesA.count("log1p_normalized_ais_circuit")!=1):
            raise ValueError("Original source model feature names/columns differ")
        hb=np.column_stack([base,new_cols])
        hba=np.column_stack([ais,new_cols])
        hbf=np.column_stack([hb,floral_jaccard])
        hbfa=np.column_stack([hba,floral_jaccard])
        designs[target]={"HRE_LAT":base,"HRE_AIS_LAT":ais,"HB":hb,"HBA":hba,
                         "HBF":hbf,"HBFA":hbfa}
        feature_names[target]={
           "HRE_LAT":names0,"HRE_AIS_LAT":namesA,
           "HB":names0+list(new),"HBA":namesA+list(new),
           "HBF":names0+list(new)+["2001_2010_general_flora_Jaccard_imputed"],
           "HBFA":namesA+list(new)+["2001_2010_general_flora_Jaccard_imputed"],
        }
        if any(x.shape[0]!=1275 or not np.isfinite(x).all() for x in designs[target].values()):
            raise ValueError("General flora source model design incomplete")
    return designs,feature_names,masks,historical

def comparisons_table(pred,fold,sample_masks):
    records=[]
    per_fold=[]
    for target in TARGETS:
        a=pred[pred.target.eq(target)]
        f=fold[fold.target.eq(target)]
        for name,reference,augmented,role in CONTRASTS:
            x=a[a.model.eq(reference)].set_index("source_pair_index",verify_integrity=True)
            y=a[a.model.eq(augmented)].set_index("source_pair_index",verify_integrity=True).reindex(x.index)
            if len(x)!=1275 or len(y)!=1275 or y.predicted.isna().any():
                raise ValueError("Original complete source pair prediction coverage missing")
            xx=x.sort_index();yy=y.sort_index()
            if not np.array_equal(xx.source_pair_index if "source_pair_index" in xx.columns else xx.index,yy.index):
                raise ValueError("Source pair comparison ID order differs")
            reference_sq=xx.squared_error.to_numpy(dtype=float)
            target_sq=yy.squared_error.to_numpy(dtype=float)
            reference_ae=xx.absolute_error.to_numpy(dtype=float)
            target_ae=yy.absolute_error.to_numpy(dtype=float)
            scope_flags={
               "all_pairs":np.ones(1275,dtype=bool),
               "within_longitude_block":xx.pair_class.eq("within_longitude_block").to_numpy(),
               "between_longitude_blocks":xx.pair_class.eq("between_longitude_blocks").to_numpy(),
               "both_early_species_indexed":sample_masks["both_early_species_indexed"],
            }
            for scope in SCOPES:
                mask=scope_flags[scope]
                n=int(mask.sum())
                if n<=0:
                    raise ValueError("Source-only botanical quality scope has zero eligible original pair comparisons")
                r_mse=float(reference_sq[mask].mean());n_mse=float(target_sq[mask].mean())
                r_mae=float(reference_ae[mask].mean());n_mae=float(target_ae[mask].mean())
                if r_mse<=0 or r_mae<=0:
                    raise ValueError("Full reference model has zero denominator")
                records.append({"target":target,"contrast":name,"role":role,"pair_class":scope,
                    "unique_evaluated_pairs":n,"MSE_reference":r_mse,"MSE_augmented":n_mse,
                    "MSE_reduction_percent":100*(r_mse-n_mse)/r_mse,
                    "MAE_reference":r_mae,"MAE_augmented":n_mae,
                    "MAE_reduction_percent":100*(r_mae-n_mae)/r_mae,
                    "pairs_with_lower_squared_error":int(np.sum(reference_sq[mask]>target_sq[mask])),
                    "pairs_with_higher_squared_error":int(np.sum(reference_sq[mask]<target_sq[mask])),
                    "pairs_with_lower_absolute_error":int(np.sum(reference_ae[mask]>target_ae[mask])),
                })
            left=f[f.model.eq(reference)].set_index("fold_id",verify_integrity=True)
            right=f[f.model.eq(augmented)].set_index("fold_id",verify_integrity=True).reindex(left.index)
            if len(left)!=15 or len(right)!=15 or right.MSE.isna().any():
                raise ValueError("Missing fixed geographic training fold")
            for foldid,row in left.iterrows():
                other=right.loc[foldid]
                if row.test_pairs!=other.test_pairs or row.train_pairs!=other.train_pairs:
                    raise ValueError("Compared models lack matching source fold support")
                per_fold.append({"target":target,"contrast":name,"fold_id":foldid,
                   "pair_class":row.pair_class,"train_pairs":int(row.train_pairs),
                   "test_pairs":int(row.test_pairs),
                   "MSE_reference":float(row.MSE),"MSE_augmented":float(other.MSE),
                   "MSE_gain":float(row.MSE-other.MSE),
                   "MAE_reference":float(row.MAE),"MAE_augmented":float(other.MAE),
                   "MAE_gain":float(row.MAE-other.MAE)})
    if len(records)!=48 or len(per_fold)!=180:
        raise ValueError("Expected 48 target/contrast/scope rows and 180 fold comparisons")
    return pd.DataFrame(records),pd.DataFrame(per_fold)

def reproduce_parent(summaries):
    if git_blob_sha(PARENT.read_bytes())!=PARENT_BLOB:
        raise ValueError("Original PR134 published AIS source model changed")
    archived=pd.read_csv(PARENT)
    old=archived[archived.contrast.eq("HRE_to_HRE_AIS_LAT")]
    here=summaries[summaries.contrast.eq("HRE_to_HRE_AIS_LAT")&
                   ~summaries.pair_class.eq("both_early_species_indexed")]
    if len(old)!=6 or len(here)!=6:
        raise ValueError("Expected original six PR134 AIS true-source contrasts")
    maximum=0.
    for record in here.itertuples(index=False):
        original=old[(old.target==record.target)&
                     (old.pair_class==record.pair_class)]
        if len(original)!=1 or int(original.iloc[0].unique_evaluated_pairs)!=record.unique_evaluated_pairs:
            raise ValueError("Original AIS baseline source pair cohort differs")
        for field in ("MSE_reference","MSE_augmented","MAE_reference","MAE_augmented"):
            error=abs(float(original.iloc[0][field])-float(getattr(record,field)))
            maximum=max(maximum,error)
            if error>1e-9:
                raise ValueError("Original AIS MSE/MAE could not be reproduced: "+field)
    return maximum

def run(source_dir,out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/p).resolve()) for p in ("data","provenance","results")):
        raise ValueError("Do not overwrite original GEB source/results")
    validate=validate_input_lock(ROOT/"results/island51_host_circuit_joint_20261106/input_lock.json")
    q,ids,_=source(PAIRS,ROSTER)
    if len(q)!=1275 or len(ids)!=51:
        raise ValueError("Original 51 island pairs source not preserved")
    flora,coverage,island_table,source_receipt=early_general_flora_source(source_dir,ids)
    preparation=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    original_ports=load_original_port_counts(ROSTER,ids,preparation)
    matrices={k:load_resistance_matrix(ROOT/INPUT_PATHS[k],ids)
              for k in ("uniform","ais","wind")}
    geographic,where=partition_geographic_blocks(GEOGRAPHY,ids)
    folds,fold_support=holdout_plan(q,ids,geographic)
    designs,feature_names,masks,historical=build_models(q,ids,original_ports,matrices,flora)
    source_pair_category=np.array(["within_longitude_block" if geographic[str(r.island_a)]==geographic[str(r.island_b)]
         else "between_longitude_blocks" for r in q.itertuples(index=False)],dtype=object)
    if int(np.sum(source_pair_category=="within_longitude_block"))!=235:
        raise ValueError("Original longitude-block classes changed")
    predictions=[];losses=[]
    a=q.island_a.to_numpy(dtype=str);b=q.island_b.to_numpy(dtype=str)
    for target in TARGETS:
        observed=q[target].to_numpy(dtype=float)
        for fid,train,test in folds:
            test_indices=np.flatnonzero(test)
            kind=str(fold_support.loc[fold_support.fold_id.eq(fid),"pair_class"].iloc[0])
            if not np.all(source_pair_category[test]==kind):
                raise ValueError("Original 15 folds do not match geography-only pair classes")
            for model in MODELS:
                estimate=train_predict(designs[target][model],observed,train,test)
                if estimate.shape!=(len(test_indices),) or not np.isfinite(estimate).all():
                    raise ValueError("Missing original source pair model predictions")
                r=(observed[test]-estimate)
                sq=r**2;ae=np.abs(r)
                losses.append({"target":target,"model":model,"fold_id":fid,
                    "pair_class":kind,"train_pairs":int(train.sum()),
                    "test_pairs":int(test.sum()),"MSE":float(sq.mean()),
                    "MAE":float(ae.mean())})
                for k,p,actual,error,absolute in zip(test_indices,estimate,observed[test],sq,ae):
                    predictions.append({"target":target,"model":model,"fold_id":fid,
                        "pair_class":kind,"source_pair_index":int(k),
                        "island_a":a[k],"island_b":b[k],
                        "early_general_flora_indexed_both_endpoints":bool(masks["both_early_species_indexed"][k]),
                        "observed":float(actual),"predicted":float(p),
                        "squared_error":float(error),"absolute_error":float(absolute)})
    results=pd.DataFrame(predictions)
    fold=pd.DataFrame(losses)
    if (len(results)!=15300 or len(fold)!=180 or
        results.duplicated(["target","model","source_pair_index"]).any()):
        raise ValueError("Original51 source prediction/model cohort incomplete")
    summary,comparison_folds=comparisons_table(results,fold,masks)
    maxdiff=reproduce_parent(summary)
    out.mkdir(parents=True,exist_ok=True)
    with gzip.open(out/"heldout_predictions.csv.gz","wt",newline="",encoding="utf-8") as f:
        results.to_csv(f,index=False,float_format="%.17g")
    file_outputs={
        "model_fold_losses.csv":fold,
        "comparison_fold_losses.csv":comparison_folds,
        "comparison_scope_summary.csv":summary,
        "original51_early_general_flora_coverage.csv":island_table,
        "source_geographic_folds.csv":fold_support,
    }
    for name,frame in file_outputs.items():
        frame.to_csv(out/name,index=False,float_format="%.17g")
    receipt={
        "schema":"chocho_original51_early_general_flora_double_unseen_prediction_v1",
        "status":"ALL_51_ISLANDS_15_DOUBLE_UNSEEN_FOLDS_PREDICTED",
        "source_input_lock_status":validate["status"],
        "source_GBIF_source_receipt_sha256":sha256(Path(source_dir)/"source_receipt.json"),
        "source_Plant_indexed_species_file_sha256":sha256(Path(source_dir)/"GBIF_original194_port_2001_2010_speciesKey.csv.gz"),
        "source_GBIF_index_mismatch_ports":source_receipt["port_count_mismatches"],
        "original_WPI_ports":194,"marine_qualified_ports":182,
        "islands":51,"source_pair_unique":1275,"source_geography_folds":15,
        "source_early_general_flora_zero_index_islands":sum(not flora[k] for k in ids),
        "source_pairs_with_both_early_taxa_indexed":int(masks["both_early_species_indexed"].sum()),
        "models":MODELS,"comparisons":[x[0] for x in CONTRASTS],
        "predictions":len(results),"model_fold_rows":len(fold),
        "comparison_fold_rows":len(comparison_folds),"scope_rows":len(summary),
        "original_AIS_source_baseline_reproduced_max_error":maxdiff,
        "feature_names":feature_names,
        "original_host_candidates_2001_2024_not_temporally_safe":True,
        "early_taxa_are_GBIF_general_Plantae_not_confirmed_host_species":True,
        "GBIF_2001_2010_live_query_not_exact_original_archive":True,
        "2001_2010_butterflies_and_general_flora_contemporaneous_not_temporal_precedence":True,
        "zero_index_imputed_and_explicitly_controlled_not_ecological_absence":True,
        "source_GBIF_plant_status_not_reproduction_or_ship_transport":True,
        "out_sha256":{p.name:sha256(p) for p in out.iterdir() if p.is_file() and p.name!="execution_receipt.json"},
    }
    (out/"execution_receipt.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("EARLY_FLORA_RESULT_JSON="+summary.to_json(orient="records",double_precision=12),flush=True)
    print("EARLY_FLORA_CHECK_JSON="+json.dumps({
       "predictions":len(results),"model_folds":len(fold),
       "comparisons":len(comparison_folds),"scopes":len(summary),
       "parent_AIS_max_error":maxdiff,"source_ports":194,
       "islands_with_early_specieskey":51-sum(not flora[k] for k in ids),
       "pair_both_recorded":int(masks["both_early_species_indexed"].sum())
    },sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    a=parser.parse_args()
    run(a.source_dir,a.out_dir)
