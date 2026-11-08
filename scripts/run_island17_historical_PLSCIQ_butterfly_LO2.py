#!/usr/bin/env python3
"""Time-ordered original 2006-10 UNCTAD port INDEX versus later butterfly community.

Source population: 24 independent WPI PortWatch-verified ports, 22 ports
with >=16 of 20 old historical index quarters, 17 original island groups,
136 unique original island pairs. Exact two-test-islands-excluded LO2 fits.
No modern AIS, 2018 liner co-service or future-indexed HOSTS in main models.
"""
from __future__ import annotations
import argparse,gzip,hashlib,itertools,json
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from exploratory_island51_circuit_network_holdout import source,train_predict
from run_island51_host_circuit_joint import INPUT_PATHS,load_original_port_counts,validate_input_lock
from run_island51_latitude_controls import ROOT,ROSTER,PAIRS,latitude_pair_controls
from run_island51_pre2011_effort_control import early_botanical_features
from run_island51_observed_ocean_circuits import sha256

WPI17=ROOT/"results/island51_PortWatch17_geoverification_20261108/WPI17_independent_PortWatch_source_geography_LOCODE.csv"
WPI13=ROOT/"results/island51_PortWatch13_geoverification_20261108/WPI13_independent_PortWatch_source_geography_LOCODE.csv"
OLD=ROOT/"results/island51_PortWatch17_geoverification_20261108/independent_WPI17_source_audit.json"
NEW=ROOT/"results/island51_PortWatch13_geoverification_20261108/independent_WPI13_source_audit.json"
UNCTAD=ROOT/"results/island51_maritime_raw_recovery_20261108/UNCTAD_PLSCIQ_originalWPI194_name_country_candidates.csv"
RECEIPT=ROOT/"results/island51_maritime_raw_recovery_20261108/UNCTAD_PLSCIQ_source_recovery_receipt.json"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
MODELS=("BASE","MED_SUM","MED_SUM_DIFF","MAX_SUM_DIFF")
COMPARISONS=(
    ("BASE_to_MED_SUM","BASE","MED_SUM","historical_median_connectedness_sum"),
    ("BASE_to_MED_SUM_DIFF","BASE","MED_SUM_DIFF","historical_median_connectedness_sum_and_difference"),
    ("MED_SUM_to_MED_SUM_DIFF","MED_SUM","MED_SUM_DIFF","historical_median_connectedness_dissimilarity_conditional_sum"),
    ("BASE_to_MAX_SUM_DIFF","BASE","MAX_SUM_DIFF","historical_maximum_port_connectedness_sum_and_difference"),
    ("MED_SUM_DIFF_to_MAX_SUM_DIFF","MED_SUM_DIFF","MAX_SUM_DIFF","alternative_island_hub_index_nonnested"),
)
SCOPES=("all17","leave_original_landmass9_out")
ORIGINAL_WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
EXPECTED_UNCTAD_HASHES={
 "0":"29a343bb6b1113b1699a9ca1821952cce1cb8f150831921da68e248e78fb3e1b",
 "500":"c82e2d2119b0f4b03e4afb18c04f9c4442f19a2bfce96efbc8d30a87aed42933",
 "1000":"7eb1eaf4349008b784d9e6d2bd834cadf29c51fa493c2cdd2f83e2d767978255",
}

def original_source():
    locked=validate_input_lock(ROOT/"results/island51_host_circuit_joint_20261106/input_lock.json")
    q,ids,_=source(PAIRS,ROSTER)
    if len(ids)!=51 or len(q)!=1275:
        raise ValueError("Original source butterflies 51-island/1275-pair panel changed")
    preparation=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    counts=load_original_port_counts(ROSTER,ids,preparation)
    if sum(counts.values())!=194:
        raise ValueError("Original WPI194 per-island source port counts changed")
    return locked,q,ids,counts

def old_port_index_source(ids):
    wpibytes=WPI.read_bytes()
    old_blob=hashlib.sha1(b"blob "+str(len(wpibytes)).encode()+b"\x00"+wpibytes).hexdigest()
    if old_blob!=ORIGINAL_WPI_BLOB:
        raise ValueError("Original WPI source geography Git blob changed")
    historical=json.loads(RECEIPT.read_text())
    if not historical["source_version"]["is_exact_2026_10_30_original_recovered"]:
        raise ValueError("Old 2026-10-30 UNCTAD historical index edition not recovered")
    if historical["series"]!=1300 or historical["total_port_quarter_rows"]!=81900 or historical["nonmissing_port_quarter_rows"]!=57312:
        raise ValueError("Original 2006-2021 PLSCIQ source coverage changed")
    pages={str(r["offset"]):r for r in historical["page_receipts"]}
    if set(pages)!=set(EXPECTED_UNCTAD_HASHES):
        raise ValueError("Historical source page provenance offset set differs")
    for off,sha in EXPECTED_UNCTAD_HASHES.items():
        if pages[off]["sha256"]!=sha or not pages[off]["matches_original_2026_10_30_sha256"]:
            raise ValueError("UNCTAD historical edition SHA256 mismatch on page "+off)
    old=json.loads(OLD.read_text())
    new=json.loads(NEW.read_text())
    if old["status"]!="ALL17_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED" or new["status"]!="ALL13_FIXED_WPI_AND_BOUNDED_REGISTRY_IDENTITY_ROWS_RECONSTRUCTED":
        raise ValueError("Full independent source port geoidentity audits absent")
    df17=pd.read_csv(WPI17,dtype={"original_WPI_port_id":str,"original_island_id":str})
    df13=pd.read_csv(WPI13,dtype={"original_WPI_port_id":str,"original_island_id":str})
    checked=pd.concat([df17,df13],ignore_index=True)
    if len(checked)!=30 or checked.original_WPI_port_id.nunique()!=30:
        raise ValueError("Frozen original 30 source candidate port cohort differs")
    flags=checked.original_WPI_independent_geographic_identity_verified.astype(str).str.lower().eq("true")
    known=checked.loc[flags].copy()
    if len(known)!=24 or known.verified_PortWatch_portid.nunique()!=24 or known.original_island_id.nunique()!=17:
        raise ValueError("24 independent WPI port geographic identities not retained")
    if known.original_island_id.astype(str).isin(ids).all() is False:
        raise ValueError("Verified historical source port outside original 51 island source")
    historical_ports=pd.read_csv(UNCTAD,dtype={"original_WPI_id":str,"original_island_id":str})
    if len(historical_ports)!=194 or historical_ports.original_WPI_id.nunique()!=194:
        raise ValueError("Original UNCTAD WPI194 source rows missing")
    joined=known.merge(historical_ports,left_on="original_WPI_port_id",right_on="original_WPI_id",
        validate="one_to_one",indicator=True,suffixes=("","_legacy"))
    if len(joined)!=24 or not joined["_merge"].eq("both").all() or not joined.original_island_id.eq(joined.original_island_id_legacy).all():
        raise ValueError("Original independently geoverified ports not aligned to historical PLSCIQ source")
    qualified=joined.loc[
       joined.status.eq("UNIQUE_COUNTRY_AND_CANONICAL_NAME_CANDIDATE_NOT_GEO_VERIFIED") &
       (joined.available_2006_2010_quarters.astype(int)>=16)].copy()
    if (len(qualified)!=22 or qualified.original_WPI_port_id.nunique()!=22 or
        qualified.original_island_id.nunique()!=17 or
        qualified.mean_2006_2010_PLCSIQ_old_index_base.isna().any() or
        (qualified.mean_2006_2010_PLCSIQ_old_index_base.astype(float)<0).any()):
        raise ValueError("Historical time-safe index supports not exactly 22 WPI ports and 17 islands")
    if set(qualified.original_island_id).issubset(set(ids)) is False:
        raise ValueError("Original source island missing")
    return checked,joined,qualified

def controls(q,ids,original_port_counts):
    l=q.i.to_numpy(dtype=int);r=q.j.to_numpy(dtype=int)
    d=np.log1p(q.port_distance_km.to_numpy(dtype=float))
    port=np.log1p(np.asarray([original_port_counts[k] for k in ids],dtype=float))
    lat=latitude_pair_controls(q,ids)["abs_latitude_separation_degrees"]
    effort,coverage,src=early_botanical_features(q,ids,original_port_counts)
    baseline={
       "log_original_port_distance":d,
       "squared_log_original_port_distance":d*d,
       "sum_log1p_original_source_WPI_port_counts":port[l]+port[r],
       "abs_difference_log1p_original_source_WPI_port_counts":np.abs(port[l]-port[r]),
       "original_latitude_separation_degrees":lat,
       "sum_log1p_2001_2010_allplant_buffer_records":effort["sum_log1p_2001_2010_allplant_buffer_records"],
       "abs_difference_log1p_2001_2010_allplant_buffer_records":effort["abs_difference_log1p_2001_2010_allplant_buffer_records"],
    }
    if any(np.shape(x)!=(1275,) or not np.isfinite(x).all() for x in baseline.values()):
        raise ValueError("Original 1275 source-pair geographic and historic recording covariates incomplete")
    return baseline,src

def designs(q,ids,old_ports,qualified):
    features,record_receipt=controls(q,ids,old_ports)
    group=qualified.groupby("original_island_id").mean_2006_2010_PLCSIQ_old_index_base.agg(["median","max","count"])
    if len(group)!=17:
        raise ValueError("Historical source island aggregation cannot omit any index-supported island")
    aggregates={k:{"median":float(group.loc[k,"median"]),"max":float(group.loc[k,"max"]),
                   "port_count":int(group.loc[k,"count"])} for k in group.index}
    i=q.i.to_numpy(dtype=int);j=q.j.to_numpy(dtype=int)
    index_arrays={}
    for name in ("median","max"):
        values=np.array([np.log1p(aggregates[k][name]) if k in aggregates else np.nan for k in ids],dtype=float)
        left,right=values[i],values[j]
        index_arrays[name]={"sum":left+right,"absdiff":np.abs(left-right)}
    model_design={}
    feature_names={}
    for target in TARGETS:
        base=dict(features)
        if target=="new_butterfly_jaccard":
            base["previous_2001_2010_butterfly_Jaccard"]=q.earlier_butterfly_jaccard.to_numpy(dtype=float)
        others={
           "BASE":dict(base),
           "MED_SUM":dict(base,historical2006_2010_PLSCIQ_median_endpoint_sum=index_arrays["median"]["sum"]),
           "MED_SUM_DIFF":dict(base,historical2006_2010_PLSCIQ_median_endpoint_sum=index_arrays["median"]["sum"],
                               historical2006_2010_PLSCIQ_median_endpoint_absdiff=index_arrays["median"]["absdiff"]),
           "MAX_SUM_DIFF":dict(base,historical2006_2010_PLSCIQ_max_endpoint_sum=index_arrays["max"]["sum"],
                               historical2006_2010_PLSCIQ_max_endpoint_absdiff=index_arrays["max"]["absdiff"]),
        }
        model_design[target]={name:np.column_stack(list(cols.values())) for name,cols in others.items()}
        feature_names[target]={name:list(cols) for name,cols in others.items()}
        for mat in model_design[target].values():
            # Outside source 17 islands the historical index is intentionally
            # unavailable. Only the frozen eligible pair mask is used.
            if np.any(~np.isfinite(mat[np.isfinite(mat).all(axis=1)])):
                raise ValueError("Invalid source model predictor values")
    return model_design,feature_names,aggregates,record_receipt,index_arrays

def write_csv(path,frame):
    frame.to_csv(path,index=False,float_format="%.17g")

def run(outdir):
    out=Path(outdir).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("data","results","provenance")):
        raise ValueError("Never overwrite original GEB/source results")
    locked,q,ids,counts=original_source()
    checked,joined,qualified=old_port_index_source(ids)
    X,feature_names,aggr,historic_effort,idx_cols=designs(q,ids,counts,qualified)
    m=q.island_a.isin(aggr)&q.island_b.isin(aggr)
    m2=m&~q.island_a.eq("9")&~q.island_b.eq("9")
    if (int(m.sum())!=136 or int(m2.sum())!=120 or
        len(aggr)!=17 or sum(k!="9" for k in aggr)!=16):
        raise ValueError("Expected exact source 17 islands/136 and excluding group9 16 islands/120 pairs")
    portLedger=joined.copy()
    portLedger["historical_index_eligible_geoverified_WPI_source_port"]=portLedger.original_WPI_port_id.isin(set(qualified.original_WPI_port_id))
    portCols=["original_WPI_port_id","original_island_id","original_WPI_port_name",
      "verified_PortWatch_portid","verified_PortWatch_LOCODE",
      "original_WPI_to_PortWatch_distance_km","available_2006_2010_quarters",
      "mean_2006_2010_PLCSIQ_old_index_base",
      "historical_index_eligible_geoverified_WPI_source_port"]
    if not set(portCols).issubset(portLedger):
        raise ValueError("Original 24 source port/legacy UNCTAD index fields missing")
    islandLedger=[]
    for k in ids:
        a=aggr.get(k)
        islandLedger.append({"original_island_id":k,"original_source_WPI_port_count":counts[k],
            "qualified_geoverified_historical_index_ports":a["port_count"] if a else 0,
            "historical_2006_2010_median_index":a["median"] if a else None,
            "historical_2006_2010_max_index":a["max"] if a else None,
            "historical_index_eligible":bool(a),
            "original_group9_multi_country_aggregation":k=="9",
            "zero_source_record_not_no_ships":True})
    selected=q.loc[m].copy()
    pairLedger=pd.DataFrame({
        "original_pair_source_index":selected.index.astype(int),
        "original_island_a":selected.island_a.astype(str).to_numpy(),
        "original_island_b":selected.island_b.astype(str).to_numpy(),
        "new_butterfly_jaccard":selected.new_butterfly_jaccard.astype(float).to_numpy(),
        "earlier_butterfly_jaccard":selected.earlier_butterfly_jaccard.astype(float).to_numpy(),
        "original_source_distance_km":selected.port_distance_km.astype(float).to_numpy(),
        "historical_index_eligible_both_islands":True,
        "sensitivity_original_group9_excluded":m2[m].to_numpy(dtype=bool),
    })
    predictions=[];folds=[]
    comparison=[]
    for scope,mask in (("all17",m),("leave_original_landmass9_out",m2)):
        ids_scope=sorted(set(q.loc[mask,"island_a"])|set(q.loc[mask,"island_b"]),key=int)
        expected_n=17 if scope=="all17" else 16
        if len(ids_scope)!=expected_n:
            raise ValueError("Unexpected original source island scope size")
        eligible_original_index=np.flatnonzero(mask.to_numpy(dtype=bool))
        folds_for_scope={}
        for z in eligible_original_index:
            a,b=str(q.iloc[z].island_a),str(q.iloc[z].island_b)
            test=np.zeros(len(q),dtype=bool);test[z]=True
            train=mask.to_numpy(dtype=bool)&~q.island_a.isin((a,b)).to_numpy()&~q.island_b.isin((a,b)).to_numpy()
            expected_train=105 if scope=="all17" else 91
            if int(train.sum())!=expected_train or int(test.sum())!=1 or np.any(train&test):
                raise ValueError("Strict double-island-unseen fold training source island endpoints leaked")
            fid=f"{scope}:{a}:{b}"
            folds_for_scope[z]=(fid,train,test)
            folds.append({"scope":scope,"fold_id":fid,"original_pair_source_index":int(z),
                "test_island_a":a,"test_island_b":b,"train_pairs":int(train.sum()),
                "test_pairs":1,"both_test_endpoints_excluded_from_training":True})
        for target in TARGETS:
            true=q[target].to_numpy(dtype=float)
            for sourceidx,(fid,train,test) in folds_for_scope.items():
                for model in MODELS:
                    arr=X[target][model]
                    if not np.isfinite(arr[train]).all() or not np.isfinite(arr[test]).all():
                        raise ValueError("Missing historic index in selected train/test cohort")
                    pred=train_predict(arr,true,train,test)
                    if pred.shape!=(1,) or not np.isfinite(pred).all():
                        raise ValueError("Invalid double-island-unseen forecast")
                    y=float(true[sourceidx]);p=float(pred[0]);error=y-p
                    predictions.append({"scope":scope,"target":target,"model":model,
                      "fold_id":fid,"original_pair_source_index":int(sourceidx),
                      "test_island_a":str(q.iloc[sourceidx].island_a),
                      "test_island_b":str(q.iloc[sourceidx].island_b),
                      "train_pairs":int(train.sum()),"test_pairs":1,
                      "observed":y,"predicted":p,"squared_error":error*error,
                      "absolute_error":abs(error)})
    preds=pd.DataFrame(predictions)
    if len(preds)!=2048 or preds.duplicated(["scope","target","model","original_pair_source_index"]).any():
        raise ValueError("Expected 2 targets x 4 models x (136+120)=2048 unique heldout forecasts")
    for scope in SCOPES:
        for target in TARGETS:
            one=preds[preds.scope.eq(scope)&preds.target.eq(target)]
            for label,ref,aug,role in COMPARISONS:
                base=one[one.model.eq(ref)].set_index("original_pair_source_index").sort_index()
                other=one[one.model.eq(aug)].set_index("original_pair_source_index").sort_index()
                if len(base)!= (136 if scope=="all17" else 120) or not base.index.equals(other.index):
                    raise ValueError("Comparison target/fold source pairs not equal")
                r_mse=float(base.squared_error.mean());a_mse=float(other.squared_error.mean())
                r_mae=float(base.absolute_error.mean());a_mae=float(other.absolute_error.mean())
                comparison.append({"scope":scope,"target":target,"contrast":label,"role":role,
                  "original_islands":17 if scope=="all17" else 16,"test_unique_island_pairs":len(base),
                  "train_pairs_per_fold":105 if scope=="all17" else 91,
                  "MSE_reference":r_mse,"MSE_augmented":a_mse,
                  "MSE_reduction_percent":100*(r_mse-a_mse)/r_mse,
                  "MAE_reference":r_mae,"MAE_augmented":a_mae,
                  "MAE_reduction_percent":100*(r_mae-a_mae)/r_mae,
                  "pairs_with_lower_squared_error":int((base.squared_error>other.squared_error).sum()),
                  "pairs_with_higher_squared_error":int((base.squared_error<other.squared_error).sum()),
                  "pairs_with_lower_absolute_error":int((base.absolute_error>other.absolute_error).sum())})
    summary=pd.DataFrame(comparison)
    if len(summary)!=20:
        raise ValueError("Expected fixed 20 contrasts")
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/"original24_geoverified_ports_historical_PLSCIQ_quality.csv",portLedger[portCols])
    write_csv(out/"original51_historical_PLSCIQ_island_eligibility.csv",pd.DataFrame(islandLedger))
    write_csv(out/"original17_136_source_pair_manifest.csv",pairLedger)
    write_csv(out/"double_unseen_source_folds.csv",pd.DataFrame(folds))
    write_csv(out/"comparison_scope_summary.csv",summary)
    with gzip.open(out/"LO2_heldout_predictions.csv.gz","wt",newline="",encoding="utf-8") as stream:
        preds.to_csv(stream,index=False,float_format="%.17g")
    receipt={
      "schema":"chocho_historical_2006_2010_port_liner_PLSCIQ_butterfly_LO2_v1",
      "status":"ALL_SOURCE_FIXED_17_ISLANDS_AND_GROUP9_EXCLUDED_LO2_FORECASTS_COMPLETE",
      "source_original_islands":51,"original_WPI_ports":194,"source_verified_ports":24,
      "source_old_UNCTAD_index_ports_min16quarters":len(qualified),
      "source_index_supported_island_groups":len(aggr),
      "source_all17_pairs":136,"source_no_group9_pairs":120,
      "targets":list(TARGETS),"models":list(MODELS),
      "contrasts":[x[0] for x in COMPARISONS],"scopes":list(SCOPES),
      "individual_predictions":len(preds),"source_double_unseen_pair_folds":len(folds),
      "comparison_rows":len(summary),"feature_names":feature_names,
      "historical_2001_2010_botanical_record_effort_receipt":historic_effort,
      "old_UNCTAD_source_receipt_sha256":sha256(RECEIPT),
      "prior_source_original24_verified_port_files_sha256":{"WPI17":sha256(WPI17),"WPI13":sha256(WPI13)},
      "source_island9_mixed_country_grouping_kept_primary_and_dropped_only_in_fixed_sensitivity":True,
      "candidate_HOSTS_2001_2024_and_modernAIS_or_modelledGLSN_not_used_as_primary_predictors":True,
      "early_butterfly_2001_2010_target_contextual_not_forward_prediction":True,
      "PLSCIQ_is_port_index_not_ship_voyage_or_actual_butterfly_colonization":True,
      "all_file_sha256":{p.name:sha256(p) for p in out.iterdir() if p.is_file()}
    }
    (out/"execution_receipt.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("HISTORICAL_PORT_PLSCIQ_BUTTERFLY_LO2_RESULTS="+summary.to_json(orient="records",double_precision=12),flush=True)
    print("HISTORICAL_PORT_PLSCIQ_SOURCE_COUNTS="+json.dumps({
       "original_WPI_historical_port_index_rows":len(qualified),
       "source_supported_original_island_groups":len(aggr),
       "source_double_unseen_pair_models":len(preds),
       "source_pair_folds":len(folds),
       "number_of_comparisons":len(summary)},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    args=p.parse_args()
    run(args.out_dir)
