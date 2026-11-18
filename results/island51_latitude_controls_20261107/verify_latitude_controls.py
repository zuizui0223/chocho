#!/usr/bin/env python3
"""Independent source-pair audit of latitude-adjusted AIS/wind source-factor comparisons."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv"
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
BASE_MODELS=("HRE","HRE_AIS","HRE_WIND","HCRE")
MODELS=BASE_MODELS+tuple(m+"_LAT" for m in BASE_MODELS)
ORIGINAL_COMPARISONS={
    "HRE_to_HRE_AIS":("HRE","HRE_AIS"),
    "HRE_to_HRE_WIND":("HRE","HRE_WIND"),
    "HRE_AIS_to_HCRE":("HRE_AIS","HCRE"),
    "HRE_WIND_to_HCRE":("HRE_WIND","HCRE"),
    "HRE_to_HCRE":("HRE","HCRE"),
}
COMPARISONS=dict(ORIGINAL_COMPARISONS)
COMPARISONS.update({name+"_LAT":(ref+"_LAT",aug+"_LAT")
                    for name,(ref,aug) in ORIGINAL_COMPARISONS.items()})
COMPARISONS.update({"HRE_to_HRE_LAT":("HRE","HRE_LAT"),
                    "HCRE_to_HCRE_LAT":("HCRE","HCRE_LAT")})

SCOPES={"all_pairs":1275,"within_longitude_block":235,"between_longitude_blocks":1040}

def load_rows(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        result=list(csv.DictReader(f))
    if not result:
        raise ValueError("Empty source/result CSV "+str(path))
    return result

def close(actual,expected,context):
    if not math.isfinite(actual) or not math.isfinite(expected) or not math.isclose(
            actual,expected,rel_tol=1e-10,abs_tol=1e-11):
        raise ValueError(context+": recorded "+str(actual)+" differs from "+str(expected))

def audit(results_dir):
    directory=Path(results_dir)
    source=load_rows(SOURCE)
    if len(source)!=1275:
        raise ValueError("Original source pair count changed")
    ids={r["island_a"] for r in source}|{r["island_b"] for r in source}
    geo=load_rows(directory/"longitude_island_blocks.csv")
    blocks={}
    for row in geo:
        island=row["island_id"]
        group=int(row["longitude_block"])
        if island in blocks or group not in range(5):
            raise ValueError("Duplicate or invalid island longitude assignment")
        blocks[island]=group
    if len(geo)!=51 or set(blocks)!=ids or Counter(blocks.values())!=Counter({0:11,1:10,2:10,3:10,4:10}):
        raise ValueError("Geographic blocks do not preserve exactly51 original islands")
    ordered=sorted(geo,key=lambda r:(float(r["representative_longitude"]),int(r["island_id"])))
    if [int(x["longitude_block"]) for x in ordered]!=[0]*11+[1]*10+[2]*10+[3]*10+[4]*10:
        raise ValueError("Longitude blocks are not source-coordinate-only quintile assignment")
    partition={}
    fold=load_rows(directory/"fold_specification.csv")
    if len(fold)!=15:
        raise ValueError("Must have 15 block-pair holdouts")
    fold_meta={}
    for row in fold:
        k,h=int(row["block_a"]),int(row["block_b"])
        name=f"B{k}_B{h}"
        if k>h or k not in range(5) or h not in range(5) or row["fold_id"]!=name or name in fold_meta:
            raise ValueError("Block-pair fold names/identities disagree")
        held={x for x,g in blocks.items() if g in (k,h)}
        saved=set(json.loads(row["heldout_island_ids"]))
        if saved!=held or int(row["heldout_islands"])!=len(held) or int(row["training_islands"])!=51-len(held):
            raise ValueError("Fold heldout island IDs incorrect")
        training=[x for x in source if x["island_a"] not in held and x["island_b"] not in held]
        tests=[i for i,x in enumerate(source) if sorted([blocks[x["island_a"]],blocks[x["island_b"]]])==[k,h]]
        expected=len(tests)
        if (len(training)!=int(row["train_pairs"]) or expected!=int(row["test_pairs"]) or
            expected!=(math.comb(11 if k==0 else 10,2) if k==h else (11 if k==0 else 10)*(11 if h==0 else 10))):
            raise ValueError("Fold source support differs from exact geographic split")
        kind="within_longitude_block" if k==h else "between_longitude_blocks"
        if row["pair_class"]!=kind or not tests:
            raise ValueError("Block fold class wrong")
        fold_meta[name]={"class":kind,"train":len(training),"test":len(tests)}
        for i in tests:
            if i in partition:
                raise ValueError("Source pair assigned to multiple test folds")
            partition[i]=name
    if len(partition)!=1275 or sum(x["test"] for x in fold_meta.values())!=1275:
        raise ValueError("Not every source pair tested once")
    loss=defaultdict(list)
    by_pair={}
    with (directory/"heldout_pair_predictions.csv").open(newline="",encoding="utf-8") as f:
        n=0
        for row in csv.DictReader(f):
            target,model=row["target"],row["model"]
            i=int(row["source_pair_index"]); fid=row["fold_id"]
            if target not in TARGETS or model not in MODELS or i<0 or i>=1275:
                raise ValueError("Unknown model or source pair")
            original=source[i]
            if fid!=partition[i] or row["pair_class"]!=fold_meta[fid]["class"]:
                raise ValueError("Prediction fold identity wrong")
            if row["island_a"]!=original["island_a"] or row["island_b"]!=original["island_b"]:
                raise ValueError("Prediction endpoint identities changed")
            observed=float(row["observed"]);pred=float(row["predicted"])
            if not (0<=observed<=1 and 0<=pred<=1):
                raise ValueError("Observed/predicted Jaccard outside [0,1]")
            close(observed,float(original[target]),"original butterfly observation")
            sq=(pred-observed)**2
            ae=abs(pred-observed)
            close(sq,float(row["squared_error"]),"original squared error")
            close(ae,float(row["absolute_error"]),"original absolute error")
            key=(target,model,i)
            if key in by_pair:
                raise ValueError("Source pair has more than one model-specific prediction")
            by_pair[key]={"sq":sq,"ae":ae,"class":row["pair_class"],"fold":fid}
            loss[(target,model,fid)].append((sq,ae))
            n+=1
    if n!=20400 or len(by_pair)!=20400 or len(loss)!=240 or any(len(v)!=fold_meta[k[2]]["test"] for k,v in loss.items()):
        raise ValueError("Incomplete 7650 predictions or 90 model-fold groupings")
    fold_rows=load_rows(directory/"model_fold_losses.csv")
    checked=set()
    for row in fold_rows:
        key=(row["target"],row["model"],row["fold_id"])
        if key in checked or key not in loss:
            raise ValueError("Duplicate model/fold")
        checked.add(key)
        vals=loss[key]
        if (int(row["test_pairs"])!=len(vals) or int(row["train_pairs"])!=fold_meta[key[2]]["train"] or
            row["pair_class"]!=fold_meta[key[2]]["class"]):
            raise ValueError("Model fold sizes wrong")
        close(math.fsum(z[0] for z in vals)/len(vals),float(row["MSE"]),"fold MSE")
        close(math.fsum(z[1] for z in vals)/len(vals),float(row["MAE"]),"fold MAE")
    if len(checked)!=240:
        raise ValueError("240 model folds not audited")
    comparison_rows=load_rows(directory/"contrast_fold_losses.csv")
    checked_contrast=set()
    for row in comparison_rows:
        name=row["comparison"];target=row["target"];fid=row["fold_id"]
        if target not in TARGETS or name not in COMPARISONS or fid not in fold_meta:
            raise ValueError("Unexpected contrast fold")
        key=(target,name,fid)
        if key in checked_contrast:
            raise ValueError("Duplicate contrast fold")
        checked_contrast.add(key)
        ref,aug=COMPARISONS[name]
        xa=loss[(target,ref,fid)]; xb=loss[(target,aug,fid)]
        r_mse=math.fsum(z[0] for z in xa)/len(xa)
        a_mse=math.fsum(z[0] for z in xb)/len(xb)
        r_mae=math.fsum(z[1] for z in xa)/len(xa)
        a_mae=math.fsum(z[1] for z in xb)/len(xb)
        close(r_mse,float(row["MSE_reference"]),"fold comparison MSE ref")
        close(a_mse,float(row["MSE_augmented"]),"fold comparison MSE aug")
        close(r_mse-a_mse,float(row["MSE_gain"]),"fold comparison MSE gain")
        close(r_mae,float(row["MAE_reference"]),"fold comparison MAE ref")
        close(a_mae,float(row["MAE_augmented"]),"fold comparison MAE aug")
        close(r_mae-a_mae,float(row["MAE_gain"]),"fold comparison MAE gain")
        if ((row["MSE_improved"].lower()=="true")!=(r_mse>a_mse) or
            row["pair_class"]!=fold_meta[fid]["class"]):
            raise ValueError("Fold improvement sign/class mismatched")
    if len(checked_contrast)!=360:
        raise ValueError("Expected all 360 contrast-fold records")
    summaries=load_rows(directory/"contrast_scope_summary.csv")
    summary_checked=set()
    for row in summaries:
        target,kind,name=row["target"],row["pair_class"],row["contrast"]
        key=(target,kind,name)
        if target not in TARGETS or kind not in SCOPES or name not in COMPARISONS or key in summary_checked:
            raise ValueError("Unexpected or duplicated summary scope")
        summary_checked.add(key)
        ref,aug=COMPARISONS[name]
        indices=[i for i in range(1275) if kind=="all_pairs" or by_pair[(target,ref,i)]["class"]==kind]
        if len(indices)!=SCOPES[kind] or int(row["unique_evaluated_pairs"])!=len(indices):
            raise ValueError("Summary pair population differs")
        x=[by_pair[(target,ref,i)] for i in indices]
        y=[by_pair[(target,aug,i)] for i in indices]
        r_mse=math.fsum(t["sq"] for t in x)/len(x)
        a_mse=math.fsum(t["sq"] for t in y)/len(y)
        r_mae=math.fsum(t["ae"] for t in x)/len(x)
        a_mae=math.fsum(t["ae"] for t in y)/len(y)
        for col,val in (
            ("MSE_reference",r_mse),("MSE_augmented",a_mse),
            ("MSE_reduction_percent",100*(r_mse-a_mse)/r_mse),
            ("MAE_reference",r_mae),("MAE_augmented",a_mae),
            ("MAE_reduction_percent",100*(r_mae-a_mae)/r_mae),
            ("pairs_with_lower_squared_error",sum(m["sq"]>n["sq"] for m,n in zip(x,y))),
            ("pairs_with_higher_squared_error",sum(m["sq"]<n["sq"] for m,n in zip(x,y))),
            ("pairs_with_lower_absolute_error",sum(m["ae"]>n["ae"] for m,n in zip(x,y))),
        ):
            close(val,float(row[col]),"summary "+str(key)+" "+col)
    if len(summary_checked)!=72:
        raise ValueError("All72 comparisons/scopes must be saved")
    result={"schema":"chocho_double_unseen_AIS_wind_latitude_independent_audit_v1",
      "status":"ALL_1275_PAIR_IDENTITIES_AND_ALL_DOUBLE_UNSEEN_SCORES_RECONSTRUCTED",
      "geographic_island_count":len(blocks),"unique_source_pair_count":len(partition),
      "longitude_fold_count":len(fold_meta),"saved_predictions_reconstructed":len(by_pair),
      "model_folds_reconstructed":len(checked),"contrast_folds_reconstructed":len(checked_contrast),
      "scope_summaries_reconstructed":len(summary_checked),
      "source_butterfly_observations_verified":2*1275,
      "no_test_pair_endpoint_in_its_training_population":True}
    (directory/"independent_audit.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("LATITUDE_CONTROL_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--results-dir",type=Path,required=True)
    audit(p.parse_args().results_dir)
