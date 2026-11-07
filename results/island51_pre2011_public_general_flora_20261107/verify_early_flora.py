#!/usr/bin/env python3
"""Independent Python-stdlib audit of GBIF194 early general flora and 51-island predictions.

No pandas/numpy, no access to the fitting or GBIF-acquisition module. Rebuild source
speciesKey unions, historic pair covariates, all fold memberships, held-out errors,
and all 48 model contrasts from saved raw source tables and original butterfly pairs.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
ROSTER=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv"
PAIR_SOURCE=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv"
GEO=ROOT/"results/island51_host_circuit_atlas_20261106/island_geography.csv"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
EARLY=ROOT/"provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_plant_availability.csv"
PARENT=ROOT/"results/island51_latitude_controls_20261107/contrast_scope_summary.csv"
MODELS=("HRE_LAT","HRE_AIS_LAT","HB","HBA","HBF","HBFA")
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
COMPARE={
    "HRE_to_HRE_AIS_LAT":("HRE_LAT","HRE_AIS_LAT"),
    "HRE_LAT_to_HB":("HRE_LAT","HB"),
    "HB_to_HBF":("HB","HBF"),
    "HBA_to_HBFA":("HBA","HBFA"),
    "HB_to_HBA":("HB","HBA"),
    "HBF_to_HBFA":("HBF","HBFA"),
}
CATEGORIES=("all_pairs","within_longitude_block","between_longitude_blocks","both_early_species_indexed")

def read(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def read_gz(path):
    with gzip.open(path,"rt",newline="",encoding="utf-8") as f:
        yield from csv.DictReader(f)

def equal(a,b,detail):
    if not math.isfinite(float(a)) or not math.isfinite(float(b)) or not math.isclose(
          float(a),float(b),rel_tol=1e-10,abs_tol=1e-10):
        raise ValueError(f"Numeric source mismatch {detail}: {a} vs {b}")

def original():
    roster=[r for r in read(ROSTER) if r["eligible"].lower()=="true"]
    ids={r["landmass_id"] for r in roster}
    ports=read(WPI)
    early=read(EARLY)
    if len(ids)!=51 or len(ports)!=194 or len(early)!=194:
        raise ValueError("Original 51-island/194-port source changed")
    pids={}
    for r in ports:
        pid=int(r["port_id"]);island=r["island_id"]
        if pid in pids or island!=r["landmass_id"] or island not in ids:
            raise ValueError("Original WPI port identity/landmass disagree")
        pids[pid]=island
    original_early={}
    for r in early:
        pid=int(r["port_id"])
        if pid in original_early or pid not in pids or pids[pid]!=r["landmass_id"]:
            raise ValueError("Original year-indexed source port identity mismatch")
        original_early[pid]=int(r["early_plant_records"])
    if len(pids)!=194 or len(original_early)!=194 or sum(original_early.values())!=155002:
        raise ValueError("Archived early total 155002 index missing")
    pairs=read(PAIR_SOURCE)
    if len(pairs)!=1275:
        raise ValueError("Original 1275 butterfly/source pairs missing")
    seen=set()
    for pair in pairs:
        a,b=pair["island_a"],pair["island_b"]
        if a not in ids or b not in ids or a==b:
            raise ValueError("Source butterfly pair unknown island")
        key=tuple(sorted((a,b),key=int))
        if key in seen:
            raise ValueError("Duplicated original unordered butterfly island pair")
        seen.add(key)
        for t in TARGETS:
            y=float(pair[t])
            if not 0<=y<=1:
                raise ValueError("Original butterfly target outside Jaccard range")
    geo=read(GEO)
    if len(geo)!=51 or {x["island_id"] for x in geo}!=ids:
        raise ValueError("Original source longitude cohort changed")
    ordered=sorted(geo,key=lambda r:(float(r["representative_longitude"]),int(r["island_id"])))
    blocks={}
    for index,r in enumerate(ordered):
        blocks[r["island_id"]]=0 if index<11 else 1+(index-11)//10
    if Counter(blocks.values())!=Counter({0:11,1:10,2:10,3:10,4:10}):
        raise ValueError("Original five longitude blocks changed")
    fids={}
    category={}
    folds={}
    for i,pair in enumerate(pairs):
        k,h=sorted((blocks[pair["island_a"]],blocks[pair["island_b"]]))
        fids[i]=f"B{k}_B{h}"
        category[i]="within_longitude_block" if k==h else "between_longitude_blocks"
    for k in range(5):
        for h in range(k,5):
            fid=f"B{k}_B{h}"
            held={sid for sid,g in blocks.items() if g in (k,h)}
            test={i for i in fids if fids[i]==fid}
            train=sum(r["island_a"] not in held and r["island_b"] not in held for r in pairs)
            if len(test)!=(math.comb(11 if k==0 else 10,2) if k==h else
                           (11 if k==0 else 10)*(11 if h==0 else 10)):
                raise ValueError("Original test pair count changed")
            if train!=math.comb(51-len(held),2):
                raise ValueError("Original geographic heldout source training island leaked")
            folds[fid]={"n_test":len(test),"n_train":train,"class":
                        "within_longitude_block" if k==h else "between_longitude_blocks"}
    if Counter(category.values())!=Counter({"within_longitude_block":235,"between_longitude_blocks":1040}):
        raise ValueError("Expected original 235/1040 source pair classes")
    return ids,pids,original_early,pairs,fids,category,folds

def source_flora(source_dir,ids,pids,original_early):
    s=Path(source_dir)
    receipt=json.loads((s/"source_receipt.json").read_text())
    if receipt.get("status")!="COMPLETE_ALL_194_SOURCE_PORTS_INDEPENDENT_GBIF_LIVE_ACQUISITION":
        raise ValueError("Acquisition incomplete")
    coverage=read(s/"GBIF_original194_port_2001_2010_coverage.csv")
    if len(coverage)!=194:
        raise ValueError("Original 194 live GBIF port coverage rows incomplete")
    live={}
    for p in coverage:
        pid=int(p["port_id"])
        if pid in live or pid not in pids or p["island_id"]!=pids[pid]:
            raise ValueError("Original194 live port/landmass source identity changed")
        if int(p["original_early_port_buffer_records"])!=original_early[pid]:
            raise ValueError("Archived GBIF count source altered")
        live[pid]={"indexed":int(p["live_GBIF_2001_2010_plant_records"]),
                   "n_facets":int(p["live_faceted_species_keys"])}
        if live[pid]["indexed"]<0 or live[pid]["n_facets"]<0:
            raise ValueError("Negative GBIF source count")
    flora={sid:set() for sid in ids}
    records=Counter()
    facet_count=Counter()
    keys=set()
    n=0
    used=set()
    for row in read_gz(s/"GBIF_original194_port_2001_2010_speciesKey.csv.gz"):
        pid=int(row["port_id"]);island=row["island_id"];key=int(row["GBIF_speciesKey"])
        val=int(row["indexed_records_2001_2010"])
        if pid not in live or island!=pids[pid] or key<=0 or val<=0 or (pid,key) in used:
            raise ValueError("GBIF source duplicate/invalid port × indexed plant speciesKey")
        used.add((pid,key))
        flora[island].add(key)
        keys.add(key)
        records[pid]+=val;facet_count[pid]+=1;n+=1
    if n!=int(receipt["total_port_speciesKey_rows"]):
        raise ValueError("Live archived source speciesKey rows incomplete")
    for pid in pids:
        if records[pid]>live[pid]["indexed"] or facet_count[pid]!=live[pid]["n_facets"]:
            raise ValueError("Live GBIF port facet counts disagree with recorded index")
    coverage_island=read(s/"GBIF_original51_early_general_flora_islands.csv")
    if len(coverage_island)!=51 or {r["island_id"] for r in coverage_island}!=ids:
        raise ValueError("Original51 island early flora coverage table missing")
    port_groups=Counter(pids.values())
    for r in coverage_island:
        sid=r["island_id"]
        if (int(r["original_port_count"])!=port_groups[sid] or
            int(r["unique_early_indexed_general_plant_speciesKeys"])!=len(flora[sid])):
            raise ValueError("Original island general plant species set union disagrees")
    if (len(keys)!=int(receipt["distinct_GBIF_speciesKeys_within_source_query"]) or
        sum(not x for x in flora.values())!=int(receipt["islands_with_no_speciesKey_not_true_absence"])):
        raise ValueError("Original plant speciesKey count/recorded-zero coverage differs")
    return flora,live,n,len(keys)

def check_pair_features(output,pairs,flora):
    features=read(output/"source_pair_general_flora_covariates.csv")
    if len(features)!=1275:
        raise ValueError("Missing all 1275 archived pair-specific general flora predictors")
    indexed={}
    for k,row in enumerate(features):
        idx=int(row["source_pair_index"])
        if idx!=k or row["island_a"]!=pairs[k]["island_a"] or row["island_b"]!=pairs[k]["island_b"]:
            raise ValueError("General-flora feature source-pair order incorrect")
        a,b=flora[row["island_a"]],flora[row["island_b"]]
        aa,bb=len(a),len(b)
        denominator=len(a|b)
        jaccard=len(a&b)/denominator if denominator else 0.
        expected={
            "early_general_plant_specieskey_jaccard_imputed":jaccard,
            "sum_log1p_2001_2010_general_plant_specieskey_richness":math.log1p(aa)+math.log1p(bb),
            "abs_difference_log1p_2001_2010_general_plant_specieskey_richness":abs(math.log1p(aa)-math.log1p(bb)),
            "only_one_endpoint_no_2001_2010_plant_specieskey":int((aa==0) != (bb==0)),
            "both_endpoints_no_2001_2010_plant_specieskey":int(aa==0 and bb==0),
        }
        for field,value in expected.items():
            equal(float(row[field]),value,"source early flora "+field)
        recorded=(aa>0 and bb>0)
        if (row["both_early_species_indexed"].lower()=="true")!=recorded:
            raise ValueError("Unrecorded historical flora island was treated as botanical absence")
        indexed[k]=recorded
    return indexed

def audit_predictions(output,pairs,fids,category,folds,indexed):
    errors={}
    losses=defaultdict(list)
    n=0
    for row in read_gz(output/"heldout_predictions.csv.gz"):
        t,m=row["target"],row["model"];idx=int(row["source_pair_index"])
        fid=row["fold_id"]
        if (t not in TARGETS or m not in MODELS or idx not in range(1275) or
            row["pair_class"]!=category[idx] or fid!=fids[idx] or
            row["island_a"]!=pairs[idx]["island_a"] or row["island_b"]!=pairs[idx]["island_b"]):
            raise ValueError("Pair prediction not in the original double-endpoint-unseen heldout fold")
        key=(t,m,idx)
        if key in errors:
            raise ValueError("Original pair represented more than once in its target/model")
        actual=float(row["observed"]);prediction=float(row["predicted"])
        if not(0<=actual<=1 and 0<=prediction<=1):
            raise ValueError("Source butterfly observed/predicted community Jaccard outside range")
        equal(actual,float(pairs[idx][t]),"original source butterfly target")
        sq=(actual-prediction)**2;ae=abs(actual-prediction)
        equal(sq,float(row["squared_error"]),"saved original pair squared error")
        equal(ae,float(row["absolute_error"]),"saved original pair absolute error")
        if (row["early_general_flora_indexed_both_endpoints"].lower()=="true")!=indexed[idx]:
            raise ValueError("Source pair recorded-early-plants flag differs")
        errors[key]=(sq,ae)
        losses[(t,m,fid)].append((sq,ae))
        n+=1
    if n!=15300 or len(errors)!=15300 or len(losses)!=180:
        raise ValueError("Expected 15,300 predictions and 180 target/model/fold combinations")
    if any(len(v)!=folds[key[2]]["n_test"] for key,v in losses.items()):
        raise ValueError("Incomplete original double-unseen fold source test pair count")
    model_folds=read(output/"model_fold_losses.csv")
    checked=set()
    for row in model_folds:
        key=(row["target"],row["model"],row["fold_id"])
        if key in checked or key not in losses:
            raise ValueError("Unknown/duplicate model fold record")
        checked.add(key)
        values=losses[key]
        if (int(row["train_pairs"])!=folds[key[2]]["n_train"] or
            int(row["test_pairs"])!=folds[key[2]]["n_test"] or
            row["pair_class"]!=folds[key[2]]["class"]):
            raise ValueError("Changed train/test island exclusion or geography class")
        equal(math.fsum(x[0] for x in values)/len(values),float(row["MSE"]),"model fold MSE")
        equal(math.fsum(x[1] for x in values)/len(values),float(row["MAE"]),"model fold MAE")
    if len(checked)!=180:
        raise ValueError("Expected all 180 model fold errors")
    return errors,losses,n

def audit_comparisons(output,errors,losses,category,indexed,folds):
    saved=read(output/"comparison_fold_losses.csv")
    seen=set()
    for row in saved:
        t,name,fid=row["target"],row["contrast"],row["fold_id"]
        key=(t,name,fid)
        if t not in TARGETS or name not in COMPARE or key in seen or fid not in folds:
            raise ValueError("Unknown/duplicate source comparison fold")
        seen.add(key)
        ref,aug=COMPARE[name]
        x=losses[(t,ref,fid)];y=losses[(t,aug,fid)]
        rx=math.fsum(r[0] for r in x)/len(x)
        ay=math.fsum(r[0] for r in y)/len(y)
        rma=math.fsum(r[1] for r in x)/len(x)
        ama=math.fsum(r[1] for r in y)/len(y)
        for field,value in (
            ("MSE_reference",rx),("MSE_augmented",ay),("MSE_gain",rx-ay),
            ("MAE_reference",rma),("MAE_augmented",ama),("MAE_gain",rma-ama)
        ):
            equal(float(row[field]),value,"comparison fold "+field)
        if (int(row["test_pairs"])!=folds[fid]["n_test"] or
            int(row["train_pairs"])!=folds[fid]["n_train"] or
            row["pair_class"]!=folds[fid]["class"]):
            raise ValueError("Source comparison training fold or test scope changed")
    if len(saved)!=180 or len(seen)!=180:
        raise ValueError("Expected all 180 source fold comparisons")
    results=read(output/"comparison_scope_summary.csv")
    audited=set()
    for row in results:
        t,c,scope=row["target"],row["contrast"],row["pair_class"]
        key=(t,c,scope)
        if t not in TARGETS or c not in COMPARE or scope not in CATEGORIES or key in audited:
            raise ValueError("Unknown/duplicate source comparison scope")
        audited.add(key)
        ref,aug=COMPARE[c]
        indices=[i for i in range(1275) if
                scope=="all_pairs" or
                scope=="both_early_species_indexed" and indexed[i] or
                scope in ("within_longitude_block","between_longitude_blocks") and category[i]==scope]
        if not indices or len(indices)!=int(row["unique_evaluated_pairs"]):
            raise ValueError("Source target/contrast/scope population changed")
        x=[errors[(t,ref,i)] for i in indices]
        y=[errors[(t,aug,i)] for i in indices]
        rx=math.fsum(v[0] for v in x)/len(x);ay=math.fsum(v[0] for v in y)/len(y)
        rma=math.fsum(v[1] for v in x)/len(x);ama=math.fsum(v[1] for v in y)/len(y)
        for field,value in (
            ("MSE_reference",rx),("MSE_augmented",ay),
            ("MSE_reduction_percent",100*(rx-ay)/rx),
            ("MAE_reference",rma),("MAE_augmented",ama),
            ("MAE_reduction_percent",100*(rma-ama)/rma),
            ("pairs_with_lower_squared_error",sum(a[0]>b[0] for a,b in zip(x,y))),
            ("pairs_with_higher_squared_error",sum(a[0]<b[0] for a,b in zip(x,y))),
            ("pairs_with_lower_absolute_error",sum(a[1]>b[1] for a,b in zip(x,y))),
        ):
            equal(float(row[field]),value,"source pair pooled/scope "+field)
    if len(results)!=48 or len(audited)!=48:
        raise ValueError("Expected all 48 source contrast/scope outcomes")
    parent=read(PARENT)
    original=[r for r in parent if r["contrast"]=="HRE_to_HRE_AIS_LAT"]
    if len(original)!=6:
        raise ValueError("Expected six immutable source parent AIS comparisons")
    maxdiff=0.
    for p in original:
        observed=next((q for q in results if q["target"]==p["target"] and
                     q["contrast"]=="HRE_to_HRE_AIS_LAT" and q["pair_class"]==p["pair_class"]),None)
        if not observed or observed["unique_evaluated_pairs"]!=p["unique_evaluated_pairs"]:
            raise ValueError("Original AIS source comparison population changed")
        for fld in ("MSE_reference","MSE_augmented","MAE_reference","MAE_augmented"):
            diff=abs(float(observed[fld])-float(p[fld]))
            maxdiff=max(maxdiff,diff)
            if diff>1e-9:
                raise ValueError("Original AIS baseline not reproduced")
    return maxdiff

def audit(source_dir,results_dir):
    source_dir=Path(source_dir);results_dir=Path(results_dir)
    ids,pids,archive,pairs,fids,category,folds=original()
    flora,ports,source_rows,taxa=source_flora(source_dir,ids,pids,archive)
    indexed=check_pair_features(results_dir,pairs,flora)
    errs,losses,n=audit_predictions(results_dir,pairs,fids,category,folds,indexed)
    delta=audit_comparisons(results_dir,errs,losses,category,indexed,folds)
    outcome={
        "schema":"chocho_original194_GBIF_early_general_flora_independent_source_and_prediction_audit_v1",
        "status":"COMPLETE_194_PORT_EARLY_FLORA_SOURCE_AND_51_ISLAND_PREDICTIONS_AUDITED",
        "source_ports_verified":len(pids),"source_islands_verified":len(ids),
        "original_butterfly_pair_target_observations":len(pairs)*len(TARGETS),
        "source_port_speciesKey_rows_reconstructed":source_rows,
        "distinct_general_Plantae_speciesKeys_reconstructed":taxa,
        "original_pair_general_flora_features_verified":len(indexed),
        "source_pair_predictions_verified":n,
        "source_model_folds_verified":len(losses),
        "source_contrast_fold_results_verified":180,
        "source_comparison_scope_rows_verified":48,
        "double_unseen_training_island_exclusions_verified":len(folds),
        "original_true_AIS_parent_max_absolute_difference":delta,
        "all194_original_ports_acquired_before_analysis":True,
        "historic_taxa_not_proven_larval_host_or_introduced":True,
        "archived_indexed_zero_not_absence":True,
    }
    (results_dir/"independent_audit.json").write_text(json.dumps(outcome,indent=2,sort_keys=True)+"\n")
    print("EARLY_FLORA_INDEPENDENT_AUDIT_JSON="+json.dumps(outcome,sort_keys=True),flush=True)
    return outcome

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--results-dir",type=Path,required=True)
    x=p.parse_args()
    audit(x.source_dir,x.results_dir)
