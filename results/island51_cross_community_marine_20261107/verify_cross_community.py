#!/usr/bin/env python3
"""Independent stdlib original194-port/990-pair three-community prediction auditor.

Reconstructs early indexed Plantae taxa and all source pair Jaccards from archived
GBIF source rows, original fixed longitude groups, held-out observations, and
every prediction/fold/scope error. No pandas, numpy, or fitting-module imports.
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
SOURCES=ROOT/"results/island51_pre2011_public_general_flora_20261107/source"
PAIR_SOURCE=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv"
ISLAND_SOURCE=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv"
GEO=ROOT/"results/island51_host_circuit_atlas_20261106/island_geography.csv"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
PANELS={
    "early_general_flora":"early_general_flora_jaccard",
    "earlier_butterflies":"earlier_butterfly_jaccard",
    "later_butterflies":"new_butterfly_jaccard",
    "later_butterflies_retro_hosts":"new_butterfly_jaccard",
}
MODELS=("BASE","AIS","WIND","BOTH")
CONTRASTS={
    "BASE_to_AIS":("BASE","AIS"),
    "BASE_to_WIND":("BASE","WIND"),
    "AIS_to_BOTH":("AIS","BOTH"),
    "WIND_to_BOTH":("WIND","BOTH"),
    "BASE_to_BOTH":("BASE","BOTH"),
}
SCOPES=("all_both_plants_indexed","within_longitude_block","between_longitude_blocks","both_plants_at_least_10_specieskeys")

def read_csv(path):
    with Path(path).open("rt",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def iter_gz(path):
    with gzip.open(path,"rt",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f)

def approx(actual,expected,context):
    a=float(actual);b=float(expected)
    if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10):
        raise ValueError(f"{context}: {a} is not equal to {b}")

def source_data(source_dir):
    original=read_csv(PAIR_SOURCE)
    eligible_islands=[r for r in read_csv(ISLAND_SOURCE) if r["eligible"].lower()=="true"]
    ids={r["landmass_id"] for r in eligible_islands}
    if len(ids)!=51 or len(original)!=1275 or len({(r["island_a"],r["island_b"]) for r in original})!=1275:
        raise ValueError("Original butterfly 51-island source pair cohort changed")
    for r in original:
        if r["island_a"] not in ids or r["island_b"] not in ids or r["island_a"]==r["island_b"]:
            raise ValueError("Original source pair endpoints missing")
        for target in ("earlier_butterfly_jaccard","new_butterfly_jaccard"):
            if not 0<=float(r[target])<=1:
                raise ValueError("Original indexed butterfly Jaccard invalid")
    ports=read_csv(WPI)
    port_ids={int(r["port_id"]):r["island_id"] for r in ports}
    if (len(ports)!=194 or len(port_ids)!=194 or
        sum(r["included"].lower()=="true" for r in ports)!=182 or
        {r["island_id"] for r in ports}!=ids):
        raise ValueError("Original 194 source WPI ports not retained")
    receipt=json.loads((Path(source_dir)/"source_receipt.json").read_text())
    if (receipt.get("status")!="COMPLETE_ALL_194_SOURCE_PORTS_INDEPENDENT_GBIF_LIVE_ACQUISITION" or
        receipt.get("source_port_original_geo_count")!=194 or
        receipt.get("distinct_GBIF_speciesKeys_within_source_query")!=4876):
        raise ValueError("Archived direct GBIF original194 plant source incomplete")
    coverage=read_csv(Path(source_dir)/"GBIF_original194_port_2001_2010_coverage.csv")
    if len(coverage)!=194 or len({int(r["port_id"]) for r in coverage})!=194:
        raise ValueError("Original 194 GBIF botanical coverage records missing")
    plant_sets={island:set() for island in ids}
    facets=set()
    species_keys=set()
    count_per_port=Counter()
    for row in iter_gz(Path(source_dir)/"GBIF_original194_port_2001_2010_speciesKey.csv.gz"):
        pid=int(row["port_id"]);island=row["island_id"];key=int(row["GBIF_speciesKey"])
        weight=int(row["indexed_records_2001_2010"])
        if (pid not in port_ids or port_ids[pid]!=island or
            key<=0 or weight<=0 or (pid,key) in facets):
            raise ValueError("Source botanical speciesKey×port invalid/duplicated")
        facets.add((pid,key))
        species_keys.add(key)
        count_per_port[pid]+=weight
        plant_sets[island].add(key)
    if len(facets)!=27598 or len(species_keys)!=4876:
        raise ValueError("Archived original 194 source speciesKey census not complete")
    for row in coverage:
        pid=int(row["port_id"])
        if row["island_id"]!=port_ids[pid] or count_per_port[pid]>int(row["live_GBIF_2001_2010_plant_records"]):
            raise ValueError("Archived plant species faceting exceeds original port counts")
    count={k:len(v) for k,v in plant_sets.items()}
    if sum(n>0 for n in count.values())!=45 or sum(n>=10 for n in count.values())!=39:
        raise ValueError("Historical source-only plant-index coverage changed")
    geo=read_csv(GEO)
    if len(geo)!=51 or {x["island_id"] for x in geo}!=ids:
        raise ValueError("Source geographic original 51 islands not covered")
    ordered=sorted(geo,key=lambda r:(float(r["representative_longitude"]),int(r["island_id"])))
    longitude={}
    for idx,row in enumerate(ordered):
        longitude[row["island_id"]]=0 if idx<11 else 1+(idx-11)//10
    if Counter(longitude.values())!=Counter({0:11,1:10,2:10,3:10,4:10}):
        raise ValueError("Frozen five longitude blocks changed")
    indexed={}
    high={}
    pair_folds={}
    pair_class={}
    plant_jaccard={}
    scopes={k:set() for k in SCOPES}
    for idx,row in enumerate(original):
        a,b=row["island_a"],row["island_b"]
        ab=len(plant_sets[a]&plant_sets[b])
        total=len(plant_sets[a]|plant_sets[b])
        plant_jaccard[idx]=ab/total if total else 0.0
        indexed[idx]=count[a]>0 and count[b]>0
        high[idx]=count[a]>=10 and count[b]>=10
        g,h=sorted([longitude[a],longitude[b]])
        fold=f"B{g}_B{h}"
        pair_folds[idx]=fold
        kind="within_longitude_block" if g==h else "between_longitude_blocks"
        pair_class[idx]=kind
        if indexed[idx]:
            scopes["all_both_plants_indexed"].add(idx)
            scopes[kind].add(idx)
            if high[idx]:
                scopes["both_plants_at_least_10_specieskeys"].add(idx)
    if (len(scopes["all_both_plants_indexed"])!=990 or
        len(scopes["both_plants_at_least_10_specieskeys"])!=741 or
        len(scopes["within_longitude_block"]|scopes["between_longitude_blocks"])!=990):
        raise ValueError("Exactly 990 original matched taxa pairs and 741 high-coverage pairs required")
    original_folds={}
    for k in range(5):
        for h in range(k,5):
            fid=f"B{k}_B{h}"
            held={sid for sid,g in longitude.items() if g in (k,h)}
            test={i for i in scopes["all_both_plants_indexed"] if pair_folds[i]==fid}
            train={i for i in scopes["all_both_plants_indexed"]
                   if original[i]["island_a"] not in held and original[i]["island_b"] not in held}
            if not test or len(train)<20 or any(original[i]["island_a"] in held or original[i]["island_b"] in held for i in train):
                raise ValueError("Source-matched double-unseen training cohort invalid")
            original_folds[fid]={"test":test,"train":train,
                 "class":"within_longitude_block" if k==h else "between_longitude_blocks",
                 "held":held}
    if len(original_folds)!=15 or sum(len(v["test"]) for v in original_folds.values())!=990:
        raise ValueError("Original double-unseen fold pair partition incomplete")
    return original,plant_sets,plant_jaccard,indexed,high,pair_folds,pair_class,scopes,original_folds

def check_ledger(out,original,plant_jaccard,indexed,high,pair_class):
    rows=read_csv(out/"source_pair_community_panel.csv")
    if len(rows)!=1275:
        raise ValueError("Original full 1275 pair source-response ledger not published")
    for idx,r in enumerate(rows):
        if (int(r["source_pair_index"])!=idx or
            r["island_a"]!=original[idx]["island_a"] or
            r["island_b"]!=original[idx]["island_b"]):
            raise ValueError("Plant/community full source pair index mismatched")
        approx(r["early_general_flora_jaccard_imputed"],plant_jaccard[idx],"old botanical Jaccard")
        approx(r["earlier_butterfly_jaccard"],original[idx]["earlier_butterfly_jaccard"],"original early butterfly target")
        approx(r["later_butterfly_jaccard"],original[idx]["new_butterfly_jaccard"],"original later butterfly target")
        approx(r["source_port_distance_km"],original[idx]["port_distance_km"],"original port-source distance")
        if (r["early_general_flora_both_islands_indexed"].lower()=="true")!=indexed[idx]:
            raise ValueError("Matched historical indexed-botanical cohort changed")
        if (r["early_general_flora_both_islands_at_least_10_specieskeys"].lower()=="true")!=high[idx]:
            raise ValueError("Matched high-specieskey source pair cohort changed")
        if r["longitude_pair_class"]!=pair_class[idx]:
            raise ValueError("Original source pair geography scope changed")

def check_folds(out,original_folds):
    rows=read_csv(out/"source_eligible_geo_folds.csv")
    if len(rows)!=15:
        raise ValueError("Expected exactly fifteen original source geography blocks")
    seen=set()
    for r in rows:
        fid=r["fold_id"]
        if fid in seen or fid not in original_folds:
            raise ValueError("Unknown or repeated source fold")
        seen.add(fid)
        actual=original_folds[fid]
        if (r["pair_class"]!=actual["class"] or
            int(r["train_source_pairs"])!=len(actual["train"]) or
            int(r["test_source_pairs"])!=len(actual["test"]) or
            set(json.loads(r["source_heldout_island_ids"]))!=actual["held"]):
            raise ValueError("Saved original heldout training/sample source support disagrees")
    if len(seen)!=15:
        raise ValueError("Incomplete 15 source folds")

def prediction_audit(out,original,jaccard,indexed,fids,kind,folds):
    errors={}
    by_fold=defaultdict(list)
    n=0
    for row in iter_gz(out/"heldout_predictions.csv.gz"):
        panel=row["panel"];model=row["model"]
        fid=row["fold_id"];i=int(row["source_pair_index"])
        if (panel not in PANELS or model not in MODELS or i not in indexed or not indexed[i] or
            fid!=fids[i] or fid not in folds or row["pair_class"]!=kind[i] or
            i not in folds[fid]["test"] or
            row["island_a"]!=original[i]["island_a"] or row["island_b"]!=original[i]["island_b"]):
            raise ValueError("Source pair forecast outside original plant-matched double-unseen sample")
        observed=float(row["observed"])
        true=jaccard[i] if panel=="early_general_flora" else float(original[i][PANELS[panel]])
        approx(observed,true,"original source panel response")
        predicted=float(row["predicted"])
        if not 0<=predicted<=1:
            raise ValueError("Forecast outside valid response Jaccard")
        sq=(true-predicted)**2
        ae=abs(true-predicted)
        approx(row["squared_error"],sq,"individual squared error")
        approx(row["absolute_error"],ae,"individual absolute error")
        key=(panel,model,i)
        if key in errors:
            raise ValueError("Original source pair repeated for model")
        errors[key]=(sq,ae)
        by_fold[(panel,model,fid)].append((sq,ae))
        n+=1
    if n!=15840 or len(errors)!=15840 or len(by_fold)!=240:
        raise ValueError("Expected 15840 source-pair forecasts and 240 model folds")
    for panel in PANELS:
        for model in MODELS:
            if any((panel,model,i) not in errors for i in indexed if indexed[i]):
                raise ValueError("Original eligible 990 source pairs not tested once")
    for key,v in by_fold.items():
        if len(v)!=len(folds[key[2]]["test"]):
            raise ValueError("Original source-fold test size changed")
    tab=read_csv(out/"model_fold_losses.csv")
    seen=set()
    for r in tab:
        key=(r["panel"],r["model"],r["fold_id"])
        if key in seen or key not in by_fold:
            raise ValueError("Unknown/duplicated model fold")
        seen.add(key)
        values=by_fold[key];support=folds[key[2]]
        if (int(r["test_pairs"])!=len(support["test"]) or
            int(r["train_pairs"])!=len(support["train"]) or
            r["pair_class"]!=support["class"]):
            raise ValueError("Source matched training/test fold support altered")
        approx(r["MSE"],math.fsum(v[0] for v in values)/len(values),"fold MSE")
        approx(r["MAE"],math.fsum(v[1] for v in values)/len(values),"fold MAE")
    if len(tab)!=240 or len(seen)!=240:
        raise ValueError("Missing all 240 source-pair model folds")
    return errors,by_fold,n

def comparisons_audit(out,errors,by_fold,scopes,folds):
    rows=read_csv(out/"comparison_fold_losses.csv")
    seen=set()
    for r in rows:
        panel,contrast,fid=r["panel"],r["contrast"],r["fold_id"]
        key=(panel,contrast,fid)
        if key in seen or panel not in PANELS or contrast not in CONTRASTS or fid not in folds:
            raise ValueError("Source comparison fold invalid/duplicate")
        seen.add(key)
        ref,aug=CONTRASTS[contrast]
        x=by_fold[(panel,ref,fid)];y=by_fold[(panel,aug,fid)]
        a=math.fsum(z[0] for z in x)/len(x);b=math.fsum(z[0] for z in y)/len(y)
        c=math.fsum(z[1] for z in x)/len(x);d=math.fsum(z[1] for z in y)/len(y)
        for name,value in (("MSE_reference",a),("MSE_augmented",b),
                           ("MSE_gain",a-b),("MAE_reference",c),
                           ("MAE_augmented",d),("MAE_gain",c-d)):
            approx(r[name],value,"source comparison fold "+name)
        if int(r["train_pairs"])!=len(folds[fid]["train"]) or int(r["test_pairs"])!=len(folds[fid]["test"]):
            raise ValueError("Source comparison fold population changed")
    if len(rows)!=300 or len(seen)!=300:
        raise ValueError("Expected 300 source panel/model contrast folds")
    summary=read_csv(out/"comparison_scope_summary.csv")
    done=set()
    for row in summary:
        panel,contrast,scope=row["panel"],row["contrast"],row["pair_class"]
        key=(panel,contrast,scope)
        if (key in done or panel not in PANELS or contrast not in CONTRASTS or
            scope not in SCOPES):
            raise ValueError("Duplicate/missing source-community summary identity")
        done.add(key)
        ref,aug=CONTRASTS[contrast]
        indices=sorted(scopes[scope])
        if int(row["unique_source_pairs"])!=len(indices):
            raise ValueError("Original indexed botanically matched pair count changed")
        x=[errors[(panel,ref,i)] for i in indices]
        y=[errors[(panel,aug,i)] for i in indices]
        r_mse=math.fsum(v[0] for v in x)/len(x)
        a_mse=math.fsum(v[0] for v in y)/len(y)
        r_mae=math.fsum(v[1] for v in x)/len(x)
        a_mae=math.fsum(v[1] for v in y)/len(y)
        for name,value in (
           ("MSE_reference",r_mse),("MSE_augmented",a_mse),
           ("MSE_gain_percent",100*(r_mse-a_mse)/r_mse),
           ("MAE_reference",r_mae),("MAE_augmented",a_mae),
           ("MAE_gain_percent",100*(r_mae-a_mae)/r_mae),
           ("pairs_with_lower_squared_error",sum(z[0]>w[0] for z,w in zip(x,y))),
           ("pairs_with_higher_squared_error",sum(z[0]<w[0] for z,w in zip(x,y))),
           ("pairs_with_lower_absolute_error",sum(z[1]>w[1] for z,w in zip(x,y))),
        ):
            approx(row[name],value,"independent community scope "+name)
    if len(summary)!=80 or len(done)!=80:
        raise ValueError("All 80 source-matched community/source/geographic comparisons required")
    return len(summary)

def audit(out):
    directory=Path(out)
    original,sets,jaccard,indexed,high,fids,kind,scopes,folds=source_data(SOURCES)
    check_ledger(directory,original,jaccard,indexed,high,kind)
    check_folds(directory,folds)
    errors,by_fold,predictions=prediction_audit(directory,original,jaccard,indexed,fids,kind,folds)
    comparisons=comparisons_audit(directory,errors,by_fold,scopes,folds)
    result={
      "schema":"chocho_990_pairs_early_plants_early_late_butterfly_cross_community_independent_audit_v1",
      "status":"ALL_990_SOURCE_PAIRS_EVERY_PANEL_MARINE_MODEL_AND_CONTRAST_AUDITED",
      "original_source_islands":51,
      "historical_general_Plant_specieskey_islands":45,
      "source_port_count":194,
      "source_plant_specieskey_rows":27598,
      "matched_original_source_pairs":990,
      "high_plant_indexed_pair_count":741,
      "original_source_pair_rows_reconstructed":len(jaccard),
      "geography_both_endpoints_unseen_folds":len(folds),
      "individual_source_pair_predictions_verified":predictions,
      "model_folds_verified":len(by_fold),
      "contrast_folds_verified":300,
      "all_scoped_comparisons_verified":comparisons,
      "original_source_plant_index_nonrecords_are_not_absence":True,
      "AIS_source_postdates_early_plants_and_early_butterflies":True,
      "no_independent_pair_level_p_values_or_actual_shipment_claims":True,
    }
    (directory/"independent_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("CROSS_COMMUNITY_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("--results-dir",type=Path,required=True)
    audit(arg.parse_args().results_dir)
