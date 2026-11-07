#!/usr/bin/env python3
"""Independent no-pandas auditor of source-driven top-port exclusions and CV predictions."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
PAIRS=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv"
ISLANDS=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv"
GEOGRAPHY=ROOT/"results/island51_host_circuit_atlas_20261106/island_geography.csv"
MODELS=("HRE_LAT","HRE_AIS_LAT","HRE_WIND_LAT","HCRE_LAT")
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
COMPARE={
    "HRE_to_HRE_AIS_LAT":("HRE_LAT","HRE_AIS_LAT"),
    "HRE_to_HRE_WIND_LAT":("HRE_LAT","HRE_WIND_LAT"),
    "HRE_AIS_to_HCRE_LAT":("HRE_AIS_LAT","HCRE_LAT"),
    "HRE_WIND_to_HCRE_LAT":("HRE_WIND_LAT","HCRE_LAT"),
    "HRE_to_HCRE_LAT":("HRE_LAT","HCRE_LAT"),
}
SCENARIOS={"FULL":(),"DROP_15":("15",),
           "DROP_9":("9",),"DROP_15_9":("15","9")}
SCOPES={"all_pairs":1275,"within_longitude_block":235,"between_longitude_blocks":1040}

def csv_rows(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        items=list(csv.DictReader(f))
    if not items:
        raise ValueError("Empty audit input: "+str(path))
    return items

def same(a,b,detail):
    if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(
        a,b,rel_tol=1e-10,abs_tol=1e-11):
        raise ValueError("Audit mismatch "+detail+": "+str(a)+" != "+str(b))

def original_source():
    src=csv_rows(PAIRS)
    eligible=[r for r in csv_rows(ISLANDS) if r["eligible"].lower()=="true"]
    ids={r["landmass_id"] for r in eligible}
    counts={r["landmass_id"]:int(r["ports"]) for r in eligible}
    if len(src)!=1275 or len(ids)!=51 or sum(counts.values())!=194:
        raise ValueError("The original island and port roster changed")
    ranking=sorted(ids,key=lambda k:(-counts[k],int(k)))
    if [(x,counts[x]) for x in ranking[:2]]!=[("15",39),("9",21)]:
        raise ValueError("Original top two port groups are not source-defined")
    src_endpoints=set()
    for row in src:
        a,b=row["island_a"],row["island_b"]
        if a not in ids or b not in ids or a==b:
            raise ValueError("Unknown pair/source island ID")
        key=tuple(sorted([a,b],key=int))
        if key in src_endpoints:
            raise ValueError("Duplicate source unordered pair")
        src_endpoints.add(key)
    if len(src_endpoints)!=1275:
        raise ValueError("Incomplete source unordered pair support")
    geo=csv_rows(GEOGRAPHY)
    if len(geo)!=51 or len({x["island_id"] for x in geo})!=51 or {
        x["island_id"] for x in geo}!=ids:
        raise ValueError("Incomplete source original geography")
    for row in geo:
        if row["island_id"]!=row["landmass_id"] or int(row["original_ports"])!=counts[row["island_id"]]:
            raise ValueError("WPI original port counts differ from original source roster")
    ordered=sorted(geo,key=lambda r:(float(r["representative_longitude"]),int(r["island_id"])))
    blocks={}
    for i,row in enumerate(ordered):
        blocks[row["island_id"]]=(0 if i<11 else 1+(i-11)//10)
    if Counter(blocks.values())!=Counter({0:11,1:10,2:10,3:10,4:10}):
        raise ValueError("Original longitude blocks are not source-defined")
    return src,counts,ranking,ordered,blocks

def verify_manifest(out,src,counts,ranking,ordered,blocks):
    ranking_saved=csv_rows(out/"original_port_roster_rank.csv")
    if len(ranking_saved)!=51:
        raise ValueError("Original port rank has missing island")
    for n,(actual,key) in enumerate(zip(ranking_saved,ranking),1):
        if (int(actual["original_port_rank"])!=n or
            actual["island_id"]!=key or
            int(actual["original_port_count"])!=counts[key] or
            (actual["excluded_by_fixed_port_training_scenarios"].lower()=="true")!=(key in ("15","9"))):
            raise ValueError("Wrong original port rank for "+str(n))
    geo=csv_rows(out/"source_longitude_island_blocks.csv")
    if len(geo)!=51:
        raise ValueError("Original latitude/longitude roster incomplete")
    for a,b in zip(geo,ordered):
        if (a["island_id"]!=b["island_id"] or
            int(a["longitude_block"])!=blocks[b["island_id"]]):
            raise ValueError("Output longitude block assignment differs")
        for key in ("representative_longitude","representative_latitude"):
            same(float(a[key]),float(b[key]),"original source "+key)
        if a["group_example_port"]!=b["example_port_name"]:
            raise ValueError("Original example-port identity changed")

    p=csv_rows(out/"scenario_fold_specification.csv")
    if len(p)!=60:
        raise ValueError("Must have 4 scenarios × 15 original block-pair folds")
    manifest={}
    expected_tests={}
    for i,source in enumerate(src):
        x,y=blocks[source["island_a"]],blocks[source["island_b"]]
        fid=f"B{min(x,y)}_B{max(x,y)}"
        expected_tests.setdefault(fid,[]).append(i)
    if len(expected_tests)!=15 or sum(map(len,expected_tests.values()))!=1275:
        raise ValueError("Original source block/pair mapping incomplete")
    ids=set(blocks)
    for row in p:
        scenario=row["scenario"];fid=row["fold_id"]
        if scenario not in SCENARIOS or fid not in expected_tests or (scenario,fid) in manifest:
            raise ValueError("Unknown or repeated source/scenario fold")
        k,h=map(int,[row["block_a"],row["block_b"]])
        if fid!=f"B{k}_B{h}" or not(0<=k<=h<=4):
            raise ValueError("Unordered block pair identity mismatch")
        held={sid for sid in ids if blocks[sid] in (k,h)}
        drop=set(SCENARIOS[scenario])
        excluded=held|drop
        if (set(json.loads(row["source_heldout_island_ids"]))!=held or
            tuple(json.loads(row["additional_training_excluded_ids"]))!=SCENARIOS[scenario] or
            set(json.loads(row["total_training_excluded_island_ids"]))!=excluded):
            raise ValueError("Training exclusions not the frozen heldouts + source-ranked port groups")
        test=expected_tests[fid]
        train=[i for i,s in enumerate(src) if s["island_a"] not in excluded and s["island_b"] not in excluded]
        original_train=[i for i,s in enumerate(src) if s["island_a"] not in held and s["island_b"] not in held]
        klass="within_longitude_block" if k==h else "between_longitude_blocks"
        if (row["pair_class"]!=klass or
            len(test)!=int(row["test_pairs"]) or
            len(train)!=int(row["train_pairs"]) or
            len(original_train)!=int(row["original_fold_train_pairs"]) or
            len(ids)-len(excluded)!=int(row["training_islands"]) or
            len(train)!=math.comb(len(ids)-len(excluded),2) or
            row["pair_class"] not in SCOPES):
            raise ValueError("Original pair test or scenario training size mismatch")
        manifest[(scenario,fid)]={"class":klass,"test_indices":set(test),
            "train_size":len(train),"test_size":len(test),"excluded":excluded}
    if len(manifest)!=60:
        raise ValueError("Incomplete source-driven training-exclusion manifests")
    return manifest,expected_tests

def audit(out_dir):
    out=Path(out_dir)
    src,counts,ranking,ordered,blocks=original_source()
    manifest,fold_tests=verify_manifest(out,src,counts,ranking,ordered,blocks)
    losses=defaultdict(list)
    by_pair={}
    count=0
    with (out/"heldout_pair_predictions.csv").open(encoding="utf-8",newline="") as f:
        for row in csv.DictReader(f):
            scenario,target,model=row["scenario"],row["target"],row["model"]
            fid=row["fold_id"]
            idx=int(row["source_pair_index"])
            if (scenario not in SCENARIOS or target not in TARGETS or model not in MODELS or
                idx<0 or idx>=1275 or (scenario,fid) not in manifest):
                raise ValueError("Saved prediction model, scenario, fold or pair outside original fixed population")
            settings=manifest[(scenario,fid)]
            original=src[idx]
            if (idx not in settings["test_indices"] or
                row["pair_class"]!=settings["class"] or
                row["island_a"]!=original["island_a"] or
                row["island_b"]!=original["island_b"]):
                raise ValueError("Saved source pair or exact both-unseen test fold wrong")
            actual,pred=float(row["observed"]),float(row["predicted"])
            if not (0<=actual<=1 and 0<=pred<=1):
                raise ValueError("Invalid butterfly similarity observation/prediction")
            same(actual,float(original[target]),"source butterfly observation")
            sq=(actual-pred)**2
            ae=abs(actual-pred)
            same(sq,float(row["squared_error"]),"saved squared error")
            same(ae,float(row["absolute_error"]),"saved absolute error")
            key=(scenario,target,model,idx)
            if key in by_pair:
                raise ValueError("Original pair appears more than once in one scenario/target/model")
            by_pair[key]={"sq":sq,"ae":ae,"class":settings["class"],"fold":fid}
            losses[(scenario,target,model,fid)].append((sq,ae))
            count+=1
    if (count!=40800 or len(by_pair)!=40800 or len(losses)!=480 or
        any(len(values)!=manifest[(key[0],key[3])]["test_size"]
            for key,values in losses.items())):
        raise ValueError("All 40,800 original-pair predictions and 480 model folds required")
    for scenario in SCENARIOS:
        for target in TARGETS:
            for model in MODELS:
                if any((scenario,target,model,i) not in by_pair for i in range(1275)):
                    raise ValueError("Incomplete original pair support in a scenario/target/model")

    model_folds=csv_rows(out/"model_fold_losses.csv")
    saved_folds=set()
    for row in model_folds:
        key=(row["scenario"],row["target"],row["model"],row["fold_id"])
        if key in saved_folds or key not in losses:
            raise ValueError("Repeated or unknown model/fold")
        saved_folds.add(key)
        values=losses[key]
        meta=manifest[(key[0],key[3])]
        if (len(values)!=int(row["test_pairs"]) or
            meta["train_size"]!=int(row["train_pairs"]) or
            meta["class"]!=row["pair_class"]):
            raise ValueError("Model fold source training or test counts differ")
        same(math.fsum(z[0] for z in values)/len(values),float(row["MSE"]),"model fold MSE")
        same(math.fsum(z[1] for z in values)/len(values),float(row["MAE"]),"model fold MAE")
    if len(model_folds)!=480 or len(saved_folds)!=480:
        raise ValueError("Required 480 complete model folds")

    comparisons=csv_rows(out/"contrast_fold_losses.csv")
    seen=set()
    for row in comparisons:
        scenario,target,comp,fid=row["scenario"],row["target"],row["comparison"],row["fold_id"]
        key=(scenario,target,comp,fid)
        if (key in seen or scenario not in SCENARIOS or
            target not in TARGETS or comp not in COMPARE or (scenario,fid) not in manifest):
            raise ValueError("Unknown comparison fold")
        seen.add(key)
        ref,aug=COMPARE[comp]
        ra=losses[(scenario,target,ref,fid)]
        aa=losses[(scenario,target,aug,fid)]
        r_mse=math.fsum(z[0] for z in ra)/len(ra);a_mse=math.fsum(z[0] for z in aa)/len(aa)
        r_mae=math.fsum(z[1] for z in ra)/len(ra);a_mae=math.fsum(z[1] for z in aa)/len(aa)
        for field,value in (("MSE_reference",r_mse),("MSE_augmented",a_mse),
                            ("MSE_gain",r_mse-a_mse),("MAE_reference",r_mae),
                            ("MAE_augmented",a_mae),("MAE_gain",r_mae-a_mae)):
            same(value,float(row[field]),"fold contrast "+field)
        meta=manifest[(scenario,fid)]
        if (row["pair_class"]!=meta["class"] or
            int(row["test_pairs"])!=meta["test_size"] or
            int(row["train_pairs"])!=meta["train_size"] or
            (row["MSE_improved"].lower()=="true")!=(r_mse>a_mse)):
            raise ValueError("Saved comparison sign, test or training support differs")
    if len(comparisons)!=600 or len(seen)!=600:
        raise ValueError("Expected exactly 600 independent comparison-fold audits")
    scopes=csv_rows(out/"contrast_scope_summary.csv")
    checked=set()
    for row in scopes:
        scenario,target,comp,cat=row["scenario"],row["target"],row["contrast"],row["pair_class"]
        key=(scenario,target,comp,cat)
        if (key in checked or scenario not in SCENARIOS or target not in TARGETS or
            comp not in COMPARE or cat not in SCOPES):
            raise ValueError("Unexpected or duplicate scenario/scope summary")
        checked.add(key)
        ref,aug=COMPARE[comp]
        ids=[i for i in range(1275)
             if cat=="all_pairs" or by_pair[(scenario,target,ref,i)]["class"]==cat]
        if len(ids)!=SCOPES[cat] or int(row["unique_evaluated_pairs"])!=len(ids):
            raise ValueError("Wrong geographical pair population")
        x=[by_pair[(scenario,target,ref,i)] for i in ids]
        y=[by_pair[(scenario,target,aug,i)] for i in ids]
        r_mse=math.fsum(v["sq"] for v in x)/len(x);a_mse=math.fsum(v["sq"] for v in y)/len(y)
        r_mae=math.fsum(v["ae"] for v in x)/len(x);a_mae=math.fsum(v["ae"] for v in y)/len(y)
        for name,val in (("MSE_reference",r_mse),("MSE_augmented",a_mse),
                         ("MSE_reduction_percent",100*(r_mse-a_mse)/r_mse),
                         ("MAE_reference",r_mae),("MAE_augmented",a_mae),
                         ("MAE_reduction_percent",100*(r_mae-a_mae)/r_mae),
                         ("pairs_with_lower_squared_error",sum(a["sq"]>b["sq"] for a,b in zip(x,y))),
                         ("pairs_with_higher_squared_error",sum(a["sq"]<b["sq"] for a,b in zip(x,y))),
                         ("pairs_with_lower_absolute_error",sum(a["ae"]>b["ae"] for a,b in zip(x,y)))):
            same(float(row[name]),float(val),"full pair-level scenario summary "+name)
    if len(scopes)!=120 or len(checked)!=120:
        raise ValueError("Expected exactly 120 scope summaries")
    receipt={"schema":"chocho_top_port_training_influence_independent_audit_v1",
       "status":"ALL_SOURCE_PAIR_TRAINING_INFLUENCE_AND_PREDICTION_AUDITS_PASS",
       "original_islands":51,"source_ports":194,"ranked_top_port_islands":["15","9"],
       "all_scenarios_checked":list(SCENARIOS),
       "fixed_original_pair_population":1275,
       "fixed_longitude_block_pair_folds":15,
       "scenario_fold_source_rows":len(manifest),
       "source_butterfly_pair_target_observations":2550,
       "saved_predictions_reconstructed":count,"model_folds_verified":len(model_folds),
       "comparison_folds_verified":len(comparisons),"scope_comparisons_verified":len(scopes),
       "no_test_pair_or_excluded_training_island_in_training":True,
       "all_four_scenarios_evaluate_the_same_1275_pairs_once":True,
       "folds_and_pair_observations_are_not_independent":True}
    (out/"independent_audit.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("TOP_PORT_INFLUENCE_AUDIT_JSON="+json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--results-dir",type=Path,required=True)
    audit(ap.parse_args().results_dir)
