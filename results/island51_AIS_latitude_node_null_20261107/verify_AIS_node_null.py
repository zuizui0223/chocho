#!/usr/bin/env python3
"""Independent stdlib source-pair, latitude-label, and 1,025,100 prediction audit.

Reconstructs saved errors, 2400 scores and 12 empirical ranks without importing
the fitting/label-shuffling module or using pandas or numpy. The original
source-locked fitting workflow separately checks 12 parent CV reference scores.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCE_PAIRS=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv"
SOURCE_ROSTER=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv"
SOURCE_GEO=ROOT/"results/island51_host_circuit_atlas_20261106/island_geography.csv"
PARENT=ROOT/"results/island51_top_port_training_influence_20261107/contrast_scope_summary.csv"
TARGETS=("new_butterfly_jaccard","earlier_butterfly_jaccard")
CONTRASTS=("AIS_without_wind","AIS_after_wind")
SCOPES={"all_pairs":1275,"within_longitude_block":235,"between_longitude_blocks":1040}
PARENT_NAMES={"AIS_without_wind":"HRE_to_HRE_AIS_LAT","AIS_after_wind":"HRE_WIND_to_HCRE_LAT"}

def rows(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        data=list(reader)
    if not data:
        raise ValueError("Source or saved file contains no rows: "+str(path))
    return data

def close(a,b,what):
    if not (math.isfinite(a) and math.isfinite(b) and
            math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10)):
        raise ValueError(what+" mismatch: "+str(a)+" versus "+str(b))

def sources():
    original=rows(SOURCE_PAIRS)
    eligible=[r for r in rows(SOURCE_ROSTER) if r["eligible"].lower()=="true"]
    ids={r["landmass_id"] for r in eligible}
    if len(ids)!=51 or len(original)!=1275 or sum(int(r["ports"]) for r in eligible)!=194:
        raise ValueError("Original butterfly or biological port source changed")
    source_pairs=set()
    for row in original:
        a,b=row["island_a"],row["island_b"]
        if a not in ids or b not in ids or a==b:
            raise ValueError("Source pair endpoints invalid")
        ordered=tuple(sorted((a,b),key=int))
        if ordered in source_pairs:
            raise ValueError("Repeated source island pair")
        source_pairs.add(ordered)
        for target in TARGETS:
            value=float(row[target])
            if not math.isfinite(value) or not 0<=value<=1:
                raise ValueError("Original observed butterfly similarity invalid")
    geo=rows(SOURCE_GEO)
    if len(geo)!=51 or len({r["island_id"] for r in geo})!=51:
        raise ValueError("Original 51 source geography rows absent")
    if {r["island_id"] for r in geo}!=ids:
        raise ValueError("Original pair and geographical island identities differ")
    latitude={}
    for r in geo:
        sid=r["island_id"]
        value=float(r["representative_latitude"])
        if r["landmass_id"]!=sid or not math.isfinite(value) or value==0 or not(-90<=value<=90):
            raise ValueError("Invalid original source signed latitude")
        latitude[sid]=value
    south=sorted([s for s in ids if latitude[s]<0],key=lambda s:(latitude[s],int(s)))
    north=sorted([s for s in ids if latitude[s]>0],key=lambda s:(latitude[s],int(s)))
    if len(south)!=10 or len(north)!=41:
        raise ValueError("Expected ten southern and forty-one northern source islands")
    band={}
    for code,original_order,sizes in (
            ("S",south,(5,5)),("N",north,(6,5,5,5,5,5,5,5))):
        cursor=0
        for k,n in enumerate(sizes):
            for sid in original_order[cursor:cursor+n]:
                band[sid]=code+str(k)
            cursor+=n
        if cursor!=len(original_order):
            raise ValueError("Latitude-band assignments incomplete")
    sorted_lon=sorted(geo,key=lambda r:(float(r["representative_longitude"]),int(r["island_id"])))
    long_group={}
    for idx,r in enumerate(sorted_lon):
        long_group[r["island_id"]]=0 if idx<11 else 1+(idx-11)//10
    if Counter(long_group.values())!=Counter({0:11,1:10,2:10,3:10,4:10}):
        raise ValueError("Longitude-block support differs")
    fid={}
    kind={}
    for idx,r in enumerate(original):
        k,h=sorted([long_group[r["island_a"]],long_group[r["island_b"]]])
        fid[idx]="B"+str(k)+"_B"+str(h)
        kind[idx]="within_longitude_block" if k==h else "between_longitude_blocks"
    if Counter(kind.values())!=Counter({"within_longitude_block":235,"between_longitude_blocks":1040}):
        raise ValueError("Original source pair geographic class differs")
    folds={}
    for k in range(5):
        for h in range(k,5):
            name="B"+str(k)+"_B"+str(h)
            held={sid for sid,group in long_group.items() if group in (k,h)}
            train=[r for r in original if r["island_a"] not in held and r["island_b"] not in held]
            test=[i for i in fid if fid[i]==name]
            if len(train)!=math.comb(51-len(held),2):
                raise ValueError("Original double-unseen training population does not match held-out islands")
            n0=math.comb(11 if k==0 else 10,2) if k==h else (
                (11 if k==0 else 10)*(11 if h==0 else 10))
            if len(test)!=n0:
                raise ValueError("Original longitude test fold count differs")
            folds[name]={"held":held,"train_count":len(train),"test_count":len(test)}
    return original,ids,latitude,band,fid,kind,folds

def validate_assignments(output,ids,latitude,band):
    source_bands=rows(output/"latitude_bands.csv")
    if len(source_bands)!=51 or {r["island_id"] for r in source_bands}!=ids:
        raise ValueError("Saved original latitude bands lack 51 islands")
    for row in source_bands:
        sid=row["island_id"]
        if (row["latitude_band"]!=band[sid] or
            row["hemisphere"]!=band[sid][0] or
            not math.isclose(float(row["representative_latitude"]),latitude[sid],abs_tol=1e-12) or
            int(row["band_size"])!=(6 if band[sid]=="N0" else 5)):
            raise ValueError("Saved target-independent source band altered")
    ledger=rows(output/"AIS_vertex_assignments.csv")
    if len(ledger)!=10200:
        raise ValueError("Expected 200 original/199 null maps times 51 islands")
    maps=defaultdict(dict)
    for row in ledger:
        replicate=int(row["replicate"])
        a,b=row["source_island_id"],row["donor_island_id"]
        if replicate not in range(200) or a not in ids or b not in ids or a in maps[replicate]:
            raise ValueError("Invalid/duplicate real or surrogate node mapping")
        if band[a]!=band[b] or row["latitude_band"]!=band[a]:
            raise ValueError("AIS node label illegally crosses hemisphere/latitude group")
        maps[replicate][a]=b
    if set(maps)!=set(range(200)):
        raise ValueError("200 complete identity/surrogate label maps missing")
    encoded=set()
    order=sorted(ids,key=int)
    for rep in range(200):
        mapping=maps[rep]
        if set(mapping)!=ids or set(mapping.values())!=ids:
            raise ValueError("Saved AIS permutation not bijective")
        flat=tuple(mapping[k] for k in order)
        if flat in encoded:
            raise ValueError("Duplicated original AIS node permutation")
        encoded.add(flat)
        if (rep==0 and any(mapping[k]!=k for k in ids) or
            rep>0 and all(mapping[k]==k for k in ids)):
            raise ValueError("Observed rep 0 not identity, or surrogate is the identity")
    return maps

def check_reference(output,original,fid):
    refs={}
    counts=Counter()
    baseline=rows(output/"reference_predictions.csv")
    if len(baseline)!=5100:
        raise ValueError("Expected exactly 5100 baseline held-out source pair predictions")
    for row in baseline:
        t,c=row["target"],row["contrast"]
        idx=int(row["source_pair_index"])
        if t not in TARGETS or c not in CONTRASTS or idx not in range(1275):
            raise ValueError("Reference model source/target outside 51-island sample")
        source=original[idx]
        if (row["fold_id"]!=fid[idx] or row["island_a"]!=source["island_a"] or
            row["island_b"]!=source["island_b"]):
            raise ValueError("Reference pair/fold source does not match original")
        observed,pred=float(row["observed"]),float(row["predicted"])
        if not(0<=observed<=1 and 0<=pred<=1):
            raise ValueError("Saved reference observation/prediction outside [0,1]")
        close(observed,float(source[t]),"Original butterfly source observation")
        sq=(observed-pred)**2
        ae=abs(observed-pred)
        close(sq,float(row["squared_error"]),"Baseline squared error")
        close(ae,float(row["absolute_error"]),"Baseline absolute error")
        key=(t,c,idx)
        if key in refs:
            raise ValueError("Duplicate source reference prediction")
        refs[key]=(sq,ae)
        counts[(t,c)]+=1
    if len(refs)!=5100 or len(counts)!=4 or set(counts.values())!={1275}:
        raise ValueError("All baseline source predictions are incomplete")
    return refs

def audit_predictions(output,original,fid,kind,refs):
    aggregates=defaultdict(lambda:[0,0.,0.,0.,0.,0,0,0])
    groups=0
    last_group=None
    seen=set()
    n=0
    with gzip.open(output/"AIS_reassigned_augmented_predictions.csv.gz","rt",
                   newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rep=int(row["replicate"]);t=row["target"];c=row["contrast"];i=int(row["source_pair_index"])
            if rep not in range(200) or t not in TARGETS or c not in CONTRASTS or i not in range(1275):
                raise ValueError("Unexpected source pair, butterfly period or AIS null replicate")
            grp=(rep,t,c)
            if last_group is None:
                last_group=grp
            if grp!=last_group:
                if len(seen)!=1275:
                    raise ValueError("Original pair coverage incomplete for a permutation/group")
                groups+=1
                last_group=grp
                seen=set()
            if i in seen:
                raise ValueError("Repeated original pair prediction in a null source map")
            seen.add(i)
            if row["fold_id"]!=fid[i]:
                raise ValueError("AIS-null prediction tested on wrong source longitude fold")
            observed=float(row["observed"]);pred=float(row["predicted"])
            if not(0<=observed<=1 and 0<=pred<=1):
                raise ValueError("AIS-null observed/predicted Jaccard outside [0,1]")
            close(observed,float(original[i][t]),"Original AIS-null target observation")
            sq=(pred-observed)**2;ae=abs(pred-observed)
            close(sq,float(row["squared_error"]),"AIS-null augmented squared error")
            close(ae,float(row["absolute_error"]),"AIS-null augmented absolute error")
            r_sq,r_abs=refs[(t,c,i)]
            for scope in ("all_pairs",kind[i]):
                s=aggregates[(rep,t,c,scope)]
                s[0]+=1;s[1]+=r_sq;s[2]+=sq;s[3]+=r_abs;s[4]+=ae
                s[5]+=int(r_sq>sq);s[6]+=int(r_sq<sq);s[7]+=int(r_abs>ae)
            n+=1
    if last_group is not None:
        if len(seen)!=1275:
            raise ValueError("Last AIS-null permutation lacked original island pairs")
        groups+=1
    if n!=1020000 or groups!=800 or len(aggregates)!=2400:
        raise ValueError("Require exactly 1,020,000 predictions / 800 groups / 2400 scores")
    return aggregates,n,groups

def reconstruct_summaries(output,aggregates):
    scores={}
    results=rows(output/"AIS_relabeling_scope_scores.csv")
    if len(results)!=2400:
        raise ValueError("2400 original and placebo scope summaries required")
    for r in results:
        rep=int(r["replicate"])
        t,c,scope=r["target"],r["contrast"],r["pair_class"]
        key=(rep,t,c,scope)
        if key in scores or key not in aggregates or scope not in SCOPES:
            raise ValueError("Unexpected/duplicate AIS-null scoped comparison")
        n,sr,sa,ar,aa,better,worse,abetter=aggregates[key]
        if n!=SCOPES[scope] or int(r["n_pairs"])!=n:
            raise ValueError("Wrong source pair scope count in AIS placebo scores")
        r_mse=sr/n;a_mse=sa/n;r_mae=ar/n;a_mae=aa/n
        actual={
          "MSE_reference":r_mse,"MSE_augmented":a_mse,
          "MSE_gain_percent":100*(r_mse-a_mse)/r_mse,
          "MAE_reference":r_mae,"MAE_augmented":a_mae,
          "MAE_gain_percent":100*(r_mae-a_mae)/r_mae,
          "pairs_MSE_improved":better,"pairs_MSE_worse":worse,
          "pairs_MAE_improved":abetter,
        }
        for field,value in actual.items():
            close(float(r[field]),float(value),"AIS-null source scope "+field)
        scores[key]=actual
    if len(scores)!=2400:
        raise ValueError("Incomplete original and null scope score set")
    return scores

def quantile_linear(numbers,p):
    s=sorted(numbers)
    if not len(s):
        raise ValueError("No null samples for percentile")
    x=(len(s)-1)*p
    low=int(math.floor(x));high=int(math.ceil(x))
    return s[low]*(high-x)+s[high]*(x-low) if high!=low else s[low]

def verify_ranks(output,scores):
    saved=rows(output/"observed_vs_latitude_node_null.csv")
    if len(saved)!=12:
        raise ValueError("Expected 12 actual-versus-199 placebo rank comparisons")
    seen=set()
    for r in saved:
        t,c,scope=r["target"],r["contrast"],r["pair_class"]
        key=(t,c,scope)
        if key in seen or t not in TARGETS or c not in CONTRASTS or scope not in SCOPES:
            raise ValueError("Duplicated/unknown AIS-null rank record")
        seen.add(key)
        if int(r["unique_test_pairs"])!=SCOPES[scope] or int(r["null_replicates"])!=199:
            raise ValueError("Original source-rank support differs")
        for measure in ("MSE","MAE"):
            column=measure+"_gain_percent"
            observed=scores[(0,t,c,scope)][column]
            null=[scores[(rep,t,c,scope)][column] for rep in range(1,200)]
            q025,q050,q975=(quantile_linear(null,q) for q in (.025,.5,.975))
            exceed=sum(x>=observed for x in null)
            expected={
                measure+"_observed_gain_percent":observed,
                measure+"_null_025_gain_percent":q025,
                measure+"_null_median_gain_percent":q050,
                measure+"_null_975_gain_percent":q975,
                measure+"_null_mean_gain_percent":math.fsum(null)/len(null),
                measure+"_null_ge_observed_count":exceed,
                measure+"_descriptive_upper_tail_rank_fraction":(1+exceed)/200,
            }
            for field,val in expected.items():
                close(float(r[field]),float(val),"AIS empirical source rank "+field)
    if len(seen)!=12:
        raise ValueError("Missing all 12 AIS-versus-placebo source comparisons")
    return len(seen)

def validate_parent_true_result(scores):
    parent=rows(PARENT)
    if len(parent)!=120:
        raise ValueError("Earlier PR135 source comparisons incomplete")
    original=[r for r in parent if r["scenario"]=="FULL" and
              r["contrast"] in PARENT_NAMES.values()]
    if len(original)!=12:
        raise ValueError("Cannot independently confirm all 12 parent AIS contrasts")
    maxdiff=0.0
    for r in original:
        t,scope=r["target"],r["pair_class"]
        c=next(k for k,v in PARENT_NAMES.items() if v==r["contrast"])
        observed=scores[(0,t,c,scope)]
        if int(r["unique_evaluated_pairs"])!=SCOPES[scope]:
            raise ValueError("Historical source model contrast population changed")
        for savedkey,expectedkey in (("MSE_reference","MSE_reference"),
                                    ("MSE_augmented","MSE_augmented"),
                                    ("MAE_reference","MAE_reference"),
                                    ("MAE_augmented","MAE_augmented")):
            val=float(r[savedkey]);other=observed[expectedkey]
            maxdiff=max(maxdiff,abs(val-other))
            if abs(val-other)>1e-9:
                raise ValueError("True AIS graph failed parent independent model reproduction")
    return maxdiff

def audit(result_dir):
    output=Path(result_dir)
    original,ids,latitude,bands,fids,kinds,folds=sources()
    maps=validate_assignments(output,ids,latitude,bands)
    refs=check_reference(output,original,fids)
    aggregates,n,groups=audit_predictions(output,original,fids,kinds,refs)
    scores=reconstruct_summaries(output,aggregates)
    rank_count=verify_ranks(output,scores)
    discrepancy=validate_parent_true_result(scores)
    result={
        "schema":"chocho_AIS_latitude_node_relabel_independent_source_audit_v1",
        "status":"ALL_MAPS_SOURCE_PAIRS_PREDICTIONS_SCOPE_SCORES_AND_RANKS_RECONSTRUCTED",
        "source_islands":51,"source_pairs":1275,"botanical_ports":194,
        "latitude_bands_verified":10,"AIS_source_node_assignment_rows":10200,
        "different_nonidentity_source_maps":199,
        "reference_predictions_verified":len(refs),
        "AIS_augmented_relabel_predictions_verified":n,
        "source_model_target_replicate_groups":groups,
        "source_scoped_loss_rows_verified":len(scores),
        "actual_vs_null_ranks_verified":rank_count,
        "original_parent_source_scores_max_abs_difference":discrepancy,
        "all_test_island_endpoints_excluded_from_train":True,
        "AIS_surrogate_maps_preserve_original_hemisphere_and_latitude_band":True,
        "descriptive_permutation_ranks_not_causal_null_p_values":True,
    }
    (output/"independent_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("AIS_NODE_NULL_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",type=Path,required=True)
    audit(p.parse_args().results_dir)
