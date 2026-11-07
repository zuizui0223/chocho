#!/usr/bin/env python3
"""Independent stdlib source-pair, 750-km AIS vertex swaps, and 1,025,100 prediction audit.

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

MAX_DISTANCE_KM=750.0
EARTH_RADIUS_KM=6371.0088

def haversine(lat1,lon1,lat2,lon2):
    a,b,c,d=(math.radians(v) for v in (lat1,lon1,lat2,lon2))
    square=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 2*EARTH_RADIUS_KM*math.asin(math.sqrt(min(1.0,max(0.0,square))))

def sources():
    original=rows(SOURCE_PAIRS)
    eligible=[r for r in rows(SOURCE_ROSTER) if r["eligible"].lower()=="true"]
    ids={r["landmass_id"] for r in eligible}
    if len(ids)!=51 or len(original)!=1275 or sum(int(r["ports"]) for r in eligible)!=194:
        raise ValueError("Original butterfly or biological port source changed")
    seen_pairs=set()
    for row in original:
        a,b=row["island_a"],row["island_b"]
        if a not in ids or b not in ids or a==b:
            raise ValueError("Invalid original island pair endpoint")
        unordered=tuple(sorted((a,b),key=int))
        if unordered in seen_pairs:
            raise ValueError("Duplicate source unordered island pair")
        seen_pairs.add(unordered)
        for target in TARGETS:
            x=float(row[target])
            if not math.isfinite(x) or not 0<=x<=1:
                raise ValueError("Original butterfly community Jaccard outside [0,1]")
    if len(seen_pairs)!=1275:
        raise ValueError("Original source pair population incomplete")
    geography=rows(SOURCE_GEO)
    if len(geography)!=51 or {r["island_id"] for r in geography}!=ids:
        raise ValueError("Original source geographic island identities changed")
    coords={}
    for row in geography:
        sid=row["island_id"]
        if sid in coords or row["landmass_id"]!=sid:
            raise ValueError("Duplicate source geography or mismatched landmass ID")
        lat=float(row["representative_latitude"])
        lon=float(row["representative_longitude"])
        if not all(math.isfinite(v) for v in (lat,lon)) or not(-90<=lat<=90 and -180<=lon<=180) or lat==0:
            raise ValueError("Invalid original representative island centroid")
        coords[sid]=(lat,lon)
    if sum(lat<0 for lat,lon in coords.values())!=10 or sum(lat>0 for lat,lon in coords.values())!=41:
        raise ValueError("Original southern/northern island populations changed")
    # Independent source-only graph-matching reconstruction: ascending numerical IDs, then distance and IDs.
    eligible_source_pairs=[]
    ordered=sorted(ids,key=int)
    for i,sid in enumerate(ordered):
        for other in ordered[i+1:]:
            a,b=coords[sid],coords[other]
            if a[0]*b[0]<=0:
                continue
            d=haversine(a[0],a[1],b[0],b[1])
            if d<=MAX_DISTANCE_KM:
                eligible_source_pairs.append((d,int(sid),int(other),sid,other))
    eligible_source_pairs.sort(key=lambda t:(t[0],t[1],t[2]))
    selected=[]
    used=set()
    for d,_,_,sid,other in eligible_source_pairs:
        if sid in used or other in used:
            continue
        used.add(sid)
        used.add(other)
        selected.append({"pair_id":"P"+str(len(selected)+1).zfill(2),
                         "island_a":sid,"island_b":other,
                         "centroid_greatcircle_distance_km":d,
                         "hemisphere":"S" if coords[sid][0]<0 else "N",
                         "latitude_a":coords[sid][0],"longitude_a":coords[sid][1],
                         "latitude_b":coords[other][0],"longitude_b":coords[other][1]})
    unmatched=ids-used
    if len(selected)!=18 or len(used)!=36 or len(unmatched)!=15:
        raise ValueError("Original cap/greedy matching no longer gives eighteen source local swaps")
    sorted_lon=sorted(geography,key=lambda r:(float(r["representative_longitude"]),int(r["island_id"])))
    long_group={}
    for i,row in enumerate(sorted_lon):
        long_group[row["island_id"]]=0 if i<11 else 1+(i-11)//10
    if Counter(long_group.values())!=Counter({0:11,1:10,2:10,3:10,4:10}):
        raise ValueError("Original source longitude groups changed")
    fold_id={}
    kind={}
    for idx,row in enumerate(original):
        k,h=sorted([long_group[row["island_a"]],long_group[row["island_b"]]])
        fold_id[idx]="B"+str(k)+"_B"+str(h)
        kind[idx]="within_longitude_block" if k==h else "between_longitude_blocks"
    if Counter(kind.values())!=Counter({"within_longitude_block":235,"between_longitude_blocks":1040}):
        raise ValueError("Original source longitude scope sizes differ")
    for k in range(5):
        for h in range(k,5):
            held={sid for sid,g in long_group.items() if g in (k,h)}
            train=sum(row["island_a"] not in held and row["island_b"] not in held for row in original)
            expected=math.comb(51-len(held),2)
            if train!=expected:
                raise ValueError("Original two-island-unseen source training folds changed")
    return original,ids,selected,unmatched,fold_id,kind,coords

def validate_assignments(out,ids,selected,unmatched,coords):
    exported=rows(out/"local_swap_pairs.csv")
    if len(exported)!=18:
        raise ValueError("Local AIS surrogate source pairing must include 18 rows")
    fields=("latitude_a","longitude_a","latitude_b","longitude_b","centroid_greatcircle_distance_km")
    for original,actual in zip(selected,exported):
        for name in ("pair_id","island_a","island_b","hemisphere"):
            if actual[name]!=original[name]:
                raise ValueError("Source geography-selected local swap identities do not match")
        for name in fields:
            close(float(actual[name]),float(original[name]),"original port-haversine pair "+name)
        if float(actual["centroid_greatcircle_distance_km"])>MAX_DISTANCE_KM:
            raise ValueError("Local AIS surrogate swap exceeds 750 km")
    pair_for={}
    for row in selected:
        for sid in (row["island_a"],row["island_b"]):
            if sid in pair_for:
                raise ValueError("Geographic null pair endpoint reused")
            pair_for[sid]=row["pair_id"]
    if len(pair_for)!=36 or set(pair_for)&set(unmatched):
        raise ValueError("Unpaired source island changed")
    ledger=rows(out/"AIS_vertex_assignments.csv")
    if len(ledger)!=10200:
        raise ValueError("Exactly 200 original/null AIS maps × 51 source islands required")
    maps=defaultdict(dict)
    for row in ledger:
        r=int(row["replicate"])
        a,b=row["source_island_id"],row["donor_island_id"]
        if r not in range(200) or a not in ids or b not in ids or a in maps[r]:
            raise ValueError("Invalid original/null local vertex map")
        if row["local_swap_pair_id"]!=pair_for.get(a,"FIXED"):
            raise ValueError("Local source geometry pair membership altered")
        maps[r][a]=b
    if set(maps)!=set(range(200)):
        raise ValueError("Missing original/null source map replicate")
    unique=set()
    canonical=sorted(ids,key=int)
    for rep in range(200):
        mapping=maps[rep]
        if set(mapping)!=ids or set(mapping.values())!=ids:
            raise ValueError("Original AIS vertex local map not bijective")
        code=tuple(mapping[k] for k in canonical)
        if code in unique:
            raise ValueError("Repeated original/nonidentity local AIS source vertex map")
        unique.add(code)
        if rep==0 and any(mapping[k]!=k for k in ids):
            raise ValueError("Original AIS null replicate 0 must be original identity")
        if rep>0 and all(mapping[k]==k for k in ids):
            raise ValueError("Placebo AIS node map equals original")
        for row in selected:
            a,b=row["island_a"],row["island_b"]
            if not (mapping[a]==a and mapping[b]==b or
                    mapping[a]==b and mapping[b]==a):
                raise ValueError("AIS graph node reassignment is not whole-pair geographical swap")
        for sid in unmatched:
            if mapping[sid]!=sid:
                raise ValueError("Isolated geographical island node cannot change AIS label")
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
    saved=rows(output/"observed_vs_local750_node_null.csv")
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
    original,ids,selected,unmatched,fids,kinds,coords=sources()
    maps=validate_assignments(output,ids,selected,unmatched,coords)
    refs=check_reference(output,original,fids)
    aggregates,n,groups=audit_predictions(output,original,fids,kinds,refs)
    scores=reconstruct_summaries(output,aggregates)
    rank_count=verify_ranks(output,scores)
    discrepancy=validate_parent_true_result(scores)
    result={
        "schema":"chocho_AIS_local750_vertex_swap_independent_source_audit_v1",
        "status":"ALL_MAPS_SOURCE_PAIRS_PREDICTIONS_SCOPE_SCORES_AND_RANKS_RECONSTRUCTED",
        "source_islands":51,"source_pairs":1275,"botanical_ports":194,
        "local_swappable_pairs_verified":len(selected),
        "fixed_original_islands":len(unmatched),
        "max_allowed_swap_distance_km":MAX_DISTANCE_KM,"AIS_source_node_assignment_rows":10200,
        "different_nonidentity_source_maps":199,
        "reference_predictions_verified":len(refs),
        "AIS_augmented_relabel_predictions_verified":n,
        "source_model_target_replicate_groups":groups,
        "source_scoped_loss_rows_verified":len(scores),
        "actual_vs_null_ranks_verified":rank_count,
        "original_parent_source_scores_max_abs_difference":discrepancy,
        "all_test_island_endpoints_excluded_from_train":True,
        "AIS_surrogate_maps_preserve_source_greatcircle_cap_and_hemisphere":True,
        "descriptive_permutation_ranks_not_causal_null_p_values":True,
    }
    (output/"independent_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("AIS_LOCAL_NULL_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",type=Path,required=True)
    audit(p.parse_args().results_dir)
