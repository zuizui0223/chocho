#!/usr/bin/env python3
"""Compare AIS network alignment to local (at most 750 km) whole-graph vertex swaps."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from exploratory_island51_circuit_network_holdout import source, train_predict
from run_island51_host_circuit_joint import INPUT_PATHS, TARGETS, load_original_port_counts, load_resistance_matrix, validate_input_lock
from run_island51_latitude_controls import (
    ROOT, PAIRS, ROSTER, GEOGRAPHY, GEOGRAPHY_GIT_BLOB,
    make_features, latitude_pair_controls, partition_geographic_blocks, holdout_plan,
)
from run_island51_pre2011_effort_control import git_blob_sha
from run_island51_observed_ocean_circuits import sha256

SEED=20261107
NULL_REPS=199
COMPARE=(
    ("AIS_without_wind","HRE","HRE_AIS","HRE_to_HRE_AIS_LAT"),
    ("AIS_after_wind","HRE_WIND","HCRE","HRE_WIND_to_HCRE_LAT"),
)
SCOPES=("all_pairs","within_longitude_block","between_longitude_blocks")
PARENT=ROOT/"results/island51_top_port_training_influence_20261107/contrast_scope_summary.csv"
PARENT_GIT_BLOB="6f752f251edcb8fc47ae164416a59bde7dc42cf1"
MAX_SOURCE_SWAP_DISTANCE_KM=750.0
EARTH_RADIUS_KM=6371.0088

def source_haversine_km(lat1,lon1,lat2,lon2):
    a,b,c,d=np.radians([lat1,lon1,lat2,lon2])
    dlat=c-a
    dlon=d-b
    square=np.sin(dlat/2)**2+np.cos(a)*np.cos(c)*np.sin(dlon/2)**2
    return float(2*EARTH_RADIUS_KM*np.arcsin(np.sqrt(min(1.0,max(0.0,float(square))))))

def freeze_local_swap_pairs(ids):
    if git_blob_sha(GEOGRAPHY.read_bytes())!=GEOGRAPHY_GIT_BLOB:
        raise ValueError("Source WPI representative island coordinates changed")
    geo=pd.read_csv(GEOGRAPHY,dtype={"island_id":str,"landmass_id":str})
    if (len(geo)!=51 or geo.island_id.nunique()!=51 or
        set(geo.island_id)!=set(ids) or not geo.island_id.eq(geo.landmass_id).all()):
        raise ValueError("Original 51 island geography identifiers changed")
    lat=dict(zip(geo.island_id,geo.representative_latitude.astype(float)))
    lon=dict(zip(geo.island_id,geo.representative_longitude.astype(float)))
    if (len(lat)!=51 or len(lon)!=51 or
        any(not np.isfinite(lat[k]) or not np.isfinite(lon[k]) or
            abs(lat[k])>90 or abs(lon[k])>180 or lat[k]==0 for k in ids)):
        raise ValueError("Invalid original signed port representative coordinates")
    if sum(lat[k]<0 for k in ids)!=10 or sum(lat[k]>0 for k in ids)!=41:
        raise ValueError("Source hemisphere cohort differs from frozen original")
    ordered=sorted(ids,key=int)
    eligible=[]
    for index,x in enumerate(ordered):
        for y in ordered[index+1:]:
            if lat[x]*lat[y]<=0:
                continue
            distance=source_haversine_km(lat[x],lon[x],lat[y],lon[y])
            if distance<=MAX_SOURCE_SWAP_DISTANCE_KM:
                eligible.append((distance,int(x),int(y),x,y))
    eligible.sort(key=lambda x:(x[0],x[1],x[2]))
    consumed=set()
    membership=[]
    table=[]
    position={v:i for i,v in enumerate(ids)}
    for distance,_,_,x,y in eligible:
        if x in consumed or y in consumed:
            continue
        consumed.update((x,y))
        pair_id="P"+str(len(table)+1).zfill(2)
        table.append({
           "pair_id":pair_id,"island_a":x,"island_b":y,
           "centroid_greatcircle_distance_km":distance,
           "hemisphere":"S" if lat[x]<0 else "N",
           "latitude_a":lat[x],"longitude_a":lon[x],
           "latitude_b":lat[y],"longitude_b":lon[y],
        })
        membership.append(np.array([position[x],position[y]],dtype=int))
    unmatched=sorted(set(ids)-consumed,key=int)
    if len(membership)!=18 or len(table)!=18 or len(consumed)!=36 or len(unmatched)!=15:
        raise ValueError("Original source geometry no longer gives exactly 18 local pairs")
    if (max(x["centroid_greatcircle_distance_km"] for x in table)>
        MAX_SOURCE_SWAP_DISTANCE_KM or
        not all(x["centroid_greatcircle_distance_km"]>0 for x in table)):
        raise ValueError("One local swap is too far or invalid")
    return membership,pd.DataFrame(table),unmatched

def generate_local_vertex_assignments(ids,memberships):
    if len(memberships)!=18 or len({i for pair in memberships for i in pair})!=36:
        raise ValueError("Need 18 nonoverlapping local node pairs")
    rng=np.random.default_rng(SEED)
    identity=np.arange(len(ids),dtype=int)
    result=[identity.copy()]
    seen={tuple(identity)}
    attempts=0
    while len(result)<NULL_REPS+1:
        attempts+=1
        if attempts>100000:
            raise ValueError("Too many repeated local AIS null maps")
        mapping=identity.copy()
        switches=rng.integers(0,2,size=len(memberships))
        for selected,members in zip(switches,memberships):
            if selected:
                i,j=int(members[0]),int(members[1])
                mapping[i],mapping[j]=mapping[j],mapping[i]
        key=tuple(mapping)
        if key in seen:
            continue
        if sorted(mapping.tolist())!=list(range(len(ids))):
            raise ValueError("AIS source vertex map is not bijective")
        if any(set(mapping[members])!=set(members) for members in memberships):
            raise ValueError("AIS source node moved outside its local pair")
        all_assigned={int(i) for pair in memberships for i in pair}
        if any(int(mapping[i])!=i for i in range(len(ids)) if i not in all_assigned):
            raise ValueError("Unpaired AIS source vertex shifted")
        seen.add(key)
        result.append(mapping)
    return result,attempts

def ais_source_matrix(q,ids,source_pair_values):
    n=len(ids)
    if len(q)!=1275 or len(source_pair_values)!=1275:
        raise ValueError("Original native AIS support must contain all original island pairs")
    result=np.full((n,n),np.nan,dtype=float)
    np.fill_diagonal(result,0.)
    a=q.i.to_numpy(dtype=int)
    b=q.j.to_numpy(dtype=int)
    if np.any(a==b) or len(set(zip(np.minimum(a,b),np.maximum(a,b))))!=1275:
        raise ValueError("Unexpected source island-pair indexing")
    if not np.isfinite(source_pair_values).all():
        raise ValueError("Nonfinite original AIS circuit features")
    result[a,b]=source_pair_values
    result[b,a]=source_pair_values
    if not np.isfinite(result).all() or not np.allclose(result,result.T,atol=0,rtol=0):
        raise ValueError("Source AIS matrix did not reconstruct as full symmetric graph")
    return result

def metrics(reference,augmented,observation,subset):
    r_err=(observation[subset]-reference[subset])
    a_err=(observation[subset]-augmented[subset])
    rmse=float(np.mean(r_err**2))
    amse=float(np.mean(a_err**2))
    rmae=float(np.mean(np.abs(r_err)))
    amae=float(np.mean(np.abs(a_err)))
    if rmse<=0 or rmae<=0:
        raise ValueError("Positive reference MSE/MAE required")
    return {"n_pairs":int(np.count_nonzero(subset)),
            "MSE_reference":rmse,"MSE_augmented":amse,
            "MSE_gain_percent":100*(rmse-amse)/rmse,
            "MAE_reference":rmae,"MAE_augmented":amae,
            "MAE_gain_percent":100*(rmae-amae)/rmae,
            "pairs_MSE_improved":int(np.sum(r_err**2>a_err**2)),
            "pairs_MSE_worse":int(np.sum(r_err**2<a_err**2)),
            "pairs_MAE_improved":int(np.sum(np.abs(r_err)>np.abs(a_err)))}

def assert_parent_result(observed):
    if git_blob_sha(PARENT.read_bytes())!=PARENT_GIT_BLOB:
        raise ValueError("Original PR135 FULL AIS comparison table changed")
    old=pd.read_csv(PARENT)
    original=old.loc[old.scenario.eq("FULL")]
    if len(old)!=120 or len(original)!=30:
        raise ValueError("Original complete original 120-row source result missing")
    if len(observed)!=12:
        raise ValueError("Expected all twelve true original source contrasts")
    maxdiff=0.0
    for row in observed:
        contrast=next(x[3] for x in COMPARE if x[0]==row["contrast"])
        match=original[(original.target==row["target"])&
                       (original.contrast==contrast)&
                       (original.pair_class==row["pair_class"])]
        if len(match)!=1 or int(match.iloc[0].unique_evaluated_pairs)!=row["n_pairs"]:
            raise ValueError("Original source comparison reference population differs")
        for key,oldkey in (("MSE_reference","MSE_reference"),
                           ("MSE_augmented","MSE_augmented"),
                           ("MAE_reference","MAE_reference"),
                           ("MAE_augmented","MAE_augmented")):
            d=abs(float(row[key])-float(match.iloc[0][oldkey]))
            maxdiff=max(maxdiff,d)
            if d>1e-9:
                raise ValueError("Original PR135 true AIS predictor does not reproduce: "+key)
    return maxdiff

def run(out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("data","provenance","results")):
        raise ValueError("Never overwrite original source or parent scientific results")
    lock=validate_input_lock(ROOT/"results/island51_host_circuit_joint_20261106/input_lock.json")
    q,ids,_=source(PAIRS,ROSTER)
    if len(q)!=1275 or len(ids)!=51:
        raise ValueError("Original 51 landmasses or 1275 source pairs changed")
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    ports=load_original_port_counts(ROSTER,ids,prep)
    if len(ports)!=51 or sum(ports.values())!=194:
        raise ValueError("Original 194 botanical port grouping changed")
    circuits={name:load_resistance_matrix(ROOT/INPUT_PATHS[name],ids)
              for name in ("uniform","ais","wind")}
    longitude_blocks,geo_spec=partition_geographic_blocks(GEOGRAPHY,ids)
    folds,original_fold_spec=holdout_plan(q,ids,longitude_blocks)
    memberships,local_pair_table,unmatched=freeze_local_swap_pairs(ids)
    mappings,attempts=generate_local_vertex_assignments(ids,memberships)
    features,feature_names,historical=make_features(q,ids,ports,circuits)
    lat=latitude_pair_controls(q,ids)
    left=q.i.to_numpy(dtype=int);right=q.j.to_numpy(dtype=int)
    initial=np.asarray(features[TARGETS[0]]["HRE_AIS"][:,feature_names[TARGETS[0]]["HRE_AIS"].index("log1p_normalized_ais_circuit")],dtype=float)
    for target in TARGETS[1:]:
        copied=features[target]["HRE_AIS"][:,feature_names[target]["HRE_AIS"].index("log1p_normalized_ais_circuit")]
        if not np.array_equal(initial,copied):
            raise ValueError("Original AIS source feature changed between butterfly responses")
    ais_matrix=ais_source_matrix(q,ids,initial)
    within=np.asarray([longitude_blocks[ids[a]]==longitude_blocks[ids[b]]
                       for a,b in zip(left,right)],dtype=bool)
    if int(within.sum())!=235 or int((~within).sum())!=1040:
        raise ValueError("Original fixed longitude-block scope support changed")
    masks={"all_pairs":np.ones(1275,dtype=bool),
           "within_longitude_block":within,"between_longitude_blocks":~within}
    true_response={t:q[t].to_numpy(dtype=float) for t in TARGETS}
    reference_design={}
    reference_prediction={}
    source_pair_assignment=np.zeros(1275,dtype=np.uint8)
    for fid,train,test in folds:
        source_pair_assignment[test]+=1
    if not np.array_equal(source_pair_assignment,np.ones(1275,dtype=np.uint8)):
        raise ValueError("Original double-unseen folds fail unique pair coverage")
    for target in TARGETS:
        for contrast,reference,augmented,_ in COMPARE:
            original=features[target][reference]
            refdesign=np.column_stack([original]+[lat[key] for key in lat])
            reference_design[(target,contrast)]=refdesign
            estimates=np.full(1275,np.nan)
            for fid,train,test in folds:
                estimates[test]=train_predict(refdesign,true_response[target],train,test)
            if not np.isfinite(estimates).all():
                raise ValueError("Missing source baseline predictions")
            reference_prediction[(target,contrast)]=estimates
    out.mkdir(parents=True,exist_ok=True)
    assigned={str(row.island_a):str(row.pair_id) for row in local_pair_table.itertuples(index=False)}
    assigned.update({str(row.island_b):str(row.pair_id) for row in local_pair_table.itertuples(index=False)})
    if len(assigned)!=36 or any(i in assigned for i in unmatched):
        raise ValueError("Source-only local pair membership invalid")
    pd.DataFrame([{"replicate":rep,"source_island_id":ids[index],
        "donor_island_id":ids[int(mapping[index])],
        "local_swap_pair_id":assigned.get(ids[index],"FIXED")}
        for rep,mapping in enumerate(mappings) for index in range(51)]
    ).to_csv(out/"AIS_vertex_assignments.csv",index=False)
    local_pair_table.to_csv(out/"local_swap_pairs.csv",index=False,float_format="%.17g")
    base_rows=[]
    for target in TARGETS:
        for contrast,*_ in COMPARE:
            observed=true_response[target]
            reference=reference_prediction[(target,contrast)]
            for idx in range(1275):
                fid=next(fid for fid,train,test in folds if test[idx])
                base_rows.append({"target":target,"contrast":contrast,
                    "fold_id":fid,"source_pair_index":idx,
                    "island_a":str(q.iloc[idx].island_a),"island_b":str(q.iloc[idx].island_b),
                    "observed":float(observed[idx]),"predicted":float(reference[idx]),
                    "squared_error":float((reference[idx]-observed[idx])**2),
                    "absolute_error":float(abs(reference[idx]-observed[idx]))})
    pd.DataFrame(base_rows).to_csv(out/"reference_predictions.csv",index=False,float_format="%.17g")

    summary_rows=[]
    all_rows=[]
    detail=out/"AIS_reassigned_augmented_predictions.csv.gz"
    with gzip.open(detail,"wt",encoding="utf-8",newline="",compresslevel=3) as file:
        writer=csv.writer(file)
        writer.writerow(["replicate","target","contrast","fold_id","source_pair_index",
                         "observed","predicted","squared_error","absolute_error"])
        for rep,mapping in enumerate(mappings):
            new_ais=ais_matrix[mapping[left],mapping[right]]
            if not np.isfinite(new_ais).all():
                raise ValueError("AIS vertex-reassigned circuit absent for a source pair")
            if rep==0 and not np.array_equal(new_ais,initial):
                raise ValueError("Original AIS identity-map circuit changed")
            for target in TARGETS:
                obs=true_response[target]
                for contrast,*_ in COMPARE:
                    original=reference_design[(target,contrast)]
                    design=np.column_stack([original,new_ais])
                    predicted=np.full(1275,np.nan)
                    for fid,train,test in folds:
                        pred=train_predict(design,obs,train,test)
                        predicted[test]=pred
                        for idx,v in zip(np.flatnonzero(test),pred):
                            o=float(obs[idx]);p=float(v)
                            writer.writerow((rep,target,contrast,fid,int(idx),
                                             format(o,".17g"),format(p,".17g"),
                                             format((o-p)**2,".17g"),format(abs(o-p),".17g")))
                    if not np.isfinite(predicted).all():
                        raise ValueError("Incomplete predictions for null relabeling "+str(rep))
                    ref=reference_prediction[(target,contrast)]
                    for kind in SCOPES:
                        score={"replicate":rep,"target":target,"contrast":contrast,
                               "pair_class":kind}
                        score.update(metrics(ref,predicted,obs,masks[kind]))
                        summary_rows.append(score)
                        if rep==0:
                            all_rows.append(score)
            if rep==0:
                original_max_diff=assert_parent_result(all_rows)
                print("SOURCE_ALIGNED_AIS_PARENT_MAX_DIFF="+format(original_max_diff,".17g"),flush=True)
            if rep in (0,49,99,149,199):
                print("AIS_LOCAL_NULL_REPLICATE_COMPLETE="+str(rep),flush=True)

    summary=pd.DataFrame(summary_rows)
    if len(summary)!=2400 or summary.duplicated(["replicate","target","contrast","pair_class"]).any():
        raise ValueError("Expected 2400 unique observed/null model/scope scores")
    comparison=[]
    for target in TARGETS:
        for contrast,*_ in COMPARE:
            for kind in SCOPES:
                subset=summary[(summary.target==target)&(summary.contrast==contrast)&
                               (summary.pair_class==kind)]
                if len(subset)!=200:
                    raise ValueError("Observed and 199 null scores incomplete")
                true=subset[subset.replicate.eq(0)].iloc[0]
                null=subset[subset.replicate.ne(0)]
                row={"target":target,"contrast":contrast,"pair_class":kind,
                     "unique_test_pairs":int(true.n_pairs),
                     "null_replicates":199}
                for short,column in (("MSE","MSE_gain_percent"),("MAE","MAE_gain_percent")):
                    actual=float(true[column]);draws=null[column].to_numpy(dtype=float)
                    quantile=np.quantile(draws,[.025,.5,.975])
                    exceed=int(np.sum(draws>=actual))
                    row.update({
                      short+"_observed_gain_percent":actual,
                      short+"_null_025_gain_percent":float(quantile[0]),
                      short+"_null_median_gain_percent":float(quantile[1]),
                      short+"_null_975_gain_percent":float(quantile[2]),
                      short+"_null_mean_gain_percent":float(np.mean(draws)),
                      short+"_null_ge_observed_count":exceed,
                      short+"_descriptive_upper_tail_rank_fraction":(1+exceed)/200,
                    })
                comparison.append(row)
    if len(comparison)!=12:
        raise ValueError("Expected twelve fixed observed/null source comparisons")
    summary.to_csv(out/"AIS_relabeling_scope_scores.csv",index=False,float_format="%.17g")
    pd.DataFrame(comparison).to_csv(out/"observed_vs_local750_node_null.csv",index=False,float_format="%.17g")
    manifest={"schema":"chocho_original51_750km_AIS_graph_vertex_swap_null_v1",
      "status":"199_MAXIMUM_750KM_SOURCE_PAIR_VERTEX_NULLS_EVALUATED",
      "source_islands":51,"biological_ports":194,"pairs":1275,
      "butterfly_targets":list(TARGETS),"contrasts":[x[0] for x in COMPARE],
      "longitude_double_unseen_folds":15,"within_pairs":235,"between_pairs":1040,
      "max_allowed_source_swap_km":MAX_SOURCE_SWAP_DISTANCE_KM,
      "earth_radius_km":EARTH_RADIUS_KM,
      "source_geographic_pairs":len(local_pair_table),
      "paired_source_vertices":36,
      "unpaired_immutable_island_ids":list(unmatched),
      "source_swap_pairs_max_geodesic_km":float(local_pair_table.centroid_greatcircle_distance_km.max()),
      "random_seed":SEED,"rng":"numpy.default_rng PCG64",
      "observed_identity_replicate":0,"distinct_nonidentity_permutations":NULL_REPS,
      "draw_attempts_including_accepted":attempts,
      "raw_reference_predictions":len(base_rows),
      "raw_augmented_predictions":200*len(TARGETS)*len(COMPARE)*1275,
      "scope_scores":len(summary),"observed_vs_null_rows":len(comparison),
      "true_source_parent_MSE_MAE_max_absolute_difference":original_max_diff,
      "source_input_lock_status":lock["status"],
      "source_geography_sha256":sha256(GEOGRAPHY),
      "historic_recording_effort_provenance":historical,
      "null_relabelings_are_not_independent_biological_or_geographic_replicates":True,
      "source_geodesic_distance_capped_for_each_swapped_pair":True,
      "null_preserves_complete_graph_topology_under_vertex_labels_not_actual_vessel_routes":True,
      "host_identity_pooled_2001_2024_can_postdate_butterfly_observations":True,
      "no_causal_null_p_values_or_butterfly_colonization_claim":True,
      "outputs_sha256":{p.name:sha256(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!="execution_receipt.json"}}
    (out/"execution_receipt.json").write_text(json.dumps(manifest,sort_keys=True,indent=2)+"\n")
    print("AIS_LOCAL_NULL_RESULT_JSON="+pd.DataFrame(comparison).to_json(orient="records",double_precision=12),flush=True)
    print("AIS_LOCAL_NULL_CHECK_JSON="+json.dumps({"maps":len(mappings),"reference_rows":len(base_rows),
      "augmented_rows":200*2*2*1275,"scope_rows":len(summary),
      "original_parent_max_difference":original_max_diff},sort_keys=True),flush=True)
    return manifest

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)
