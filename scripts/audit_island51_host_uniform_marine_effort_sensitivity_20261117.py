#!/usr/bin/env python3
"""Audit whether a same-51-island marine+candidate-host link survives source-verified
candidate-host richness and 2001–2010 general-plant recording-effort controls.

This is an exploratory sensitivity, not a causal mediation test. Source:
  Original uniform graph and biology: CI 42193019385 / artifact 13490383689.
  Original independent longitude fold plan: CI 42325176725 / artifact 13554313981.
  Original historical all-plant indexing and host richness: CI 42323799520 /
     artifact 13553828995.
No biological records or source island assignments are recomputed.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
import audit_island51_host_vs_uniform_marine_circuit_nested_20261117 as orig

EARLY_COVERAGE_SHA="be1e14536574dd99ebf70b9bd026a0624cf5aac103bea182b68181a3c8264bdb"
EARLY_RECEIPT_SHA="4cdfd41b85d3a17c07f0e85b13c3b5d8e24ccb8999880b22df1e9aca2e4d1ccb"

def main():
    p=argparse.ArgumentParser()
    for label in ("pairs","prep","longitude-blocks","early-coverage","early-receipt","output"):
        p.add_argument("--"+label,type=Path,required=True)
    a=p.parse_args()
    orig.verify_sha(a.pairs,orig.PIN_PAIRS)
    orig.verify_sha(a.prep,orig.PIN_PREP)
    orig.verify_sha(a.longitude_blocks,orig.PIN_BLOCK)
    orig.verify_sha(a.early_coverage,EARLY_COVERAGE_SHA)
    orig.verify_sha(a.early_receipt,EARLY_RECEIPT_SHA)
    q=pd.read_csv(a.pairs,dtype={"island_a":str,"island_b":str})
    original_prep=json.loads(a.prep.read_text())
    island_geo=pd.read_csv(a.longitude_blocks,dtype={"island_id":str})
    coverage=pd.read_csv(a.early_coverage,dtype={"island_id":str})
    receipt=json.loads(a.early_receipt.read_text())
    ids=sorted(original_prep["by_island"],key=int)
    if (len(ids),len(q),len(island_geo),len(coverage))!=(51,1275,51,51):
        raise ValueError("Source original 51-island identity or 1275 pair support changed")
    if set(island_geo.island_id)!=set(ids) or set(coverage.island_id)!=set(ids):
        raise ValueError("Original early botanical / geography island IDs disagree")
    if receipt["early_plant_record_scope"]["sum_port_buffer_early_plant_records"]!=155002:
        raise ValueError("Original historical plant recording effort population drifted")
    ports={k:int(v["original_ports"]) for k,v in original_prep["by_island"].items()}
    richness={k:int(v) for k,v in receipt["candidate_host_counts"].items()}
    records={str(x.island_id):int(x.early_plant_buffer_record_count_sum) for x in coverage.itertuples()}
    if sum(ports.values())!=194 or set(richness)!=set(ids) or set(records)!=set(ids) or sum(records.values())!=155002:
        raise ValueError("Original archived 194 port and 155002 plant effort source changed")
    def margins(m):
        v=np.log1p(q.island_a.map(m).to_numpy(float))
        w=np.log1p(q.island_b.map(m).to_numpy(float))
        return [v+w,np.abs(v-w)]
    d=np.log1p(q.port_distance_km.to_numpy(float))
    u=q.log1p_normalized_uniform_circuit.to_numpy(float)
    h=q.host_jaccard.to_numpy(float)
    classes=dict(zip(island_geo.island_id,island_geo.longitude_block))
    ga=q.island_a.map(classes).to_numpy(int)
    gb=q.island_b.map(classes).to_numpy(int)
    double_folds=[]
    for j in range(5):
        for k in range(j,5):
            la=((ga==j)|(ga==k))
            lb=((gb==j)|(gb==k))
            train=(~la)&(~lb)
            test=((ga==j)&(gb==k))|((ga==k)&(gb==j))
            double_folds.append((f"B{j}_B{k}",train,test))
    results={};foldrows=[];predictions=[]
    for target in ("new_butterfly_jaccard","earlier_butterfly_jaccard"):
        base=[d,d*d]
        if target=="new_butterfly_jaccard":
            base.append(q.earlier_butterfly_jaccard.to_numpy(float))
        # Original WPI port-count controls PLUS the only new, source-verified
        # source richness and early general botanical recording-effort margins.
        base.extend(margins(ports)+margins(richness)+margins(records))
        B=np.column_stack(base)
        designs={"B":B,"BU":np.column_stack((B,u)),
                 "BH":np.column_stack((B,h)),
                 "BUH":np.column_stack((B,u,h))}
        folds=[]
        for island in ids:
            test=(q.island_a.eq(island)|q.island_b.eq(island)).to_numpy()
            folds.append((island,~test,test))
        results[target]={}
        for mode,fs in (("LOIO",folds),("double_unseen",double_folds)):
            r,f,pred=orig.scores(q,target,designs,fs,mode)
            results[target][mode]=r
            foldrows.extend(f);predictions.extend(pred)
    after=results["new_butterfly_jaccard"]
    assert abs(after["LOIO"]["relative_MSE_reduction_pct"]["BH_to_BUH"]-0.7115161610139026)<1e-9
    assert abs(after["double_unseen"]["relative_MSE_reduction_pct"]["BH_to_BUH"]-(-0.5219475448315041))<1e-9
    assert abs(after["LOIO"]["relative_MSE_reduction_pct"]["B_to_BU"]-1.501679897068296)<1e-9
    assert abs(after["double_unseen"]["relative_MSE_reduction_pct"]["B_to_BU"]-0.7635727127138323)<1e-9
    a.output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(foldrows).to_csv(a.output/"adjusted_island51_model_fold_losses.csv",index=False)
    pd.DataFrame(predictions).to_csv(a.output/"adjusted_island51_predictions.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    report={
      "schema":"chocho_source_verified_51island_host_uniform_richness_and_pre2011_allplant_effort_20261117_v1",
      "status":"SOURCE_VERIFIED_POSTHOC_RICHNESS_EFFORT_FALSIFICATION_NOT_CAUSAL_MEDIATION",
      "source_original_uniform_CI":42193019385,
      "source_original_longitude_blocks_CI":42325176725,
      "source_original_historical_plant_effort_CI":42323799520,
      "input_SHA256":{"original_pairs":orig.PIN_PAIRS,"original_prep":orig.PIN_PREP,
                      "longitude_blocks":orig.PIN_BLOCK,"historical_plant_coverage":EARLY_COVERAGE_SHA,
                      "host_counts_and_effort_receipt":EARLY_RECEIPT_SHA},
      "n_islands":51,"n_original_pairs":1275,"total_2001_2010_plant_buffer_records":155002,
      "models":{
        "B":"original port-distance log+square, sum/difference original port counts, earlier butterfly Jaccard only for later response, sum/difference marginal candidate-host richness and original 2001–2010 general plant record counts",
        "BU":"B plus GSHHG undirected uniform effective resistance",
        "BH":"B plus original pooled 2001–2024 candidate introduced HOSTS Jaccard",
        "BUH":"B plus both uniform resistance and pooled candidate hosts"},
      "results":results,
      "limits":[
         "Historically indexed all-plant record COUNTS are observation proxies, not historical source-verified host plant identities",
         "Candidate host species composition pooled 2001–2024 may postdate 2011–2022 butterfly records",
         "This is sensitivity on same original 51 selected landmasses, not independent cohort",
         "Positive or negative MSE changes describe conditional predictions, never plant mediated causal dispersal",
         "Geographically blocked pairs share island identities, folds and outcomes, do not use independent-pair p-values",
         "Uniform circuit models sea topology not observed direction or vessel cargo",
         "Original WCVP 104 botanical-region alien flora Circuitscape response is DISTINCT from this butterfly response"
      ]
    }
    (a.output/"adjusted_51island_host_uniform_marine_source_verified.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("SOURCE_VERIFIED_ADJUSTED_51_ISLAND_HOST_UNIFORM",json.dumps({
       x:{y:r["relative_MSE_reduction_pct"] for y,r in z.items()} for x,z in results.items()
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
