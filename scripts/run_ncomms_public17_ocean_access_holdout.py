#!/usr/bin/env python3
"""Unselected 225 butterflies × 17 GSHHG islands: prior records and sea-access pilot.

This is a retrospective 2026 GBIF index-replay, NOT actual colonisation,
breeding, shipment, historical host arrival, or NC-ready confirmation.
It tests if uniform-ocean effective resistance improves prediction *beyond*
prior same-butterfly source indexing, simple geodesic distance, geography,
area, port counts, and dated GBIF record effort.

Fixed pre-analysis controls and splits:
- index windows 2011–2017 and 2018–2022, each separate;
- no use of future-period records when constructing source distance or effort;
- B = background source history/geography/effort;
  D = B+closest recorded source great-circle port-to-port distance;
  C = B+closest recorded source uniform-ocean resistance;
  DC = B+both, tests only incremental circuit information beyond distance;
- 17 island-, 6 source-geographic-component-, 5 species-blocked fold tests;
- compare *out-of-fold* Brier and logloss, not significance of individual
  species×island cells and not the same OOF as old 51-island similarity models.

No historical plant/host predictor is used in this sea-only pilot.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pathlib
from collections import Counter, defaultdict
from datetime import datetime, timezone

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT=pathlib.Path(__file__).resolve().parents[1]
SOURCES={
 "butterflies":("data/exploratory/ncomms_public17_access_20261108/full17_butterfly_island_225_atrisk_source.csv",
                 "60e38e31d1d1b8896fb925aaa5ee0a405fbef2df"),
 "receipt":("data/exploratory/ncomms_public17_access_20261108/full17_query_receipt_source.json",
            "6ef70869578e1bde3bca13bb49e44ef3b4fb4360"),
 "ocean":("results/island51_observed_ocean_circuit_20261106/island51_observed_circuit_pairs.csv",
          "7832b1d9d8b1047027fd5197d7bc0d31c22cad3f"),
 "islands":("provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json",
            "54c4515947d2a3ee7545e4007985ef82a2af6868"),
 "geography":("results/island51_host_circuit_atlas_20261106/island_geography.csv",
              "f0c5b1139e54b4872d28682f0991f89f2368ab43"),
 "early_general_botany":("results/island51_pre2011_public_general_flora_20261107/source/GBIF_original51_early_general_flora_islands.csv",
                         "048fa6d644f5579879620deb394ce74ace653e26"),
 "wpi_islands":("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv",
                "18cecde5bf3066f276380167f451f1fa5cc1caa1"),
 "global_area":("data/exploratory/ncomms_public17_access_20261108/original_global_GSHHG_area20_frame_source.csv",
                "f0fabefd7b5e15d23c1ce7ad2127c8a03357f2a8"),
}
FEATURES={
 "B":["log_target_before_window_record_effort","log_prior_species_record_count_elsewhere",
      "log_prior_species_indexed_island_count","has_prior_indexed_source_island",
      "log_island_area","log_island_ports","latitude","absolute_latitude",
      "longitude_sine","longitude_cosine","log_early_general_plant_index_richness"],
 "D":["min_log_geodesic_km"],
 "C":["min_log_uniform_effective_resistance"],
}
TARGETS={
 "2011_2017":("eligible_2011_2017","first_indexed_in_2011_2017"),
 "2018_2022":("eligible_2018_2022","first_indexed_in_2018_2022"),
}
SPLITS=["island","component","species"]
MODELS={"B":FEATURES["B"],"D":FEATURES["B"]+FEATURES["D"],
        "C":FEATURES["B"]+FEATURES["C"],
        "DC":FEATURES["B"]+FEATURES["D"]+FEATURES["C"]}
CONTRASTS=[("B","D"),("B","C"),("D","DC"),("C","DC")]


def git_blob_sha(x:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(x)).encode("ascii")+b"\0"+x).hexdigest()


def load(root: pathlib.Path=ROOT)->tuple[dict,dict]:
    data,pins={},{}
    for key,(path,pin) in SOURCES.items():
        raw=(root/path).read_bytes()
        digest=git_blob_sha(raw)
        if digest!=pin:raise ValueError("Input Git blob changed: "+path+" ["+digest+"]")
        if path.endswith(".csv"):
            data[key]=list(csv.DictReader(raw.decode("utf8-sig").splitlines()))
        else:data[key]=json.loads(raw)
        pins[key]={"path":path,"git_blob_sha1":pin,"sha256":hashlib.sha256(raw).hexdigest()}
    return data,pins


def source_components(doc:dict)->dict[str,int]:
    groups={str(z["island_id"]) for z in doc["islands"]}
    if len(groups)!=17 or len(doc["source_pairs"])!=14:raise ValueError("Original 17 islands/14 near-pairs changed")
    p={x:x for x in groups}
    def root(x):
        if p[x]!=x:p[x]=root(p[x])
        return p[x]
    for row in doc["source_pairs"]:
        a,b=str(row["island_a"]),str(row["island_b"])
        if a not in groups or b not in groups:raise ValueError("Foreign island pair")
        ra,rb=root(a),root(b)
        if ra!=rb:p[rb]=ra
    dis=defaultdict(list)
    for x in groups:dis[root(x)].append(x)
    sets=sorted([tuple(sorted(a,key=int)) for a in dis.values()],key=lambda z:int(z[0]))
    if sorted(map(len,sets))!=[2,2,2,2,4,5]:
        raise ValueError("Unexpected original 17-island near-pair components")
    return {x:i for i,group in enumerate(sets) for x in group}


def prepare(data:dict, period:str)->tuple[list[dict],dict]:
    if period not in TARGETS:raise ValueError("Unexpected outcome period")
    risk=data["butterflies"]
    receipt=data["receipt"]
    if (receipt.get("status")!="PUBLIC_REINDEXED_OCCURRENCE_HISTORY_NOT_TRUE_COLONIZATION_OR_HOST_MEDIATION"
        or int(receipt["risk_cells"])!=3825):
        raise ValueError("Wrong public original GBIF 225×17 risk receipt")
    if len(risk)!=225*17:raise ValueError("Risk panel not exactly 225×17")
    keyed={(r["speciesKey"],r["landmass_id"]):r for r in risk}
    if len(keyed)!=len(risk):raise ValueError("Duplicate GBIF taxonKey×landmass")
    island_ids={r["landmass_id"] for r in risk}
    sp_ids={r["speciesKey"] for r in risk}
    if len(island_ids)!=17 or len(sp_ids)!=225:raise ValueError("Source population changed")
    if island_ids!={str(r["island_id"]) for r in data["islands"]["islands"]}:
        raise ValueError("Historical GSHHG source islands are inconsistent")
    # Parent table includes 51 islands and precisely one of each original pair.
    pair={}
    for r in data["ocean"]:
        key=tuple(sorted([str(r["island_a"]),str(r["island_b"])],key=int))
        if key in pair:raise ValueError("Repeated island pair")
        d=float(r["distance_km_verified"])
        resist=float(r["uniform_distributed_port_resistance"])
        if not(0<d<30000 and 0<resist<1e5):
            raise ValueError("Invalid source geom/ocean resistance")
        pair[key]=(d,resist)
    if len(pair)!=1275:raise ValueError("Original source 51-island circuit graph changed")
    for a in island_ids:
        for b in island_ids:
            if a!=b and tuple(sorted([a,b],key=int)) not in pair:
                raise ValueError("Missing source 17-island uniform resistance")
    geo={r["landmass_id"]:r for r in data["geography"]}
    early={r["island_id"]:r for r in data["early_general_botany"]}
    port={r["landmass_id"]:r for r in data["wpi_islands"]}
    areamap={r["landmass_id"]:r for r in data["global_area"]}
    if not island_ids.issubset(set(geo)&set(early)&set(port)&set(areamap)):
        raise ValueError("Missing original island geography/plant index/port controls")
    comp=source_components(data["islands"])
    raw_period=TARGETS[period]
    for r in risk:
        if r["source_not_future_selected"]!="True" or r["index_not_true_absence_or_colonization"]!="True":
            raise ValueError("Future-selection/inference lock was removed")
    source_by_species=defaultdict(list)
    target_effort=Counter()
    for r in risk:
        old=int(r["pre2001_2010_indexed_records"])
        mid=int(r["future2011_2017_indexed_records"])
        if min(old,mid,int(r["future2018_2022_indexed_records"]))<0:
            raise ValueError("Negative source indexed record count")
        if old>0 or (period=="2018_2022" and mid>0):
            source_by_species[r["speciesKey"]].append(r)
        target_effort[r["landmass_id"]]+=old+(mid if period=="2018_2022" else 0)
        if not ((int(r["eligible_2011_2017"]) == int(old==0)) and
                (int(r["first_indexed_in_2011_2017"]) == int(old==0 and mid>0))):
            raise ValueError("Source 2011–17 risk definition inconsistent")
        if (int(r["eligible_2018_2022"]) != int(old==0 and mid==0) or
            int(r["first_indexed_in_2018_2022"]) != int(old==0 and mid==0 and
                                                          int(r["future2018_2022_indexed_records"])>0)):
            raise ValueError("Source 2018–22 risk definition inconsistent")
    data_rows=[]
    for r in risk:
        if int(r[raw_period[0]])!=1:continue
        s=r["speciesKey"];id=r["landmass_id"]
        options=[x for x in source_by_species[s] if x["landmass_id"]!=id]
        ds=[pair[tuple(sorted([id,z["landmass_id"]],key=int))][0] for z in options]
        rs=[pair[tuple(sorted([id,z["landmass_id"]],key=int))][1] for z in options]
        count=sum(int(z["pre2001_2010_indexed_records"])+(
            int(z["future2011_2017_indexed_records"]) if period=="2018_2022" else 0)
            for z in options)
        lat=float(geo[id]["representative_latitude"])
        lon=float(geo[id]["representative_longitude"])
        area=float(areamap[id]["area_km2_equal_area"])
        island_ports=int(port[id]["ports"])
        plant_index=int(early[id]["unique_early_indexed_general_plant_speciesKeys"])
        if not(-90<=lat<=90 and -180<=lon<=180 and area>0 and island_ports>0 and plant_index>=0):
            raise ValueError("Source island geography/effort invalid")
        ff={
            "log_target_before_window_record_effort":math.log1p(target_effort[id]),
            "log_prior_species_record_count_elsewhere":math.log1p(count),
            "log_prior_species_indexed_island_count":math.log1p(len(options)),
            "has_prior_indexed_source_island":int(len(options)>0),
            "log_island_area":math.log(area),
            "log_island_ports":math.log1p(island_ports),
            "latitude":lat/90,
            "absolute_latitude":abs(lat)/90,
            "longitude_sine":math.sin(math.radians(lon)),
            "longitude_cosine":math.cos(math.radians(lon)),
            "log_early_general_plant_index_richness":math.log1p(plant_index),
            "min_log_geodesic_km":math.log1p(min(ds)) if ds else 0.,
            "min_log_uniform_effective_resistance":math.log1p(min(rs)) if rs else 0.,
        }
        data_rows.append({"species":r["species"],"speciesKey":s,"landmass_id":id,
                          "component":str(comp[id]),"record_index_case":int(r[raw_period[1]]),
                          "prior_source_islands":len(options),"features":ff})
    needed={tuple((z["speciesKey"],z["landmass_id"])) for z in data_rows}
    if len(needed)!=len(data_rows):raise ValueError("Duplicate evaluable risk observation")
    assert_num=(3584,80) if period=="2011_2017" else (3504,18)
    if (len(data_rows),sum(z["record_index_case"] for z in data_rows))!=assert_num:
        raise ValueError("Reindexed source case/denominator drift")
    return data_rows,{"period":period,"risk_rows":len(data_rows),"first_index_cases":assert_num[1],
                       "prior_source_rows":sum(z["prior_source_islands"]>0 for z in data_rows),
                       "prior_source_case_count":sum(z["record_index_case"] for z in data_rows
                                                     if z["prior_source_islands"]>0),
                       "source_unselected_butterflies":225,"source_original_islands":17,
                       "original_geographic_components":6}


def loss(y:np.ndarray,p:np.ndarray)->dict:
    p=np.clip(np.asarray(p,dtype=float),1e-6,1-1e-6)
    return {"brier":float(np.mean((y-p)**2)),
            "logloss":float(-np.mean(y*np.log(p)+(1-y)*np.log1p(-p)))}


def fit_cv(rows:list[dict],split:str)->tuple[list[dict],dict]:
    if split not in SPLITS:raise ValueError("Bad prespecified holdout kind")
    n=len(rows); y=np.array([r["record_index_case"] for r in rows],dtype=int)
    mtx={k:np.array([[r["features"][f] for f in cols] for r in rows],dtype=float)
         for k,cols in MODELS.items()}
    if any(not np.isfinite(v).all() for v in mtx.values()):
        raise ValueError("Nonfinite original source feature")
    if split=="island":keys=[r["landmass_id"] for r in rows]
    elif split=="component":keys=[r["component"] for r in rows]
    else:
        # Partition on frozen original source numeric taxonKey, not outcomes.
        keys=[str(int(r["speciesKey"])%5) for r in rows]
    unique=sorted(set(keys),key=lambda x:int(x))
    if len(unique)!={"island":17,"component":6,"species":5}[split]:
        raise ValueError("Unexpected blocked split size")
    pred={k:np.full(n,np.nan,dtype=float) for k in MODELS}
    seen=np.zeros(n,dtype=int)
    folds=[]
    key=np.asarray(keys,dtype=str)
    for group in unique:
        test=(key==group);train=~test
        if not (np.any(test) and np.any(train) and len(set(y[train]))==2):
            raise ValueError("One source holdout fold lacks training positives")
        seen[test]+=1
        for model,X in mtx.items():
            modeler=make_pipeline(StandardScaler(),LogisticRegression(
                C=1,solver="lbfgs",max_iter=2000,class_weight=None,random_state=20261108))
            modeler.fit(X[train],y[train])
            pred[model][test]=np.clip(modeler.predict_proba(X[test])[:,1],1e-6,1-1e-6)
        folds.append({"holdout_group":group,"test_cells":int(np.sum(test)),
                      "test_positive_index_cases":int(np.sum(y[test]))})
    if not np.all(seen==1) or any(not np.isfinite(p).all() for p in pred.values()):
        raise ValueError("Predictions are incomplete or holdout cells repeated")
    outcomes={m:loss(y,p) for m,p in pred.items()}
    contrasts={}
    for ref,add in CONTRASTS:
        k=ref+"_to_"+add
        contrasts[k]={l:{"reference":outcomes[ref][l],
                         "augmented":outcomes[add][l],
                         "percentage_reduction":100*(outcomes[ref][l]-outcomes[add][l])/outcomes[ref][l]}
                      for l in ("brier","logloss")}
    ids={(z["landmass_id"],z["speciesKey"]) for z in rows}
    if len(ids)!=n:raise ValueError("Holdout source pairs duplicated")
    oof=[{"species":r["species"],"speciesKey":r["speciesKey"],
          "landmass_id":r["landmass_id"],"component":r["component"],
          "first_index_case":r["record_index_case"],
          "prior_source_islands":r["prior_source_islands"],
          "pred_B":pred["B"][i],"pred_D":pred["D"][i],
          "pred_C":pred["C"][i],"pred_DC":pred["DC"][i]}
         for i,r in enumerate(rows)]
    return oof,{"split":split,"n_folds":len(folds),"n_test_cells":n,
                "n_positive_index_cases":int(y.sum()),"fold_source_audit":folds,
                "out_of_fold_scores":outcomes,"prediction_contrasts":contrasts,
                "nothing_is_fit_on_test_labels":True,
                "not_true_absence_or_establishment":True}


def run(output:pathlib.Path,root:pathlib.Path=ROOT)->dict:
    data,pins=load(root)
    result={
      "schema":"chocho_NC_17island_225butterfly_public_record_index_heldout_access_v1",
      "status":"INDEPENDENT_PUBLIC_GBIF_REINDEX_PREDICTION_NOT_CAUSAL",
      "source_pins":pins,
      "fixed_candidate_butterflies":225,
      "original_landmass_groups":17,
      "source_population":"2026-10-31 reindexed GBIF human-observation query with 2001–2010, 2011–2017, 2018–2022 windows",
      "source_access":"51 original source-defined uniform-ocean Circuitscape pairwise resistance",
      "covariates":["pre-window recorded butterfly effort by island",
                   "earlier same-species indexed source count and recorded source records",
                   "island area and indexed modern WPI port count",
                   "lat/lon proxy", "2001–2010 indexed general flora richness",
                   "nearest earlier indexed source geodesic distance",
                   "nearest earlier indexed source uniform-ocean effective resistance"],
      "hypotheses":"Does uniform-ocean resistance add beyond simple geographic distance and pre-window source/indexing controls?",
      "not_tested":["actual butterfly colonisation","historical shipping",
                    "species-specific pre2011 larval host causation",
                    "host x sea resistance interaction",
                    "climate-adjusted causal marine access effect"],
      "model_spec":"StandardScaler(train-only) + LogisticRegression L2 C=1, fixed features and fixed island/component/species folds",
      "periods":{}
    }
    output.mkdir(parents=True,exist_ok=True)
    for period in TARGETS:
        rows,meta=prepare(data,period)
        split_results={}
        for split in SPLITS:
            oof,score=fit_cv(rows,split)
            split_results[split]=score
            with (output/(period+"_"+split+"_heldout_predictions.csv")).open("w",newline="",encoding="utf8") as f:
                writer=csv.DictWriter(f,fieldnames=list(oof[0]))
                writer.writeheader();writer.writerows(oof)
        result["periods"][period]={"source":meta,"heldout":split_results}
    before=output/"source_prediction_receipt.json"
    receipt=json.dumps(result,indent=2,ensure_ascii=False)+"\n"
    before.write_text(receipt,encoding="utf8")
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--out-dir",type=pathlib.Path,required=True)
    args=cli.parse_args()
    r=run(args.out_dir)
    print(json.dumps({"status":r["status"],
           "periods":{t:{"data":q["source"],
             "contrasts":{k:{contrast:vals["brier"]["percentage_reduction"] for contrast,vals
                                in z["prediction_contrasts"].items()}
                          for k,z in q["heldout"].items()}}
                      for t,q in r["periods"].items()}},indent=2))
if __name__=="__main__":
    main()
