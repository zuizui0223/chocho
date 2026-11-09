#!/usr/bin/env python3
"""Butterfly-first comparison on independently reindexed 225 taxa × 51 landmasses.

Scores describe GBIF's historical index recording outcome only; they DO NOT
measure true biological arrival, successful colonisation, breeding, larval
host use, measured dispersal ability, shipping or historical port presence.
The 51-island panel is a NEW public 2026 GBIF snapshot, not the original
restricted 83,573-row Parquet nor directly interchangeable with old 17 data.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/"scripts/run_ncomms_public17_ocean_access_holdout.py"
TRAIT=ROOT/"scripts/run_butterfly_trait_marine_access_20261108.py"
TABLE="full51_225butterfly_period_index.csv"
RECEIPT="full51_source_validation.json"
PERIODS=("2011_2017","2018_2022")
SPLITS=("island","geographic_cluster","species")
FIXED_SEED=20261109
FEATURES_BASE=[
 "log_target_before_window_record_effort",
 "log_prior_species_record_count_elsewhere",
 "log_prior_species_indexed_island_count",
 "log_island_area","log_island_ports",
 "latitude","absolute_latitude","longitude_sine","longitude_cosine",
 "log_early_general_plant_index_richness",
 "wing_size_z","disturbance_associated","disturbance_avoidant",
]
MODELS={
 "distance":FEATURES_BASE+["min_log_geodesic_km"],
 "distance_wn":FEATURES_BASE+["min_log_geodesic_km","wing_by_geodesic"],
 "distance_sea":FEATURES_BASE+["min_log_geodesic_km","min_log_uniform_effective_resistance"],
 "distance_sea_wn":FEATURES_BASE+["min_log_geodesic_km","min_log_uniform_effective_resistance","wing_by_sea"],
 "distance_wn_sea":FEATURES_BASE+["min_log_geodesic_km","wing_by_geodesic","min_log_uniform_effective_resistance"],
 "distance_wn_sea_wn":FEATURES_BASE+["min_log_geodesic_km","wing_by_geodesic","min_log_uniform_effective_resistance","wing_by_sea"],
 "distance_sea_disturb":FEATURES_BASE+["min_log_geodesic_km","min_log_uniform_effective_resistance","disturbance_by_sea"],
 "distance_sea_wport":FEATURES_BASE+["min_log_geodesic_km","min_log_uniform_effective_resistance","wing_by_port"],
}
CONTRASTS=[
 ("distance","distance_wn","wing_by_geodesic"),
 ("distance","distance_sea","sea_over_geodesic"),
 ("distance_sea","distance_sea_wn","wing_by_sea"),
 ("distance_wn_sea","distance_wn_sea_wn","wing_by_sea_after_wing_distance"),
 ("distance_sea","distance_sea_disturb","disturbance_by_sea"),
 ("distance_sea","distance_sea_wport","wing_by_port"),
]

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise RuntimeError("Source script unavailable")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def sha_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()

def load_source(source_dir):
    path=source_dir/TABLE
    receipt_path=source_dir/RECEIPT
    record=json.loads(receipt_path.read_text(encoding="utf8"))
    if record["schema"]!="chocho_publicGBIF_225butterfly_original51_source_merged_20261109_v1":
        raise ValueError("Wrong public source cohort")
    if record["merged_table_sha256"]!=sha_file(path):
        raise ValueError("Consolidated 51-island temporal source changed")
    with path.open(encoding="utf8",newline="") as f:
        rows=list(csv.DictReader(f))
    if (len(rows),len({z["species"] for z in rows}),len({z["landmass_id"] for z in rows}))!=(11475,225,51):
        raise ValueError("Not the original 225 × 51 indexed cohort")
    if len({(r["speciesKey"],r["landmass_id"]) for r in rows})!=len(rows):
        raise ValueError("Duplicate original-source risk ID")
    return rows,record

def data_bridge():
    source=load_module(PARENT,"original_public17_frozen_source")
    trait=load_module(TRAIT,"original_butterfly_trait")
    data,pins=source.load(ROOT)
    bridge,stats=trait.build_bridge(data,trait.read_pinned_csv(*trait.TRAITS),
                                    trait.read_pinned_csv(*trait.RESOURCES))
    return source,trait,data,pins,bridge,stats

def geographic_blocks(data,ids):
    geography={r["landmass_id"]:r for r in data["geography"]}
    if not ids<=set(geography):raise ValueError("Missing original island geography")
    v=[]
    for i in sorted(ids,key=int):
        lat=math.radians(float(geography[i]["representative_latitude"]))
        lon=math.radians(float(geography[i]["representative_longitude"]))
        v.append([math.cos(lat)*math.cos(lon),math.cos(lat)*math.sin(lon),math.sin(lat)])
    fit=KMeans(n_clusters=7,n_init=12,random_state=FIXED_SEED)
    clusters=fit.fit_predict(np.array(v))
    groups={i:int(k) for i,k in zip(sorted(ids,key=int),clusters)}
    if len(set(groups.values()))!=7:
        raise ValueError("Source-only seven geographic-cluster identities invalid")
    # A singleton distant group is allowed: it is an honest island-level test,
    # not a reason to reselect a different K after observing the outcomes.
    return groups

def make_rows(panel,period,data,bridge,clusters):
    geo={r["landmass_id"]:r for r in data["geography"]}
    early={r["island_id"]:r for r in data["early_general_botany"]}
    ports={r["landmass_id"]:r for r in data["wpi_islands"]}
    areas={r["landmass_id"]:r for r in data["global_area"]}
    ids={r["landmass_id"] for r in panel}
    if not ids<=set(geo)&set(early)&set(ports)&set(areas):
        raise ValueError("Original 51 geography, area, port or flora control missing")
    pairs={}
    for x in data["ocean"]:
        pair=tuple(sorted((x["island_a"],x["island_b"]),key=int))
        if pair in pairs:raise ValueError("Duplicate source sea-distance pair")
        pairs[pair]=(float(x["distance_km_verified"]),
                     float(x["uniform_distributed_port_resistance"]))
    if len(pairs)!=1275:raise ValueError("Source-locked 51-island resistance changed")
    source_by_species=defaultdict(list)
    prior_effort=Counter()
    for r in panel:
        old=int(r["pre2001_2010_indexed_records"])
        mid=int(r["future2011_2017_indexed_records"])
        prev=old+(mid if period=="2018_2022" else 0)
        if prev>0:source_by_species[r["speciesKey"]].append(r)
        prior_effort[r["landmass_id"]]+=prev
    cov=Counter()
    out=[]
    for r in panel:
        i=r["landmass_id"];s=r["species"];k=r["speciesKey"]
        old=int(r["pre2001_2010_indexed_records"])
        mid=int(r["future2011_2017_indexed_records"])
        later=int(r["future2018_2022_indexed_records"])
        eligible=(old==0) if period=="2011_2017" else (old==0 and mid==0)
        event=(mid>0) if period=="2011_2017" else (later>0)
        if not eligible:continue
        cov["index_eligible"]+=1
        cov["first_index_events"]+=int(event)
        z=bridge[s]
        if not z.get("trait_available"):continue
        cov["trait_eligible"]+=1
        cov["trait_eligible_events"]+=int(event)
        source_options=[q for q in source_by_species[k] if q["landmass_id"]!=i]
        if not source_options:continue
        distances=[pairs[tuple(sorted((i,q["landmass_id"]),key=int))][0] for q in source_options]
        resist=[pairs[tuple(sorted((i,q["landmass_id"]),key=int))][1] for q in source_options]
        if not all(d>0 for d in distances+resist):raise ValueError("Nonpositive original circuit inputs")
        prev=sum(int(q["pre2001_2010_indexed_records"])+(
                 int(q["future2011_2017_indexed_records"]) if period=="2018_2022" else 0)
                 for q in source_options)
        lat=float(geo[i]["representative_latitude"])
        lon=float(geo[i]["representative_longitude"])
        area=float(areas[i]["area_km2_equal_area"])
        n_ports=int(ports[i]["ports"])
        botanical_index=int(early[i]["unique_early_indexed_general_plant_speciesKeys"])
        if area<=0 or n_ports<=0 or botanical_index<0 or not (-90<=lat<=90 and -180<=lon<=180):
            raise ValueError("Invalid non-temporal static original island covariate")
        f={
          "log_target_before_window_record_effort":math.log1p(prior_effort[i]),
          "log_prior_species_record_count_elsewhere":math.log1p(prev),
          "log_prior_species_indexed_island_count":math.log1p(len(source_options)),
          "log_island_area":math.log(area),
          "log_island_ports":math.log1p(n_ports),
          "latitude":lat/90,"absolute_latitude":abs(lat)/90,
          "longitude_sine":math.sin(math.radians(lon)),
          "longitude_cosine":math.cos(math.radians(lon)),
          "log_early_general_plant_index_richness":math.log1p(botanical_index),
          "wing_size_z":float(z["wing_size_z"]),
          "disturbance_associated":int(z["disturbance_associated"]),
          "disturbance_avoidant":int(z["disturbance_avoidant"]),
          "min_log_geodesic_km":math.log1p(min(distances)),
          "min_log_uniform_effective_resistance":math.log1p(min(resist)),
        }
        f["wing_by_geodesic"]=f["wing_size_z"]*f["min_log_geodesic_km"]
        f["wing_by_sea"]=f["wing_size_z"]*f["min_log_uniform_effective_resistance"]
        f["disturbance_by_sea"]=(f["disturbance_associated"]-f["disturbance_avoidant"])*f["min_log_uniform_effective_resistance"]
        f["wing_by_port"]=f["wing_size_z"]*f["log_island_ports"]
        out.append({"species":s,"speciesKey":k,"landmass_id":i,
                    "record_index_case":int(event),"wing_size_z":float(z["wing_size_z"]),
                    "disturbance":int(z["disturbance_category"]),
                    "prior_source_count":len(source_options),
                    "geographic_cluster":clusters[i],"features":f})
    if len(out)<400 or sum(r["record_index_case"] for r in out)<10:
        raise ValueError("Insufficient observed source events to validate 51-island butterfly trait pilot")
    cov.update({"prior_source_and_trait_rows":len(out),
                "prior_source_and_trait_events":sum(r["record_index_case"] for r in out),
                "species_in_risk_set":len({r["species"] for r in out}),
                "islands_in_risk_set":len({r["landmass_id"] for r in out})})
    return out,dict(cov)

def loss(y,p):
    pp=np.clip(p,1e-6,1-1e-6)
    return {"brier":float(np.mean((y-pp)**2)),
            "logloss":float(-np.mean(y*np.log(pp)+(1-y)*np.log1p(-pp)))}

def fit_oof(rows,split):
    n=len(rows);y=np.asarray([r["record_index_case"] for r in rows])
    groups=np.array([
        r["landmass_id"] if split=="island" else str(r["geographic_cluster"]) if split=="geographic_cluster"
        else str(int(r["speciesKey"])%5) for r in rows
    ])
    uniq=sorted(set(groups),key=int)
    expected={"island":51,"geographic_cluster":7,"species":5}[split]
    # Some islands have no source history for these selected butterfly taxa,
    # so leave-island-out counts can be <51 even when all 51 input islands exist.
    if split=="island":
        if len(uniq)<35:raise ValueError("Too few island groups remain at risk")
    elif len(uniq)!=expected:raise ValueError("Fixed spatial/species holdout drift")
    mats={k:np.asarray([[r["features"][field] for field in fields] for r in rows],dtype=float)
          for k,fields in MODELS.items()}
    if any(not np.isfinite(a).all() for a in mats.values()):
        raise ValueError("Non-finite source trait or access feature")
    pred={k:np.full(n,np.nan) for k in MODELS}
    for group in uniq:
        mask=groups==group
        train=~mask
        if len(set(y[train]))<2 or not mask.any():
            raise ValueError("Fold not estimable with enough training index events")
        for name,X in mats.items():
            model=make_pipeline(StandardScaler(),LogisticRegression(
                C=1,solver="lbfgs",max_iter=1200,random_state=FIXED_SEED))
            model.fit(X[train],y[train])
            pred[name][mask]=np.clip(model.predict_proba(X[mask])[:,1],1e-6,1-1e-6)
    if any(not np.isfinite(a).all() for a in pred.values()):
        raise ValueError("Non-exhaustive out-of-fold prediction")
    metrics={k:loss(y,p) for k,p in pred.items()}
    tests={}
    for a,b,desc in CONTRASTS:
        tests[desc]={metric:{
            "baseline":metrics[a][metric],
            "extended":metrics[b][metric],
            "brier_or_logloss_reduction_pct":
                100*(metrics[a][metric]-metrics[b][metric])/metrics[a][metric]}
            for metric in ("brier","logloss")}
    pred_rows=[{
        "species":r["species"],"speciesKey":r["speciesKey"],
        "landmass_id":r["landmass_id"],"geographic_cluster":r["geographic_cluster"],
        "event":r["record_index_case"],"wing_size_z":r["wing_size_z"],
        "disturbance":r["disturbance"],"prior_source_count":r["prior_source_count"],
        **{"pred_"+name:float(p[i]) for name,p in pred.items()}
    } for i,r in enumerate(rows)]
    return {"source_cells":n,"first_index_events":int(y.sum()),
            "n_holdout_folds":len(uniq),"scores":metrics,"contrasts":tests},pred_rows

def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run(source_dir,out_dir):
    panel,receipt=load_source(source_dir)
    src,traits,data,pins,bridge,bmeta=data_bridge()
    ids={r["landmass_id"] for r in panel}
    clusters=geographic_blocks(data,ids)
    result={
      "schema":"chocho_public51_butterfly_wing_disturbance_access_v1",
      "analysis_class":"EXPLORATORY_RETROSPECTIVE_GBIF_RECORD_INDEX_TRAIT_MODERATION",
      "source_panel_SHA256":receipt["merged_table_sha256"],
      "frozen_17_control_source_Git_blobs":pins,
      "traits_source_Git_blob":traits.TRAITS[1],
      "original51_island_groups":len(ids),
      "geographic_blocking":"KMeans k=7 on original geography xyz unit-sphere coordinates, source-only, seed 20261109",
      "geographic_groups":clusters,
      "trait_coverage":bmeta,"periods":{},
      "limits":[
        "2026 reindex source not the private 83573-row historical author source",
        "51 landmass groups selected via port+host source eligibility, not representative of all islands",
        "Contemporary source-index nonobservations do not prove biological absence",
        "Nearest prior-indexed same-species source only among 51 sampled landmasses",
        "Modern WPI port counts are not verified historical trajectories",
        "Wing span/forewing proxy averages are not a validated flight or movement rate",
        "Uniform ocean resistance is not a historical route, and islands do not act as stopover nodes",
        "Early general plant richness captures indexing opportunity, not dated species-specific HOSTS availability",
        "No evidence of butterfly reproduction, plant/larva carriage, successful crossing or causal stepping stones"
      ]
    }
    all_preds=[]
    out_dir.mkdir(parents=True,exist_ok=True)
    for period in PERIODS:
        rows,meta=make_rows(panel,period,data,bridge,clusters)
        parts={}
        for split in SPLITS:
            stats,pp=fit_oof(rows,split)
            parts[split]=stats
            for z in pp:
                z["period"]=period;z["holdout"]=split
            all_preds.extend(pp)
        result["periods"][period]={"cohort":meta,"heldout":parts}
    audit=out_dir/"public51_trait_access_audit.json"
    audit.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    write_csv(out_dir/"public51_trait_oof_predictions.csv",all_preds)
    summary={
      "source_cells":11475,
      "periods":{p:{"cohort":v["cohort"],
         "heldout":{s:{k:{metric:x[metric]["brier_or_logloss_reduction_pct"] for metric in ("brier","logloss")}
                    for k,x in st["contrasts"].items()} for s,st in v["heldout"].items()}}
        for p,v in result["periods"].items()}}
    (out_dir/"public51_trait_summary.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    run(args.source_dir,args.out_dir)
