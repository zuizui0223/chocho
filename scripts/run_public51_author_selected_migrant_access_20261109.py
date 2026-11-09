#!/usr/bin/env python3
"""Does observed island access transfer within independently author-listed migratory butterflies?

The all-225 proxy is author-GBIF-request-selected from an indirectly linked
migration review, NOT its original Table S1. Unlisted=UNKNOWN, not sedentary.
This is retrospective GBIF first-INDEX, not verified arrival/dispersal/colonization.
Model only the LISTED subset as a within-positive-evidence sensitivity; report
UNKNOWN group only as coverage, never as 'nonmigrant' negative control.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/"scripts/run_public51_butterfly_trait_access_20261109.py"
SOURCE_CROSSWALK=ROOT/"data/exploratory/public51_author_migrant_20261109/migration_gbif_author2025_crosswalk.csv"
CROSSWALK_BLOB="439dd1034ca16f96ff94222695d98dc39747c2e4"
KNOWN="REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"
UNKNOWN="NOT_LISTED_OR_TAXON_UNRESOLVED_UNKNOWN"
SEED=20261109
PERIODS=("2011_2017","2018_2022")
SPLITS=("island","geographic_cluster","species")
COVARIATES=[
 "log_target_before_window_record_effort",
 "log_prior_species_record_count_elsewhere",
 "log_prior_species_indexed_island_count",
 "log_island_area","log_island_ports",
 "latitude","absolute_latitude","longitude_sine","longitude_cosine",
 "log_early_general_plant_index_richness",
 "wing_size_z","disturbance_associated","disturbance_avoidant",
]
MODELS={
 "geo":COVARIATES+["min_log_geodesic_km"],
 "geo_wing":COVARIATES+["min_log_geodesic_km","wing_by_geodesic"],
 "geo_sea":COVARIATES+["min_log_geodesic_km","min_log_uniform_effective_resistance"],
 "geo_sea_wing":COVARIATES+["min_log_geodesic_km","min_log_uniform_effective_resistance","wing_by_sea"],
 "geo_wing_sea":COVARIATES+["min_log_geodesic_km","wing_by_geodesic","min_log_uniform_effective_resistance"],
 "geo_wing_sea_wing":COVARIATES+["min_log_geodesic_km","wing_by_geodesic","min_log_uniform_effective_resistance","wing_by_sea"],
}
COMPARISONS=[
 ("geo","geo_wing","wing_x_geodesic"),
 ("geo","geo_sea","sea_beyond_geodesic"),
 ("geo_sea","geo_sea_wing","wing_x_sea_given_distance"),
 ("geo_wing_sea","geo_wing_sea_wing","wing_x_sea_after_wing_x_distance"),
]

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise RuntimeError("Missing original source")
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj)
    return obj

def git_blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def migration_crosswalk(panel):
    if git_blob(SOURCE_CROSSWALK)!=CROSSWALK_BLOB:
        raise ValueError("Independent author-selected migration proxy changed")
    with SOURCE_CROSSWALK.open(encoding="utf8",newline="") as f:
        data=list(csv.DictReader(f))
    if len(data)!=225:raise ValueError("Expected exactly 225 original butterfly statuses")
    keys={(z["speciesKey"],z["species"]) for z in panel}
    if len(keys)!=225 or len({z["speciesKey"] for z in data})!=225:
        raise ValueError("Original butterfly taxon-key resolution changed")
    if {(z["speciesKey"],z["species"]) for z in data}!=keys:
        raise ValueError("Author-selected butterfly taxonomy differs from source 225")
    labels={x["speciesKey"]:x["author_2025_GBIF_requested_taxon_key_proxy"] for x in data}
    if set(labels.values())!={KNOWN,UNKNOWN}:
        raise ValueError("Migration roster must be documented vs unknown only")
    if Counter(labels.values())!={KNOWN:111,UNKNOWN:114}:
        raise ValueError("Unexpected author-GBIF-request status coverage")
    return labels

def brier_logloss(y,p):
    pp=np.clip(p,1e-6,1-1e-6)
    return {"brier":float(np.mean((y-pp)**2)),
            "logloss":float(-np.mean(y*np.log(pp)+(1-y)*np.log1p(-pp)))}

def fold_scores(rows):
    y=np.asarray([int(x["record_index_case"]) for x in rows],int)
    if int(y.sum())<12:raise ValueError("Too few independent migration proxy candidate record indices")
    matrices={
        name:np.asarray([[r["features"][k] for k in cols] for r in rows],float)
        for name,cols in MODELS.items()
    }
    if not all(np.isfinite(m).all() for m in matrices.values()):
        raise ValueError("Nonfinite source covariates")
    results={}
    stored=[]
    for split in SPLITS:
        groups=np.array([
            z["landmass_id"] if split=="island" else
            str(z["geographic_cluster"]) if split=="geographic_cluster"
            else str(int(z["speciesKey"])%5)
            for z in rows])
        unique=sorted(set(groups),key=int)
        if split=="geographic_cluster" and len(unique)!=7:raise ValueError("Unreproducible seven-region-source splits")
        if split=="species" and len(unique)!=5:raise ValueError("Unreproducible species-blocks")
        if split=="island" and len(unique)<35:raise ValueError("Too few candidate source islands")
        pred={k:np.full(len(rows),np.nan,float) for k in MODELS}
        for key in unique:
            te=groups==key;tr=~te
            if len(set(y[tr]))<2:raise ValueError("A source heldout fold lacks a positive training record")
            for name,X in matrices.items():
                pipeline=make_pipeline(StandardScaler(),
                    LogisticRegression(C=1,solver="lbfgs",random_state=SEED,max_iter=1600))
                pipeline.fit(X[tr],y[tr])
                pred[name][te]=np.clip(pipeline.predict_proba(X[te])[:,1],1e-6,1-1e-6)
        if any(not np.isfinite(p).all() for p in pred.values()):
            raise ValueError("Unscored candidate migration proxy heldout source cell")
        metrics={name:brier_logloss(y,pp) for name,pp in pred.items()}
        compares={}
        for base,augment,label in COMPARISONS:
            compares[label]={sc:{"baseline":metrics[base][sc],
                "extended":metrics[augment][sc],
                "percent_reduction":float(100*(metrics[base][sc]-metrics[augment][sc])/metrics[base][sc])}
                for sc in ("brier","logloss")}
        results[split]={"n_folds":len(unique),"risk_cells":len(rows),"first_indexed_cases":int(y.sum()),
                        "model_losses":metrics,"comparisons":compares}
        for idx,r in enumerate(rows):
            stored.append({"period":r["period"],"holdout":split,"species":r["species"],
                "speciesKey":r["speciesKey"],"landmass_id":r["landmass_id"],
                "geographic_cluster":r["geographic_cluster"],
                "event":int(y[idx]),"wing_size_z":r["wing_size_z"],
                "disturbance":r["disturbance"],
                "source_documented_migrant_proxy":"REVIEW_SELECTED_GBIF_REQUEST",
                "prior_conspecific_island_sources":r["prior_source_count"],
                **{"pred_"+name:float(vals[idx]) for name,vals in pred.items()}})
    return results,stored

def write_csv(path,rows):
    if not rows:raise ValueError("No evidence table")
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run(out):
    parent=module(PARENT,"source51_butterfly")
    panel,receipt=parent.load_source(
        ROOT/"data/exploratory/public51_butterfly_index_20261109")
    statuses=migration_crosswalk(panel)
    src,trait,model_data,pins,traits,bmeta=parent.data_bridge()
    blocks=parent.geographic_blocks(model_data,{r["landmass_id"] for r in panel})
    summary={
      "schema":"chocho_225x51_indirect_documented_migrant_subset_access_oof_20261109_v1",
      "status":"EXPLORATORY_WITHIN_REVIEW_SELECTED_MIGRANT_PROXY_NOT_NONMIGRANT_COMPARISON",
      "source_sha256":receipt["merged_table_sha256"],
      "source_author_migrant_proxy_git_blob":CROSSWALK_BLOB,
      "proxy_originally_prepared_from":"2025 author GBIF request derived indirectly from 2021 migration review; NOT original supplemental Table S1",
      "reconstructed_225_species_roster":{"author_selected_candidates":111,"unlisted_migration_status_unknown":114},
      "LepTraits_complete_species":bmeta["matched_complete_trait_species"],
      "source_selection":"previously indexed conspecific elsewhere among original 51; original source risk and complete LepTraits only",
      "models":"Fixed L2 logistic, source effort, port, geography, wing, disturbance and original Circuitscape; island/geographic-area/species out-of-fold",
      "cohorts":{},"heldout":{},
      "limits":[
       "The selected category is a positive author-selected taxon proxy, not a validated independent observation of the target-island butterflies migrating.",
       "The NOT_LISTED category is UNKNOWN, not a sedentary or non-migrant control.",
       "The source-selected migration roster likely depends on reporting and species detectability.",
       "All GBIF first indexing events are retrospective observations and not confirmed biological arrivals or establishment.",
       "Wingspan/forewing composite proxy is not flight performance.",
       "Circuitscape uniform-sea resistance is not shipping density or true movement route.",
       "No documented historical larval host at each target island or observed ship transport.",
       "Only original 51 port-eligible landmasses are observed; unseen sources elsewhere could be closer.",
       "Source-only 7-region KMeans and species modulo-5 folds test geographic/source transport of a retrospective model, not out-of-time forecast.",
      ],
    }
    all_preds=[]
    for period in PERIODS:
        panelrows,meta=parent.make_rows(panel,period,model_data,traits,blocks)
        for z in panelrows:z["period"]=period
        by_status=defaultdict(list)
        for z in panelrows:by_status[statuses[z["speciesKey"]]].append(z)
        subset=by_status[KNOWN]
        unknown=by_status[UNKNOWN]
        if not subset:raise ValueError("No independently documented author-reviewed candidate migrant rows")
        for label,rows in ((KNOWN,subset),(UNKNOWN,unknown)):
            summary["cohorts"].setdefault(period,{})[label]={
              "at_risk_rows":len(rows),
              "first_index_cases":sum(z["record_index_case"] for z in rows),
              "butterfly_species":len({z["speciesKey"] for z in rows}),
              "landmass_groups":len({z["landmass_id"] for z in rows}),
              "not_nonmigrant_control":label==UNKNOWN,
            }
        fitted,pp=fold_scores(subset)
        summary["heldout"][period]=fitted
        all_preds.extend(pp)
    out.mkdir(parents=True,exist_ok=True)
    (out/"migrant_author_selected_access_audit.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
    write_csv(out/"migrant_author_selected_OOF_predictions.csv",all_preds)
    print(json.dumps({"status":summary["status"],"source":summary["source_sha256"],
        "cohorts":summary["cohorts"],"improvements":{
          period:{split:{comp:{m:round(v[m]["percent_reduction"],5) for m in ("brier","logloss")}
                     for comp,v in q["comparisons"].items()}
                     for split,q in spl.items()} for period,spl in summary["heldout"].items()
         }},ensure_ascii=False,indent=2))
    return summary

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    args=parser.parse_args()
    run(args.out_dir)
