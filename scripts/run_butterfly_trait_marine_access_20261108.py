#!/usr/bin/env python3
"""Butterflies first: source-verified trait x marine access and pairwise fauna analysis.

All outcomes are GBIF INDEX EVENTS / indexed fauna, not arrivals, breeding or
shipborne introductions. Wing-size proxy mixes LepTraits WS and FW bounds and
IS NOT a directly measured flight/migration ability. Host Jaccard is pooled
across years and not a historical larval-host availability indicator.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from itertools import combinations

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
BASE_SCRIPT=ROOT/"scripts/run_ncomms_public17_ocean_access_holdout.py"
TRAITS=("data/frozen/figure_sources/s1_resource_descriptors.csv",
        "09541b5894529cc8629c7515ea9022ded9edead2")
RESOURCES=("data/frozen/figure_sources/anthropogenic_species_metrics.csv",
           "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd")
MODELS={
  "B":[],
  "W_DIST":["wing_by_geodesic"],
  "D_DIST":["disturbance_by_geodesic"],
  "B_SEA":["min_log_uniform_effective_resistance"],
  "W_SEA":["min_log_uniform_effective_resistance","wing_by_sea"],
  "D_SEA":["min_log_uniform_effective_resistance","disturbance_by_sea"],
  "W_PORT":["wing_by_port"],
  "W_LAT":["wing_by_latitude"],
}
CONTRASTS=[("B","W_DIST"),("B","D_DIST"),
           ("B_SEA","W_SEA"),("B_SEA","D_SEA"),
           ("B","W_PORT"),("B","W_LAT")]
PERIODS=("2011_2017","2018_2022")
SPLITS=("island","component","species")
SMALL="small_50pct"; LARGE="large_50pct"
SEED=20261108

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def read_pinned_csv(path,pin):
    raw=(ROOT/path).read_bytes()
    if git_blob(raw)!=pin:
        raise ValueError("Historical frozen input changed: "+path)
    return list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))

def load():
    spec=importlib.util.spec_from_file_location("fixed_public17", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Original fixed public17 script not available")
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    d,pins=m.load(ROOT)
    return m,d,pins,read_pinned_csv(*TRAITS),read_pinned_csv(*RESOURCES)

def disturbance_class(s):
    v=str(s).strip().lower()
    if v.startswith("disturbance-associated"):
        return 1
    if v.startswith("disturbance-avoidant"):
        return -1
    if v.startswith("no evidence") or v.startswith("mixed"):
        return 0
    return None

def build_bridge(data,raw_traits,raw_resources):
    known=sorted({r["species"] for r in data["butterflies"]})
    traits={}
    duplicate=[]
    for r in raw_traits:
        s=r["species"]
        if s in traits: duplicate.append(s)
        traits[s]=r
    if duplicate:raise ValueError("Nonunique frozen S1 species")
    resources={r["species"]:r for r in raw_resources}
    bridge={}
    for s in known:
        t=traits.get(s)
        if t is None:
            bridge[s]={"species":s,"trait_available":False}
            continue
        try:
            wing=float(t["wing_size_proxy"])
        except (ValueError,TypeError):
            wing=float("nan")
        dis=disturbance_class(t.get("disturbance_affinity",""))
        if not math.isfinite(wing) or wing<=0:
            bridge[s]={"species":s,"trait_available":False,"reason":"no_numeric_wing"}
            continue
        res=resources.get(s)
        opportunity=float(res["log_resource_expansion"]) if res and res.get("log_resource_expansion") else None
        bridge[s]={"species":s,"trait_available":dis is not None,
                   "wing_size_proxy":wing,"disturbance_category":dis,
                   "disturbance_raw":t["disturbance_affinity"],
                   "host_family_count":t.get("host_family_count",""),
                   "global_log_resource_expansion":opportunity,
                   "global_resource_scope":"contemporary_pooled_NOT_historical_island_host"}
    present=[x for x in bridge.values() if x.get("trait_available")]
    if len(present)<100:
        raise ValueError("Trait bridge unexpectedly small; must not model sparse selection")
    wing_vals=np.array([z["wing_size_proxy"] for z in present])
    center=float(wing_vals.mean());std=float(wing_vals.std())
    med=float(np.median(wing_vals))
    if std<=0:raise ValueError("No wing variation")
    for z in present:
        z["wing_size_z"]=(z["wing_size_proxy"]-center)/std
        z["wing_median_group"]=SMALL if z["wing_size_proxy"]<=med else LARGE
        z["disturbance_associated"]=int(z["disturbance_category"]==1)
        z["disturbance_avoidant"]=int(z["disturbance_category"]==-1)
    summary={
        "source_candidate_species":len(known),
        "matched_complete_trait_species":len(present),
        "missing_or_unclassified_trait_species":len(known)-len(present),
        "wing_size_proxy_mean":center,
        "wing_size_proxy_sd":std,
        "wing_size_proxy_median":med,
        "wing_proxy_mixes_ws_fw_bounds_not_dispersal_rate":True,
        "trait_groups":dict(Counter(z["wing_median_group"] for z in present)),
        "disturbance_groups":dict(Counter(str(z["disturbance_category"]) for z in present)),
        "resource_descriptor_linked_species":sum(z.get("global_log_resource_expansion") is not None for z in present),
    }
    return bridge,summary

def traits_for_row(r,bridge):
    t=bridge[r["species"]]
    f=r["features"].copy()
    geo=f["min_log_geodesic_km"]
    sea=f["min_log_uniform_effective_resistance"]
    wing=t["wing_size_z"]
    disturbance=t["disturbance_associated"]-t["disturbance_avoidant"]
    f.update({
        "wing_size_z":wing,
        "disturbance_associated":t["disturbance_associated"],
        "disturbance_avoidant":t["disturbance_avoidant"],
        "wing_by_geodesic":wing*geo,
        "disturbance_by_geodesic":disturbance*geo,
        "wing_by_sea":wing*sea,
        "disturbance_by_sea":disturbance*sea,
        "wing_by_port":wing*f["log_island_ports"],
        "wing_by_latitude":wing*f["absolute_latitude"],
    })
    return f

def probability_loss(y,prob):
    p=np.clip(np.asarray(prob),1e-6,1-1e-6)
    return {"brier":float(np.mean((y-p)**2)),
            "logloss":float(-np.mean(y*np.log(p)+(1-y)*np.log1p(-p)))}

def cluster_ci(rows,y,p0,p1,seed,n=3000):
    islands=sorted(set(r["landmass_id"] for r in rows),key=int)
    target=np.asarray([r["landmass_id"] for r in rows])
    improvement=(y-p0)**2-(y-p1)**2
    sums=np.array([improvement[target==i].sum() for i in islands])
    counts=np.array([(target==i).sum() for i in islands])
    rng=np.random.default_rng(seed)
    idx=rng.integers(0,len(islands),size=(n,len(islands)))
    vals=sums[idx].sum(axis=1)/counts[idx].sum(axis=1)
    return {
        "gain_absolute":float(improvement.mean()),
        "island_cluster_fixed_pred_95pct":[float(x) for x in np.quantile(vals,[0.025,0.975])],
        "n_island_clusters":len(islands),
        "bootstrap_draws":n,"not_model_refitted":True,
    }

def model_event_traits(source,data,bridge,period):
    original,meta=source.prepare(data,period)
    original_source=[r for r in original if r["prior_source_islands"]>0]
    rows=[r for r in original_source if bridge[r["species"]].get("trait_available")]
    for r in rows:r["features"]=traits_for_row(r,bridge)
    if len(rows)<100:raise ValueError("No eligible source-index rows with traits")
    y=np.asarray([r["record_index_case"] for r in rows],dtype=int)
    n_pos=int(y.sum())
    base=source.FEATURES["B"]+source.FEATURES["D"]+[
        "wing_size_z","disturbance_associated","disturbance_avoidant"]
    feats={k:base+extra for k,extra in MODELS.items()}
    matrices={k:np.array([[r["features"][f] for f in fields] for r in rows],dtype=float)
              for k,fields in feats.items()}
    if any(not np.isfinite(m).all() for m in matrices.values()):
        raise ValueError("Nonfinite trait, location, or resource covariate")
    scores={}
    pred_rows=[]
    for split in SPLITS:
        keys=[r["landmass_id"] if split=="island" else r["component"] if split=="component"
              else str(int(r["speciesKey"])%5) for r in rows]
        uniq=sorted(set(keys),key=int)
        arr=np.array(keys)
        expected={"island":17,"component":6,"species":5}[split]
        if len(uniq)!=expected:raise ValueError(f"{period}/{split}: changed fold count")
        pred={name:np.full(len(rows),np.nan) for name in MODELS}
        for k in uniq:
            tr=arr!=k;te=arr==k
            if len(set(y[tr]))<2:
                raise ValueError(f"Empty positive training fold: {period}/{split}/{k}")
            for name,X in matrices.items():
                fit=make_pipeline(StandardScaler(),LogisticRegression(
                    C=1,solver="lbfgs",max_iter=1500,random_state=SEED))
                fit.fit(X[tr],y[tr])
                pred[name][te]=np.clip(fit.predict_proba(X[te])[:,1],1e-6,1-1e-6)
        if any(not np.isfinite(x).all() for x in pred.values()):
            raise ValueError("Incomplete OOF predictions")
        ls={k:probability_loss(y,v) for k,v in pred.items()}
        tests={}
        for p0,p1 in CONTRASTS:
            tests[p0+"_to_"+p1]={metric:{
                "baseline":ls[p0][metric],
                "augmented":ls[p1][metric],
                "reduction_pct":100*(ls[p0][metric]-ls[p1][metric])/ls[p0][metric],
            } for metric in ("brier","logloss")}
            if split=="island":
                tests[p0+"_to_"+p1]["conditional_island_resampling"]=cluster_ci(
                    rows,y,pred[p0],pred[p1],SEED+len(tests)+100*len(period))
        scores[split]={
            "rows":len(rows),"events":n_pos,"folds":len(uniq),
            "out_of_fold_losses":ls,"contrasts":tests,
        }
        for i,r in enumerate(rows):
            pred_rows.append({
                "period":period,"fold_type":split,
                "species":r["species"],"speciesKey":r["speciesKey"],
                "landmass_id":r["landmass_id"],
                "first_index_event":int(y[i]),
                "prior_source_islands":r["prior_source_islands"],
                "wing_size_proxy":bridge[r["species"]]["wing_size_proxy"],
                "wing_group":bridge[r["species"]]["wing_median_group"],
                "disturbance_class":bridge[r["species"]]["disturbance_category"],
                **{"pred_"+k:float(v[i]) for k,v in pred.items()},
            })
    return scores,pred_rows,{
      "period":period,"full_source_rows":len(original),"full_source_events":meta["first_index_cases"],
      "prior_source_eligible":len(original_source),
      "prior_source_events":sum(r["record_index_case"] for r in original_source),
      "trait_eligible":len(rows),"trait_eligible_events":n_pos,
      "trait_species":len(set(r["species"] for r in rows)),
    }

def jac(a,b):
    union=len(a|b)
    return float(len(a&b)/union) if union>0 else None

def pairwise_index_fauna(data,bridge):
    valid={s:z for s,z in bridge.items() if z.get("trait_available")}
    groups={
        SMALL:{s for s,z in valid.items() if z["wing_median_group"]==SMALL},
        LARGE:{s for s,z in valid.items() if z["wing_median_group"]==LARGE},
        "disturbance_associated":{s for s,z in valid.items() if z["disturbance_category"]==1},
        "disturbance_avoidant":{s for s,z in valid.items() if z["disturbance_category"]==-1},
    }
    ids=sorted({r["landmass_id"] for r in data["butterflies"]},key=int)
    if len(ids)!=17:raise ValueError("Not the original 17 islands")
    pairs={tuple(sorted([r["island_a"],r["island_b"]],key=int)):r for r in data["ocean"]}
    earlier=defaultdict(set);by2022=defaultdict(set)
    for r in data["butterflies"]:
        s,i=r["species"],r["landmass_id"]
        if s not in valid:continue
        old=int(r["pre2001_2010_indexed_records"])>0
        late=int(r["future2011_2017_indexed_records"])+int(r["future2018_2022_indexed_records"])>0
        if old:earlier[i].add(s)
        if old or late:by2022[i].add(s)
    records=[]
    for g,species in groups.items():
        if len(species)<8:
            continue
        for a,b in combinations(ids,2):
            p=pairs[tuple(sorted([a,b],key=int))]
            j0=jac(earlier[a]&species,earlier[b]&species)
            j1=jac(by2022[a]&species,by2022[b]&species)
            records.append({
                "stratum":g,"group_species":len(species),"island_a":a,"island_b":b,
                "early_J":j0,"by2022_J":j1,
                "host_J_pooled_NOT_historical":float(p["host_jaccard"]),
                "distance_km":float(p["distance_km_verified"]),
                "sea_effective_resistance":float(p["uniform_distributed_port_resistance"]),
            })
    return records,groups

def pairwise_two_island_unseen(records):
    results={}
    for group in sorted({r["stratum"] for r in records}):
        r=[x for x in records if x["stratum"]==group and x["early_J"] is not None and x["by2022_J"] is not None]
        if len(r)<60:
            results[group]={"eligible_pairs":len(r),"status":"INSUFFICIENT_PAIR_UNIONS"}
            continue
        y=np.array([z["by2022_J"] for z in r])
        X=np.array([[z["early_J"],math.log1p(z["distance_km"]),
                      z["host_J_pooled_NOT_historical"],math.log1p(z["sea_effective_resistance"])]
                    for z in r])
        if not np.isfinite(X).all():
            raise ValueError("Invalid index composition predictors")
        schemas={
            "B":[0,1],"B_H":[0,1,2],"B_S":[0,1,3],"B_HS":[0,1,2,3],
        }
        scores={}
        preds={}
        for model,columns in schemas.items():
            p=np.full(len(r),np.nan)
            for idx,z in enumerate(r):
                tr=np.array([v["island_a"] not in {z["island_a"],z["island_b"]}
                             and v["island_b"] not in {z["island_a"],z["island_b"]} for v in r])
                if int(tr.sum())<50:raise ValueError("Not enough disjoint training island-pairs")
                m=make_pipeline(StandardScaler(),Ridge(alpha=1.0))
                m.fit(X[tr][:,columns],y[tr])
                p[idx]=m.predict(X[idx:idx+1][:,columns])[0]
            preds[model]=p
            scores[model]=float(np.mean((y-p)**2))
        results[group]={
            "status":"SOURCE_INDEX_OBSERVATION_NOT_REAL_COLONISATION",
            "species_in_stratum":r[0]["group_species"],
            "original_island_pair_count":136,
            "eligible_pairs_with_source_union":len(r),
            "held_two_islands_out_training_pair_formula":"training pairs exclude both tested island IDs",
            "mse":scores,
            "gain_B_to_BH_pct":100*(scores["B"]-scores["B_H"])/scores["B"],
            "gain_B_to_BS_pct":100*(scores["B"]-scores["B_S"])/scores["B"],
            "gain_BH_to_BHS_pct":100*(scores["B_H"]-scores["B_HS"])/scores["B_H"],
            "target":"cumulative GBIF-indexed trait-specific species Jaccard by 2022",
            "host_similarity_is_all_candidate_hosts_NOT_trait_partitioned":True,
        }
    return results

def write_csv(path,rows):
    if not rows:return
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def run(out):
    src,data,pins,traits,resources=load()
    bridge,bridge_meta=build_bridge(data,traits,resources)
    model_scores={}
    predictions=[]
    coverage={}
    for period in PERIODS:
        s,p,m=model_event_traits(src,data,bridge,period)
        model_scores[period]=s;predictions+=p;coverage[period]=m
    pair_rows,groups=pairwise_index_fauna(data,bridge)
    pair_scores=pairwise_two_island_unseen(pair_rows)
    result={
        "schema":"chocho_butterfly_trait_access_public17_pilot_v1",
        "analysis_class":"EXPLORATORY_TRAIT_MODERATION_OF_INDEXED_BUTTERFLY_RECORDS",
        "source_git_blobs":pins,
        "trait_git_blob":TRAITS[1],
        "global_resource_git_blob":RESOURCES[1],
        "bridge":bridge_meta,"event_coverage":coverage,
        "event_oof":model_scores,"trait_partitioned_pairwise":pair_scores,
        "limitations":[
          "wing_size_proxy combines LepTraits WS and FW numerical bounds; not a direct measure of mobility",
          "2019 WPI port counts do not establish historical voyages",
          "modern pooled host Jaccard is not a dated species-specific island resource",
          "non-indexed butterflies are not confirmed ecological absences",
          "2018–2022 prior-source risk contains few positive index events",
          "17 islands selected for preexisting geographic source analysis, not all 51 islands",
          "no causal proof of shipborne plants, butterflies, inter-island relay or colonisation",
          "out-of-fold model comparisons are exploratory and multiple tests are unadjusted",
        ],
    }
    out.mkdir(parents=True,exist_ok=True)
    (out/"SOURCE_TRAIT_AUDIT.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    write_csv(out/"original225_species_traits_join.csv",list(bridge.values()))
    write_csv(out/"trait_event_oof_predictions.csv",predictions)
    write_csv(out/"trait_stratified_island_pair_jaccard.csv",pair_rows)
    print(json.dumps({
        "trait_coverage":bridge_meta,"risk_coverage":coverage,
        "event_2011_2017_held_island":{
             k:round(v["brier"]["reduction_pct"],4) for k,v in model_scores["2011_2017"]["island"]["contrasts"].items()},
        "event_2018_2022_held_island":{
             k:round(v["brier"]["reduction_pct"],4) for k,v in model_scores["2018_2022"]["island"]["contrasts"].items()},
        "pairs":{k:{"n":v.get("eligible_pairs"),
                    "host_gain":v.get("gain_B_to_BH_pct"),
                    "marine_gain":v.get("gain_BH_to_BHS_pct")} for k,v in pair_scores.items()}
    },indent=2,ensure_ascii=False))
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    run(args.out_dir)

if __name__=="__main__":
    main()
