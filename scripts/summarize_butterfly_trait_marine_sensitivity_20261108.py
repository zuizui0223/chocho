#!/usr/bin/env python3
"""Source-exact checks: trait interaction direction, OOF influence, matched islands.

Reads prior stage's public-source outputs. LOO diagnostics NEVER refit models.
Matched island-pair sample re-fits ridge only on TRAINING pairs whose two island
IDs differ from both held-out islands. No private raw Parquet is read.
"""
from __future__ import annotations
import argparse,csv,json,importlib.util
from collections import Counter
from pathlib import Path
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge,LogisticRegression

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT/"scripts/run_butterfly_trait_marine_access_20261108.py"

def get_original():
    spec=importlib.util.spec_from_file_location("trait_original",ORIGINAL)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def read_rows(path):
    with path.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def eval_gain(rows,base,new):
    if not rows:return None
    y=np.array([int(z["first_index_event"]) for z in rows])
    a=np.array([float(z["pred_"+base]) for z in rows])
    b=np.array([float(z["pred_"+new]) for z in rows])
    loss0=np.mean((y-a)**2);loss1=np.mean((y-b)**2)
    return float(100*(loss0-loss1)/loss0)

def leave_one_summary(rows,base,new):
    pos=[r for r in rows if r["first_index_event"]=="1"]
    remove_species=sorted({x["species"] for x in pos})
    remove_islands=sorted({x["landmass_id"] for x in rows},key=int)
    by_sp={s:eval_gain([r for r in rows if r["species"]!=s],base,new) for s in remove_species}
    by_is={i:eval_gain([r for r in rows if r["landmass_id"]!=i],base,new) for i in remove_islands}
    by_evt={r["species"]+"/"+r["landmass_id"]:
             eval_gain([z for z in rows if not(z["species"]==r["species"] and
                                               z["landmass_id"]==r["landmass_id"])],base,new)
            for r in pos}
    return {
      "n":len(rows),"events":len(pos),"initial_gain_pct":eval_gain(rows,base,new),
      "source_positive_species_count":len(by_sp),
      "positive_species_deleted_gain_range_pct":[min(by_sp.values()),max(by_sp.values())],
      "positive_species_deleted_all_positive":all(q>0 for q in by_sp.values()),
      "island_deleted_gain_range_pct":[min(by_is.values()),max(by_is.values())],
      "island_deleted_all_positive":all(q>0 for q in by_is.values()),
      "positive_cell_deleted_gain_range_pct":[min(by_evt.values()),max(by_evt.values())],
      "positive_cell_deleted_all_positive":all(q>0 for q in by_evt.values()),
      "NO_MODEL_REFITTING":True,
    }

def matched_pair_tests(rows):
    found={g:[] for g in sorted({r["stratum"] for r in rows})}
    for r in rows:
        if r["early_J"]!="" and r["by2022_J"]!="":
            found[r["stratum"]].append(r)
    common=set.intersection(*({(r["island_a"],r["island_b"]) for r in grp} for grp in found.values()))
    if len(common)<60:raise ValueError("Insufficient common original island pairs")
    out={}
    for group,rr in found.items():
        same=sorted([r for r in rr if (r["island_a"],r["island_b"]) in common],
                    key=lambda r:(int(r["island_a"]),int(r["island_b"])))
        y=np.array([float(r["by2022_J"]) for r in same])
        x=np.array([[float(r["early_J"]),np.log1p(float(r["distance_km"])),
            float(r["host_J_pooled_NOT_historical"]),
            np.log1p(float(r["sea_effective_resistance"]))] for r in same])
        preds={}
        for name,cols in {"B":[0,1],"BH":[0,1,2],"BS":[0,1,3],"BHS":[0,1,2,3]}.items():
            p=np.zeros(len(same))
            for j,r in enumerate(same):
                ids={r["island_a"],r["island_b"]}
                tr=np.array([z["island_a"] not in ids and z["island_b"] not in ids for z in same])
                if tr.sum()<30:raise ValueError("Too few disjoint pair training rows")
                fit=make_pipeline(StandardScaler(),Ridge(alpha=1.0))
                fit.fit(x[tr][:,cols],y[tr])
                p[j]=fit.predict(x[j:j+1][:,cols])[0]
            preds[name]=float(np.mean((y-p)**2))
        out[group]={
          "same_original_island_pairs":len(same),
          "baseline_MSE":preds["B"],
          "host_increment_MSE_reduction_pct":100*(preds["B"]-preds["BH"])/preds["B"],
          "sea_increment_over_distance_MSE_reduction_pct":100*(preds["B"]-preds["BS"])/preds["B"],
          "sea_increment_over_distance_and_pooled_hosts_MSE_reduction_pct":
              100*(preds["BH"]-preds["BHS"])/preds["BH"],
          "host_similarity_not_trait_specific_or_historical":True,
        }
    return {"matched_pair_count":len(common),"strata":out}

def directional_fold_signs(m):
    source,data,pins,traits,resources=m.load()
    bridge,_=m.build_bridge(data,traits,resources)
    parent=importlib.util.spec_from_file_location(
        "public17_original",ROOT/"scripts/run_ncomms_public17_ocean_access_holdout.py")
    p=importlib.util.module_from_spec(parent);parent.loader.exec_module(p)
    results={}
    for period in m.PERIODS:
        rows,_=p.prepare(data,period)
        rows=[r for r in rows if r["prior_source_islands"]>0 and
              bridge[r["species"]].get("trait_available")]
        for r in rows:r["features"]=m.traits_for_row(r,bridge)
        base=p.FEATURES["B"]+p.FEATURES["D"]+[
          "wing_size_z","disturbance_associated","disturbance_avoidant"]
        y=np.array([r["record_index_case"] for r in rows])
        by_split={}
        for split in m.SPLITS:
            keys=[r["landmass_id"] if split=="island" else r["component"] if split=="component"
                  else str(int(r["speciesKey"])%5) for r in rows]
            k=np.asarray(keys)
            by_model={}
            for model,extra,interaction in [
                ("W_SEA",m.MODELS["W_SEA"],"wing_by_sea"),
                ("W_DIST",m.MODELS["W_DIST"],"wing_by_geodesic"),
                ("W_PORT",m.MODELS["W_PORT"],"wing_by_port")]:
                fields=base+extra
                X=np.array([[r["features"][f] for f in fields] for r in rows])
                coefs=[]
                for grp in sorted(set(k)):
                    train=k!=grp
                    clf=make_pipeline(StandardScaler(),LogisticRegression(
                        C=1,solver="lbfgs",max_iter=1500,random_state=m.SEED))
                    clf.fit(X[train],y[train])
                    coefs.append(float(clf.steps[-1][1].coef_[0][fields.index(interaction)]))
                by_model[model]={
                  "interaction":interaction,
                  "coef_train_standardized_min":float(min(coefs)),
                  "coef_train_standardized_max":float(max(coefs)),
                  "coef_train_standardized_mean":float(np.mean(coefs)),
                  "positive_fold_coefficient_count":sum(z>0 for z in coefs),
                  "fold_count":len(coefs),
                  "coefficient_is_conditional_not_migration_rate":True,
                }
            by_split[split]=by_model
        results[period]=by_split
    return results

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-dir",type=Path,required=True)
    args=ap.parse_args()
    s=args.source_dir
    records=read_rows(s/"trait_event_oof_predictions.csv")
    pair=read_rows(s/"trait_stratified_island_pair_jaccard.csv")
    m=get_original()
    output={
      "schema":"chocho_trait_marine_posthoc_influence_and_matchedpair_20261108",
      "source":"Original frozen 225 species x 17 islands, frozen S1 traits, GitHub Actions generated predictions",
      "warnings":["post hoc sensitivity not prior hypothesis confirmation",
                  "leave one out uses frozen OOF predictions, no model refitting",
                  "small 2018–2022 positive events",
                  "same-pair host similarity pooled across all hosts, not prior species-specific plants",
                  "not biology of verified arrival or shipping"],
      "leave_one":{},"matched_pairs":matched_pair_tests(pair),
      "directional_coefficients":directional_fold_signs(m),
    }
    for period in m.PERIODS:
        for holdout in m.SPLITS:
            sub=[r for r in records if r["period"]==period and r["fold_type"]==holdout]
            for name,ref,added in [
               ("wing_by_sea","B_SEA","W_SEA"),
               ("wing_by_geodesic","B","W_DIST"),
               ("wing_by_port","B","W_PORT")]:
                output["leave_one"]["_".join([period,holdout,name])]=leave_one_summary(sub,ref,added)
    target=s/"TRAIT_SENSITIVITY_AND_MATCHED_PAIR_AUDIT.json"
    target.write_text(json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({"same_pairs":output["matched_pairs"]["matched_pair_count"],
      "late_island":output["leave_one"]["2018_2022_island_wing_by_sea"],
      "early_island":output["leave_one"]["2011_2017_island_wing_by_sea"],
      "coef_signs":output["directional_coefficients"]["2018_2022"]["island"]},
      indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
