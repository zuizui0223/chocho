#!/usr/bin/env python3
"""Within-address repeated monarch milkweed surveys: held-out host-identity predictions.

Exploratory; archived DOI 10.6084/m9.figshare.25648644.v1.
A positive prediction is not a measured cohort survival effect.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss,log_loss
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder,StandardScaler

SOURCE_SPECIES=["Asclepias curassavica","Asclepias fascicularis",
                "Asclepias speciosa","Gomphocarpus physocarpus"]
RAW_COLS=["Eggs","Instar_4","Instar_5","Num_plants"]
NUM_COLS=["log_eggs","log_prior_late","log_plants0","log_plants1",
          "days_gap","month_sin","month_cos"]
ORIGIN_COL=["Native_status"]
SPECIES_COL=["Milkweed_sp"]

def load(path:Path):
    d=pd.read_csv(path,low_memory=False)
    d["day"]=pd.to_datetime(d["Date"],format="%m/%d/%y",errors="coerce")
    for c in RAW_COLS:d[c]=pd.to_numeric(d[c],errors="coerce")
    valid=d["day"].notna()&d[RAW_COLS].notna().all(axis=1)&(d["Num_plants"]>0)
    valid=valid&(d[RAW_COLS]>=0).all(axis=1)
    d=d.loc[valid].copy()
    if len(d)!=4480:raise RuntimeError(f"frozen valid patch row count drift {len(d)} != 4480")
    d["late"]=d["Instar_4"]+d["Instar_5"]
    g=(d.groupby(["Route","Address_comp","Milkweed_sp","day"],dropna=False)
       .agg(eggs=("Eggs","sum"),late=("late","sum"),
            plants=("Num_plants","sum"),origin_n=("Native_status","nunique"),
            Native_status=("Native_status","first"),
            source_rows=("Eggs","size")).reset_index())
    if (g.origin_n!=1).any():raise RuntimeError("same patch/date source has conflicting native origin")
    transitions=[]
    for (route,address,species),site in g.groupby(["Route","Address_comp","Milkweed_sp"],dropna=False):
        visits=site.sort_values("day").to_dict("records")
        for a,b in zip(visits,visits[1:]):
            if a["day"].year!=b["day"].year:continue
            gap=(b["day"]-a["day"]).days
            if gap<10 or gap>35:continue
            month=a["day"].month
            transitions.append({
                "Route":str(route),"Address_comp":str(address),
                "Milkweed_sp":species,"Native_status":a["Native_status"],
                "year":a["day"].year,"index_month":month,"days_gap":gap,
                "eggs0":a["eggs"],"late0":a["late"],
                "plants0":a["plants"],"plants1":b["plants"],
                "late1":b["late"],"y":int(b["late"]>0),
                "log_eggs":math.log1p(a["eggs"]),
                "log_prior_late":math.log1p(a["late"]),
                "log_plants0":math.log1p(a["plants"]),
                "log_plants1":math.log1p(b["plants"]),
                "month_sin":math.sin(2*math.pi*month/12),
                "month_cos":math.cos(2*math.pi*month/12),
            })
    t=pd.DataFrame(transitions)
    eligible=t[t.year==2022].Milkweed_sp.value_counts()
    eligible=list(eligible[eligible>=100].index)
    if set(eligible)!=set(SOURCE_SPECIES):
        raise RuntimeError(f"species data-only selection drift {eligible}")
    t=t[t.Milkweed_sp.isin(SOURCE_SPECIES)].copy().reset_index(drop=True)
    if len(t[t.year==2022])!=1901:
        raise RuntimeError(f"2022 transitions source drift {len(t[t.year==2022])}!=1901")
    if t.year.isna().any():raise RuntimeError("bad visit year")
    return t,{"complete_original_source_rows":len(d),"unique_patch_day":len(g),
        "same_date_multiple_row_groups":int((g.source_rows>1).sum()),
        "selected_species":SOURCE_SPECIES,"selection_min_2022_transitions":100,
        "eligible_transition_counts":{
            str(year):{sp:int((t.loc[t.year==year,"Milkweed_sp"]==sp).sum())
                       for sp in SOURCE_SPECIES} for year in (2022,2023,2024)}}

def pipeline(feature):
    pre=ColumnTransformer(
        transformers=[("numeric",StandardScaler(),NUM_COLS),
                      ("categorical",OneHotEncoder(handle_unknown="ignore"),feature)])
    return make_pipeline(pre,LogisticRegression(C=1,penalty="l2",solver="liblinear",
                                                random_state=20261008,max_iter=1000))
def fit_predict(train,test,feature):
    if len(train)<50 or train.y.nunique()!=2:raise RuntimeError("training set inadequate")
    model=pipeline(feature)
    model.fit(train[NUM_COLS+feature],train.y)
    return np.clip(model.predict_proba(test[NUM_COLS+feature])[:,1],1e-12,1-1e-12)

def metrics(y,pr):
    if not len(y):return None
    y=np.asarray(y,dtype=float);pr=np.asarray(pr,dtype=float)
    return {"n":len(y),"positive_events":int(y.sum()),"log_loss":float(log_loss(y,pr,labels=[0,1])),
        "brier":float(brier_score_loss(y,pr)),"pred_mean":float(pr.mean()),"observed_prevalence":float(y.mean())}

def bootstrap_route_delta(p, reps=4999,seed=20261008):
    # Draw entire routes (not individual patch observations) to respect temporal repeats.
    rng=np.random.default_rng(seed)
    groups={r:np.flatnonzero(p.Route.to_numpy()==r) for r in sorted(p.Route.unique())}
    routes=list(groups)
    obs=float((p.loss_origin-p.loss_species).mean())
    b=[]
    for _ in range(reps):
        rr=rng.choice(routes,size=len(routes),replace=True)
        ind=np.concatenate([groups[x] for x in rr])
        b.append(float((p.loss_origin-p.loss_species).to_numpy()[ind].mean()))
    q=np.quantile(b,[.025,.975])
    return {"origin_minus_species_mean_log_loss":obs,
       "route_bootstrap_ci95":[float(q[0]),float(q[1])],
       "route_bootstrap_reps":reps,"number_independent_routes":len(routes),
       "route_differences":[{"route":str(r),"n":len(groups[r]),
         "mean_origin_minus_species_loss":float((p.loss_origin-p.loss_species).to_numpy()[ix].mean())}
           for r,ix in groups.items()]}

def evaluate(t):
    tr=t.loc[t.year==2022].copy()
    if tr.Route.nunique()<10:raise RuntimeError("insufficient cross-validated survey routes")
    y=np.asarray(tr.y)
    cv_a=np.zeros(len(tr));cv_b=np.zeros(len(tr))
    route=tr.Route.to_numpy()
    for r in sorted(np.unique(route)):
        test=(route==r);train=~test
        cv_a[test]=fit_predict(tr.iloc[np.flatnonzero(train)],tr.iloc[np.flatnonzero(test)],ORIGIN_COL)
        cv_b[test]=fit_predict(tr.iloc[np.flatnonzero(train)],tr.iloc[np.flatnonzero(test)],SPECIES_COL)
    tr["pr_origin"]=cv_a
    tr["pr_species"]=cv_b
    tr["loss_origin"]=-(y*np.log(cv_a)+(1-y)*np.log(1-cv_a))
    tr["loss_species"]=-(y*np.log(cv_b)+(1-y)*np.log(1-cv_b))
    b=bootstrap_route_delta(tr)
    out={"route_holdout_2022":{
         "origin_model":metrics(tr.y,cv_a),"identity_model":metrics(tr.y,cv_b),
         "difference":b,
         "species_diagnostics":[{
            "species":sp,"origin_model":metrics(tr.loc[tr.Milkweed_sp==sp,"y"],tr.loc[tr.Milkweed_sp==sp,"pr_origin"]),
            "identity_model":metrics(tr.loc[tr.Milkweed_sp==sp,"y"],tr.loc[tr.Milkweed_sp==sp,"pr_species"])
        } for sp in SOURCE_SPECIES]}}
    for year in (2023,2024):
        te=t.loc[t.year==year].copy()
        if len(te)==0:continue
        a=fit_predict(tr,te,ORIGIN_COL);c=fit_predict(tr,te,SPECIES_COL)
        out[f"temporal_holdout_{year}"]={
          "origin_model":metrics(te.y,a),"identity_model":metrics(te.y,c),
          "origin_minus_identity_log_loss":float(log_loss(te.y,a,labels=[0,1])-log_loss(te.y,c,labels=[0,1])),
          "species":[{"species":sp,"n":int((te.Milkweed_sp==sp).sum()),
             "late_positive":int(te.loc[te.Milkweed_sp==sp,"y"].sum())}
            for sp in SOURCE_SPECIES]
        }
    delta=b["origin_minus_species_mean_log_loss"];lower=b["route_bootstrap_ci95"][0]
    temporal=out["temporal_holdout_2023"]["origin_minus_identity_log_loss"]
    out["decision"]={
        "gate_passed":bool(delta>0 and lower>0 and temporal>0),
        "if_not_passed":"No reliable host identity benefit beyond coarse origin in prospective 4th/5th instar detection; do not rescue with cutoff or selected-season tests.",
        "interpretation_boundary":"Next-visit late larva probability is detection/occupancy of different instars, NOT tracked-individual egg-to-larva survival, fitness, or demographic impact. Model has potential shared-visit and season confounding.",
        "posthoc_status":"Species identities and purpose were motivated by earlier source stage odds, so this test is exploratory despite protocol before predictive fit."
    }
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--raw-csv",type=Path,required=True)
    ap.add_argument("--protocol-json",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args()
    p=json.loads(a.protocol_json.read_text())
    if p["expected_four_species"]!=SOURCE_SPECIES:raise RuntimeError("protocol drift")
    t,source_audit=load(a.raw_csv)
    results={"schema":"chocho_monarch_revisit_host_identity_predictive_v0.1",
      "protocol":"docs/exploratory/URBAN_MILKWEED_FORWARD_HOST_IDENTITY_PROTOCOL_V01.json",
      "source":"Erickson, Schultz and Crone (2025); Figshare DOI 10.6084/m9.figshare.25648644.v1",
      "status":"EXPLORATORY_HELDOUT_DYNAMIC_STAGE_DETECTION",
      "source_audit":source_audit,**evaluate(t)}
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(results,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(results,indent=2),flush=True)
if __name__=="__main__":main()
