#!/usr/bin/env python3
"""Out-of-year validation of an explicit host-decline x host-backup interaction.

Observational nonoccupancy only. Source and model gate frozen before fitted results.
No ecological causality, invasion history, demographic or evolutionary mechanism.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline

import audit_melitaea_19year_fallback_support as src

URL=src.SOURCE_URL
MD5=src.MD5
NET="data/patch_network.tsv"
PRIMARY_FIELDS=["log_nests","plantago_t","veronica_t","plantago_prev","veronica_prev",
                "veronica_decline","plantago_high","year_linear","year_squared",
                "log_area","x_10km","y_10km","grazing_presence","grazing_intensity",
                "plantago_dry","veronica_dry"]

def get_source():
    req=Request(URL,headers={"User-Agent":"chocho-audited-research/1.0"})
    with urlopen(req,timeout=45) as res:
        raw=res.read(src.MAX_BYTES+1)
    if len(raw)>src.MAX_BYTES or hashlib.md5(raw).hexdigest()!=MD5:
        raise ValueError("source MD5 not exact")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.namelist().count(src.MEMBER)!=1 or z.namelist().count(NET)!=1:
            raise ValueError("source missing exact TSVs")
        return z.read(src.MEMBER),z.read(NET),hashlib.sha256(raw).hexdigest()

def read_static(tsv):
    reader=csv.DictReader(io.StringIO(tsv.decode("utf-8-sig")),delimiter="\t")
    if not {"patch","x","y","area"}.issubset(reader.fieldnames or []):
        raise ValueError("missing static position or area fields")
    out={}
    for r in reader:
        p=r["patch"].strip()
        if p in out:
            raise ValueError("duplicated patch geography")
        out[p]={k:src.parse_num(r.get(k)) for k in ["x","y","area"]}
        if out[p]["area"] is not None and out[p]["area"]<0:
            raise ValueError("negative site area")
    return out

def build_panel(rows,net):
    panel=[]
    counts=Counter()
    for (patch,year),r in sorted(rows.items()):
        before=rows.get((patch,year-1))
        after=rows.get((patch,year+1))
        if before is None or after is None:
            continue
        if r["population"] is None or r["population"]<=0 or after["population"] is None:
            continue
        if any(x is None for x in (before["veronica"],r["veronica"],r["plantago"])):
            continue
        counts["total"]+=1
        last_pl=before["plantago"]
        decline=int(before["veronica"]>r["veronica"] and before["veronica"]>=1)
        high=int(r["plantago"]>=2)
        geo=net.get(patch,{})
        def div(v, scale):
            return None if v is None else v/scale
        area=geo.get("area")
        base=[
            math.log1p(r["population"]),r["plantago"],r["veronica"],
            last_pl,before["veronica"],decline,high,year-2009,(year-2009)**2,
            math.log1p(area) if area is not None else None,
            div(geo.get("x"),10000),div(geo.get("y"),10000),
            r["grazing_presence"],r["grazing_intensity"],r["plantago_dry"],r["veronica_dry"],
        ]
        panel.append({"patch":patch,"year":year,"y":int(after["population"]==0),
                      "decline":decline,"backup_high":high,"base":base,
                      "interaction":decline*high})
    return panel

def predictions(panel, increment=False):
    years=sorted(set(r["year"] for r in panel))
    arr=np.array([[np.nan if x is None else x for x in r["base"]] +
                  ([r["interaction"]] if increment else []) for r in panel],dtype=float)
    ys=np.array([r["y"] for r in panel],dtype=int)
    pred=np.full(len(panel),np.nan)
    for yr in years:
        test=np.array([r["year"]==yr for r in panel])
        train=~test
        if not test.any() or len(set(ys[train]))<2:
            raise ValueError("invalid temporal holdout")
        mod=make_pipeline(SimpleImputer(strategy="median",add_indicator=True,keep_empty_features=True),
                          StandardScaler(),
                          LogisticRegression(C=1.0,solver="lbfgs",max_iter=1000))
        mod.fit(arr[train],ys[train])
        pred[test]=mod.predict_proba(arr[test])[:,1]
    if not np.isfinite(pred).all():
        raise ValueError("nonfinite heldout predictions")
    return np.clip(pred,1e-9,1-1e-9)

def cluster_interval(diff, groups, seed, iterations=2000):
    u=sorted(set(groups))
    idx={g:np.array([i for i,k in enumerate(groups) if k==g]) for g in u}
    rng=np.random.default_rng(seed)
    result=np.empty(iterations)
    for j in range(iterations):
        ids=rng.choice(u,size=len(u),replace=True)
        arr=np.concatenate([idx[g] for g in ids])
        result[j]=diff[arr].mean()
    return [float(x) for x in np.quantile(result,[0.025,0.975])]

def assess(panel):
    if len(panel)<100:
        raise ValueError("missing source support")
    a=predictions(panel,False)
    b=predictions(panel,True)
    y=np.array([r["y"] for r in panel],dtype=int)
    ll=lambda p: -(y*np.log(p)+(1-y)*np.log(1-p))
    losses=ll(a)-ll(b)
    brier=(a-y)**2-(b-y)**2
    obs={}
    years=sorted(set(r["year"] for r in panel))
    for yr in years:
        ix=np.array([r["year"]==yr for r in panel])
        obs[str(yr)]={"records":int(ix.sum()),"logloss_gain":float(losses[ix].mean())}
    patch_ci=cluster_interval(losses,[r["patch"] for r in panel],20261010)
    year_ci=cluster_interval(losses,[r["year"] for r in panel],20261011)
    cells={}
    for decline in [0,1]:
        for high in [0,1]:
            i=np.array([r["decline"]==decline and r["backup_high"]==high for r in panel])
            cells[f"decline_{decline}_backup_{high}"]={
                "n":int(i.sum()),"next_year_zero":int(y[i].sum()),
                "raw_nonoccupancy":float(y[i].mean()) if i.any() else None,
                "model_base_predicted_nonoccupancy":float(a[i].mean()) if i.any() else None,
                "model_extra_predicted_nonoccupancy":float(b[i].mean()) if i.any() else None
            }
    k=sum(v["logloss_gain"]>0 for v in obs.values())
    gain=float(losses.mean())
    passed=gain>0 and patch_ci[0]>0 and year_ci[0]>0 and k>=12
    return {"n":len(panel),"patches":len(set(r["patch"] for r in panel)),
            "years":len(years),"event_cells":cells,
            "baseline_logloss":float(ll(a).mean()),"increment_logloss":float(ll(b).mean()),
            "mean_heldout_logloss_gain":gain,"mean_heldout_brier_gain":float(brier.mean()),
            "patch_cluster_95ci":patch_ci,"year_cluster_95ci":year_ci,
            "positive_heldout_years":k,"per_year":obs,
            "frozen_predictive_gate_pass":bool(passed),
            "decision":("OBSERVATIONAL_PREDICTIVE_SIGNAL_ONLY" if passed
                        else "NO_ROBUST_HELDOUT_HOST_BACKUP_GAIN"),
            "evolutionary_or_causal_hysteresis_estimated":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    survey,network,sha=get_source()
    rows,meta=src.ingest(survey)
    support=src.support(rows)
    if not support["data_support_pass"]:
        raise ValueError("predeclared source support gate did not pass")
    net=read_static(network)
    panel=build_panel(rows,net)
    if len(panel)!=support["three_year_consecutive_occupied_risk_set"]:
        raise ValueError("source support comparison mismatch")
    outcome=assess(panel)
    receipt={"schema":"chocho_19year_dynamic_host_backup_v01",
             "model_protocol":"docs/exploratory/MELITAEA_19YEAR_FALLBACK_FIT_ADDENDUM_V01.json",
             "data_doi":"10.5061/dryad.ksn02v707","source_MD5":MD5,"source_sha256":sha,
             "source_distinct_from_DiLeo":True,
             "observational_not_host_removal":True,"source_support":support,
             "model":outcome,"GEB_pr38_untouched":True}
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
