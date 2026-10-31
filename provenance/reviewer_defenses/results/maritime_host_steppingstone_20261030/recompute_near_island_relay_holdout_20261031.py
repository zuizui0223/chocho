#!/usr/bin/env python3
"""Recompute early botanical support and time-held-out island record predictions.

Raw species x landmass x year Parquet stays private. This analyzes *recording*,
not actual butterfly arrival, reproduction, ship transport, or causal host help.
Needs pandas, numpy, scikit-learn, pyarrow (for reading the private Parquet).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
SRC=Path("provenance/reviewer_defenses/results")
PAIR=SRC/"island_host_network_20261030/island_host_network_pairs.csv"
PLANT=SRC/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
OUT=SRC/"maritime_host_steppingstone_20261030/near_island_host_opportunity_holdout_20261031.json"
PRIVATE="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
PAIR_BLOB="5663d2ffe42419f5c4496b92c9f3bbf76175907b"
PLANT_BLOB="0a01831056c2cd59a4bd41bd1dbd4c84d77cf922"
NEAR=[(79,110),(81,214),(169,276),(169,715),(169,1577),(206,397),
      (214,242),(242,257),(260,753),(276,321),(276,715),(276,1577),
      (715,1577),(821,951)]
FOCAL={(169,276):{5137612,5137706},(276,715):{5130682}}
BASE=["log_target_effort","log_source_effort","log_global_islands",
      "log_area","log_min_dist","year_ctr","n_source","first_seen_age",
      "last_seen_age","n_recent_source"]

def file_sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()

def blob_sha(p):
    raw=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def analyze(path:Path, root:Path=ROOT, decoded_pickle:Path|None=None):
    import numpy as np
    import pandas as pd
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss
    if file_sha256(path)!=PRIVATE:
        raise ValueError("Original-year raw SHA-256 mismatch")
    if blob_sha(root/PAIR)!=PAIR_BLOB or blob_sha(root/PLANT)!=PLANT_BLOB:
        raise ValueError("Public source Git blob mismatch")
    D=pd.read_pickle(decoded_pickle) if decoded_pickle else pd.read_parquet(path)
    if (len(D),D.speciesKey.nunique(),D.landmass_id.nunique(),
        int(D.year.min()),int(D.year.max()),int(D.records.sum()))!=(83573,224,3384,2001,2024,14425536):
        raise ValueError("Original private Parquet cohort changed")
    if D.duplicated(["speciesKey","landmass_id","year"]).any() or not D.records.eq(D.precise1km_records).all():
        raise ValueError("Original source duplicates/precision mismatch")
    D=D.copy()
    D["landmass_id"]=D.landmass_id.astype(int)
    D["speciesKey"]=D.speciesKey.astype(int)
    pairs={}
    for row in csv.DictReader((root/PAIR).open(newline="",encoding="utf8")):
        a,b=sorted((int(row["island_a"]),int(row["island_b"])))
        pairs[a,b]=(float(row["port_distance_km"]),int(row["shared_new_species"]))
    near={k:v for k,v in pairs.items() if v[0]<200}
    if len(pairs)!=1275 or set(near)!=set(NEAR):
        raise ValueError("14 near port-representative island pairs changed")
    plants=defaultdict(list)
    for row in csv.DictReader((root/PLANT).open(newline="",encoding="utf8")):
        plants[int(row["landmass_id"])].append(int(row["early_plant_records"]))
    plants={k:max(v) for k,v in plants.items()}
    graph=defaultdict(list)
    for a,b in NEAR:
        dist=near[a,b][0]
        graph[a].append((b,dist))
        graph[b].append((a,dist))
    species=sorted(int(k) for k in D.speciesKey.unique())
    years={(int(k),int(i)):sorted(int(t) for t in g.year)
           for (k,i),g in D.groupby(["speciesKey","landmass_id"])}
    local=D[D.landmass_id.isin(graph)]
    size=local.groupby("landmass_id").landmass_area_km2.first().to_dict()
    obs=local.groupby(["landmass_id","year"]).agg(taxa=("speciesKey","nunique"),
                  rec=("records","sum"))
    effort={(int(i),int(t)):(int(r.taxa),int(r.rec)) for (i,t),r in obs.iterrows()}
    seen=defaultdict(lambda:defaultdict(int))
    for (k,i),v in years.items():
        for t in range(max(v[0]+1,2001),2023):
            seen[k][t]+=1
    groups=defaultdict(lambda:[0,0,set()])
    all_new=0
    all_chain=0
    for (a,b),(dist,published_count) in near.items():
        n=0
        for k in species:
            ay=years.get((k,a),[])
            by=years.get((k,b),[])
            if not ay or not by or not(2011<=ay[0]<=2022 and 2011<=by[0]<=2022):
                continue
            n+=1
            if ay[0]<by[0]:up,down=ay,by
            elif by[0]<ay[0]:up,down=by,ay
            else:up,down=[],[]
            chain=bool(up and any(up[0]<t<down[0] and t<=up[0]+5 for t in up))
            focal=k in FOCAL.get((a,b),set())
            bot=min(plants[a],plants[b])>=100
            key=("focal_" if focal else "other_")+("high" if bot else "low")
            g=groups[key]
            g[0]+=1;g[1]+=int(chain);g[2].add((a,b))
            all_chain+=int(chain)
        if n!=published_count:
            raise ValueError(f"Public count and original year-table disagree: {(a,b)}")
        all_new+=n
    if (all_new,all_chain)!=(41,21):
        raise ValueError("Full near-pair annual chronology changed")
    cohort={key:{"n":v[0],"chains":v[1],"pairs":len(v[2])} for key,v in groups.items()}
    if (cohort["focal_high"]["n"],cohort["focal_high"]["chains"],
        cohort["other_high"]["n"],cohort["other_high"]["chains"])!=(3,3,9,3):
        raise ValueError("Fixed earlier botanical opportunity comparison changed")
    # As-of-time source status. Each risk row is unique species/destination/year.
    rows=[]
    for target,neighbors in graph.items():
        for k in species:
            dest=years.get((k,target),[])
            fy=dest[0] if dest else 9999
            for t in range(2011,min(2023,fy+1)):
                prev=[effort.get((target,y),(0,0)) for y in range(t-5,t)]
                if sum(v[0]>0 for v in prev)<2:
                    continue
                src=[]
                for i,d in neighbors:
                    ys=[y for y in years.get((k,i),[]) if y<t]
                    if ys:src.append((i,d,ys))
                if not src:continue
                src_eff=sum(effort.get((i,y),(0,0))[1]
                            for i,_,_ in src for y in range(t-5,t))
                rows.append(dict(target=target,speciesKey=k,year=t,first=int(t==fy),
                    log_target_effort=math.log1p(sum(z[1] for z in prev)),
                    log_source_effort=math.log1p(src_eff),
                    log_global_islands=math.log1p(seen[k][t]),
                    log_area=math.log1p(size[target]),
                    log_min_dist=math.log1p(min(d for _,d,_ in src)),
                    year_ctr=t-2014,
                    n_source=len(src),
                    first_seen_age=t-min(ys[0] for _,_,ys in src),
                    last_seen_age=t-max(ys[-1] for _,_,ys in src),
                    n_recent_source=sum(any(y>=t-5 for y in ys) for _,_,ys in src),
                    source_repeat=int(any(len(ys)>=2 for _,_,ys in src)),
                    source_recentrepeat=int(any(sum(y>=t-5 for y in ys)>=2 for _,_,ys in src))))
    P=pd.DataFrame(rows).sort_values(["target","speciesKey","year"]).reset_index(drop=True)
    if (len(P),int(P["first"].sum()),P.speciesKey.nunique(),P.target.nunique())!=(865,91,70,17):
        raise ValueError("Annual source/destination risk cohort changed")
    if P.duplicated(["target","speciesKey","year"]).any():
        raise ValueError("Pseudoreplicated species/destination/year risk row")
    scores=[]
    for train_last,test_first,test_last in [(2015,2016,2018),
                                             (2017,2018,2020),
                                             (2019,2020,2022),
                                             (2017,2018,2022)]:
        tr=P[P.year<=train_last]
        te=P[P.year.between(test_first,test_last)]
        y=te["first"].to_numpy()
        predictions={}
        ll={}
        for label,extra in [("baseline",[]),("ever_repeat",["source_repeat"]),
                            ("recent_repeat",["source_recentrepeat"])]:
            cols=BASE+extra
            mod=make_pipeline(StandardScaler(),
                LogisticRegression(C=1.0,max_iter=1000,solver="lbfgs"))
            mod.fit(tr[cols],tr["first"])
            preds=mod.predict_proba(te[cols])[:,1]
            predictions[label]=preds
            ll[label]=float(log_loss(y,preds,labels=[0,1]))
        r=dict(train_through=train_last,test=[test_first,test_last],
               train_rows=len(tr),train_events=int(tr["first"].sum()),
               test_rows=len(te),test_events=int(y.sum()),logloss=ll)
        if (train_last,test_first,test_last)==(2017,2018,2022):
            target=te.target.to_numpy()
            ids=sorted(set(target))
            idx={i:np.flatnonzero(target==i) for i in ids}
            loss=lambda p:-(y*np.log(p)+(1-y)*np.log1p(-p))
            for label in ("ever_repeat","recent_repeat"):
                rng=np.random.default_rng(20261031)
                gain=loss(predictions["baseline"])-loss(predictions[label])
                distribution=[]
                for _ in range(1500):
                    sampled=rng.choice(ids,len(ids),replace=True)
                    inds=np.concatenate([idx[i] for i in sampled])
                    distribution.append(gain[inds].mean())
                r.setdefault("island_bootstrap",{})[label]={
                    "gain":float(gain.mean()),
                    "ci95":[float(v) for v in np.quantile(distribution,[0.025,0.975])]}
        scores.append(r)
    return {"cohort":cohort,
            "annual_risk":{"rows":len(P),"events":int(P["first"].sum()),
                "source_unrepeated":[int(sum(P.source_repeat==0)),
                    int(P[P.source_repeat==0]["first"].sum())],
                "source_repeated":[int(sum(P.source_repeat==1)),
                    int(P[P.source_repeat==1]["first"].sum())]},
            "prediction":scores,"source_sha256":PRIVATE,
            "NO_CAUSAL_HOST_SHIPPING_OR_RELAY_INFERENCE":True}

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--original-year-parquet",type=Path,required=True)
    a.add_argument("--source-root",type=Path,default=ROOT)
    a.add_argument("--decoded-local-pickle",type=Path,default=None,
       help="Optional previously SHA-verified *local* decoding, not to publish")
    a.add_argument("--output",type=Path,required=True)
    opts=a.parse_args()
    result=analyze(opts.original_year_parquet,opts.source_root,opts.decoded_local_pickle)
    opts.output.parent.mkdir(parents=True,exist_ok=True)
    opts.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(result["annual_risk"]))

if __name__=="__main__":
    main()
