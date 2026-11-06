#!/usr/bin/env python3
"""Source-locked island-pair circuit-resistance *geometric baseline*.

Uses the existing 51 island landmasses / 1,275 island pairs, with representative
WPI port great-circle distances and Jaccard host/butterfly similarities.

IMPORTANT: Distance-only conductance is NOT a map of currents, AIS shipping,
wind direction, oceanographic currents or historical transport. This is a
discrete undirected NETWORK analogue of Circuitscape pairwise mode, calculated
exactly from a graph Laplacian pseudoinverse. The topology is predeclared as
k nearest islands + a minimum spanning tree for global connectivity.

Validation is truly island-held-out: for each of 51 island nodes, fit ONLY
pairs among the other 50 and predict the 50 held-out pairs. Graph is built
from outcome-independent geography and may include held-out node positions.
No pair-independent confidence intervals or p-values are computed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

N_ISLANDS=51
N_PAIRS=N_ISLANDS*(N_ISLANDS-1)//2
PREDECLARED=[
    ("primary_knn3_linear_distance",3,1.0),
    ("sensitivity_knn2_linear_distance",2,1.0),
    ("sensitivity_knn5_linear_distance",5,1.0),
    ("sensitivity_knn8_linear_distance",8,1.0),
    ("sensitivity_knn3_squared_distance",3,2.0),
]


def source(pairs_path:Path,islands_path:Path):
    p=pd.read_csv(pairs_path,dtype={"island_a":str,"island_b":str})
    ids=pd.read_csv(islands_path,dtype={"landmass_id":str})
    ids=ids[ids.eligible.astype(str).str.lower().eq("true")].copy()
    if (len(ids),ids.landmass_id.nunique())!=(N_ISLANDS,N_ISLANDS):
        raise ValueError("51 source-eligible island endpoints not retained")
    if len(p)!=N_PAIRS:raise ValueError("Expected exact 1,275 original island pairs")
    fixed=sorted(ids.landmass_id.tolist(),key=lambda s:int(s))
    names={s:i for i,s in enumerate(fixed)}
    distances=np.zeros((N_ISLANDS,N_ISLANDS),dtype=float)
    seen=set()
    for row in p.itertuples(index=False):
        a,b=str(row.island_a),str(row.island_b)
        if a not in names or b not in names or a==b:
            raise ValueError("Unknown source island IDs / self-distance")
        i,j=names[a],names[b]
        pair=(min(i,j),max(i,j))
        if pair in seen:raise ValueError("Repeated island pair")
        seen.add(pair)
        d=float(row.port_distance_km)
        if not np.isfinite(d) or d<=0:raise ValueError("Missing or nonpositive source distance")
        distances[i,j]=d;distances[j,i]=d
    if len(seen)!=N_PAIRS or (distances[np.triu_indices(N_ISLANDS,1)]<=0).any():
        raise ValueError("Unfilled source great-circle distance matrix")
    for col in ("host_jaccard","earlier_butterfly_jaccard","new_butterfly_jaccard"):
        if not np.isfinite(p[col].to_numpy()).all() or not p[col].between(0,1).all():
            raise ValueError("Observed Jaccard source changed "+col)
    q=p.copy()
    q["i"]=q.island_a.map(names).astype(int)
    q["j"]=q.island_b.map(names).astype(int)
    q["distance_km_verified"]=distances[q.i.to_numpy(),q.j.to_numpy()]
    if not np.allclose(q.distance_km_verified,q.port_distance_km,atol=1e-9):
        raise ValueError("Geographic baseline mismatched")
    return q,fixed,distances


def mst_edges(dist):
    """Deterministic, no outcomes: sorted Kruskal spanning tree."""
    n=len(dist)
    parent=list(range(n))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    pairs=sorted((float(dist[i,j]),i,j) for i in range(n) for j in range(i+1,n))
    edges=[]
    for _,i,j in pairs:
        a,b=root(i),root(j)
        if a!=b:
            parent[a]=b
            edges.append((i,j))
            if len(edges)==n-1:break
    if len(edges)!=n-1:raise ValueError("MST incomplete")
    return set(edges)


def adjacency(dist,k=3,power=1.0):
    n=len(dist)
    if n<3 or not (1<=k<n) or power<=0:
        raise ValueError("Invalid circuit parameters")
    edges=set()
    for i in range(n):
        order=np.argsort(np.where(np.arange(n)==i,np.inf,dist[i]))
        for j in order[:k]:
            edges.add((min(i,int(j)),max(i,int(j))))
    original=len(edges)
    mst=mst_edges(dist)
    added=mst-edges
    edges|=mst
    lap=np.zeros((n,n),dtype=float)
    all_d=[]
    for i,j in edges:
        d=float(dist[i,j])
        all_d.append(d)
        conductance=(d/1000.0)**(-power)
        lap[i,i]+=conductance
        lap[j,j]+=conductance
        lap[i,j]-=conductance
        lap[j,i]-=conductance
    diag={
        "n_nodes":n,"knn":k,"edge_cost_distance_power":power,
        "n_knn_undirected_edges":original,"mst_bridge_edges_added":len(added),
        "n_total_edges":len(edges),
        "mean_edge_geodesic_km":float(np.mean(all_d)),
        "median_edge_geodesic_km":float(np.median(all_d)),
        "q90_edge_geodesic_km":float(np.quantile(all_d,.90)),
        "max_edge_geodesic_km":float(np.max(all_d)),
        "mst_bridge_geodesic_km":[float(dist[i,j]) for i,j in sorted(added)],
        "note":"MST bridges can cross whole oceans; these are hypothetical geographic links, not ship routes."
    }
    return lap,diag


def resistance_from_laplacian(lap):
    if lap.shape[0]!=lap.shape[1]:raise ValueError("Laplacian not square")
    if not np.allclose(lap,lap.T):raise ValueError("Undirected circuit must be symmetric")
    eig,vec=np.linalg.eigh(lap)
    if eig[-1]<=0 or sum(eig<max(eig[-1]*1e-11,1e-12))!=1:
        raise ValueError("Graph disconnected or Laplacian invalid")
    lp=(vec[:,1:]/eig[1:])@vec[:,1:].T
    d=np.diag(lp)
    r=d[:,None]+d[None,:]-2*lp
    r=np.maximum((r+r.T)/2,0)
    np.fill_diagonal(r,0)
    if not np.allclose(r,r.T) or np.min(r)<0 or not np.isfinite(r).all():
        raise ValueError("Invalid effective resistance")
    return r


def rank_corr(x,y):
    # Ranked source summary only, no invalid independent pair p-values.
    a=pd.Series(x).rank(method="average").to_numpy()
    b=pd.Series(y).rank(method="average").to_numpy()
    return float(np.corrcoef(a,b)[0,1])


def make_features(q,eff,*,response):
    d=np.log1p(q.port_distance_km.to_numpy(dtype=float))
    features={
      "log_great_circle_distance":d,
      "squared_log_great_circle_distance":d*d
    }
    if response=="new_butterfly_jaccard":
        features["earlier_butterfly_Jaccard"]=q.earlier_butterfly_jaccard.to_numpy(dtype=float)
    if eff is not None:
        r=eff[q.i.to_numpy(),q.j.to_numpy()]
        features["log1p_effective_resistance"]=np.log1p(r)
    return features


def train_predict(X,y,train,test):
    # Scaling uses ONLY training island-pair rows, no test outcomes.
    a=X[train];b=X[test]
    mu=np.mean(a,axis=0)
    sd=np.std(a,axis=0)
    sd=np.where(sd<1e-12,1.0,sd)
    A=np.column_stack([np.ones(len(a)),(a-mu)/sd])
    B=np.column_stack([np.ones(len(b)),(b-mu)/sd])
    beta=np.linalg.lstsq(A,y[train],rcond=None)[0]
    return np.clip(B@beta,0,1)


def holdout(q,base,aug,target,ids):
    a=np.column_stack(list(base.values()))
    b=np.column_stack(list(aug.values()))
    y=q[target].to_numpy(dtype=float)
    errs=[]
    for i,name in enumerate(ids):
        test=(q.i.to_numpy()==i)|(q.j.to_numpy()==i)
        train=~test
        if (int(train.sum()),int(test.sum()))!=(1225,50):
            raise ValueError("Not a full 51-island LOIO split")
        p0=train_predict(a,y,train,test)
        p1=train_predict(b,y,train,test)
        obs=y[test]
        l0=float(np.mean((obs-p0)**2))
        l1=float(np.mean((obs-p1)**2))
        errs.append({"heldout_island_id":name,"test_pairs":int(test.sum()),
             "train_pairs":int(train.sum()),
             "baseline_mse":l0,"circuit_augmented_mse":l1,
             "delta_mse_circuit_minus_baseline":l1-l0,
             "circuit_improves_this_island":l1<l0,
             "MAE_base":float(np.mean(np.abs(obs-p0))),
             "MAE_circuit":float(np.mean(np.abs(obs-p1)))})
    m0=float(np.mean([x["baseline_mse"] for x in errs]))
    m1=float(np.mean([x["circuit_augmented_mse"] for x in errs]))
    return {
      "target":target,"baseline_features":list(base),
      "augmented_features":list(aug),
      "n_islands_left_out":len(ids),
      "n_predictions_heldout_island_pairs":len(ids)*50,
      "n_unique_observed_pairs":len(q),
      "pairs_each_predicted_twice_one_per_endpoint":True,
      "MSE_baseline":m0,"MSE_augmented":m1,
      "relative_MSE_reduction_percent":float((m0-m1)/m0*100) if m0>0 else None,
      "absolute_MSE_gain":m0-m1,
      "islands_with_better_augmented_prediction":sum(x["circuit_improves_this_island"] for x in errs),
      "islands_with_worse_augmented_prediction":sum(not x["circuit_improves_this_island"] for x in errs),
      "MSE_ratio_aug_over_base":m1/m0 if m0>0 else None,
      "MAE_baseline":float(np.mean([x["MAE_base"] for x in errs])),
      "MAE_augmented":float(np.mean([x["MAE_circuit"] for x in errs])),
      "per_island_LOIO":errs,
      "no_pairwise_p_values_or_CI":True,
      "no_butterfly_colonization_or_ship_transport_claim":True
    }


def run(pairpath:Path,islandpath:Path):
    q,nodes,dist=source(pairpath,islandpath)
    out={
     "schema":"chocho_island_51_discrete_circuit_distance_steppingstone_holdout_v1",
     "status":"SOURCE_VERIFIED_GEOMETRIC_CIRCUIT_NETWORK_NOT_AIS_WIND_OR_RASTER_CIRCUITSCAPE",
     "sources":{"island_host_pairs_SHA256":hashlib.sha256(pairpath.read_bytes()).hexdigest(),
                "island_host_islands_SHA256":hashlib.sha256(islandpath.read_bytes()).hexdigest()},
     "n_eligible_island_landmasses":len(nodes),"n_dependent_pairs":len(q),
     "graph_design":"undirected k nearest islands plus global MST connectivity, edge conductance = (km/1000)^(-power)",
     "predeclared_primary":{"graph":"primary_knn3_linear_distance",
        "selection":"independent island-out MSE relative to log1p distance and its square",
        "butterfly_comparator":"earlier 2001–2010 butterfly Jaccard included in both butterfly models",
        "plant_comparator":"direct host Jaccard and geographic distance only"},
     "ecological_limits":[
        "This is a mathematical NETWORK circuit calculated from a source distance matrix, not Circuitscape.jl raster resistance model.",
        "No AIS shipping lanes, marine wind field, coastline/land-ocean raster, currents or observed island pair traffic were used.",
        "Arbitrary kNN/MST hypothetical edge configuration contains long inter-ocean bridges; robustness across 2,3,5,8 neighbors is mandatory.",
        "Great-circle distance from average port locations is not actual sea distance, port distance, water flow or historical transport exposure.",
        "The old and later butterfly Jaccard source indices are retrospective index sets, not first biological arrivals or verified breeding.",
        "All 1,275 pairwise observations share 51 island endpoints; only 51 island-heldout folds are meaningful, but folds are dependent and have no independent p-values.",
        "Graph topology including held-out locations is permissible geographic transduction but is not forecasting unseen geography without its distance matrix.",
        "Host candidate plants were pooled in 2001–2024 and could postdate butterfly records, so the models cannot be prospective transport validation.",
        "Do not choose graph k or power based on this exploratory sample and claim independent optimal transport pathways."
     ],
     "models":{}
    }
    for label,k,power in PREDECLARED:
        lap,graphmeta=adjacency(dist,k=k,power=power)
        r=resistance_from_laplacian(lap)
        dv=q.port_distance_km.to_numpy()
        rv=r[q.i.to_numpy(),q.j.to_numpy()]
        graphmeta["rank_rho_resistance_vs_geodesic_km"]=rank_corr(rv,dv)
        graphmeta["rank_rho_resistance_vs_host_Jaccard"]=rank_corr(rv,q.host_jaccard)
        graphmeta["rank_rho_resistance_vs_new_butterfly_Jaccard"]=rank_corr(rv,q.new_butterfly_jaccard)
        results={"graph":graphmeta,"targets":{}}
        for target in ("host_jaccard","new_butterfly_jaccard","earlier_butterfly_jaccard"):
            features=make_features(q,r,response=target)
            original=features.copy()
            original.pop("log1p_effective_resistance")
            if target=="earlier_butterfly_jaccard":
                # Earlier butterfly assemblage as negative-comparison target,
                # no earlier butterfly composition used as its own predictor.
                pass
            model=holdout(q,original,features,target,nodes)
            results["targets"][target]=model
        out["models"][label]=results
        print("ISLAND_CIRCUIT_GEOMETRY_LOIO",label,json.dumps({
            "graph":{x:y for x,y in graphmeta.items() if x!="mst_bridge_geodesic_km"},
            "targets":{k:{"base":v["MSE_baseline"],
                         "augmented":v["MSE_augmented"],
                         "incremental_MSE_pct":v["relative_MSE_reduction_percent"],
                         "improved_islands":v["islands_with_better_augmented_prediction"]}
                       for k,v in results["targets"].items()}
           },sort_keys=True),flush=True)
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairs",type=Path,required=True)
    p.add_argument("--islands",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=run(a.pairs,a.islands)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    m=z["models"]["primary_knn3_linear_distance"]["targets"]
    print("ISLAND_CIRCUIT_51_PRIMARY_RESULT",json.dumps({
        "host_MSE_reduction_pct":m["host_jaccard"]["relative_MSE_reduction_percent"],
        "new_butterfly_MSE_reduction_pct":m["new_butterfly_jaccard"]["relative_MSE_reduction_percent"],
        "old_butterfly_negative_control_pct":m["earlier_butterfly_jaccard"]["relative_MSE_reduction_percent"]
    },sort_keys=True),flush=True)


if __name__=="__main__":main()
