#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def average_ranks(values):
    x=np.asarray(values,dtype=float)
    order=np.argsort(x,kind="mergesort")
    out=np.empty(len(x),dtype=float)
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and x[order[j]]==x[order[i]]:
            j+=1
        out[order[i:j]]=0.5*((i+1)+j)
        i=j
    return out


def corr(a,b):
    a=np.asarray(a,dtype=float); b=np.asarray(b,dtype=float)
    if len(a)<3 or np.std(a)==0 or np.std(b)==0:
        return None
    return float(np.corrcoef(a,b)[0,1])


def spearman(a,b):
    return corr(average_ranks(a),average_ranks(b))


def residualize(values, controls):
    y=np.asarray(values,dtype=float)
    X=np.column_stack([np.ones(len(y))]+[np.asarray(x,dtype=float) for x in controls])
    return y-X@np.linalg.lstsq(X,y,rcond=None)[0]


def group_residualize(values, groups):
    x=np.asarray(values,dtype=float)
    out=x.copy()
    buckets=defaultdict(list)
    for i,g in enumerate(groups):
        buckets[g].append(i)
    for inds in buckets.values():
        idx=np.asarray(inds,dtype=int)
        out[idx]-=np.mean(x[idx])
    return out


def load_pairs(path):
    consumers=defaultdict(set)
    family={}
    with path.open(newline="",encoding="utf-8") as handle:
        reader=csv.DictReader(handle)
        for row in reader:
            insect=str(row["insect_species"]).strip()
            host=str(row["accepted_plant_name_id"]).strip()
            fam=str(row["family"]).strip()
            if insect and host and fam:
                consumers[host].add(insect)
                family[host]=fam
    return consumers,family


def load_units(path):
    out=defaultdict(set)
    with path.open(newline="",encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            h=str(row["accepted_plant_name_id"]).strip()
            u=str(row["area_code_l3"]).strip()
            if h and u: out[h].add(u)
    return out


def quantile_bins(values,n=5):
    x=np.asarray(values,dtype=float)
    qs=np.quantile(x,np.linspace(0,1,n+1))
    # protect duplicated cut points
    return np.searchsorted(qs[1:-1],x,side="right")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--native-distribution-csv",type=Path,required=True)
    ap.add_argument("--contemporary-distribution-csv",type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=4999)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    consumers,family=load_pairs(args.insect_host_csv)
    native=load_units(args.native_distribution_csv)
    contemporary=load_units(args.contemporary_distribution_csv)

    rows=[]
    for host,users in consumers.items():
        n=set(native.get(host,set()))
        c=set(contemporary.get(host,set()))
        fam=family.get(host,"")
        if not fam or not n or not n<=c:
            continue
        rows.append({
            "host":host,
            "family":fam,
            "degree":len(users),
            "native_units":len(n),
            "contemporary_units":len(c),
            "added_units":len(c-n),
            "log_expansion":math.log1p(len(c))-math.log1p(len(n)),
        })
    if len(rows)<100:
        raise RuntimeError("too few eligible host plants")

    degree=np.asarray([r["degree"] for r in rows],float)
    native_units=np.asarray([r["native_units"] for r in rows],float)
    added=np.asarray([r["added_units"] for r in rows],float)
    logexp=np.asarray([r["log_expansion"] for r in rows],float)
    fam=np.asarray([r["family"] for r in rows],dtype=object)

    rd=average_ranks(np.log1p(degree))
    rn=average_ranks(np.log1p(native_units))
    ry=average_ranks(logexp)
    ra=average_ranks(added)

    partial_logexp=corr(residualize(rd,[rn]),residualize(ry,[rn]))
    partial_added=corr(residualize(rd,[rn]),residualize(ra,[rn]))

    family_adj_logexp=corr(group_residualize(rd,fam),group_residualize(ry,fam))

    bins=quantile_bins(np.log1p(native_units),5)
    strata=np.asarray([f"{f}|q{b}" for f,b in zip(fam,bins)],dtype=object)
    xd=group_residualize(rd,strata)
    yy=group_residualize(ry,strata)
    observed_stratified=corr(xd,yy)

    buckets=defaultdict(list)
    for i,g in enumerate(strata):
        buckets[g].append(i)
    permutable=[np.asarray(v,dtype=int) for v in buckets.values() if len(v)>1]
    rng=np.random.default_rng(20260928)
    null=np.empty(args.permutations,float)
    # Because permutation is within strata, each stratum mean is invariant.
    # Permuting the already residualized degree ranks is exactly equivalent
    # and avoids rebuilding group means on every replicate.
    for k in range(args.permutations):
        rp=xd.copy()
        for idx in permutable:
            rp[idx]=rng.permutation(xd[idx])
        null[k]=corr(rp,yy)

    # Descriptive degree classes avoid empty quantile bins caused by the large
    # mass of plants with degree == 1.
    degree_classes=[]
    masks=[
        ("1 consumer", degree == 1),
        ("2 consumers", degree == 2),
        ("3–5 consumers", (degree >= 3) & (degree <= 5)),
        ("6+ consumers", degree >= 6),
    ]
    for label,mask in masks:
        idx=np.where(mask)[0]
        degree_classes.append({
            "class":label,
            "plants":int(len(idx)),
            "median_degree":float(np.median(degree[idx])),
            "expanded_fraction":float(np.mean(added[idx]>0)),
            "median_log_expansion":float(np.median(logexp[idx])),
            "median_added_units":float(np.median(added[idx])),
        })

    payload={
        "schema":"chocho_host_plant_prominence_expansion_v0.1",
        "status":"POSTHOC_PLANT_LEVEL_NETWORK_PROMINENCE_ANALYSIS",
        "plants":len(rows),
        "plant_families":len(set(fam)),
        "degree_definition":"Number of distinct Lepidoptera species in the fixed HOSTS-WCVP interaction reconstruction using the plant.",
        "associations":{
            "spearman_log_degree_vs_log_expansion":spearman(np.log1p(degree),logexp),
            "partial_rank_controlling_log_native_breadth":partial_logexp,
            "family_adjusted_rank_correlation":family_adj_logexp,
            "partial_rank_degree_vs_absolute_added_controlling_log_native_breadth":partial_added,
            "family_and_native_breadth_quintile_stratified_rank_correlation":observed_stratified,
            "stratified_permutation_p_two_sided":float((1+np.sum(np.abs(null)>=abs(observed_stratified)))/(len(null)+1)),
            "stratified_null_q025":float(np.quantile(null,0.025)),
            "stratified_null_median":float(np.median(null)),
            "stratified_null_q975":float(np.quantile(null,0.975)),
        },
        "degree_classes":degree_classes,
        "interpretation_rule":"A positive within-family/native-breadth association would show that plant species prominent across the Lepidoptera-host network are disproportionately anthropogenically redistributed, providing a plant-level explanation for why usage-weighted host nulls absorb the apparent butterfly host-identity excess.",
        "claim_boundary":"HOSTS consumer degree conflates ecological host breadth/commonness with study and recording intensity. This analysis cannot partition those components.",
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
