#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def average_ranks(values):
    x=np.asarray(values,dtype=float)
    order=np.argsort(x,kind="mergesort")
    ranks=np.empty(len(x),dtype=float)
    start=0
    while start<len(order):
        stop=start+1
        while stop<len(order) and x[order[stop]]==x[order[start]]:
            stop+=1
        ranks[order[start:stop]]=0.5*((start+1)+stop)
        start=stop
    return ranks


def spearman(a,b):
    x=average_ranks(a); y=average_ranks(b)
    if len(x)<3 or np.std(x)==0 or np.std(y)==0:
        return None
    return float(np.corrcoef(x,y)[0,1])


def bootstrap_ci(rows, reps=49999, seed=20260929):
    if len(rows)<3:
        return None
    x=np.asarray([r["host_family_count"] for r in rows],float)
    y=np.asarray([r["log_resource_expansion"] for r in rows],float)
    rng=np.random.default_rng(seed)
    vals=np.empty(reps,float)
    for i in range(reps):
        idx=rng.integers(0,len(rows),size=len(rows))
        vals[i]=spearman(x[idx],y[idx])
    return {
        "replicates":reps,
        "seed":seed,
        "ci90":np.quantile(vals,[0.05,0.95]).tolist(),
        "ci95":np.quantile(vals,[0.025,0.975]).tolist(),
    }


def summary(rows):
    rho=spearman(
        [r["host_family_count"] for r in rows],
        [r["log_resource_expansion"] for r in rows],
    )
    return {
        "species":len(rows),
        "expanded_species":sum(r["introduced_added_units"]>0 for r in rows),
        "median_log_resource_expansion":float(np.median([r["log_resource_expansion"] for r in rows])),
        "median_added_units":float(np.median([r["introduced_added_units"] for r in rows])),
        "spearman_host_family_vs_log_expansion":rho,
        "bootstrap":bootstrap_ci(rows),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--anthropogenic-csv",type=Path,required=True)
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    rows=[]
    with args.anthropogenic_csv.open(newline="",encoding="utf-8") as handle:
        reader=csv.DictReader(handle)
        for r in reader:
            rows.append({
                "species":r["species"].strip(),
                "host_family_count":float(r["host_family_count"]),
                "log_resource_expansion":float(r["log_resource_expansion"]),
                "introduced_added_units":int(r["introduced_added_units"]),
            })
    if len(rows)!=239:
        raise RuntimeError(f"expected 239 resource-eligible species; got {len(rows)}")

    families=defaultdict(set)
    with args.insect_host_csv.open(newline="",encoding="utf-8") as handle:
        reader=csv.DictReader(handle)
        required={"insect_species","family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for r in reader:
            insect=r["insect_species"].strip()
            family=r["family"].strip()
            if insect and family:
                families[insect].add(family)

    annotated=[]
    for r in rows:
        fams=families.get(r["species"],set())
        if not fams:
            raise RuntimeError(f"no resolved host family for {r['species']}")
        q=dict(r)
        q["uses_poaceae"]="Poaceae" in fams
        q["resolved_poaceae_only"]=(fams=={"Poaceae"})
        q["family_level_specialist_poaceae"]=(
            q["resolved_poaceae_only"] and q["host_family_count"]==1
        )
        annotated.append(q)

    any_poaceae=[r for r in annotated if r["uses_poaceae"]]
    non_poaceae=[r for r in annotated if not r["uses_poaceae"]]
    poaceae_only=[r for r in annotated if r["family_level_specialist_poaceae"]]
    without_poaceae_only=[r for r in annotated if not r["family_level_specialist_poaceae"]]

    payload={
        "schema":"chocho_butterfly_poaceae_sensitivity_v0.1",
        "status":"SUCCESS_POSTHOC_POACEAE_HOST_GUILD_SENSITIVITY",
        "definition":{
            "any_poaceae_user":"At least one resolved HOSTS-WCVP host species belongs to Poaceae.",
            "family_level_specialist_poaceae":"All resolved host families are Poaceae and LepTraits host_family_count equals 1.",
        },
        "full_panel":summary(annotated),
        "poaceae_users":summary(any_poaceae),
        "non_poaceae_users":summary(non_poaceae),
        "excluding_any_poaceae_user":summary(non_poaceae),
        "family_level_poaceae_specialists":summary(poaceae_only),
        "excluding_family_level_poaceae_specialists":summary(without_poaceae_only),
        "claim_boundary":"This sensitivity asks whether the near-zero full-panel diet-breadth association is driven by butterflies using Poaceae, especially one-family Poaceae specialists. It does not isolate plant-use causality or infer that Poaceae redistribution caused butterfly range expansion.",
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
