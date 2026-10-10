#!/usr/bin/env python3
"""Post-hoc plot-level uncertainty audit of ORIGINAL published monarch experiment.

NOT a new biological discovery, not independent study replication, and not
experimental evidence for heterospecific pathogen transfer. Original study is
Ragonese et al. (2025), DOI 10.1111/een.70010.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.request import Request, urlopen

COMMIT="74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2"
BLOB="e7a625911d8dc9b871a52d7657314a5d05472216"
URL=f"https://raw.githubusercontent.com/IRagonese/MilkweedWarming2024/{COMMIT}/MWwarming_comp_May16.csv"
ITERATIONS=9999
SEED=20261010
TEMP=("ambient","elevated")
HOST=("tropical","swamp")
OE=("infected","control")

def git_blob_sha(content:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()

def mean(a):
    return sum(a)/len(a)

def parse(raw):
    if git_blob_sha(raw)!=BLOB:
        raise ValueError("pinned original publication-source Git blob SHA mismatch")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(rows)!=240:
        raise ValueError("expected 240 original larvae")
    required={"ID","Plant_ID","PlotNum","Temp","Milkweed","OE_treatment","Surv_adult","Notes"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("missing source columns")
    ids=set()
    plants={}
    plots={}
    missing=[]
    for row in rows:
        id_=row["ID"]
        if id_ in ids or not id_:
            raise ValueError("duplicate larva ID")
        ids.add(id_)
        code=(row["Temp"],row["Milkweed"],row["OE_treatment"])
        if code[0] not in TEMP or code[1] not in HOST or code[2] not in OE:
            raise ValueError("invalid treatment")
        p=row["PlotNum"]
        plant=row["Plant_ID"]
        if not p or not plant:
            raise ValueError("missing plant or plot")
        value=row["Surv_adult"].strip()
        if value in ("0","1"):
            survived=int(value)
        elif value in ("NA",""):
            survived=None
            missing.append({"ID":id_,"Plant_ID":plant,"PlotNum":p,"Notes":row["Notes"]})
        else:
            raise ValueError("invalid adult survival")
        plots.setdefault(p,{"temperature":code[0],"hosts":{"tropical":[],"swamp":[]},"plants":set()})
        item=plots[p]
        if item["temperature"]!=code[0]:
            raise ValueError("plot assigned multiple temperatures")
        item["hosts"][code[1]].append(survived)
        item["plants"].add(plant)
        if plant not in plants:plants[plant]=[]
        plants[plant].append((p,code))
    if len(plots)!=30 or len(plants)!=120:
        raise ValueError("not 30 plots and 120 plants")
    if any(len(v)!=2 or len(set(v))!=1 for v in plants.values()):
        raise ValueError("plants not exactly 2 larvae with one treatment each")
    if any(len(p["plants"])!=4 or len(p["hosts"]["tropical"])!=4 or
           len(p["hosts"]["swamp"])!=4 for p in plots.values()):
        raise ValueError("incorrect plot structure")
    if len(missing)!=1 or missing[0]["ID"]!="22b":
        raise ValueError("original accidental-injury fate count/identity changed")
    return plots,missing

def lcg(seed):
    state=seed
    while True:
        state=(1664525*state+1013904223)&0xFFFFFFFF
        yield state/4294967296

def plot_differences(plots,injury_as_mortality=False):
    strata={t:[] for t in TEMP}
    for p in plots.values():
        def prop(host):
            vals=p["hosts"][host]
            denom=([0 if x is None else x for x in vals] if injury_as_mortality
                   else [x for x in vals if x is not None])
            if not denom:
                raise ValueError("no observed outcomes in host")
            return mean(denom)
        strata[p["temperature"]].append(prop("tropical")-prop("swamp"))
    if any(len(strata[t])!=15 for t in TEMP):
        raise ValueError("not 15 independent temperature plots per arm")
    return strata

def bootstrap_ci(strata,seed=SEED,iterations=ITERATIONS):
    random=lcg(seed)
    samples=[]
    for _ in range(iterations):
        values={}
        for t in TEMP:
            vals=strata[t]
            values[t]=mean([vals[int(next(random)*len(vals))] for _ in vals])
        samples.append(values["elevated"]-values["ambient"])
    samples.sort()
    return [samples[int(math.floor(iterations*.025))],
            samples[int(math.floor(iterations*.975))]]

def audit(raw):
    plots,missing=parse(raw)
    observations=plot_differences(plots,False)
    worst=plot_differences(plots,True)
    did=mean(observations["elevated"])-mean(observations["ambient"])
    worst_did=mean(worst["elevated"])-mean(worst["ambient"])
    counts={t:{"independent_plots":len(observations[t]),
               "mean_plot_level_tropical_minus_swamp":mean(observations[t])} for t in TEMP}
    return {
        "schema":"chocho_posthoc_monarch_plot_level_robustness_v01",
        "original_source":URL,
        "source_blob_sha1":BLOB,
        "source_sha256":hashlib.sha256(raw).hexdigest(),
        "published_paper":"Ragonese et al. (2025) DOI 10.1111/een.70010",
        "design":"30 plots, 15 per temperature; each plot contains 4 plants and 8 larvae; exactly 2 larvae per plant",
        "adult_outcome":"Surv_adult adult eclosion; not mate-ready/flight-capable recruitment",
        "original_assigned_larvae":240,
        "plots":30,"plants":120,
        "unknown_fates":missing,
        "known_fates":239,
        "known_adult_eclosions":sum(x for p in plots.values() for v in p["hosts"].values() for x in v if x is not None),
        "comparison":"(tropical-swamp elevated) - (tropical-swamp ambient), collapsing original OE randomization within plots",
        "plot_level":counts,
        "difference_in_differences":did,
        "bootstrap_plot_stratified_95ci":bootstrap_ci(observations),
        "bootstrap_draws":ITERATIONS,
        "bootstrap_seed":SEED,
        "injury_mortality_conservative_difference_in_differences":worst_did,
        "original_publication_already_analyzed_host_temperature_survival":True,
        "posthoc_descriptive_method":True,
        "can_prove_heterospecific_pathogen_transmission":False,
        "can_claim_generalized_plant_origin_causal_effect":False,
        "GEB_submission_PR38_unchanged":True
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--receipt",type=Path,required=True)
    args=parser.parse_args()
    with urlopen(Request(URL,headers={"User-Agent":"chocho-audit/1.0"}),timeout=30) as response:
        raw=response.read(500001)
    if len(raw)>500000:
        raise ValueError("source too large")
    result=audit(raw)
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
