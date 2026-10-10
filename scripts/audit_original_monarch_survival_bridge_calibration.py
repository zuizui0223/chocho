#!/usr/bin/env python3
"""Audit original monarch host×temperature×OE adult-emergence data.

Descriptive/source-quality ONLY. Original effects were already published in
Ragonese et al. 2025, doi:10.1111/een.70010. Neither butterfly cross-species
transmission nor native/exotic randomized host origins are tested.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

ORIGINAL_COMMIT="74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2"
ORIGINAL_GIT_BLOB_SHA1="e7a625911d8dc9b871a52d7657314a5d05472216"
URL=f"https://raw.githubusercontent.com/IRagonese/MilkweedWarming2024/{ORIGINAL_COMMIT}/MWwarming_comp_May16.csv"
GROUPS=["ambient","elevated"]
PLANTS=["tropical","swamp"]
TREATMENTS=["infected","control"]

def git_sha1(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def audit(raw):
    sha=git_sha1(raw)
    if sha!=ORIGINAL_GIT_BLOB_SHA1:
        raise ValueError("Original scientific source blob mismatch")
    records=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(records)!=240:
        raise ValueError("Source has unexpectedly changed size")
    required={"ID","Plant_ID","Temp","OE_treatment","Milkweed","Surv_adult","Surv_pupa","Notes"}
    if not required.issubset(records[0]):
        raise ValueError("Missing expected source variables")
    ids=set()
    groups=defaultdict(lambda:{"assigned":0,"known_outcomes":0,"adult":0,"pupa":0,"unknown":0,"plants":set()})
    plants=defaultdict(list)
    missing=[]
    for r in records:
        rid=r["ID"]
        if not rid or rid in ids:
            raise ValueError("Duplicate/missing monarch individual ID")
        ids.add(rid)
        key=(r["Temp"],r["Milkweed"],r["OE_treatment"])
        if key[0] not in GROUPS or key[1] not in PLANTS or key[2] not in TREATMENTS:
            raise ValueError("Unknown treatment code")
        cl=groups[key]
        cl["assigned"]+=1
        plant_id=r["Plant_ID"]
        plants[plant_id].append(key)
        cl["plants"].add(plant_id)
        adult=r["Surv_adult"].strip()
        pupa=r["Surv_pupa"].strip()
        if adult not in {"0","1"} or pupa not in {"0","1"}:
            missing.append({"ID":rid,"Plant_ID":plant_id,"condition":" / ".join(key),
                            "Surv_adult":adult,"Surv_pupa":pupa,"Notes":r["Notes"]})
            cl["unknown"]+=1
            continue
        if adult=="1" and pupa!="1":
            raise ValueError("Adult without pupation")
        cl["known_outcomes"]+=1
        cl["adult"]+=int(adult)
        cl["pupa"]+=int(pupa)
    if len(ids)!=240 or len(plants)!=120 or any(len(x)!=2 or len(set(x))!=1 for x in plants.values()):
        raise ValueError("Experimental units not exactly 2 larvae per plant")
    output={}
    for temp in GROUPS:
        for plant in PLANTS:
            for treatment in TREATMENTS:
                k=(temp,plant,treatment)
                x=groups[k]
                if x["assigned"]!=30 or len(x["plants"])!=15:
                    raise ValueError("Factorial cell unbalanced from original")
                label=" / ".join(k)
                output[label]={name:value for name,value in x.items() if name!="plants"}
                output[label]["distinct_plants"]=len(x["plants"])
                output[label]["known_adult_fraction"]=x["adult"]/x["known_outcomes"] if x["known_outcomes"] else None
    return {
        "schema":"chocho_original_monarch_survival_calibration_v01",
        "source_url":URL,
        "source_git_sha1":sha,
        "source_sha256":hashlib.sha256(raw).hexdigest(),
        "source_publication":"Ragonese et al. 2025 Ecological Entomology doi:10.1111/een.70010",
        "n_assigned":240,
        "unique_individuals":len(ids),
        "unique_plants":len(plants),
        "larvae_per_plant":2,
        "n_known_survival":240-len(missing),
        "adult_eclosions":sum(x["adult"] for x in groups.values()),
        "uncertain_fates":missing,
        "groups":output,
        "observational_summaries_only":True,
        "cross_species_pathogen_bridge_tested":False,
        "plant_origin_heritable_hysteresis_tested":False,
        "GEB_submission_PR38_untouched":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    req=Request(URL,headers={"User-Agent":"chocho-original-monarch-source-check/1.0"})
    with urlopen(req,timeout=25) as response:
        raw=response.read(500000)
    result=audit(raw)
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
