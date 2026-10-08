#!/usr/bin/env python3
"""Source-verifiable already-published monarch adult-fate factorial audit.

This is an 8-cell descriptive extraction, not a new biological mechanism
or a population-level causal test of host globalization.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

REQUIRED = (
    "ID", "Plant_ID", "PlotNum", "Temp", "Milkweed",
    "OE_treatment", "Surv_adult", "Surv_pupa"
)


def original_source(protocol):
    doc=protocol["study"]
    sha=doc["pinned_commit"]
    path=doc["path"]
    url=f"https://raw.githubusercontent.com/{doc['original_repo']}/{sha}/{path}"
    req=urllib.request.Request(url,headers={
        "User-Agent":"chocho-butterfly-ecology-verification/0.1",
        "Accept":"text/csv,application/octet-stream,*/*",
    })
    with urllib.request.urlopen(req,timeout=50) as response:
        raw=response.read()
    if len(raw)<1000:raise RuntimeError("original source unexpectedly small")
    return raw,{"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)}


def binary(text,field):
    v=str(text or "").strip()
    if v in ("0","1"):return int(v)
    if v.upper() in ("NA","N/A",""):return None
    raise RuntimeError(f"Unrecognized {field} binary value: {v}")


def audit(raw,protocol):
    rs=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if not rs:raise RuntimeError("empty cohort")
    missing=set(REQUIRED)-set(rs[0])
    if missing:raise RuntimeError("missing source fields "+repr(sorted(missing)))
    expected=protocol["fixed_level_values"]
    all_ids=set()
    strata=defaultdict(list)
    missing_fates=Counter()
    for r in rs:
        uid=r["ID"].strip()
        if not uid or uid in all_ids:raise RuntimeError("Non-unique individual cohort ID")
        all_ids.add(uid)
        for k in ("Temp","Milkweed","OE_treatment"):
            if r[k].strip() not in expected[k]:
                raise RuntimeError(f"Unexpected experimental {k} level: {r[k]}")
        adult=binary(r["Surv_adult"],"Surv_adult")
        pupa=binary(r["Surv_pupa"],"Surv_pupa")
        if adult is None:missing_fates["adult"]+=1
        if pupa is None:missing_fates["pupa"]+=1
        if adult == 1 and pupa == 0:raise RuntimeError("adult survived without pupation")
        group=(r["Temp"].strip(),r["Milkweed"].strip(),r["OE_treatment"].strip())
        strata[group].append({"adult":adult,"pupa":pupa,
            "plant":r["Plant_ID"].strip(),"plot":r["PlotNum"].strip()})
    cells=[]
    lookup={}
    for temp in expected["Temp"]:
        for weed in expected["Milkweed"]:
            for parasite in expected["OE_treatment"]:
                key=(temp,weed,parasite)
                members=strata[key]
                if not members:raise RuntimeError("Empty factorial cell "+repr(key))
                valid=[v["adult"] for v in members if v["adult"] is not None]
                pupa=[v["pupa"] for v in members if v["pupa"] is not None]
                item={"Temp":temp,"Milkweed":weed,"OE_treatment":parasite,
                      "all_individuals":len(members),"measured_adult_fates":len(valid),
                      "adult_survivors":sum(valid),
                      "adult_eclosion_fraction":sum(valid)/len(valid) if valid else None,
                      "measured_pupation_fates":len(pupa),
                      "pupation_fraction":sum(pupa)/len(pupa) if pupa else None,
                      "distinct_plants":len({v["plant"] for v in members}),
                      "distinct_plots":len({v["plot"] for v in members})}
                cells.append(item)
                lookup[key]=item
    contrasts=[]
    for parasite in expected["OE_treatment"]:
        estimates=[]
        for temp in expected["Temp"]:
            tropical=lookup[(temp,"tropical",parasite)]
            swamp=lookup[(temp,"swamp",parasite)]
            x=tropical["adult_eclosion_fraction"]
            y=swamp["adult_eclosion_fraction"]
            if x is None or y is None:
                raise RuntimeError("Cannot evaluate observed cell contrast")
            estimates.append(x-y)
        contrasts.append({
            "parasite_exposure":parasite,
            "tropical_minus_swamp_adult_eclosion_ambient":estimates[0],
            "tropical_minus_swamp_adult_eclosion_elevated":estimates[1],
            "warming_shift_of_host_difference":estimates[1]-estimates[0],
            "inference":"SOURCE_DESCRIPTIVE_ALREADY_PUBLISHED_FACTORIAL_NOT_NEW_HYPOTHESIS_TEST"
        })
    return {"status":"VERIFIED_ORIGINAL_INDIVIDUAL_ADULT_FATE_DESCRIPTIVE_NOT_NEW_DISCOVERY",
            "original_records":len(rs),"distinct_ID":len(all_ids),
            "missing_fates":dict(missing_fates),
            "factorial_cells":cells,"descriptive_host_by_temperature_contrasts":contrasts,
            "limits":["Original experiment already analyzed host, thermal and parasite effects in 2025.",
              "Fates are individual larvae with repeated Plant_ID and PlotNum; no individual-row p values.",
              "Do not infer success on host introduced botanical regions or butterfly establishment.",
              "Adult eclosion is not multigeneration population growth or field colonization.",
              "Results are a transportability/calibration example, not a newly discovered demographic mechanism."
            ]}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--protocol-json",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    a=p.parse_args()
    protocol=json.loads(a.protocol_json.read_text(encoding="utf-8"))
    if protocol["schema"]!="chocho_monarch_warming_host_survival_calibration_v0.1":
        raise RuntimeError("Incorrect protocol")
    raw,prov=original_source(protocol)
    result=audit(raw,protocol)
    result["original_source"]=prov
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)


if __name__=="__main__":
    main()
