#!/usr/bin/env python3
"""Count predeclared real host-decline events in original Schulz et al 2020 data.

Pure feasibility gate: DOES NOT fit response models or estimate ecological effects.
Source is NOT the distinct DiLeo 2024 Dryad archive.
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

SOURCE_URL="https://zenodo.org/records/4987060/files/ECOG-04799.zip?download=1"
MD5="69122a1d82b1fb970fb6638b02da3db4"
MEMBER="data/survey_data.tsv"
MAX_BYTES=12_000_000
NEED=("year","patch","population","plantago","veronica")

def parse_num(value):
    v=str(value).strip()
    if v.lower() in ("","na","nan","null",".","none"):
        return None
    x=float(v)
    if not math.isfinite(x):
        raise ValueError("nonfinite score")
    return x

def ingest(tsv:bytes):
    rows={}
    all_years=Counter()
    missing=Counter()
    host_value_counts={"plantago":Counter(),"veronica":Counter()}
    duplicate=0
    invalid=Counter()
    raw=csv.DictReader(io.StringIO(tsv.decode("utf-8-sig")),delimiter="\t")
    if not set(NEED).issubset(raw.fieldnames or []):
        raise ValueError(f"missing required headers: {sorted(set(NEED)-set(raw.fieldnames or []))}")
    for n,r in enumerate(raw,2):
        patch=(r.get("patch") or "").strip()
        yr=parse_num(r["year"])
        if not patch or yr is None or yr!=int(yr):
            invalid["patch_or_year"]+=1
            continue
        year=int(yr)
        key=(patch,year)
        if key in rows:
            duplicate+=1
            continue
        row={"year":year,"patch":patch}
        for field in ["population","plantago","veronica","previous_population","grazing_presence","grazing_intensity","plantago_dry","veronica_dry"]:
            try:
                val=parse_num(r.get(field,""))
            except ValueError:
                invalid[field]+=1
                val=None
            row[field]=val
            if val is None:
                missing[field]+=1
        rows[key]=row
        all_years[year]+=1
        for name in host_value_counts:
            if row[name] is not None:host_value_counts[name][str(row[name])]+=1
    return rows,{"row_count":sum(all_years.values())+duplicate+invalid["patch_or_year"],
                 "unique_patch_years":len(rows),"unique_patches":len({k[0] for k in rows}),
                 "years":dict(sorted(all_years.items())),"duplicates":duplicate,
                 "missing":dict(missing),"invalid":dict(invalid),
                 "plant_value_counts":{k:dict(v) for k,v in host_value_counts.items()}}

def support(rows):
    scores=[r[z] for r in rows.values() for z in ("plantago","veronica") if r[z] is not None]
    valid_coding=bool(scores) and all(v in (0.,1.,2.,3.) for v in scores)
    status={
        "host_scores_ordinal_0to3":valid_coding,
        "three_year_consecutive_occupied_risk_set":0,
        "decline_event_counts":{},
        "data_support_pass":False,
        "ecological_effect_estimated":False
    }
    if not valid_coding:
        status["stop_reason"]="host abundance not confirmed 0-3 ordinal; no event grouping"
        return status
    groups=defaultdict(lambda:{"events":0,"patches":set(),"years":set(),"next_year_zero":0})
    for (patch,yr),r in rows.items():
        before=rows.get((patch,yr-1))
        after=rows.get((patch,yr+1))
        if before is None or after is None:continue
        if r["population"] is None or r["population"]<=0:continue
        if after["population"] is None:continue
        if any(q is None for q in (before["veronica"],r["veronica"],r["plantago"])):
            continue
        status["three_year_consecutive_occupied_risk_set"]+=1
        decline=before["veronica"]>r["veronica"] and before["veronica"]>=1
        h=r["plantago"]>=2
        key=("decline_" if decline else "nondecline_")+("high" if h else "low")
        g=groups[key]
        g["events"]+=1
        g["patches"].add(patch)
        g["years"].add(yr)
        g["next_year_zero"]+=int(after["population"]==0)
    for name,v in groups.items():
        status["decline_event_counts"][name]={"patch_years":v["events"],
                                              "patches":len(v["patches"]),
                                              "years":len(v["years"]),
                                              "next_year_zero_count":v["next_year_zero"]}
    a=status["decline_event_counts"].get("decline_high",{})
    b=status["decline_event_counts"].get("decline_low",{})
    gates={
        "40_events_each":min(a.get("patch_years",0),b.get("patch_years",0))>=40,
        "30_patches_each":min(a.get("patches",0),b.get("patches",0))>=30,
        "3_calendar_years_each":min(a.get("years",0),b.get("years",0))>=3,
        "20_next_year_zeros":sum(g.get("next_year_zero_count",0) for g in [a,b])>=20,
        "nondeclining_comparators":all(k in status["decline_event_counts"] for k in ("nondecline_high","nondecline_low")),
    }
    status["support_checks"]=gates
    status["data_support_pass"]=all(gates.values())
    return status

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--receipt",type=Path,required=True)
    args=parser.parse_args()
    req=Request(SOURCE_URL,headers={"User-Agent":"chocho-original-source-v02"})
    with urlopen(req,timeout=45) as response:
        raw=response.read(MAX_BYTES+1)
    if len(raw)>MAX_BYTES or hashlib.md5(raw).hexdigest()!=MD5:
        raise ValueError("original Zenodo ZIP missing or MD5 mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        if archive.namelist().count(MEMBER)!=1:
            raise ValueError("exact survey-data file missing")
        tsv=archive.read(MEMBER)
    rows,meta=ingest(tsv)
    counts=support(rows)
    result={
        "schema":"chocho_melitaea_19year_fallback_source_support_v01",
        "source":"Schulz et al 2020, independent Zenodo 4987060, not DiLeo 2024",
        "source_md5_verified":True,
        "source_sha256":hashlib.sha256(raw).hexdigest(),
        "survey_table_sha256":hashlib.sha256(tsv).hexdigest(),
        "raw_table":MEMBER,
        "metadata":meta,"support":counts,
        "model_run":False,"ecological_fallback_effect_estimated":False,
        "main_GEB_PR38_unchanged":True
    }
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
