#!/usr/bin/env python3
"""Exploratory early/peak/late matched caterpillar branch presence odds ratios.

Pinned original Whitener et al 2026 Caterpillars Count source with proper
site-year-family-method matched surveys, 2000 site-cluster bootstrap draws.
No native/alien "rescue" or butterfly population causal conclusions allowed.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path
import random

import numpy as np

from audit_caterpillars_count_original_branch_season_source_20261026 import (
    RAW_BLOB,ORIGIN_BLOB,EARLY,PEAK,LATE,EXCLUDE,
    git_blob,read_inferred,read_origin,safe_int,season,audit
)

ROOT=Path(__file__).resolve().parents[1]
GATE=ROOT/"provenance/reviewer_defenses/results/caterpillars_count_native_alien_season_source_gate_20261026.json"
LABELS=("early","peak","late")
B=2000
SEED=20261026


def source_events(raw_path,inferred,origin):
    surveys={}
    with raw_path.open(newline="",encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            ident=str(row.get("ID") or "").strip()
            if not ident:
                continue
            site=str(row.get("Name") or "").strip()
            year=safe_int(row.get("Year"))
            day=safe_int(row.get("julianday"))
            method=str(row.get("ObservationMethod") or "").strip()
            family=str(row.get("Family") or "").strip()
            sci=str(row.get("sciName") or "").strip()
            if str(row.get("Species") or "").strip()=="N/A":
                sci=inferred.get(str(row.get("PlantFK") or "").strip(),sci)
            if sci.lower() in {"","na","n/a","nan","none"}:sci=""
            try:lon=float(row.get("Longitude") or "nan")
            except ValueError:lon=float("nan")
            wet=str(row.get("WetLeaves") or "").lower().strip() in {
                "true","t","yes","y","1"}
            eligible=(all([site,family,sci,method])
                      and year is not None and 2010<=year<=2024
                      and day is not None and 1<=day<=366
                      and math.isfinite(lon) and lon>=-100
                      and not wet and site not in EXCLUDE)
            if not eligible:
                continue
            key=(site,year,family,method,season(day),origin.get(sci),sci)
            if ident not in surveys:
                surveys[ident]={"key":key,"positive":False}
            elif surveys[ident]["key"]!=key:
                raise ValueError("Original branch survey ID conflicts in botanical/site metadata")
            if str(row.get("Group") or "")=="caterpillar":
                try:amount=float(row.get("Quantity") or "nan")
                except ValueError:amount=float("nan")
                if not math.isfinite(amount) or amount<0:
                    raise ValueError("Original caterpillar quantity invalid")
                surveys[ident]["positive"] |= amount>0
    return surveys


def survey_to_matched(surveys):
    source=defaultdict(lambda:Counter())
    for x in surveys.values():
        site,year,family,method,period,ori,sci=x["key"]
        if period not in LABELS or ori not in {"alien","native"}:
            continue
        k=(site,year,family,method,period)
        source[k][ori+"_n"]+=1
        source[k][ori+"_pos"]+=int(x["positive"])
    strata=defaultdict(list)
    for (site,year,fam,method,period),table in source.items():
        if table["native_n"] and table["alien_n"]:
            n=table["native_n"]
            a=table["alien_n"]
            ap=table["alien_pos"]
            np_=table["native_pos"]
            if not (0<=ap<=a and 0<=np_<=n):
                raise ValueError("Malformed matched survey table")
            # a=alien+, b=alien-, c=native+, d=native-.
            strata[period].append({
                "site":site,
                "method":method,
                "family":fam,
                "cells":(ap,a-ap,np_,n-np_),
            })
    return strata


def counts(source):
    arr={}
    for period,records in source.items():
        x=Counter()
        for z in records:
            a,b,c,d=z["cells"]
            x["matched_strata"]+=1
            x["alien_n"]+=a+b
            x["native_n"]+=c+d
            x["alien_positive"]+=a
            x["native_positive"]+=c
        x["sites"]=len({s["site"] for s in records})
        x["families"]=len({s["family"] for s in records})
        arr[period]=dict(x)
    return arr


def mh_logor(records):
    ab=bc=0.0
    for rec in records:
        a,b,c,d=rec["cells"]
        n=a+b+c+d
        if n<=0:
            raise ValueError("Impossible empty matched survey table")
        ab+=a*d/n
        bc+=b*c/n
    if ab<=0 or bc<=0:
        return None
    return math.log(ab/bc)


def calculate(source,method=None):
    original={}
    for period in LABELS:
        sample=[s for s in source[period]
                if method is None or s["method"]==method]
        nsites=len({s["site"] for s in sample})
        a=sum(s["cells"][0]+s["cells"][1] for s in sample)
        n=sum(s["cells"][2]+s["cells"][3] for s in sample)
        pos=sum(s["cells"][0]+s["cells"][2] for s in sample)
        passes=len(sample)>=30 and nsites>=10 and a>=200 and n>=200 and pos>=10
        original[period]={
            "matched_strata":len(sample),"matched_sites":nsites,
            "native_surveys":n,"alien_surveys":a,
            "positive_larva_surveys":pos,"source_gate":passes,
            "log_oddsratio_alien_to_native":mh_logor(sample) if passes else None,
        }
    if not all(v["source_gate"] and v["log_oddsratio_alien_to_native"] is not None
               for v in original.values()):
        return {"status":"NOT_IDENTIFIABLE_SOME_SEASONS","observed":original}
    data={period:[r for r in source[period]
                  if method is None or r["method"]==method] for period in LABELS}
    logOR={s:original[s]["log_oddsratio_alien_to_native"] for s in LABELS}
    allsites=sorted({x["site"] for records in data.values() for x in records})
    bysite={p:defaultdict(list) for p in LABELS}
    for p in LABELS:
        for x in data[p]:
            bysite[p][x["site"]].append(x)
    rng=random.Random(SEED)
    b_draw={p:[] for p in LABELS}
    late_delta=[]
    early_delta=[]
    failed=Counter()
    for run in range(B):
        sampled=[allsites[rng.randrange(len(allsites))] for _ in allsites]
        values={}
        for period in LABELS:
            merged=[r for st in sampled for r in bysite[period][st]]
            values[period]=mh_logor(merged)
        if any(v is None for v in values.values()):
            failed["MH_any_zero_cell_crossproducts"]+=1
            continue
        for p in LABELS:b_draw[p].append(values[p])
        late_delta.append(values["late"]-values["peak"])
        early_delta.append(values["early"]-values["peak"])
    if len(late_delta)<1900:
        return {
            "status":"NOT_IDENTIFIABLE_SITE_CLUSTER_BOOTSTRAP_DEGENERACY",
            "observed":original,"bootstrap_success":len(late_delta),"failed":dict(failed)
        }
    def qq(x):
        return [float(v) for v in np.quantile(x,[.025,.5,.975],method="linear")]
    observed={
        p:{**original[p],
           "OR_alien_vs_native":float(math.exp(logOR[p])),
           "95pct_site_cluster_bootstrap_logOR":qq(b_draw[p]),
           "95pct_site_cluster_bootstrap_OR":[float(math.exp(k))
                                                for k in (qq(b_draw[p])[0],qq(b_draw[p])[2])]}
        for p in LABELS
    }
    ld=float(logOR["late"]-logOR["peak"])
    ed=float(logOR["early"]-logOR["peak"])
    return {
        "status":"ORIGINAL_SOURCE_MATCHED_SEASONAL_LARVAL_OCCUPANCY_ESTIMATE_QUALIFIED",
        "observed":observed,
        "late_minus_peak_logOR":ld,
        "early_minus_peak_logOR":ed,
        "late_minus_peak_OR_ratio":math.exp(ld),
        "early_minus_peak_OR_ratio":math.exp(ed),
        "late_minus_peak_site_cluster_95_logOR":qq(late_delta),
        "early_minus_peak_site_cluster_95_logOR":qq(early_delta),
        "late_minus_peak_positive_bootstrap_fraction":sum(v>0 for v in late_delta)/len(late_delta),
        "bootstrap_draws":B,
        "bootstrap_successes":len(late_delta),
        "bootstrap_source_cluster_sites":len(allsites),
        "bootstrap_degenerate":dict(failed),
        "source_frozen_years":"2010–2024",
        "field_data_biological_unit":"branch survey; positive unidentified caterpillar",
        "not_population_or_butterfly_species_abundance":True,
    }


def run(p):
    root=p
    source=root/"data/fullDataset_2025-04-18.csv"
    origin=root/"data/plant_origin_status.csv"
    inferred=root/"data/inferredPlantNames_2024-08-20.csv"
    source_g=json.loads(GATE.read_text(encoding="utf-8"))
    if source_g["status"]!="ALL_3_SEASONS_SOURCE_GATE_QUALIFIED":
        raise ValueError("Original independent source failed predeclared 3 season gate")
    git_blob(source,RAW_BLOB)
    git_blob(origin,ORIGIN_BLOB)
    surveys=source_events(source,read_inferred(inferred),read_origin(origin))
    matched=survey_to_matched(surveys)
    matched_count=counts(matched)
    for s in LABELS:
        r=source_g["matched_strata_by_season"][s]
        q=matched_count[s]
        for a,b in [
            (q["matched_strata"],r["matched_strata"]),
            (q["sites"],r["matched_distinct_sites"]),
            (q["native_n"],r["matched_native_surveys"]),
            (q["alien_n"],r["matched_alien_surveys"]),
            (q["native_positive"],r["matched_native_positive"]),
            (q["alien_positive"],r["matched_alien_positive"]),
        ]:
            if a!=b:
                raise ValueError(f"Independent source gate {s} does not reproduce: {a} != {b}")
    return {
        "schema":"chocho_CaterpillarsCount_native_alien_season_MH_site_bootstrap_v01",
        "status":"EXPLORATORY_SEASONAL_ACTUAL_CATERPILLAR_BRANCH_OCCUPANCY",
        "source_commit":"2c5690ac7c6dba7eca506da886c49781c37617bd",
        "source_study_doi":"10.1002/ecs2.70615",
        "original_source_verified_survey_peak":source_g["original_author_peak_DOY_152_212_survey_count"],
        "source_original_matched_strata":matched_count,
        "main_all_methods":calculate(matched),
        "method_visual_only_sensitivity":calculate(matched,"Visual"),
        "limitations":[
            "Branch arthropods classified as caterpillars include moths; not butterfly-specific host utilization",
            "Observational native/alien plant species identity, actual host quality and seasonal growth confounded",
            "Same site-year plant FAMILY and observation method matched; individual plant species not controlled",
            "Not an oviposition, survival, recruitment, population abundance or invasive host arrival estimate",
            "Raw native-vs-alien June–July difference already documented in Whitener et al. (2026)",
            "No source raw survey/site/location/observer/plant individual rows uploaded to repository",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=run(a.original_repo)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "main":result["main_all_methods"],
        "visual":result["method_visual_only_sensitivity"],
        "source_counts":result["source_original_matched_strata"],
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
