#!/usr/bin/env python3
"""Audit exact 2023 Kyoto butterfly competition ORIGINAL experimental units.

Determines which variables are truly longitudinal and what causal contrasts
cannot be estimated. No ecological response models, no p-values, no mediation.
DOIs: 10.1002/ece3.10164, 10.6084/m9.figshare.23170898
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

FILES={
 "larval":{"url":"https://ndownloader.figshare.com/files/41146991","md5":"0e9ff924e4ae38cb78529ec540887c42"},
 "plant":{"url":"https://ndownloader.figshare.com/files/41146994","md5":"744bc591502233354ecf40dc9d7dc561"}
}

def fetch(name):
    source=FILES[name]
    req=Request(source["url"],headers={"User-Agent":"chocho-butterfly-experimental-unit-audit/1.0"})
    with urlopen(req,timeout=30) as response:
        content=response.read(500001)
    if len(content)>500000 or hashlib.md5(content).hexdigest()!=source["md5"]:
        raise ValueError(f"{name}: not MD5-identical to published original")
    return content

def num(value):
    v=str(value or "").strip()
    if v.lower() in ("","na","nan"):
        raise ValueError("missing numeric value in original design keys")
    return float(v)

def source_rows(raw):
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))

def audit(larval,plant):
    L=source_rows(larval)
    P=source_rows(plant)
    if len(L)!=690 or len(P)!=30:
        raise ValueError("original row count / paper-source mismatch")
    requiredL={"plot.no","treat.no","s.density","a.density","date","day","total.s","total.a","cumul.sp","cumul.ap"}
    requiredP={"plot.no","treat.no","s.density","a.density",*[f"defoliation.{i}" for i in range(1,5)]}
    if not requiredL.issubset(L[0]) or not requiredP.issubset(P[0]):
        raise ValueError("wrong original experimental columns")

    larva_plots=defaultdict(list)
    plant_plots={}
    for r in L:
        plot=str(r["plot.no"]).strip()
        if not plot:
            raise ValueError("empty cage identity")
        larva_plots[plot].append(r)
    for r in P:
        plot=str(r["plot.no"]).strip()
        if not plot or plot in plant_plots:
            raise ValueError("missing / duplicate plant cage identity")
        plant_plots[plot]=r
    if set(larva_plots)!=set(plant_plots):
        raise ValueError("plant and butterfly cages don't align")
    if len(larva_plots)!=30:
        raise ValueError("expected 30 independent cages")
    treatment_counter=Counter()
    time_count=Counter()
    density_count=Counter()
    pooled_initial=Counter()
    for cage,rows in larva_plots.items():
        ref=rows[0]
        sr,ar=int(num(ref["s.density"])),int(num(ref["a.density"]))
        if any(int(num(x["s.density"]))!=sr or int(num(x["a.density"]))!=ar for x in rows):
            raise ValueError("initial treatment changed within cage")
        pr=plant_plots[cage]
        if (int(num(pr["s.density"])),int(num(pr["a.density"])))!=(sr,ar):
            raise ValueError("plant and larval treatment disagrees")
        days=[num(r["day"]) for r in rows]
        if len(days)!=len(set(days)):
            raise ValueError("duplicated day within cage")
        if len(rows)<3:
            raise ValueError("no larval longitudinal support")
        time_count[len(rows)]+=1
        treat=(sr,ar)
        treatment_counter[treat]+=1
        density_count[sr+ar]+=1
        pooled_initial["Sericinus"]+=sr
        pooled_initial["Atrophaneura"]+=ar
        for i in range(1,5):
            # Four different potted plants measured at a shared FINAL endpoint.
            # Neither column suffix nor column order is a time/date index.
            num(pr[f"defoliation.{i}"])

    if len(treatment_counter)!=15 or any(c!=2 for c in treatment_counter.values()):
        raise ValueError("not exactly 15 density combinations x two cage replicates")
    if set(density_count)!={4,8,12}:
        raise ValueError("incorrect original assigned density totals")
    return {
       "schema":"chocho_kyoto_butterfly_identifiability_v01",
       "source":"Hashimoto & Ohgushi 2023, Ecology and Evolution DOI 10.1002/ece3.10164",
       "figshare_article":23170898,
       "source_files":{name:{"url":f["url"],"md5":f["md5"]} for name,f in FILES.items()},
       "n_larval_longitudinal_rows":len(L),
       "n_plant_final_rows":len(P),
       "n_independent_cages":len(larva_plots),
       "n_density_combinations":len(treatment_counter),
       "replicates_per_density_combination":2,
       "cage_larval_timepoint_distribution":dict(time_count),
       "cage_count_by_total_initial_larval_density":dict(density_count),
       "treatment_count_by_initial_species_counts":{
           f"S{sr}_A{ar}":n for (sr,ar),n in sorted(treatment_counter.items())
       },
       "cumulative_initial_larvae":dict(pooled_initial),
       "plant_variables":{
          "names":[f"defoliation.{i}" for i in range(1,5)],
          "meaning":"four separate potted plants within EACH cage measured at final endpoint; NOT FOUR TEMPORAL ASSESSMENTS",
          "n_temporal_plant_biomass_measurements":0,
          "independently_manipulated_host_regrowth_timing":False,
          "temporal_plant_quality_or_induced_defence_measured":False
       },
       "outcome_variables":{
          "larval_survival_and_stages":"Repeated larval stage counts, not individually followed adult offspring",
          "pupation":"cumul.sp and cumul.ap are stage cumulative pupation, not adult eclosion",
          "competitor_contact_or_behavior_measured":False
       },
       "identified_effects_from_original_randomization":{
          "randomized_initial_two_species_density":True,
          "randomized_plant_regrowth_timing":False,
          "randomized_donor_pretreatment_at_equal_defoliation":False,
          "randomized_competitor_physical_contact_vs_plant_tissue_exposure":False
       },
       "new_plant_quality_mechanism_estimated":False,
       "GEB_PR38_unchanged":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--receipt",type=Path,required=True)
    args=ap.parse_args()
    receipt=audit(fetch("larval"),fetch("plant"))
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
