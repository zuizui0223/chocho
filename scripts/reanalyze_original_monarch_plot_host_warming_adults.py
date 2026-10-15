#!/usr/bin/env python3
"""Original monarch host-by-temperature adult eclosion, correctly paired by PLOT.

POSTHOC replication of a published 2025 finding, not a new ecological effect.
Original experimental date is 2021, 30 plots (15 ambient/15 elevated) with
4 plants/plot and 2 monarch larvae/plant: 240 original larvae. Infected and
control plants of each host occur in every plot, so neither plant, larva nor
OE-infected survivors constitute 240 independent temperature replicates.

Uncertainty: resample whole plots *within each temperature*, never larvae;
retain missing original adult outcomes as explicit worst-case successes/fails.
No ecological-trap, native-vs-exotic family, site-colonization or adult flight
inference, and no bootstrap "p-value" from already published host differences.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

ORIGINAL_COMMIT = "74e3e2cd8079401c9ccd4359f9e39e2d52d1f7d2"
ORIGINAL_BLOB = "e7a625911d8dc9b871a52d7657314a5d05472216"
ORIGINAL_URL = (
    "https://raw.githubusercontent.com/IRagonese/MilkweedWarming2024/"
    f"{ORIGINAL_COMMIT}/MWwarming_comp_May16.csv"
)
TEMP = ("ambient", "elevated")
HOST = ("tropical", "swamp")
INOC = ("infected", "control")
STAGE = ("Surv_pupa", "Surv_adult")
N_PLOTS_PER_TEMP = 15
LARVAE_PER_PLANT = 2
PLANTS_PER_PLOT = 4
N_ASSIGNED = 240


def _git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load_original(raw: bytes, *, require_blob: bool = True) -> tuple[list[dict], dict]:
    blob = _git_blob(raw)
    if require_blob and blob != ORIGINAL_BLOB:
        raise ValueError("Original author repository Git blob changed")
    records = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(records) != N_ASSIGNED or any(None in r for r in records):
        raise ValueError("Exactly 240 complete original source rows required")
    required = {
        "ID", "Plant_ID", "PlotNum", "Temp", "Milkweed", "OE_treatment",
        "Surv_adult", "Surv_pupa", "Lineage", "Notes",
    }
    if not required.issubset(records[0]) or any(not required.issubset(r) for r in records):
        raise ValueError("Author source outcome/plant/plot factor schema changed")
    seen_larvae = set()
    source_plants: dict[str, list[dict]] = defaultdict(list)
    source_plots: dict[str, list[dict]] = defaultdict(list)
    eight_cell = Counter()
    unknown = []
    by_stage = Counter()
    for row in records:
        cid, plant, plot = (
            row["ID"].strip(), row["Plant_ID"].strip(), row["PlotNum"].strip()
        )
        if not cid or cid in seen_larvae or not plant or not plot:
            raise ValueError("Missing/duplicate original individual, plant or plot identity")
        seen_larvae.add(cid)
        if (
            row["Temp"] not in TEMP or row["Milkweed"] not in HOST
            or row["OE_treatment"] not in INOC
        ):
            raise ValueError("Not an original balanced 2x2x2 treatment assignment")
        for field in STAGE:
            val = row[field].strip()
            if val not in {"0", "1", "NA", ""}:
                raise ValueError("Source pupil/adult indicator not documented binary or missing")
        a, p = row["Surv_adult"].strip(), row["Surv_pupa"].strip()
        if a in {"0","1"} and p in {"0","1"}:
            if a=="1" and p=="0":
                raise ValueError("Recorded original adult with no original pupation")
            by_stage["both_observed"] += 1
            by_stage["pupa_"+p] += 1
            by_stage["adult_"+a] += 1
            by_stage["pupa_recorded_adult_missing"] += 0
        else:
            by_stage["ambiguous_stage_fate"] += 1
            unknown.append({
                "original_larva_ID":cid, "original_Plant_ID":plant,
                "original_PlotNum":plot, "Temp":row["Temp"],
                "Milkweed":row["Milkweed"], "OE_treatment":row["OE_treatment"],
                "Surv_adult_original":a, "Surv_pupa_original":p,
                "Notes_original":row.get("Notes",""),
            })
        eight_cell[(row["Temp"],row["Milkweed"],row["OE_treatment"])] += 1
        source_plants[plant].append(row)
        source_plots[plot].append(row)
    if (
        len(source_plots)!=30 or len(source_plants)!=120 or len(seen_larvae)!=N_ASSIGNED
        or any(x!=30 for x in eight_cell.values()) or len(eight_cell)!=8
    ):
        raise ValueError("Original 30 plots, 120 plants or eight 30-larva cells mismatch")
    temp_count = Counter()
    for plot, rr in source_plots.items():
        if len(rr)!=8 or len({row["Plant_ID"] for row in rr})!=PLANTS_PER_PLOT:
            raise ValueError("Original plot is not four whole plants x two larvae")
        t={row["Temp"] for row in rr}
        if len(t)!=1:
            raise ValueError("Warming intervention recorded at PLOT level, not mixed within plot")
        t=next(iter(t))
        temp_count[t]+=1
        treatments=Counter((row["Milkweed"],row["OE_treatment"]) for row in rr)
        if set(treatments)!={(h,e) for h in HOST for e in INOC} or set(treatments.values())!={2}:
            raise ValueError("Plot missing one host x OE assignment of two larvae")
        if any(
            len(x)!=2 or len({y["PlotNum"] for y in x})!=1
            or len({(y["Temp"],y["Milkweed"],y["OE_treatment"]) for y in x})!=1
            for x in (source_plants[row["Plant_ID"]] for row in rr)
        ):
            raise ValueError("Original plant-level two-larva cluster not intact")
    if any(temp_count[t]!=N_PLOTS_PER_TEMP for t in TEMP):
        raise ValueError("Temperature requires exactly 15 original independent plots each")
    return records, {
        "author_original_git_blob_sha1":blob,
        "original_sha256":hashlib.sha256(raw).hexdigest(),
        "n_original_larvae":len(records),
        "n_distinct_experimental_plots":len(source_plots),
        "n_distinct_independent_host_plants":len(source_plants),
        "larvae_per_host_plant":LARVAE_PER_PLANT,
        "plants_per_plot":PLANTS_PER_PLOT,
        "plots_by_temperature":dict(temp_count),
        "original_status":dict(by_stage),
        "uncertain_original_larvae":unknown,
        "original_temperature_is_plot_assignment":True,
        "original_temperature_plot_randomization_documented":False,
        "source_mutation_or_posttreatment_selection":False,
    }


def plot_observed_adults(records: list[dict], *, unknown_success: bool) -> dict[str,dict]:
    by_plot: dict[str, dict] = defaultdict(lambda: {
        "temp":None, "counts":Counter(), "pupae":Counter(),
        "hosts":defaultdict(set),
    })
    for row in records:
        t,h,e=row["Temp"],row["Milkweed"],row["OE_treatment"]
        plot=row["PlotNum"]
        entry=by_plot[plot]
        if entry["temp"] is None:
            entry["temp"]=t
        if entry["temp"]!=t:
            raise ValueError("Contradictory plot-level warming labels")
        status=row["Surv_adult"].strip()
        if status=="1":
            is_adult=1
        elif status=="0":
            is_adult=0
        elif status in ("","NA"):
            is_adult=1 if unknown_success else 0
        else:
            raise ValueError("Undocumented adult source status")
        entry["counts"][h]+=is_adult
        entry["counts"][f"{h}|{e}"]+=is_adult
        entry["hosts"][h].add(row["Plant_ID"])
        if row["Surv_pupa"]=="1":
            entry["pupae"][h]+=1
    out={}
    for plot,v in by_plot.items():
        if set(v["counts"]) and not all(len(v["hosts"][h])==2 for h in HOST):
            raise ValueError("Plot missing original two host plants per host species")
        for h in HOST:
            if v["counts"][h] > 4:
                raise ValueError("More adults than originally assigned four larvae per host")
        out[plot] = {
            "plot_ID":plot, "temp":v["temp"],
            "tropical_adult_risk":v["counts"]["tropical"]/4,
            "swamp_adult_risk":v["counts"]["swamp"]/4,
            "tropical_minus_swamp":(v["counts"]["tropical"]-v["counts"]["swamp"])/4,
            "tropical_pupae":v["pupae"]["tropical"],
            "swamp_pupae":v["pupae"]["swamp"],
        }
    if len(out)!=30:
        raise ValueError("Need all 30 original plots without outcome restriction")
    return out


def _percentile(x: list[float], probability: float) -> float:
    if not x or not 0<=probability<=1:
        raise ValueError("Empty bootstrap or invalid percentile")
    ordered=sorted(x)
    index=probability*(len(ordered)-1)
    k=int(index)
    if k>=len(ordered)-1:
        return ordered[-1]
    w=index-k
    return ordered[k]*(1-w)+ordered[k+1]*w


def contrast_by_plot(plots: dict, *, draws:int=9999, seed:int=20261015) -> dict:
    if type(draws) is not int or not 256<=draws<=100000:
        raise ValueError("Bootstrap draws require 256..100000")
    rows=defaultdict(list)
    for x in plots.values():
        rows[x["temp"]].append(x["tropical_minus_swamp"])
    if any(len(rows[temp])!=15 for temp in TEMP):
        raise ValueError("Two actual plot-level temperature strata with 15 plots required")
    mean={key:statistics.mean(rows[key]) for key in TEMP}
    original=mean["elevated"]-mean["ambient"]
    rng=random.Random(seed)
    intervals={"interaction":[], "ambient_host_difference":[], "elevated_host_difference":[]}
    for _ in range(draws):
        amb=statistics.mean(rng.choices(rows["ambient"],k=15))
        elev=statistics.mean(rng.choices(rows["elevated"],k=15))
        intervals["interaction"].append(elev-amb)
        intervals["ambient_host_difference"].append(amb)
        intervals["elevated_host_difference"].append(elev)
    leave={}
    for plot,v in sorted(plots.items()):
        group=v["temp"]
        remaining=[o["tropical_minus_swamp"] for k,o in plots.items() if k!=plot and o["temp"]==group]
        if len(remaining)!=14:
            raise RuntimeError("One nonindependent plot incorrectly omitted")
        changed=statistics.mean(remaining)
        leave[plot]=(mean["elevated"]-changed if group=="ambient"
                     else changed-mean["ambient"])
    return {
        "observed_plot_paired_host_difference_ambient":mean["ambient"],
        "observed_plot_paired_host_difference_elevated":mean["elevated"],
        "observed_temperature_by_species_adult_emergence_risk_interaction":original,
        "stratified_whole_plot_bootstrap_95_intervals":{
            key:[_percentile(v,.025),_percentile(v,.975)]
            for key,v in intervals.items()
        },
        "bootstrap_resampling_unit":"original field plot, stratified by experimental ambient/elevated assignment (original temperature plot randomization not explicitly documented)",
        "source_temperature_plot_replicates_each":15,
        "draws":draws,"seed":seed,
        "leave_one_whole_plot_out_delta_min":min(leave.values()),
        "leave_one_whole_plot_out_delta_max":max(leave.values()),
        "leave_one_plot_out_values":leave,
        "bootstrap_interval_is_posthoc_not_original_publication_model":True,
    }


def run(raw:bytes, *, draws:int=9999, seed:int=20261015, require_blob=True) -> dict:
    records, provenance=load_original(raw,require_blob=require_blob)
    cell=defaultdict(lambda:{"n":0,"known":0,"adult":0,"pupa":0,"unknown":0,"plants":set(),"plots":set()})
    for r in records:
        k=f"{r['Temp']}|{r['Milkweed']}|{r['OE_treatment']}"
        x=cell[k]
        x["n"]+=1
        a,p=r["Surv_adult"].strip(),r["Surv_pupa"].strip()
        x["adult"]+=int(a=="1")
        x["pupa"]+=int(p=="1")
        x["known"]+=int(a in {"0","1"})
        x["unknown"]+=int(a not in {"0","1"})
        x["plants"].add(r["Plant_ID"])
        x["plots"].add(r["PlotNum"])
    groups={
        key:{
            "larvae_originally_assigned":v["n"],
            "known_adult_outcomes":v["known"],
            "adult_eclosions":v["adult"],
            "source_pupation_successes":v["pupa"],
            "uncertain_original_fates":v["unknown"],
            "whole_plants":len(v["plants"]),
            "warming_plots":len(v["plots"]),
        }
        for key,v in sorted(cell.items())
    }
    sensitivity={}
    for unknown_success,label in ((False,"missing_as_no_recorded_adult"),
                                  (True,"missing_as_one_possible_adult")):
        plots=plot_observed_adults(records,unknown_success=unknown_success)
        sensitivity[label]=contrast_by_plot(plots,draws=draws,seed=seed)
    return {
        "schema":"chocho_monarch_original_plotlevel_adult_eclosion_calibration_v01",
        "status":"PUBLISHED_PRIOR_ART_POSTHOC_DESCRIPTIVE_NO_NEW_BUTTERFLY_FINDING",
        "source":{
            "author_repository":"IRagonese/MilkweedWarming2024",
            "author_pinned_commit":ORIGINAL_COMMIT,
            "original_csv":"MWwarming_comp_May16.csv",
            "original_git_blob_sha1":provenance["author_original_git_blob_sha1"],
            "original_csv_sha256":provenance["original_sha256"],
            "paper_doi":"10.1111/een.70010",
            "original_field_year":2021,
            "original_study_region":"Athens, Georgia, USA",
        },
        "source_units":provenance,
        "source_original_eight_factorial_cells":groups,
        "pupae_vs_adult_recorded_count_difference":(
            provenance["original_status"].get("pupa_1",0) -
            provenance["original_status"].get("adult_1",0)
        ),
        "unknown_original_larvae_never_imputed_as_observed":True,
        "primary_original_cohort_plot_stratified_missing_fate_bounds":sensitivity,
        "botanical_origin_not_randomized_across_species":True,
        "species_by_temp_prior_art_original_2025_publication":True,
        "not_larval_survival_probability_by_planted_introduced_global_region":True,
        "not_flight_capable_adult_or_reproductive_recruitment":True,
        "no_go_for_independent_new_mechanism_paper_from_single_published_source":True,
        "GEB_main_and_PR38_unmodified":True,
    }


def download():
    with urlopen(
        Request(ORIGINAL_URL,headers={"User-Agent":"chocho-author-monarch-plot-eclosion-source-audit/1.0"}),
        timeout=40,
    ) as f:
        raw=f.read(200000)
    if not raw or len(raw)>=200000:
        raise ValueError("Original author monarch CSV empty or too large")
    return raw


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receipt",required=True,type=Path)
    ap.add_argument("--bootstrap-draws",type=int,default=9999)
    args=ap.parse_args()
    result=run(download(),draws=args.bootstrap_draws)
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
