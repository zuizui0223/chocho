#!/usr/bin/env python3
"""Source-locked year-aggregated GBIF butterfly observation-effort thinning.

This is a source-conditional observation-intensity sensitivity experiment, NOT a
detection-corrected occupancy analysis. Neither organism-level presences nor
new host arrivals/colonisations are inferred from GBIF records.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path

SOURCE_SHA = {
    "cells": "f1bd777e2faffdeb0d689b93dca193e8612385fa264dd3bbdb6e104c075f6525",
    "effort": "310f8f58ff1b4bedc9bf3f9b39b85c4fb11248916f2b13793795af7ec8925e82",
}
PIN = "dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0"
N_DRAW = 499
SEED = 20261119
RAW_J = {"prior": 0.11421295035628568, "later": 0.16788018941797136}
N_REGIONS, N_PAIRS, N_SPECIES = 169, 14196, 224
PERIODS = {"prior": (2001, 2010), "later": (2018, 2024)}

def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda:f.read(1<<20), b""):
            h.update(c)
    return h.hexdigest()

def qtile(values,p):
    s=sorted(values)
    pos=(len(s)-1)*p
    lo=math.floor(pos)
    hi=math.ceil(pos)
    return s[lo]+(pos-lo)*(s[hi]-s[lo])

def describe(x):
    if not x:
        raise RuntimeError("Empty distribution")
    return {
        "min":float(min(x)),
        "median":float(qtile(x,.5)),
        "q025":float(qtile(x,.025)),
        "q975":float(qtile(x,.975)),
        "max":float(max(x)),
        "mean":float(sum(x)/len(x)),
    }

def period(year):
    for p,(a,b) in PERIODS.items():
        if a<=year<=b:
            return p
    return None

def source_regions(path):
    data={}
    with path.open(newline="",encoding="utf8") as f:
        for r in csv.DictReader(f):
            if r["model"]!="M0_all_flora_effort_prior":
                continue
            code=r["wgsrpd3_code"]
            if code in data:
                raise RuntimeError(f"Duplicate original held-region source: {code}")
            data[code]=r["world_division"]
    if len(data)!=N_REGIONS:
        raise RuntimeError("Original 169-region matched GBIF sample not found")
    return data

def load_sources(a):
    for k,p in (("cells",a.cells),("effort",a.effort)):
        if sha256(p)!=SOURCE_SHA[k]:
            raise RuntimeError(f"Frozen original GBIF annual {k} SHA256 differs")
    original=json.loads(a.audit.read_text(encoding="utf8"))
    if (original["records"]["rows"]!=14637152 or
        original["coverage"]["positive_species_region_year_cells"]!=91882 or
        original["coverage"]["region_year_cells_with_any_candidate_record"]!=5390 or
        original["coverage"]["candidate_species_matched_to_a_369_region"]!=224):
        raise RuntimeError("Original 14.6M GBIF source extraction audit changed")
    for key,name in (("cells",a.cells.name),("effort",a.effort.name)):
        if original["files"][name]["sha256"]!=SOURCE_SHA[key]:
            raise RuntimeError("Source audit filename hash mismatch")
    codes=source_regions(a.region_scores)
    effort={}
    annual_count={}
    n_effort=0
    with a.effort.open(encoding="utf8",newline="") as f:
        for r in csv.DictReader(f):
            n_effort+=1
            code=r["wgsrpd3_code"]
            year=int(r["year"])
            if not 2001<=year<=2024:
                raise RuntimeError("Source effort year outside 2001–2024")
            k=(code,year)
            if k in effort:
                raise RuntimeError("Duplicate source annual effort record")
            d=int(r["distinct_gbif_datasetkey_calendar_days"])
            n=int(r["indexed_candidate_butterfly_records"])
            sp=int(r["indexed_candidate_butterfly_species"])
            if min(d,n,sp)<0:
                raise RuntimeError("Negative source observation")
            effort[k]={"days":d,"records":n,"taxa":sp}
    if n_effort!=5390:
        raise RuntimeError("Original 5,390 source regional effort cells changed")

    all_species=set()
    original_row_count=0
    original_mapped_count=0
    seen_cell=set()
    group_record_count=defaultdict(int)
    group_species_count=defaultdict(int)
    species_period_counts=defaultdict(int)
    for r in _iter_gz_rows(a.cells):
        original_row_count+=1
        code=r["wgsrpd3_code"]
        year=int(r["year"])
        insect=r["butterfly_species"].strip()
        n=int(r["indexed_gbif_records"])
        if n<=0 or not insect:
            raise RuntimeError("Original source has zero or unnamed positive species cell")
        k=(insect,code,year)
        if k in seen_cell:
            raise RuntimeError("Duplicated source butterfly-region-year cell")
        seen_cell.add(k)
        original_mapped_count+=n
        all_species.add(insect)
        group_record_count[(code,year)]+=n
        group_species_count[(code,year)]+=1
        p=period(year)
        if p is not None and code in codes:
            species_period_counts[(code,p,insect)]+=n
    if (original_row_count!=91882 or original_mapped_count!=13587452 or
        len(all_species)!=N_SPECIES or len(group_record_count)!=n_effort):
        raise RuntimeError("Annual species data cannot reproduce original GBIF source records")
    for k,n in group_record_count.items():
        if k not in effort or effort[k]["records"]!=n or effort[k]["taxa"]!=group_species_count[k]:
            raise RuntimeError(f"Original source yearly effort/species aggregate mismatch {k}")

    period_budget={c:{"prior":{"days":0,"records":0},"later":{"days":0,"records":0}}
                   for c in codes}
    for (c,y),r in effort.items():
        p=period(y)
        if c in codes and p is not None:
            period_budget[c][p]["days"]+=r["days"]
            period_budget[c][p]["records"]+=r["records"]
    if any(period_budget[c]["prior"]["days"]<10 or
           period_budget[c]["later"]["days"]<20 for c in codes):
        raise RuntimeError("Original primary 169-region recording-eligibility gate drift")
    return codes,sorted(all_species),period_budget,species_period_counts,{
        "original_yearly_species_rows":original_row_count,
        "original_year_effort_rows":n_effort,
        "original_positive_butterfly_species":len(all_species),
        "original_mapped_candidate_records":original_mapped_count,
        "cells_sha256":SOURCE_SHA["cells"],"effort_sha256":SOURCE_SHA["effort"],
        "original_GBIF_download":original["source"]["gbif_doi"],
        "source_PIN":PIN,
    }

def _iter_gz_rows(p):
    with gzip.open(p,"rt",encoding="utf8",newline="") as f:
        yield from csv.DictReader(f)

def construct(codes,plants,budgets,sp):
    region_order=sorted(codes)
    names=sorted(plants)
    ix={s:i for i,s in enumerate(names)}
    per_region={c:{"prior":[],"later":[]} for c in region_order}
    for (c,p,s),n in sp.items():
        per_region[c][p].append((ix[s],n))
    # Multiple years already summed at species×region×period grain.
    world=[codes[c] for c in region_order]
    if len(ix)!=224:
        raise RuntimeError("Source 224 names not recovered")
    for c in region_order:
        for p in PERIODS:
            if not per_region[c][p]:
                raise RuntimeError(f"Unexpected empty originally observed species list {c}:{p}")
    return region_order,per_region,world

def matrices(region_order,per_region,world):
    pair_info=[]
    for i in range(len(region_order)):
        for j in range(i+1,len(region_order)):
            pair_info.append((i,j,world[i]==world[j]))
    if len(pair_info)!=N_PAIRS:
        raise RuntimeError("Original 14,196 pairs differ")
    return pair_info

def original_masks(region_order,per_region):
    masks={}
    for p in PERIODS:
        a=[]
        for c in region_order:
            bits=0
            for i,n in per_region[c][p]:
                bits|=1<<i
            a.append(bits)
        masks[p]=a
    return masks

def jaccard_summary(masks, pair_info):
    all_sum,all_nonempty=0.,0
    within_sum,within_nonempty,within_total=0.,0,0
    between_sum,between_nonempty,between_total=0.,0,0
    undefined=0
    for i,j,same in pair_info:
        a,b=masks[i],masks[j]
        union=(a|b).bit_count()
        if same:
            within_total+=1
        else:
            between_total+=1
        if not union:
            undefined+=1
            continue
        value=(a&b).bit_count()/union
        all_sum+=value
        all_nonempty+=1
        if same:
            within_sum+=value
            within_nonempty+=1
        else:
            between_sum+=value
            between_nonempty+=1
    if not all_nonempty or not within_nonempty or not between_nonempty:
        raise RuntimeError("Original regional Jaccard is undefined throughout a group")
    return {
        "all_defined_Jaccard_mean":all_sum/all_nonempty,
        "all_fixed_zero_fill_Jaccard_mean":all_sum/len(pair_info),
        "within_world_Jaccard_mean":within_sum/within_nonempty,
        "between_world_Jaccard_mean":between_sum/between_nonempty,
        "defined_pairs":all_nonempty,
        "both_empty_pairs":undefined,
        "within_world_pairs":within_total,
        "between_world_pairs":between_total,
        "empty_regions":sum(x==0 for x in masks),
        "mean_indexed_species_regions":sum(m.bit_count() for m in masks)/len(masks),
    }

def keep_probability(q,n):
    if q<0 or q>1 or n<1:
        raise RuntimeError("Invalid source intensity/sample rate")
    if q==1:
        return 1.
    if q==0:
        return 0.
    return -math.expm1(n*math.log1p(-q))

def make_sampling_vectors(region_order,per_region,budgets,proxy):
    rates={}
    vectors={}
    for c in region_order:
        e0=budgets[c]["prior"][proxy]
        e1=budgets[c]["later"][proxy]
        if min(e0,e1)<=0:
            raise RuntimeError(f"Source has nonpositive effort for {c}/{proxy}")
        target=min(e0,e1)
        for p,value in (("prior",e0),("later",e1)):
            q=target/value
            if not (0<q<=1):
                raise RuntimeError("Impossible thinning rate")
            rates[(c,p)]=q
            vectors[(c,p)]=[(i,keep_probability(q,n)) for i,n in per_region[c][p]]
    return rates,vectors

def simulate(region_order,vec,rng):
    masks={"prior":[],"later":[]}
    for p in PERIODS:
        for c in region_order:
            b=0
            for species,prob in vec[(c,p)]:
                if prob==1 or rng.random()<prob:
                    b|=1<<species
            masks[p].append(b)
    return masks

def summarize_proxy(kind,region_order,per_region,budgets,pairs,seed):
    rates,vec=make_sampling_vectors(region_order,per_region,budgets,kind)
    rng=random.Random(seed)
    raw={}
    for p in PERIODS:
        raw[p]=jaccard_summary(original_masks(region_order,per_region)[p],pairs)
    simulated={key:[] for key in (
        "prior_mean","later_mean","late_minus_prior","prior_fixedzero","later_fixedzero",
        "late_minus_prior_fixedzero",
        "late_minus_prior_within_world","late_minus_prior_between_world",
        "prior_empty_pairs","later_empty_pairs","prior_empty_regions","later_empty_regions",
        "prior_mean_richness","later_mean_richness")}
    for _ in range(N_DRAW):
        b=simulate(region_order,vec,rng)
        stats={p:jaccard_summary(b[p],pairs) for p in PERIODS}
        for p in PERIODS:
            simulated[p+"_mean"].append(stats[p]["all_defined_Jaccard_mean"])
            simulated[p+"_fixedzero"].append(stats[p]["all_fixed_zero_fill_Jaccard_mean"])
            simulated[p+"_empty_pairs"].append(stats[p]["both_empty_pairs"])
            simulated[p+"_empty_regions"].append(stats[p]["empty_regions"])
            simulated[p+"_mean_richness"].append(stats[p]["mean_indexed_species_regions"])
        simulated["late_minus_prior"].append(stats["later"]["all_defined_Jaccard_mean"]-stats["prior"]["all_defined_Jaccard_mean"])
        simulated["late_minus_prior_fixedzero"].append(stats["later"]["all_fixed_zero_fill_Jaccard_mean"]-stats["prior"]["all_fixed_zero_fill_Jaccard_mean"])
        for part in ("within_world","between_world"):
            simulated["late_minus_prior_"+part].append(
                stats["later"][part+"_Jaccard_mean"]-stats["prior"][part+"_Jaccard_mean"])
    q0=[rates[(c,"prior")] for c in region_order]
    q1=[rates[(c,"later")] for c in region_order]
    ratio=[budgets[c]["later"][kind]/budgets[c]["prior"][kind] for c in region_order]
    return {
        "proxy":kind,
        "common_region_budget":"minimum of observed prior and later period proxy, fixed within a region",
        "n_monte_carlo_draws":N_DRAW,
        "fixed_rng_seed":seed,
        "n_prior_period_proxy_downscaled_regions":sum(q<1 for q in q0),
        "n_later_period_proxy_downscaled_regions":sum(q<1 for q in q1),
        "prior_record_retention_q":describe(q0),
        "later_record_retention_q":describe(q1),
        "later_over_prior_total_proxy_ratio":describe(ratio),
        "period_total_proxy":{
            p:{"sum":int(sum(budgets[c][p][kind] for c in region_order)),
               "median_per_region":float(qtile([budgets[c][p][kind] for c in region_order],.5))}
            for p in PERIODS
        },
        "conditional_monte_carlo_distributions":{k:describe(v) for k,v in simulated.items()},
        "fraction_draws_late_minus_prior_positive":sum(x>0 for x in simulated["late_minus_prior"])/N_DRAW,
        "limitation":"Period-specific per-species total indexed occurrence counts were independently thinned with probability derived from a region-level effort proxy, not a measured sampling visit."
    }

def main():
    ap=argparse.ArgumentParser()
    for name in ("cells","effort","audit","region-scores","original-readout","out"):
        ap.add_argument("--"+name,required=True,type=Path)
    a=ap.parse_args()
    codes,all_species,budgets,spec,audit=load_sources(a)
    reg,per_region,world=construct(codes,all_species,budgets,spec)
    pairs=matrices(reg,per_region,world)
    original={p:jaccard_summary(m,pairs) for p,m in original_masks(reg,per_region).items()}
    prior,late=original["prior"],original["later"]
    target=json.loads(a.original_readout.read_text(encoding="utf8"))
    if target["cohort"]["region_count"]!=169 or target["cohort"]["global_source_mapped_candidate_species"]!=224:
        raise RuntimeError("Independent prior/later source cohort identity drift")
    for p,field in (("prior","prior_2001_2010_mean"),("later","late_2018_2024_mean")):
        target_val=target["outcome_Jaccard_summary"][field]
        if abs(original[p]["all_defined_Jaccard_mean"]-RAW_J[p])>1e-10 or abs(RAW_J[p]-target_val)>1e-10:
            raise RuntimeError(f"Original GBIF observed butterfly source Jaccard did not reproduce: {p}")
        if original[p]["both_empty_pairs"]!=0:
            raise RuntimeError("Original sample has unexpected undefined pair Jaccard")
    observed_delta=late["all_defined_Jaccard_mean"]-prior["all_defined_Jaccard_mean"]
    tests={
        "dataset_key_calendar_day":summarize_proxy("days",reg,per_region,budgets,pairs,SEED),
        "indexed_candidate_record_count":summarize_proxy("records",reg,per_region,budgets,pairs,SEED+1),
    }
    result={
        "schema":"chocho_224_GBIF_original_annual_region_record_intensity_matched_thinning_20261119_v1",
        "status":"ORIGINAL_SOURCE_MATCHED_EXPLORATORY_RECORDING_INTENSITY_SENSITIVITY_NOT_DETECTION_CORRECTION",
        "source":audit,
        "cohort":{"species":N_SPECIES,"regions":N_REGIONS,"region_pairs":N_PAIRS,
                  "period_years":PERIODS,
                  "regions_by_world_division":{w:world.count(w) for w in sorted(set(world))}},
        "original_observed_recorded_community":original,
        "original_late_minus_prior_Jaccard":observed_delta,
        "tests":tests,
        "interpretation_limits":[
            "All source records are positively indexed HUMAN_OBSERVATION, not standardized visits; zero counts do not establish absences.",
            "DatasetKey×calendar-day is response-dependent recorded activity, not an independent complete survey-effort measure.",
            "Retaining individual indexed records with period-specific probability assumes spatial/temporal exchangeability within region/period, and species-specific observer bias may violate it.",
            "Prior (10 calendar years) and later (7 calendar years) durations differ; proxy total equalisation does not equalize years or contemporaneous ecology.",
            "Simulated variability quantifies randomized thinning conditional on indexed records, not ecological sampling uncertainty or iid 95% confidence intervals.",
            "Source 169 regions selected by prior/later GBIF effort, so this is descriptive conditional on the original panel.",
            "Taxonomy, host candidate availability and timing of alien introduction were not modeled; no butterfly recruitment, genuine homogenization or colonization mechanism is identified."
        ],
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    print("SOURCE_AND_RAW_VERIFIED",json.dumps({"regions":len(reg),
        "raw_prior":prior["all_defined_Jaccard_mean"],
        "raw_later":late["all_defined_Jaccard_mean"],
        "raw_gap":observed_delta,
        "cells_sha256":audit["cells_sha256"],"effort_sha256":audit["effort_sha256"]}))
    for key,v in tests.items():
        print("EFFORT_THINNING",json.dumps({
            "proxy":key,
            "prior_downscaled":v["n_prior_period_proxy_downscaled_regions"],
            "later_downscaled":v["n_later_period_proxy_downscaled_regions"],
            "median_q_later":v["later_record_retention_q"]["median"],
            "delta":v["conditional_monte_carlo_distributions"]["late_minus_prior"],
            "zero_filled_delta":v["conditional_monte_carlo_distributions"]["late_minus_prior_fixedzero"],
            "positive_fraction":v["fraction_draws_late_minus_prior_positive"]
        },ensure_ascii=False))

if __name__=="__main__":
    main()
