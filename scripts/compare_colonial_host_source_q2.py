#!/usr/bin/env python3
"""Q2-B historical colony × native-metropole host-source alignment.

Joins independently audited WCVP host-source islands and Lenzner 2022/GloNAF
historical empire polygons. No inference of vessel transport, first arrival,
butterfly occupancy, or causal colonization effect is permitted.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

EMPIRES=("Great_Britain","Spain","Portugal","Netherlands")
SEED=20261102


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as h:
        yield from csv.DictReader(h)


def write_csv(path, rows):
    if not rows: return
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def unique_perms(values):
    # Only exact-enumerate blocks of <=8 observations, else caller uses Monte-Carlo.
    return sorted(set(itertools.permutations(values)))


def shuffle_null(rows, groups, permutations=19999):
    """Score matched metropole native-larval-host support, equal island weight."""
    if len(rows)<3:return None
    values=[r["support_fractions"] for r in rows]
    lab=[r["empire"] for r in rows]
    bygroup=defaultdict(list)
    for i,x in enumerate(groups):bygroup[x].append(i)
    mutable=[inds for inds in bygroup.values() if len({lab[i] for i in inds})>1]
    effective=len({i for inds in mutable for i in inds})
    empirical=sum(values[i][e] for i,e in enumerate(lab))/len(lab)
    if not mutable:
        return {"observed_score":empirical,"null_status":"NO_EXCHANGEABLE_EMPIRE_LABELS",
                "effective_islands":0,"independent_exchangeable_blocks":0}
    permutationsets=[]
    complexity=1
    exact=True
    for inds in mutable:
        p=unique_perms([lab[i] for i in inds]) if len(inds)<=8 else None
        if p is None:
            exact=False
            break
        complexity*=len(p)
        if complexity>30000:
            exact=False
            break
        permutationsets.append((inds,p))
    draws=[]
    if exact:
        for assignment in itertools.product(*[p for _,p in permutationsets]):
            shuffled=lab.copy()
            for (inds,_),vals in zip(permutationsets,assignment):
                for i,v in zip(inds,vals):shuffled[i]=v
            draws.append(sum(values[i][e] for i,e in enumerate(shuffled))/len(lab))
    else:
        rng=random.Random(SEED)
        for _ in range(permutations):
            shuffled=lab.copy()
            for inds in mutable:
                labels=[shuffled[i] for i in inds]
                rng.shuffle(labels)
                for i,e in zip(inds,labels):shuffled[i]=e
            draws.append(sum(values[i][e] for i,e in enumerate(shuffled))/len(lab))
    mean_null=statistics.mean(draws)
    delta=empirical-mean_null
    greater=sum(v>=empirical-1e-12 for v in draws)
    two=sum(abs(v-mean_null)>=abs(delta)-1e-12 for v in draws)
    if exact:
        p_high=greater/len(draws)
        p_two=two/len(draws)
    else:
        p_high=(greater+1)/(len(draws)+1)
        p_two=(two+1)/(len(draws)+1)
    return {
        "observed_score":empirical,
        "null_mean":mean_null,
        "observed_minus_null_percentage_points":100*delta,
        "null_q025":sorted(draws)[int(.025*(len(draws)-1))],
        "null_q975":sorted(draws)[int(.975*(len(draws)-1))],
        "null_p_high":p_high,
        "null_p_two_sided":p_two,
        "null_type":"EXACT" if exact else "MONTE_CARLO",
        "null_assignments":len(draws),
        "effective_islands":effective,
        "independent_exchangeable_blocks":len(mutable),
        "assumptions":"Empire label permutations preserve exactly the observed source-metropole prevalence matrix and the label frequency within each stratum; p values are descriptive and do not correct nonrandom historic locations or host data bias.",
    }


def rankr(values):
    order=sorted(range(len(values)),key=lambda i:values[i])
    result=[0.]*len(values)
    j=0
    while j<len(order):
        k=j+1
        while k<len(order) and values[order[k]]==values[order[j]]:
            k+=1
        r=0.5*(j+1+k)
        for pos in order[j:k]:result[pos]=r
        j=k
    return result


def spearman(a,b):
    if len(a)<4:return None
    x,y=rankr(a),rankr(b)
    mx,my=statistics.mean(x),statistics.mean(y)
    sx=sum((v-mx)**2 for v in x)
    sy=sum((v-my)**2 for v in y)
    if sx==0 or sy==0:return None
    return sum((v-mx)*(w-my) for v,w in zip(x,y))/math.sqrt(sx*sy)


def main():
    p=argparse.ArgumentParser()
    for k in ("host-source-islands","host-source-audit","empire-islands",
              "empire-overlay-audit","area-metrics","output-dir"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    host_audit=json.loads(a.host_source_audit.read_text(encoding="utf-8"))
    empires_audit=json.loads(a.empire_overlay_audit.read_text(encoding="utf-8"))
    island_sources={r["island_code"]:r for r in read_csv(a.host_source_islands)}
    geographic={r["area_code_l3"]:r for r in read_csv(a.empire_islands)}
    area_metrics={r["code"]:r for r in read_csv(a.area_metrics)}
    if len(geographic)!=125 or empires_audit["frozen_island_regions"]!=125:
        raise RuntimeError("Historical source panel mismatch")
    if host_audit["scope"]["added_strict_island_cells"]!=5674:
        raise RuntimeError("Donor native source panel mismatch")
    rows=[]
    all_strict=Counter()
    for code,hist in sorted(geographic.items()):
        if hist["strict_historical_empire"]!="1":continue
        label=hist["classified_historical_empire"]
        if label not in EMPIRES:continue
        all_strict[label]+=1
        island=island_sources.get(code)
        if not island:continue
        n=int(island["new_butterfly_resource_cells"])
        if n==0:continue
        donor_fraction={e:int(island[e+"_native_supported_cells"])/n for e in EMPIRES}
        years=float(hist[label+"_weighted_exposure_years"])
        area_km2=float(area_metrics[code]["area_km2"])
        rows.append({
            "area_code_l3":code,"region_name":hist["region_name"],"level1":hist["level1"],
            "empire":label,"area_km2":area_km2,"years":years,
            "new_resource_cells":n,
            "support_fractions":donor_fraction,
            "actual_support_share":donor_fraction[label],
            "alternate_mean_support_share":statistics.mean(donor_fraction[e] for e in EMPIRES if e!=label),
            "within_island_actual_minus_other_pp":100*(donor_fraction[label]-statistics.mean(donor_fraction[e] for e in EMPIRES if e!=label)),
        })
    if len(rows)<8:raise RuntimeError("Too few strict historical island assignments with resource opportunity")

    def evaluate(sample):
        if len(sample)<3:return {"included":len(sample),"status":"LOW_SUPPORT"}
        region_counts=Counter(r["level1"] for r in sample)
        groups={"global": ["all"]*len(sample),
                "level1": [r["level1"] for r in sample]}
        area_cut=statistics.median(r["area_km2"] for r in sample)
        groups["level1_area_half"]=[r["level1"]+"|"+("L" if r["area_km2"]<=area_cut else "H") for r in sample]
        raw={
            "included":len(sample),
            "number_by_empire":dict(Counter(r["empire"] for r in sample)),
            "total_resource_cells":sum(r["new_resource_cells"] for r in sample),
            "actual_empire_supported_butterfly_region_cells":sum(
                r["new_resource_cells"]*r["actual_support_share"] for r in sample),
            "mean_own_metropole_support_share":statistics.mean(r["actual_support_share"] for r in sample),
            "mean_other_metropoles_support_share":statistics.mean(r["alternate_mean_support_share"] for r in sample),
            "mean_raw_actual_minus_other_pp":statistics.mean(r["within_island_actual_minus_other_pp"] for r in sample),
            "own_exposure_years_vs_actual_host_source_spearman":spearman(
                [r["years"] for r in sample],[r["actual_support_share"] for r in sample]),
            "permutation":{k:shuffle_null(sample,v) for k,v in groups.items()},
            "min_added_cells":min(r["new_resource_cells"] for r in sample),
            "area_median_km2":area_cut,
            "region_sample":dict(region_counts),
        }
        return raw
    cohorts={
      "strict_geographic_all":rows,
      "strict_and_ge_50_years":[r for r in rows if r["years"]>=50],
      "strict_and_ge_100_years":[r for r in rows if r["years"]>=100],
      "strict_and_ge_50_years_min25_cells":[r for r in rows if r["years"]>=50 and r["new_resource_cells"]>=25],
      "strict_and_ge_50_years_non_european":[r for r in rows if r["years"]>=50 and r["level1"]!="1"],
    }
    outcome={
       "schema":"chocho_colonial_host_source_match_q2_v0.1",
       "status":"OBSERVATIONAL_EMPIRE_SOURCE_ALIGNMENT_SCREEN",
       "sources":{"botanical":host_audit["source_commits"],
                  "empire_history":empires_audit["source"]},
       "source_overlay_strict_islands":sum(all_strict.values()),
       "source_overlay_strict_empire_counts":dict(all_strict),
       "source_native_island_total":host_audit["scope"]["islands_with_added_resource"],
       "n_with_strict_hist_and_resource":len(rows),
       "cohorts":{name:evaluate(sample) for name,sample in cohorts.items()},
       "interpretation_guardrails":[
         "The source-native indicator measures plants native in an empire's European metropole; it does not identify butterfly native ranges or import routes.",
         "The regression of plantation history and plant source match is descriptive. No conditional independent effect of colonization is identified.",
         "The published Lenzner mapping does not cover all colonial empires, including French, and may flag brief occupations; min50years sensitivity matters.",
         "Small island and empire cohort, geographic overlap restrictions, and host-taxonomy resolution must be retained in the manuscript.",
         "Within Level1 stratum permutation retains historical empire frequency but cannot remove unobserved climate/human history confounding.",
         "Record co-occurrence does not prove larvae actually utilized a plant, species established, or port introduction."
       ]
    }
    (a.output_dir/"colonial_host_source_alignment_audit.json").write_text(
      json.dumps(outcome,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    printable=[{k:v for k,v in r.items() if k!="support_fractions"}|
      {e+"_native_host_opportunity_share":r["support_fractions"][e] for e in EMPIRES}
      for r in rows]
    write_csv(a.output_dir/"colonial_host_islands.csv",printable)
    print(json.dumps(outcome,ensure_ascii=False,indent=2))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
