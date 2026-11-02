#!/usr/bin/env python3
"""Q4 botanical observation-effort sensitivity of port ↔ introduced plants.

Joins frozen 126-island WCVP all-flora + HOSTS-selected introduced taxa + modern
WPI ports to source-audited GBIF all-plant indexed occurrence counts (one
WGSRPD3 island region, same polygon, 1700–2026). This is deliberately not
causal transport, colonization, or historical shipping-frequency inference.

Only successful botanical GBIF effort queries are comparable: failures and
complex/antimeridian geometry are missing, NEVER zero effort. The 2023 WPI
port count is cross-sectional; the contemporary WCVP source has no arrival time.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

RNG_SEED = 20261102
PERMUTATIONS = 4999

def rows(path):
    with path.open(encoding="utf-8", newline="") as f:
        yield from csv.DictReader(f)

def float_value(s):
    try: return float(s)
    except (ValueError, TypeError): return None

def coeff(y, x, nuisance):
    # Frisch-Waugh-Lovell partial effect (avoid repeated full-model inversion).
    if len(y)!=len(x) or len(nuisance)!=len(x):
        raise ValueError("Mismatched input lengths")
    b_y=np.linalg.lstsq(nuisance,y,rcond=None)[0]
    b_x=np.linalg.lstsq(nuisance,x,rcond=None)[0]
    yr=y-nuisance@b_y
    xr=x-nuisance@b_x
    denom=float(xr@xr)
    if denom<1.e-10: return None
    return float(xr@yr/denom)

def model(data, outcome, extra, *,
          stratify="within_level1_area_half", permutations=PERMUTATIONS):
    n=len(data)
    if n<50: raise RuntimeError(f"Too few complete island effort observations: {n}")
    y=np.array([v[outcome] for v in data],dtype=float)
    x=np.array([v["log2ports"] for v in data],dtype=float)
    levels=sorted({r["level1"] for r in data})
    controls=[np.ones(n),np.array([v["log_area"] for v in data])]
    labels=["intercept","ln_WGSRPD3_region_area_km2"]
    for field in extra:
        controls.append(np.array([v[field] for v in data],dtype=float))
        labels.append(field)
    for level in levels[1:]:
        controls.append(np.array([1. if r["level1"]==level else 0. for r in data]))
        labels.append(f"Level1_{level}")
    Z=np.column_stack(controls)
    rank=int(np.linalg.matrix_rank(Z))
    if rank!=Z.shape[1]:
        raise RuntimeError(f"Singular effort-controlled model: rank {rank}/{Z.shape[1]} {labels}")
    observed=coeff(y,x,Z)
    if observed is None:raise RuntimeError("Port effect has no remaining variability")
    idx_by_level=defaultdict(list)
    for i,row in enumerate(data):
        idx_by_level[row["level1"]].append(i)
    # Stratification uses within-biogeographic-area median among COMPLETE effort data;
    # unrelated to response, ports, or wing size and held fixed under permutations.
    area_median_by_level={
        level:float(np.median([data[i]["log_area"] for i in indices]))
        for level,indices in idx_by_level.items()
    }
    strata=defaultdict(list)
    for i,row in enumerate(data):
        if stratify=="within_level1":
            stratum=row["level1"]
        else:
            stratum=row["level1"]+"|"+("small" if row["log_area"]<=area_median_by_level[row["level1"]] else "large")
        strata[stratum].append(i)
    groups=[ind for ind in strata.values() if len(ind)>1 and len({float(x[i]) for i in ind})>1]
    movable={i for indices in groups for i in indices}
    rng=random.Random(RNG_SEED)
    null=[]
    for _ in range(permutations):
        xp=x.copy()
        for ids in groups:
            order=ids.copy()
            rng.shuffle(order)
            xp[ids]=x[order]
        value=coeff(y,xp,Z)
        if value is not None: null.append(value)
    if not null:raise RuntimeError("No valid exposure permutations")
    mu=float(np.mean(null))
    tail=(1+sum(abs(z-mu)>=abs(observed-mu)-1e-12 for z in null))/(len(null)+1)
    positive_tail=(1+sum(z>=observed-1e-12 for z in null))/(len(null)+1)
    return {
        "n":n,
        "outcome":outcome,
        "exposure":"log2(1+2023_WPI_indexed_port_count)",
        "nuisance_covariates":labels,
        "regression_beta_per_log2_port_unit":observed,
        "multiplicative_change_in_1plus_outcome" :
            math.exp(observed)-1 if outcome.startswith("ln1p_") else None,
        "permutation":{
            "stratify":stratify,
            "n_permutations":len(null),
            "exchangeable_islands":len(movable),
            "exchangeable_blocks":len(groups),
            "null_mean":mu,
            "null_q025":float(np.quantile(null,.025)),
            "null_q975":float(np.quantile(null,.975)),
            "observed_minus_null_mean":observed-mu,
            "p_two_sided_centered":tail,
            "p_one_sided_positive":positive_tail,
        }
    }

def main():
    p=argparse.ArgumentParser()
    for name in ("all-flora","host-taxa","gbif-effort","canary-effort-audit","out-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    flora=list(rows(a.all_flora))
    hosts={r["code"]:r for r in rows(a.host_taxa)}
    effort=list(rows(a.gbif_effort))
    canary=json.loads(a.canary_effort_audit.read_text(encoding="utf-8"))
    if len(flora)!=126 or len(hosts)!=123 or len(effort)!=122:
        raise RuntimeError("Source population drift: need 126 WCVP islands, 123 HOSTS, 122 Q4 effort islands")
    if canary["source"]["q4_v02_population"]["introduced_host_plant_x_island_pairs"]!=5611:
        raise RuntimeError("Corrected Canary Q2/Q4 source drift")
    e={}
    statuses=Counter()
    for row in effort:
        k=row["island_code"];status=row["gbif_all_plants_index_status"]
        statuses[status]+=1
        if status!="DATED_OCCURRENCE_FOUND":continue
        n=float_value(row["gbif_all_plants_total_records"])
        pre=float_value(row["gbif_all_plants_pre1900_count"])
        pre2001=float_value(row["gbif_all_plants_pre2001_count"])
        if n is None or n<=0 or pre is None or pre2001 is None:
            raise RuntimeError(f"Invalid GBIF total records for success {k}")
        e[k]={"total":n,"pre1900":pre,"pre2001":pre2001,
              "first_year":int(row["gbif_all_plants_first_record_year"])}
    c=canary["all_botanical_indexed_effort"]
    if c["status"]!="DATED_OCCURRENCE_FOUND" or int(c["total_indexed_records"])<=0:
        raise RuntimeError("Canary all-botanical effort not verified")
    if "CNY" in e:
        raise RuntimeError("CNY has been in the original uncorrected pilot unexpectedly")
    e["CNY"]={"total":float(c["total_indexed_records"]),
              "pre1900":float(c["pre1900"]),
              "pre2001":None,
              "first_year":int(c["first_year"])}
    fcodes={r["code"] for r in flora}
    if len(fcodes)!=126 or len(set(e)-fcodes):
        raise RuntimeError("GBIF WCVP WGSRPD regions misaligned")
    enriched=[]
    excluded=[]
    for row in flora:
        code=row["code"]
        if code not in e:
            excluded.append({"code":code,"reason":"GBIF_EFFORT_NOT_INDEXED_OR_FAILED"})
            continue
        n_native=int(row["native"])
        n_intro=int(row["introduced"])
        if n_native<0 or n_intro<0:raise RuntimeError("Negative native/introduced WCVP counts")
        if code not in hosts:
            # Host taxa network has zero opportunity on 3 islands; all flora
            # analysis keeps those islands, subset analyses drop them.
            host_taxa=None; native_host_baseline=None
        else:
            host_taxa=int(hosts[code]["introduced_larval_host_taxa"])
            native_host_baseline=int(hosts[code]["native"])
        enriched.append({
            "code":code,"level1":row["level1"],
            "n_native_flora":n_native,"n_introduced_flora":n_intro,
            "host_taxa":host_taxa,"native_host_baseline":native_host_baseline,
            "log_area":math.log(float(row["area_km2"])),
            "WPI_ports":int(row["wpi_port_count"]),
            "log2ports":math.log2(1+int(row["wpi_port_count"])),
            "gbif_all_plant_indexed_records":e[code]["total"],
            "gbif_pre1900_plant_records":e[code]["pre1900"],
            "first_botanical_indexed_year":e[code]["first_year"],
            "ln1p_gbif_total":math.log1p(e[code]["total"]),
            "ln1p_gbif_pre1900":math.log1p(e[code]["pre1900"]),
            "ln1p_native_flora":math.log1p(n_native),
            "ln1p_introduced_flora":math.log1p(n_intro),
            "log_ratio_introduced_to_native":math.log1p(n_intro)-math.log1p(n_native),
            "ln1p_introduced_host_taxa":math.log1p(host_taxa) if host_taxa is not None else None,
            "ln1p_native_host_baseline":math.log1p(native_host_baseline) if native_host_baseline is not None else None,
        })
    if len(enriched)<80:raise RuntimeError(f"GBIF survey-effort coverage unexpectedly small n={len(enriched)}")
    present_host=[x for x in enriched if x["host_taxa"] is not None]
    m={}
    comparisons=[
      ("all_flora_unnormalized_no_effort",enriched,"ln1p_introduced_flora",["ln1p_native_flora"]),
      ("all_flora_with_gbif_total_effort",enriched,"ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_total"]),
      ("all_flora_with_pre1900_effort",enriched,"ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_pre1900"]),
      ("all_flora_intro_native_logratio_with_effort",enriched,"log_ratio_introduced_to_native",["ln1p_gbif_total"]),
      ("host_taxa_no_effort",present_host,"ln1p_introduced_host_taxa",["ln1p_native_host_baseline"]),
      ("host_taxa_with_effort",present_host,"ln1p_introduced_host_taxa",["ln1p_native_host_baseline","ln1p_gbif_total"]),
    ]
    for name,values,response,controls in comparisons:
        m[name]=model(values,response,controls)
    # Additional null: exact Level1 placement, not restricted to area half.
    m["all_flora_with_effort_within_Level1_only"]=model(
        enriched,"ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_total"],
        stratify="within_level1")
    # Scientific sample exclusion sensitivity: omit region with most modern indexed ports
    max_port=max(row["WPI_ports"] for row in enriched)
    m["all_flora_effort_exclude_largest_port"]=model(
        [row for row in enriched if row["WPI_ports"]<max_port],
        "ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_total"])
    m["all_flora_effort_exclude_pacific_antarctic"]=model(
        [row for row in enriched if row["level1"] not in {"6","9"}],
        "ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_total"])
    # Check whether ports predict the reporting process itself; association
    # here is not proof that all port–flora results are measurement artifacts.
    m["botanical_reporting_effort_vs_ports"]=model(
        enriched,"ln1p_gbif_total",["ln1p_native_flora"])
    # Port-intensive islands with at least one WPI entry provide a separate
    # exposure gradient, distinct from the binary catalogued-port difference.
    m["all_flora_with_effort_only_port_positive"]=model(
        [row for row in enriched if row["WPI_ports"]>0],
        "ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_total"])
    m["all_flora_with_effort_high_botanical_reporting"]=model(
        [row for row in enriched if row["gbif_all_plant_indexed_records"]>=1000],
        "ln1p_introduced_flora",["ln1p_native_flora","ln1p_gbif_total"])
    all_flora_map={r["code"]:r for r in flora}
    excluded_full=[all_flora_map[x["code"]] for x in excluded]
    def median_float(arr,field):
        return float(np.median([float(x[field]) for x in arr])) if arr else None
    selection_audit={
       "n_observed":len(enriched),
       "n_excluded":len(excluded_full),
       "median_area_km2_observed":float(np.median([math.exp(x["log_area"]) for x in enriched])),
       "median_area_km2_excluded":median_float(excluded_full,"area_km2"),
       "median_WPI_ports_observed":float(np.median([x["WPI_ports"] for x in enriched])),
       "median_WPI_ports_excluded":median_float(excluded_full,"wpi_port_count"),
       "median_introduced_flora_observed":float(np.median([x["n_introduced_flora"] for x in enriched])),
       "median_introduced_flora_excluded":median_float(excluded_full,"introduced"),
       "excluded_high_port_regions":sorted(
           [{"code":x["code"],"ports":int(x["wpi_port_count"]),"introduced_flora":int(x["introduced"])}
            for x in excluded_full],key=lambda z:-z["ports"])[:12],
       "interpretation":"Effort availability is size- and region-selective; models do not generalize automatically to all 126 islands"
    }
    summary={
      "schema":"chocho_WPI_introduced_flora_GBIF_plant_index_effort_sensitivity_v0.1",
      "status":"DESCRIPTIVE_CROSS_SECTIONAL_SURVEY_EFFORT_AUDIT_NOT_TRANSPORT_CAUSATION",
      "sources":{
        "wcvp_all_flora_126":"results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv",
        "butterfly_host_taxa_123":"results/three-gate-port-introduced-host-taxa/distinct_introduced_host_taxa_and_port_by_island.csv",
        "gbif_botanical_effort_122":"results/q4-gbif-first-year-pilot/island_plant_effort_year_facets.csv",
        "canary_botanical_effort":"results/q4-canary-plant-firstyear-v02/CNY_three_plant_firstyear_audit.json",
        "historical_index_definition":"GBIF all Plantae taxonKey=6; 1700–2026, WGSRPD3 polygon, valid coordinates, all GBIF providers",
      },
      "scope":{
         "total_WCVP_island_regions":126,
         "original_Q4_botanical_effort_status":dict(statuses),
         "CNY_effort_added":True,
         "island_regions_with_successful_botanical_index_effort":len(enriched),
         "islands_with_HOSTS_and_botanical_effort":len(present_host),
         "excluded_regions":excluded,
         "cohort_counts_by_level1":dict(Counter(x["level1"] for x in enriched)),
         "cohort_islands_with_WPI_indexed_ports":sum(x["WPI_ports"]>0 for x in enriched),
         "cohort_islands_without_WPI_indexed_ports":sum(x["WPI_ports"]==0 for x in enriched),
      },
      "source_effort_selection":selection_audit,
      "models":m,
      "interpretation_boundaries":[
        "WPI port count comes from a modern indexed 2023 snapshot; it is not historical shipping intensity, port accessibility in all years, or a treatment.",
        "WCVP floristic status and GBIF botanical effort may share collecting/digitization biases and are not independent data collection systems in a causal sense.",
        "Total indexed plant records is a post-exposure variable: adjustment is a bias sensitivity only, and may block real human-linked botanical introductions as well as collection effort.",
        "Strict missing GBIF effort sources are omitted, not treated as botanical zero; the complete-case cohort is geographically nonrandom.",
        "Matched geography still aggregates physical islands, not ports or breeding sites.",
        "Area and Level1 strata do not eliminate population-density, land-use, gardening/nursery plants, trade or spatial dependence.",
        "Permutations are descriptive against exchangeable port-label baselines, not evidence of independent transport causality.",
      ]
    }
    out_json=a.out_dir/"botanical_effort_adjusted_port_flora_audit.json"
    out_json.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    with (a.out_dir/"botanical_effort_flora_island_inputs.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["code","level1","n_native_flora","n_introduced_flora","host_taxa","native_host_baseline","WPI_ports",
                "log_area","gbif_all_plant_indexed_records","gbif_pre1900_plant_records","first_botanical_indexed_year"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in enriched:w.writerow({k:r.get(k) for k in fields})
    print(json.dumps({
        "status":summary["status"],
        "scope":summary["scope"],
        "model_summaries":{k:{
             "n":v["n"],
             "beta":v["regression_beta_per_log2_port_unit"],
             "p_two":v["permutation"]["p_two_sided_centered"],
             "percent_change":v["multiplicative_change_in_1plus_outcome"],
        } for k,v in m.items()}
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
