#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


def exact_one_sided_binom(successes: int, trials: int) -> float:
    if trials <= 0:
        return 1.0
    return sum(math.comb(trials, k) for k in range(successes, trials + 1)) / (2 ** trials)


def quantile(values, p):
    xs=sorted(float(x) for x in values)
    if not xs: return None
    pos=p*(len(xs)-1); lo=int(math.floor(pos)); hi=int(math.ceil(pos))
    if lo==hi: return xs[lo]
    f=pos-lo
    return xs[lo]*(1-f)+xs[hi]*f


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--queried-csv",type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=99999)
    ap.add_argument("--seed",type=int,default=20261007)
    ap.add_argument("--output-json",type=Path,required=True)
    ap.add_argument("--output-cell-csv",type=Path,required=True)
    a=ap.parse_args()

    rows=list(csv.DictReader(a.queried_csv.open(newline="",encoding="utf-8")))
    by_cell=defaultdict(list)
    for r in rows:
        by_cell[(str(r["species"]).strip(),str(r["wgsrpd3_code"]).strip())].append(r)

    cells=[]
    for (sp,region), vals in sorted(by_cell.items()):
        actual=[r for r in vals if r["assignment_type"]=="actual"]
        pseudo1=[r for r in vals if r["assignment_type"]=="pseudo" and int(float(r["match_rank"]))==1]
        # Primary design is defined only for cells with exactly one actual introduced known host.
        if len(actual)!=1 or len(pseudo1)!=1:
            continue
        ar=actual[0]; pr=pseudo1[0]
        ay=str(ar.get("local_host_first_record_year") or "").strip()
        py=str(pr.get("local_host_first_record_year") or "").strip()
        by=str(ar.get("butterfly_first_record_year_full") or "").strip()
        if not ay or not py or not by:
            complete=0
            actual_year=pseudo_year=butterfly_year=None
            diff=None
        else:
            complete=1
            actual_year=int(float(ay)); pseudo_year=int(float(py)); butterfly_year=int(float(by))
            diff=pseudo_year-actual_year
        cells.append({
            "species":sp,
            "wgsrpd3_code":region,
            "actual_host":ar["candidate_name"],
            "paired_pseudoresource":pr["candidate_name"],
            "match_distance":float(pr["match_distance"]),
            "actual_host_first_record_year":actual_year,
            "pseudo_first_record_year":pseudo_year,
            "butterfly_first_record_year":butterfly_year,
            "complete_primary_chronology":complete,
            "pseudo_minus_actual_resource_year":diff,
            "actual_resource_precedes_butterfly":None if not complete else int(actual_year<=butterfly_year),
            "pseudo_resource_precedes_butterfly":None if not complete else int(pseudo_year<=butterfly_year),
            "actual_strictly_earlier_than_pseudo":None if not complete else int(actual_year<pseudo_year),
            "pseudo_strictly_earlier_than_actual":None if not complete else int(pseudo_year<actual_year),
        })

    complete=[r for r in cells if int(r["complete_primary_chronology"])==1]
    if not complete:
        raise RuntimeError("no complete single-host paired chronology cells")

    diffs=[int(r["pseudo_minus_actual_resource_year"]) for r in complete]
    actual_earlier=sum(d>0 for d in diffs)
    pseudo_earlier=sum(d<0 for d in diffs)
    ties=sum(d==0 for d in diffs)
    non_ties=actual_earlier+pseudo_earlier
    sign_p=exact_one_sided_binom(actual_earlier,non_ties)

    # Precedence relative to the same butterfly clock: McNemar-style exact sign test
    # among discordant pairs.
    actual_only=0; pseudo_only=0; both=0; neither=0
    for r in complete:
        aa=int(r["actual_resource_precedes_butterfly"])
        pp=int(r["pseudo_resource_precedes_butterfly"])
        if aa and pp: both+=1
        elif aa and not pp: actual_only+=1
        elif pp and not aa: pseudo_only+=1
        else: neither+=1
    precedence_p=exact_one_sided_binom(actual_only,actual_only+pseudo_only)

    observed_mean=sum(diffs)/len(diffs)
    observed_median=float(statistics.median(diffs))
    rng=random.Random(a.seed)
    null_means=[]; null_medians=[]
    for _ in range(a.permutations):
        x=[d if rng.random()<0.5 else -d for d in diffs]
        null_means.append(sum(x)/len(x))
        null_medians.append(float(statistics.median(x)))
    p_mean=(1+sum(x>=observed_mean for x in null_means))/(a.permutations+1)
    p_median=(1+sum(x>=observed_median for x in null_medians))/(a.permutations+1)

    # Cluster-level robustness: cells are the frozen primary unit, but repeated
    # species, actual hosts and regions can induce dependence. Aggregate paired
    # differences within each cluster and test only the direction of cluster medians.
    def cluster_direction(field):
        grouped=defaultdict(list)
        for r in complete:
            grouped[str(r[field])].append(int(r["pseudo_minus_actual_resource_year"]))
        med={k:float(statistics.median(v)) for k,v in sorted(grouped.items())}
        pos=sum(v>0 for v in med.values()); neg=sum(v<0 for v in med.values()); zero=sum(v==0 for v in med.values())
        return {
            "clusters":len(med),
            "positive_median":pos,
            "negative_median":neg,
            "zero_median":zero,
            "exact_one_sided_sign_p_positive":exact_one_sided_binom(pos,pos+neg),
            "median_difference_by_cluster":med,
        }

    species_robust=cluster_direction("species")
    host_robust=cluster_direction("actual_host")
    region_robust=cluster_direction("wgsrpd3_code")

    a.output_cell_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_cell_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(cells[0])); w.writeheader(); w.writerows(cells)

    payload={
        "schema":"chocho_resource_consumer_single_host_timing_v0.2",
        "status":"IDENTITY_SPECIFIC_SINGLE_HOST_PAIRED_TEMPORAL_TEST",
        "seed":a.seed,
        "permutations":a.permutations,
        "cells_in_input_with_single_actual_and_rank1_pseudo":len(cells),
        "cells_with_complete_actual_pseudo_and_butterfly_chronology":len(complete),
        "paired_local_resource_timing":{
            "median_pseudo_minus_actual_years":observed_median,
            "mean_pseudo_minus_actual_years":observed_mean,
            "actual_host_strictly_earlier_cells":actual_earlier,
            "pseudo_strictly_earlier_cells":pseudo_earlier,
            "ties":ties,
            "exact_one_sided_sign_p_actual_earlier":sign_p,
            "sign_flip_mean_p":p_mean,
            "sign_flip_median_p":p_median,
            "null_mean_ci95":[quantile(null_means,0.025),quantile(null_means,0.975)],
            "null_median_ci95":[quantile(null_medians,0.025),quantile(null_medians,0.975)]
        },
        "resource_precedence_relative_to_butterfly":{
            "both_actual_and_pseudo_precede_butterfly":both,
            "actual_only_precedes_butterfly":actual_only,
            "pseudo_only_precedes_butterfly":pseudo_only,
            "neither_precedes_butterfly":neither,
            "exact_one_sided_discordant_p_actual_precedence":precedence_p
        },
        "cluster_robustness":{
            "species":species_robust,
            "actual_host":host_robust,
            "wgsrpd3_region":region_robust
        },
        "decision":{
            "identity_specific_resource_timing_supported":bool(sign_p<=0.05 and observed_median>0),
            "identity_specific_resource_precedence_supported":bool(precedence_p<=0.05 and actual_only>pseudo_only),
            "broad_generality_supported":bool(
                sign_p<=0.05 and observed_median>0
                and species_robust["exact_one_sided_sign_p_positive"]<=0.05
                and host_robust["exact_one_sided_sign_p_positive"]<=0.05
            ),
            "interpretation_if_cluster_robustness_fails":"Retain a cell-level identity-specific timing result only as a concentrated pattern; do not present it as a general butterfly resource-tracking principle."
        },
        "claim_boundary":[
            "First records are detection/digitization dates rather than establishment dates.",
            "The paired pseudo-resource is an alien non-host in the same WGSRPD3 region matched before local chronology is inspected.",
            "Primary inference is restricted to exactly one actual known host versus one paired pseudo-resource per cell.",
            "Known-host presence does not prove local larval use or demographic benefit."
        ]
    }
    a.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))

if __name__=="__main__":
    main()
