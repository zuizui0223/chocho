#!/usr/bin/env python3
"""Source-locked cohort bridge: 239-resource taxa vs 155 recording-forecast taxa.

This AUDITS SELECTION, not prediction gain by species. Frozen public forecast
outputs contain global/split MSE and a 225-row selection ledger, but not the
original per-butterfly M3 and M5 out-of-fold losses or prediction Parquet.
Do not derive a species-level forecasting gain from source-aggregate MSE.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

DATA = {
    "traits": "data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "eligible": "provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_all_225_final_coverage.csv",
    "gap": "results/island_resource_gap_filling/by_species.csv",
    "forecast": "provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_validation_all_metrics.csv",
}
SEED = 20261115
REPS = 4999
MODELED = "MODELED_ALL_PLANNED_VALIDATIONS"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def source_tables(root):
    result, digests = {}, {}
    for key, relative in DATA.items():
        path = root / relative
        result[key] = read_csv(path)
        digests[key] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result, digests


def describe(rows):
    if not rows:
        return {"n_species": 0}
    n = len(rows)
    def mean(key):
        return sum(z[key] for z in rows)/n
    ns_i = sum(z["island_native"] for z in rows)
    as_i = sum(z["island_added"] for z in rows)
    ns_m = sum(z["mainland_native"] for z in rows)
    as_m = sum(z["mainland_added"] for z in rows)
    return {
        "n_species": n,
        "host_catalogue_mean_known_species": mean("known_hosts"),
        "global_mean_added_resource_regions": mean("added_global"),
        "global_positive_added_butterflies": sum(z["added_global"] > 0 for z in rows),
        "global_native_resource_region_mean": mean("native_global"),
        "species_equal_island_gap_fill_pct": mean("fill_i"),
        "species_equal_mainland_gap_fill_pct": mean("fill_m"),
        "species_equal_island_minus_mainland_fill_pp": mean("delta"),
        "island_more_fill_species": sum(z["delta"] > 0 for z in rows),
        "mainland_more_fill_species": sum(z["delta"] < 0 for z in rows),
        "fill_ties": sum(z["delta"] == 0 for z in rows),
        "island_native_resource_cells": ns_i,
        "island_introduced_added_resource_cells": as_i,
        "island_resource_gain_pct_native_denominator": 100*as_i/ns_i if ns_i else None,
        "mainland_native_resource_cells": ns_m,
        "mainland_introduced_added_resource_cells": as_m,
        "mainland_resource_gain_pct_native_denominator": 100*as_m/ns_m if ns_m else None,
        "diet_family_classes": dict(Counter(z["family"] for z in rows)),
        "taxonomy_adequate_species": sum(z["taxonomy_adequate"] for z in rows),
    }


def diff(rows, key):
    a = [z[key] for z in rows if z["modelled"]]
    b = [z[key] for z in rows if not z["modelled"]]
    if not a or not b:
        raise ValueError("Cannot compare a missing source cohort")
    return sum(a)/len(a)-sum(b)/len(b)


def bootstrap(rows, key, reps=REPS, seed=SEED):
    rng = random.Random(seed)
    a = [z[key] for z in rows if z["modelled"]]
    b = [z[key] for z in rows if not z["modelled"]]
    draw=[]
    for _ in range(reps):
        ma=sum(a[rng.randrange(len(a))] for _ in a)/len(a)
        mb=sum(b[rng.randrange(len(b))] for _ in b)/len(b)
        draw.append(ma-mb)
    draw.sort()
    return [draw[int((reps-1)*.025)],draw[int((reps-1)*.975)]]


def preserve_strata(rows, key, by, reps=REPS, seed=SEED):
    groups=defaultdict(list)
    for idx,row in enumerate(rows):
        group=tuple(row[k] for k in by)
        groups[group].append(idx)
    supported={}
    excluded=[]
    for name,indices in groups.items():
        n=sum(rows[i]["modelled"] for i in indices)
        if 0<n<len(indices):
            supported[name] = (indices,n)
        else:
            excluded.extend(indices)
    keep=sorted(i for indices,_ in supported.values() for i in indices)
    if len(keep)<50:
        raise ValueError("Too few common-support species to permute")
    sub=[rows[i] for i in keep]
    observed=diff(sub,key)
    # Convert original row indices to indices of the supported subset.
    idmap={i:j for j,i in enumerate(keep)}
    strata=[([idmap[i] for i in ids],n) for ids,n in supported.values()]
    n_m=sum(row["modelled"] for row in sub)
    rng=random.Random(seed)
    null=[]
    for _ in range(reps):
        mask=[False]*len(sub)
        for ids,n in strata:
            for ix in rng.sample(ids,n):
                mask[ix]=True
        a=sum(z[key] for i,z in enumerate(sub) if mask[i])/n_m
        b=sum(z[key] for i,z in enumerate(sub) if not mask[i])/(len(sub)-n_m)
        null.append(a-b)
    mu=sum(null)/len(null)
    sorted_null=sorted(null)
    return {
        "strata": [str(k) for k in sorted(supported)],
        "n_support_strata":len(supported),
        "n_excluded_one_group_only":len(excluded),
        "n_matched_support_species":len(sub),
        "n_modeled_in_support":n_m,
        "n_not_modeled_in_support":len(sub)-n_m,
        "observed_difference":observed,
        "null_mean":mu,
        "null_95pct_interval":[sorted_null[int((reps-1)*.025)],sorted_null[int((reps-1)*.975)]],
        "permutation_two_sided_centered_p":(1+sum(abs(z-mu)>=abs(observed-mu)-1e-12 for z in null))/(reps+1),
        "permutation_replicates":reps,
        "null_def":"Relabel forecast eligible within each original host-family-breadth and/or global-native-resource quartile stratum, keeping selection counts fixed",
        "not_causal_or_phylogeny_adjusted":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument("--out",type=Path,default=Path("results/resource_forecast_cohort_selection_20261115"))
    a=p.parse_args()
    src,hashes=source_tables(a.root)
    traits={x["species"]:x for x in src["traits"]}
    eligibility={x["species"]:x for x in src["eligible"]}
    paired=defaultdict(dict)
    for r in src["gap"]:
        if r["group"] in paired[r["species"]]:
            raise ValueError("Duplicated species × geography resource rows")
        paired[r["species"]][r["group"]]=r
    if len(traits)!=239 or len(eligibility)!=225 or len(paired)!=239:
        raise ValueError("Original 239 / 225 source cohort sizes changed")
    if not set(eligibility).issubset(traits) or set(paired)!=set(traits):
        raise ValueError("Butterfly scientific-name identity does not join")
    sample=[]
    for sp in sorted(traits):
        t=traits[sp];g=paired[sp]
        if set(g)!={"island","mainland"}:raise ValueError("Missing original geography for "+sp)
        i,m=g["island"],g["mainland"]
        if (+i["empty"]<=0 or +m["empty"]<=0):
            raise ValueError("No native-host-gaps denominator for "+sp)
        inframe=sp in eligibility
        eligible=inframe and eligibility[sp]["final_disposition"]==MODELED
        sample.append({
            "species":sp,
            "forecast_category":"modeled155" if eligible else ("excluded70" if inframe else "outsideOriginal225"),
            "modelled":eligible,
            "family":t["host_breadth_stratum"],
            "taxonomy_adequate":t["host_taxonomy_lower_bound_adequate"]=="True",
            "known_hosts":int(t["resolved_host_species"]),
            "native_global":int(t["native_resource_units"]),
            "added_global":int(t["introduced_added_units"]),
            "island_native":int(i["native"]),
            "island_added":int(i["added"]),
            "mainland_native":int(m["native"]),
            "mainland_added":int(m["added"]),
            "fill_i":100*float(i["fill"]),
            "fill_m":100*float(m["fill"]),
            "delta":100*(float(i["fill"])-float(m["fill"])),
            "forecast_test_records":int(eligibility[sp]["temporal_test_pairs"]) if inframe else 0,
            "original_source_records":int(float(eligibility[sp]["source_records"])) if inframe else 0,
        })
    cohorts={c:[z for z in sample if z["forecast_category"]==c] for c in
             ("modeled155","excluded70","outsideOriginal225")}
    if [len(cohorts[x]) for x in cohorts] != [155,70,14]:
        raise ValueError("Original model-status cohort sizes changed")
    modeled_test=sum(z["forecast_test_records"] for z in cohorts["modeled155"])
    if modeled_test!=208482:
        raise ValueError("155 forecast test row sum changed")

    overlap=cohorts["modeled155"]+cohorts["excluded70"]
    # Fixed quartiles of the original 225 source taxa's *global native* host
    # opportunities, pre-existing before the phenotype response test.
    ns=sorted(z["native_global"] for z in overlap)
    quantiles=[ns[int((len(ns)-1)*q)] for q in (.25,.5,.75)]
    for z in overlap:
        z["native_quartile"]=sum(z["native_global"]>cut for cut in quantiles)
    audit={
        "study":"Forecast model eligibility is a selected subset of global HOSTS resource butterflies",
        "frozen_input_sha256":hashes,
        "counts":{"original_resource_species":239,"forecast_candidate_species":225,
                  "forecast_modeled_species":155,"forecast_excluded_species":70,
                  "resource_species_not_in_forecast_candidates":14},
        "forecast_record_test_row_sum":modeled_test,
        "cohort_metrics":{k:describe(v) for k,v in cohorts.items()},
        "original_resource_239":describe(sample),
        "selected_minus_nonselected_within_original225":{
            key:{
                "difference":diff(overlap,key),
                "descriptive_species_bootstrap_95pct":bootstrap(overlap,key,seed=SEED+j),
            } for j,key in enumerate(("known_hosts","added_global","native_global",
                                      "fill_i","fill_m","delta"))
        },
        "conditional_selection_diagnostics":{
            "frozen_global_native_quartile_cutoffs":quantiles,
            "host_family_breadth_margin_preserved":preserve_strata(overlap,"delta",["family"],seed=SEED),
            "family_and_native_resource_quartile_preserved":preserve_strata(
                overlap,"delta",["family","native_quartile"],seed=SEED+1000
            ),
        },
        "forecast_source_lock":{
            "modeled155_temporal_model_scores_only_aggregate":True,
            "modeled155_species_specific_M3_M5_MSE_or_fitted_predictions_in_this_committed_source":False,
            "impossible_from_committed_source":"Correlation of individual butterfly host expansion with individual forecast MSE gain",
            "prior_recording_results":"Known species future and held spatial MSE improved by 4.12% and 3.00% in different evaluation scopes; no unseen-species gain",
        },
        "limits":[
            "Selection is fixed by original recording support and complete environmental/host inputs; groups are not randomized",
            "The selected and excluded groups share source diet data and phylogeny; bootstrap/permutation are descriptive and not causal/evolutionary hypothesis tests",
            "225 forecast candidates and 239 resource butterflies differ by 14 taxa, which must not be silently counted as unsuccessful forecasts",
            "Source forecast metrics are aggregated over species, and cannot be converted into individual species MSE gains",
            "Regional potential plant resources are not actual larval feeding, butterfly establishment or abundance",
            "The nonselected 70 have insufficient initial recording support, so unobserved forecasting outcomes cannot be imputed to them",
            "Host-breadth and global-native-resource stratification reduce one source of selection but do not control geography or observation effort fully",
        ],
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    with (a.out/"239_species_source_crosswalk.csv").open("w",newline="",encoding="utf-8") as file:
        head=[x for x in sample[0] if x!="modelled"]
        writer=csv.DictWriter(file,fieldnames=head);writer.writeheader()
        writer.writerows([{k:v for k,v in row.items() if k!="modelled"} for row in sample])
    print(json.dumps({
        "counts":audit["counts"],
        "cohorts":{k:{key:v[key] for key in (
            "n_species","host_catalogue_mean_known_species","global_mean_added_resource_regions",
            "species_equal_island_minus_mainland_fill_pp","island_more_fill_species",
            "mainland_more_fill_species","island_resource_gain_pct_native_denominator")}
            for k,v in audit["cohort_metrics"].items()},
        "selection_difference":audit["selected_minus_nonselected_within_original225"]["delta"],
        "conditional_null":{k:{kk:v[kk] for kk in (
            "n_matched_support_species","observed_difference","null_mean",
            "permutation_two_sided_centered_p")}
            for k,v in audit["conditional_selection_diagnostics"].items() if isinstance(v,dict)},
    },ensure_ascii=False))


if __name__=="__main__":
    main()
