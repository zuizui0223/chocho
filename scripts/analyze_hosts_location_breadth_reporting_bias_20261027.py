#!/usr/bin/env python3
"""Catalogue reporting-bias audit: original 239 HOSTS Location-origin evidences.

Statistical response is source documented introduction-area host reports,
NOT adopted host, consumer introduction, larval survival, causal niche change.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SCHEMA = "chocho_239_hosts_location_breadth_reporting_bias_20261027_v01"
MODEL_SOURCE = "chocho_original_hosts_location_wcvp_plant_origin_20261027_v01"
SEED = 20261027
BOOT = 499
RIDGE = 1e-3
OUTCOMES = {
    "introduced_plant_region_host_report": "source_mapped_introduced_triples",
    "native_plant_region_host_report": "source_mapped_native_triples",
    "WCVP_unlisted_plant_region_host_report": "source_WCVP_unlisted_triples",
}


def logit_inverse(x):
    x = np.asarray(x, dtype=float)
    return 1 / (1 + np.exp(-np.clip(x, -35, 35)))


def regression(y, x, *, ridge=RIDGE, tolerance=1e-8, maxiter=100):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or len(y) != len(x) or x.shape[1] != 4:
        raise ValueError("Unexpected original host-recording logistic design")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Nonfinite original regression source")
    if not set(np.unique(y)) <= {0.0, 1.0} or y.sum() < 3 or (
            len(y) - y.sum() < 3):
        raise ValueError("Binary source-index outcome has too few positives or negatives")
    if ridge <= 0:
        raise ValueError("Penalization must be prespecified positive")
    penalty = np.diag([0, ridge, ridge, ridge])
    beta = np.zeros(x.shape[1])
    p0 = min(.999, max(.001, np.mean(y)))
    beta[0] = math.log(p0/(1-p0))
    for it in range(maxiter):
        p = logit_inverse(x @ beta)
        g = x.T @ (y-p) - penalty @ beta
        w = p*(1-p)
        H = (x.T * w) @ x + penalty
        step = np.linalg.solve(H, g)
        original_obj = float(np.sum(y*np.log(np.maximum(p, 1e-15)) +
                                    (1-y)*np.log(np.maximum(1-p, 1e-15))) -
                             .5*beta @ penalty @ beta)
        scale = 1.
        accepted = False
        for _ in range(25):
            candidate = beta + scale*step
            pp = logit_inverse(x @ candidate)
            obj = float(np.sum(y*np.log(np.maximum(pp, 1e-15)) +
                               (1-y)*np.log(np.maximum(1-pp, 1e-15))) -
                        .5*candidate @ penalty @ candidate)
            if obj >= original_obj - 1e-8:
                accepted = True
                break
            scale *= .5
        if not accepted:
            raise ValueError("Logistic Newton descent failed")
        beta = candidate
        if np.linalg.norm(scale*step) < tolerance:
            break
    else:
        raise ValueError("Logistic IRLS failed to converge")
    return beta


def get_source_rows(source):
    if source.get("schema") != MODEL_SOURCE:
        raise ValueError("Not exact WCVP-host accepted report geography source")
    primary = source["result"]["primary_and_sensitivity"]["0.95"]
    secondary = source["result"]["primary_and_sensitivity"]["0.8"]
    records = primary["per_butterfly_aggregate_only"]
    if len(records) != 239 or len({r["species"] for r in records}) != 239:
        raise ValueError("Not original complete 239 source butterfly sample")
    if len(secondary["per_butterfly_aggregate_only"]) != 239:
        raise ValueError("Original geographic area sensitivity was lost")
    if sum(r["source_mapped_introduced_triples"] for r in records) != 179 or \
       sum(r["source_mapped_native_triples"] for r in records) != 399 or \
       sum(r["source_WCVP_unlisted_triples"] for r in records) != 303:
        raise ValueError("Source per-species reporting totals do not match original 95% source margins")
    if sum(r["has_any_introduced_area_host_report"] for r in records) != 62:
        raise ValueError("Source introduced-area documented butterfly count changed")
    if {r["species"] for r in records} != {
            r["species"] for r in secondary["per_butterfly_aggregate_only"]}:
        raise ValueError("Threshold sensitivity source species changed")
    return sorted(records, key=lambda r: r["species"])


def design(rows, *, feature="family", standard=None):
    if feature not in ("family", "species_richness"):
        raise ValueError("Unfrozen dietary breadth feature")
    feature_label = ("original_source_host_family_count" if feature=="family"
                     else "original_source_resolved_host_species")
    columns = np.array([
        [np.log1p(float(r[feature_label])),
         np.log1p(float(r["source_unique_accepted_host_location_triples"])),
         np.log1p(float(r["original_source_introduced_added_resource_units"]))]
        for r in rows
    ], dtype=float)
    if standard is None:
        means = columns.mean(axis=0)
        sds = columns.std(axis=0, ddof=1)
    else:
        means, sds = standard
    if np.any(sds < 1e-10):
        raise ValueError("Source predictor unidentifiable due zero cohort variance")
    scaled = (columns-means) / sds
    x = np.column_stack([np.ones(len(rows)), scaled])
    return x, (means, sds)


def bootstrap_by_genus(rows, *, field="family", indices=None, standard=None,
                       draws=BOOT, seed=SEED):
    if indices is None:
        indices = np.arange(len(rows), dtype=int)
    indices = np.array(indices, dtype=int)
    x, params = design(rows,feature=field,standard=standard)
    genera=defaultdict(list)
    for k in indices:
        genera[rows[k]["species"].split()[0]].append(int(k))
    names=sorted(genera)
    if len(names)<10:
        raise ValueError("Insufficient genus-cluster replication")
    yy = {
        name: np.array([int(r[field]>0) for r in rows], dtype=float)
        for name, field in OUTCOMES.items()
    }
    original = {
        name: regression(y[indices], x[indices]) for name, y in yy.items()
    }
    rng = np.random.default_rng(seed)
    estimates = defaultdict(list)
    convergence = Counter()
    for b in range(draws):
        which=rng.integers(0, len(names), len(names))
        sampled = np.array([k for pos in which for k in genera[names[int(pos)]]],
                           dtype=int)
        bb = {}
        for name, y in yy.items():
            try:
                bb[name]=float(regression(y[sampled], x[sampled])[1])
            except (np.linalg.LinAlgError, ValueError):
                convergence[name]+=1
                bb[name]=None
            if bb[name] is not None and math.isfinite(bb[name]):
                estimates[name].append(bb[name])
        for a,bname in (
            ("introduced_plant_region_host_report",
             "native_plant_region_host_report"),
            ("introduced_plant_region_host_report",
             "WCVP_unlisted_plant_region_host_report"),
        ):
            if bb.get(a) is not None and bb.get(bname) is not None:
                estimates[a+"__minus__"+bname].append(bb[a]-bb[bname])
    report = {}
    for key, values in estimates.items():
        if len(values)<.90*draws:
            report[key]={"status":"BOOTSTRAP_CONVERGENCE_BELOW90PCT",
                         "valid_draws":len(values)}
            continue
        lo,hi=np.quantile(values,[.025,.975])
        report[key]={
            "valid_draws":len(values),
            "log_odds_coefficient_percentile_95ci":[float(lo),float(hi)],
            "odds_ratio_interval":[float(np.exp(np.clip(lo,-50,50))),
                                   float(np.exp(np.clip(hi,-50,50)))],
        }
    labels = {
        name: {
            "source_positive_butterflies":int(np.sum(yy[name][indices])),
            "source_negative_or_unreported_butterflies":int(
                len(indices)-np.sum(yy[name][indices])),
            "log_odds_breadth_coefficient_per_original_SD":
                float(original[name][1]),
            "odds_ratio_per_original_SD":float(np.exp(original[name][1])),
            "bootstrap":report.get(name),
        }
        for name in OUTCOMES
    }
    for a,bname in (
        ("introduced_plant_region_host_report","native_plant_region_host_report"),
        ("introduced_plant_region_host_report","WCVP_unlisted_plant_region_host_report"),
    ):
        diff = float(original[a][1]-original[bname][1])
        labels[a+"__minus__"+bname]={
            "log_odds_coeff_difference":diff,
            "bootstrap":report.get(a+"__minus__"+bname)
        }
    return {
        "sample_size":len(indices),
        "source_genus_clusters":len(names),
        "source_breadth_predictor":field,
        "bootstrap_draws":draws,
        "penalty":RIDGE,
        "failed_genus_bootstrap_attempts":dict(convergence),
        "original_full239_standardization":{
            "log1p_predictor_mean":[float(v) for v in params[0]],
            "log1p_predictor_sd":[float(v) for v in params[1]]
        },
        "design_condition_number":float(np.linalg.cond(x[indices])),
        "source_predictor_pairwise_correlation":np.corrcoef(
            x[indices,1:],rowvar=False).round(6).tolist(),
        "outcome_models":labels,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-aggregate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args=parser.parse_args()
    source=json.loads(args.source_aggregate.read_text())
    rows=get_source_rows(source)
    original_mean_sd = design(rows)[1]
    primary=bootstrap_by_genus(rows,field="family",standard=original_mean_sd)
    rich=bootstrap_by_genus(rows,field="species_richness",seed=SEED+12)
    native_evidence=[i for i,r in enumerate(rows)
                     if r["source_mapped_native_triples"]>0]
    native_sensitivity=bootstrap_by_genus(
        rows,field="family",indices=native_evidence,
        standard=original_mean_sd,seed=SEED+101)
    inf_status = (
        "EXPLORATORY_DOCUMENTATION_BIAS_ONLY" if (
            primary["source_genus_clusters"]>=25 and
            primary["outcome_models"]["introduced_plant_region_host_report"][
                "source_positive_butterflies"]>=30
        ) else "UNIDENTIFIED_SOURCE_EVIDENCE_COARSE"
    )
    r={
        "schema":SCHEMA,
        "status":inf_status,
        "source_239_original_butterflies":239,
        "source_accepted_origin_archive_sha256":hashlib.sha256(
            args.source_aggregate.read_bytes()).hexdigest(),
        "source_original_introduced_plant_location_reports":62,
        "primary_host_family_breadth":primary,
        "host_species_richness_breadth_sensitivity":rich,
        "source_native_plant_area_report_positive_only":native_sensitivity,
        "biological_limits":[
            "These models predict whether an interaction is DOCUMENTED in plant introduced-area HOSTS Location reports, not a real newly adopted host.",
            "No source report is interpreted as a biological negative or absence.",
            "Original diet breadth is derived from the same globally compiled HOSTS source and therefore is not an independent ecological predictor.",
            "Geographic sampling and plant WCVP region availability are incompletely controlled by raw source record count and host opportunity.",
            "Bootstrap genera are not phylogenetically independent randomized host introductions.",
            "A similar breadth coefficient for WCVP native/unlisted records falsifies introduced-host-specific interpretation.",
            "Crop use, true butterfly invasiveness, feeding trials, climate niche colonization and larval fitness are not measured.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    pp=primary["outcome_models"]
    print(json.dumps({
        "status":inf_status,
        "introduced_documentation":pp["introduced_plant_region_host_report"],
        "native_documentation":pp["native_plant_region_host_report"],
        "unlisted_documentation":pp["WCVP_unlisted_plant_region_host_report"],
        "intro_minus_native":pp[
            "introduced_plant_region_host_report__minus__native_plant_region_host_report"],
        "source_native_positive_subset_n":len(native_evidence),
        "condition_number":primary["design_condition_number"],
    },sort_keys=True))


if __name__=="__main__":
    main()
