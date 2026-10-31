#!/usr/bin/env python3
"""Audit WHICH documented GloNAF HOSTS taxa drive butterfly-resource geography.

The previous source-matched 284-region analysis (v2) yielded two associations:
  - all source-identified documented hosts, Jaccard partial rho about +0.156,
  - prevalence-matchable documented hosts, partial rho about +0.032.
This script decomposes the observed association WITHOUT picking favorable taxa.
A complete prevalence bin is removed at a time; common-support host-richness
and non-documented-plant controls stay frozen.

Scope: catalogued potential larval resources only; no measured butterfly
occupancy, no colonial transport or host usage, and the outcome partly
shares WCVP/HOSTS source knowledge with the predictor.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from analyze_documented_host_specificity import (
    BINS,GLO_NAF,TDWG,HOSTS,PANEL,focal_hosts,glonaf_index,prevalence_bins,
    pair_jaccard,
)
from analyze_documented_host_specificity_v2 import (
    association,overlap_supported_mask,
)
from analyze_matched_glonaf_consumer_history import produce_matched_regions
from analyze_colonial_legacy_resource_beta import (
    load as load_colonial_butterfly,orthogonal_controls,residual,correlation,
)

ROOT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
OUT=ROOT/"documented_host_prevalence_drivers_v1.json"
V2=ROOT/"documented_host_specificity_v2.json"
ALREADY_REPLAYED_BINS=ROOT/"documented_host_prevalence_leave_bin_out_v1.json"


def split_by_prevalence(true_mask,frequency,bins=BINS):
    """Complete taxon bins, with no response-informed species selection."""
    groups=defaultdict(int)
    indices=defaultdict(list)
    t=true_mask
    while t:
        bit=t&-t
        i=bit.bit_length()-1
        t^=bit
        if frequency[i]<=0:
            continue
        key=bisect_left(bins,int(frequency[i]))
        groups[key] |= bit
        indices[key].append(i)
    if not groups:
        raise RuntimeError("No active documented hosts in the tested regions")
    masks=list(groups.values())
    if sum(mask.bit_count() for mask in masks)!=sum(len(v) for v in indices.values()):
        raise RuntimeError("Host prevalence bins not mutually exclusive")
    if any(a&b for i,a in enumerate(masks) for b in masks[i+1:]):
        raise RuntimeError("Host source IDs assigned to multiple prevalence bins")
    return dict(groups),{k:len(v) for k,v in indices.items()}


def reindex_response_pairs(data,keep):
    """Keep original 284 region indices; filter only dependent pair endpoints."""
    if keep.dtype!=bool or len(keep)!=len(data["ia"]):
        raise ValueError("Pair mask must match fixed 284-region comparisons")
    return {
        **data,
        "ia":data["ia"][keep],
        "ib":data["ib"][keep],
        "exposure":data["exposure"][keep],
        "resource_delta":data["resource_delta"][keep],
        "control_sets":{
            k:[c[keep] for c in val] for k,val in data["control_sets"].items()
        },
    }


def pair_predictor(masks,ia,ib):
    return pair_jaccard(masks,ia,ib,undefined_zero=True)


def unique_partial(x,y,control):
    basis=orthogonal_controls(control)
    result=correlation(residual(x,basis),residual(y,basis))
    if result is None:
        raise RuntimeError("Predictor became constant on fixed support")
    return float(result)


def analyse():
    v2=json.loads(V2.read_text(encoding="utf8"))
    if v2["schema"]!="chocho_documented_HOSTS_focal_host_GloNAF_beta_common_support_v2":
        raise RuntimeError("Source-matched v2 receipt changed")
    data=load_colonial_butterfly()
    if len(data["ids"])!=284 or len(data["ia"])!=40171:
        raise RuntimeError("Original WGSRPD3 common source drift")
    source,_,audit=produce_matched_regions(GLO_NAF,TDWG,(.05,))
    full=[source[.05][code] for code in data["ids"]]
    hostnames,host_stats=focal_hosts(HOSTS,PANEL)
    canonical,names=glonaf_index(GLO_NAF)
    original={idx for host in hostnames & canonical.keys() for idx in canonical[host]}
    rawmask=sum(1<<idx for idx in sorted(original))
    frequency,buckets=prevalence_bins(full,names["distinct_glonaf_names"])
    # 58 of the 1312 source-matched bit IDs never appear inside the
    # fixed 284-region subset. They are NOT the 86 widespread host IDs
    # dropped because frequency-matched alternative taxa are insufficient.
    active_mask=sum(1<<idx for idx in original if frequency[idx]>0)
    masks,bin_counts=split_by_prevalence(active_mask,frequency)
    source_active=sum(bin_counts.values())
    supportmask,support=overlap_supported_mask(active_mask,frequency,buckets)
    excluded=active_mask & ~supportmask
    if (source_active!=v2["joined_names"]["source_region_active_known_taxa"]
            or excluded.bit_count()!=v2["overlap_only_prevalence_support"]["excluded_taxa_due_to_no_pseudohost_support"]):
        raise RuntimeError("Documented host source identity no longer reproduces v2")
    others=[x&~rawmask for x in full]

    all_case,internals=association(rawmask,full,others,data)
    supported_case,_=association(supportmask,full,others,data)
    if abs(all_case["observed_documented_host_rho_conditioning_on_other_aliens"]-
           v2["all_documented_host_case_without_strict_matched_null"]["observed_documented_host_rho_conditioning_on_other_aliens"])>1e-8:
        raise RuntimeError("Exact all-host partial rho failed fixed-source consistency")
    if abs(supported_case["observed_documented_host_rho_conditioning_on_other_aliens"]-
           v2["matched_frequency_control"]["observed_documented_host_rho_conditioning_on_other_aliens"])>1e-8:
        raise RuntimeError("Matched support source v2 rho changed")

    # Previously source-recomputed full-bin deletion is an independent
    # branch receipt. Validate that result instead of repeating the identical
    # ~8 partial-rank fits and risking analysis version drift.
    influence=json.loads(ALREADY_REPLAYED_BINS.read_text(encoding="utf8"))
    if influence["schema"]!="chocho_posthoc_documented_host_prevalence_influence_v1":
        raise RuntimeError("Prior response-blind all-bin diagnostic missing")
    if abs(influence["baseline_host_partial_rank_rho"]-
           all_case["observed_documented_host_rho_conditioning_on_other_aliens"])>1e-9:
        raise RuntimeError("Previous drop-bin results do not match source")
    leave_bin={}
    by_bin={}
    for row in influence["influence_by_bin"]:
        key=int(row["removed_prevalence_bin"])
        if key not in masks or int(row["excluded_original_matched_focal_host_taxa"])!=bin_counts[key]:
            raise RuntimeError("Previously computed complete prevalence bin changed")
        leave_bin[str(key)]=row
        by_bin[str(key)]={"prevalence_upper_bound":BINS[key],"source_taxa":bin_counts[key]}
    if set(map(int,leave_bin))!=set(masks):
        raise RuntimeError("Full source prevalence class coverage changed")

    # Unique high-prevalence and retained (low/mid) host information:
    # exactly the same 40,171 source region pairs, *symmetric* zero handling.
    # This avoids using a different response sample just because one class
    # has both-empty plant sets.
    ia,ib=data["ia"],data["ib"]
    all_host=[m&rawmask for m in full]
    low_host=[m&supportmask for m in full]
    high_host=[m&excluded for m in full]
    all_j,missing_all=pair_predictor(all_host,ia,ib)
    low_j,missing_low=pair_predictor(low_host,ia,ib)
    high_j,missing_high=pair_predictor(high_host,ia,ib)
    other_j,missing_other=pair_predictor(others,ia,ib)
    full_rich=np.asarray([m.bit_count() for m in full],float)
    richness=[
        np.log1p(np.minimum(full_rich[ia],full_rich[ib])),
        np.log1p(np.maximum(full_rich[ia],full_rich[ib])),
    ]
    base=(data["control_sets"]["colonial_vs_resource_beta_after_flora"]+
          [data["exposure"],*richness,other_j])
    unique={
        "high_prevalence_unique_after_retained_known_hosts":unique_partial(
            high_j,data["resource_delta"],base+[low_j]),
        "low_mid_prevalence_unique_after_excluded_widespread_hosts":unique_partial(
            low_j,data["resource_delta"],base+[high_j]),
        "all_documented_hosts_fixed_pair_support":unique_partial(
            all_j,data["resource_delta"],base),
        "non_documented_GloNAF_alien_plants_given_both_host_classes":unique_partial(
            other_j,data["resource_delta"],
            data["control_sets"]["colonial_vs_resource_beta_after_flora"]+
            [data["exposure"],*richness,low_j,high_j]),
    }

    # Repeated-pair dependency: sensitivity to removing entire WGSRPD Level1
    # world division, not a binomial 'n = 40,171 independent pairs'.
    world=np.asarray(data["continent"])
    if len(world)!=len(data["ids"]):
        raise RuntimeError("Missing world-division labels")
    divisions=sorted(set(map(str,world)))
    leave_div={}
    for level in divisions:
        included=np.array([str(world[i])!=level and str(world[j])!=level
                           for i,j in zip(ia,ib)],dtype=bool)
        subset=reindex_response_pairs(data,included)
        full_fold,_=association(rawmask,full,others,subset)
        low_fold,_=association(supportmask,full,others,subset)
        leave_div[level]={
            "excluded_wgsrpd_level1":level,
            "remaining_region_pairs":int(included.sum()),
            "remaining_region_nodes":int(np.sum(world.astype(str)!=level)),
            "all_documented_host_partial":full_fold[
                "observed_documented_host_rho_conditioning_on_other_aliens"],
            "frequency_matchable_host_partial":low_fold[
                "observed_documented_host_rho_conditioning_on_other_aliens"],
        }

    return {
        "schema":"chocho_documented_focal_HOSTS_284region_prevalence_driver_v1",
        "status":"SOURCE_LOCKED_POSTHOC_PREVALENCE_BIN_SENSITIVITY_NOT_CAUSAL",
        "source":{
            "v2_receipt_sha256":hashlib.sha256(V2.read_bytes()).hexdigest(),
            "full_bin_drop_diagnostic_sha256":hashlib.sha256(
                ALREADY_REPLAYED_BINS.read_bytes()).hexdigest(),
            "glonaf_zenodo_md5":"16f5585d2500676951a4b664bbf0ffde",
            "wgsrpd3_region_nodes":284,
            "non_independent_region_pairs":len(ia),
            "butterfly_panel_species":239,
            "source_host_names":host_stats["focal_host_names"],
            "source_glonaf_original_taxa":names["distinct_glonaf_names"],
            "geographic_source_test":"same TDWG 5%/50% polygon overlap",
        },
        "source_coverage":{
            "source_active_documented_host_taxa":source_active,
            "all_source_matched_host_taxa_including_outside_284_regions":len(original),
            "source_matched_host_taxa_not_active_in_284_regions":len(original)-source_active,
            "prevalence_matchable_hosts":supportmask.bit_count(),
            "excluded_due_to_unmatched_alternative_prevalence":excluded.bit_count(),
            "prevalence_strata":by_bin,
        },
        "all_documented_host_case":all_case,
        "prevalence_matchable_case":supported_case,
        "drop_entire_prevalence_bin_one_at_a_time":leave_bin,
        "same_fixed_pairs_unique_host_strata":{
            "metric":"rank partial correlation of resource Jaccard delta",
            "n_non_independent_region_pairs":len(ia),
            "all_host_both_empty_zero_coded":missing_all,
            "retained_host_both_empty_zero_coded":missing_low,
            "excluded_widespread_host_both_empty_zero_coded":missing_high,
            "other_alien_both_empty_zero_coded":missing_other,
            **unique,
        },
        "leave_one_world_division_out":leave_div,
        "limits":[
            "The butterfly-resource outcome uses the same documented HOSTS network as the known-host predictor; source coupling does not independently establish host-specific biological effects.",
            "Non-documented GloNAF alien plant means not a known focal HOSTS host, NOT a demonstrated nonhost.",
            "High-prevalence known host class had too few non-documented alternatives for one-to-one matching; its association is descriptive rather than a prevalence-matched causal test.",
            "Leave-one-bin-out changes total host richness and Jaccard denominator; rho differences are descriptive, not additive contributions or a mediation decomposition.",
            "Both-empty pair assignments are zero only in the explicitly labeled same-fixed-pairs comparison.",
            "The same 284 region nodes enter 40,171 correlated region-pairs; no pairwise iid p-values are reported.",
            "GloNAF mapped current alien flora and WCVP native/current reconstruction are not historical observations made in the colonial period.",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=OUT)
    args=ap.parse_args()
    result=analyse()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    def cleaner(v):
        if isinstance(v,np.integer):return int(v)
        if isinstance(v,np.floating):return float(v)
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError("Nonserializable output "+str(type(v)))
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                      default=cleaner)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k!="limits"},
                     indent=2,ensure_ascii=False,default=cleaner))


if __name__=="__main__":
    main()
