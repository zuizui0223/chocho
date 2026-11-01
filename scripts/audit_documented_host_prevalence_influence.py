#!/usr/bin/env python3
"""Post-hoc diagnostic: which GloNAF source-prevalence band drives host beta?

The diagnosed unsupported bin 7 corresponds to 101-200 of 284 TDWG3 regions.
For each bin drop ALL documented 239-butterfly HOSTS raw plant IDs in the bin,
never pick species by their response, and compute the same region-pair partial
rank statistic. Complete source 284x284 cohort and covariates are retained;
if a changed host portfolio has some pairwise double-empty opportunities, all
matching support changes are reported. Not a biological p-value or independent
HOSTS outcome.
"""
from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
import argparse
import hashlib
import json
from pathlib import Path

from analyze_documented_host_specificity import (
    GLO_NAF,TDWG,HOSTS,PANEL,BINS,focal_hosts,glonaf_index,prevalence_bins
)
from analyze_documented_host_specificity_v2 import association
from analyze_matched_glonaf_consumer_history import produce_matched_regions
from analyze_colonial_legacy_resource_beta import load as load_colonial_butterfly

BASE=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
V2=BASE/"documented_host_specificity_v2.json"
OUT=BASE/"documented_host_prevalence_leave_bin_out_v1.json"


def compute():
    documented=json.loads(V2.read_text(encoding="utf8"))
    if documented.get("schema")!="chocho_documented_HOSTS_focal_host_GloNAF_beta_common_support_v2":
        raise RuntimeError("v2 original receipt not source-consistent")
    data=load_colonial_butterfly()
    source_masks,_,origin=produce_matched_regions(GLO_NAF,TDWG,(.05,))
    full=[source_masks[.05][code] for code in data["ids"]]
    host_names,_=focal_hosts(HOSTS,PANEL)
    aliases,proof=glonaf_index(GLO_NAF)
    if proof["distinct_glonaf_names"]!=origin["glonaf"]["distinct_species"]:
        raise RuntimeError("Bit indices no longer preserved across original GloNAF")
    ids={j for n in host_names.intersection(aliases) for j in aliases[n]}
    mask=0
    for j in ids:
        mask |= (1<<j)
    if len(ids)!=1312:
        raise RuntimeError("HOSTS/GloNAF raw matched plant IDs drift")
    freq,buckets=prevalence_bins(full,proof["distinct_glonaf_names"])
    other=[v&~mask for v in full]
    baseline,_=association(mask,full,other,data)
    orig= documented["all_documented_host_case_without_strict_matched_null"]
    if abs(baseline["observed_documented_host_rho_conditioning_on_other_aliens"]
           -orig["observed_documented_host_rho_conditioning_on_other_aliens"])>1e-9:
        raise RuntimeError("Baseline host correlation differs from source-pinned v2")
    active={j for j in ids if freq[j]>0}
    groups=defaultdict(list)
    for j in sorted(active):
        groups[bisect_left(BINS,int(freq[j]))].append(j)
    stats=[]
    for group,plantids in sorted(groups.items()):
        removal=0
        for j in plantids:
            removal |= (1<<j)
        kept=mask&~removal
        try:
            result,_=association(kept,full,other,data)
            outcome={
                "state":"IDENTIFIED_ON_POSTHOC_REMAINING_SOURCE",
                "partial_host_rho":result[
                    "observed_documented_host_rho_conditioning_on_other_aliens"],
                "partial_other_alien_rho":result[
                    "other_alien_rho_conditioning_on_documented_hosts"],
                "n_supported_region_pairs":result["tested_region_pairs"],
                "pairs_both_no_documented_host":result["omitted_pairs_host_both_empty"],
                "delta_host_rho_relative_full":result[
                    "observed_documented_host_rho_conditioning_on_other_aliens"]
                     - baseline["observed_documented_host_rho_conditioning_on_other_aliens"],
            }
        except RuntimeError as e:
            outcome={"state":"UNIDENTIFIABLE","reason":str(e)[:180]}
        stats.append({
            "removed_prevalence_bin":int(group),
            "source_284region_prevalence_range_inclusive":[
                1 if group==0 else BINS[group-1]+1,
                BINS[group]],
            "excluded_original_matched_focal_host_taxa":len(plantids),
            **outcome
        })
    return {
        "schema":"chocho_posthoc_documented_host_prevalence_influence_v1",
        "status":"LEAVE_ENTIRE_SOURCE_PREVALENCE_BINS_OUT_NOT_INDEPENDENT_ECOLOGY",
        "source_fixed_239_host_aliases":len(host_names.intersection(aliases)),
        "source_active_documented_host_IDs":len(active),
        "GloNAF_284region_prevalence_bins":list(BINS),
        "baseline_host_partial_rank_rho":baseline[
            "observed_documented_host_rho_conditioning_on_other_aliens"],
        "baseline_other_alien_partial_rank_rho":baseline[
            "other_alien_rho_conditioning_on_documented_hosts"],
        "baseline_region_pairs":baseline["tested_region_pairs"],
        "influence_by_bin":stats,
        "v2_source_receipt_sha256":hashlib.sha256(V2.read_bytes()).hexdigest(),
        "limits":[
            "Post-hoc structural influence diagnostic; several hypotheses inspected.",
            "Dropping one frequency class changes number of source host taxa and sometimes 284-region pair support.",
            "Original 284-region HOSTS-derived resource beta is not independent of HOSTS host identities.",
            "The 101-200 geography band has fewer unrecorded-as-host alternative taxa, so no prevalence-matched null for that whole band.",
            "GloNAF source names are normalized to Genus epithet aliases but do not resolve all WCVP synonym ambiguities.",
            "No plant chronology, butterfly feeding, colonial shipping or butterfly occupancy is measured.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    result=compute()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                           encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k!="limits"},indent=2))


if __name__=="__main__":
    main()
