#!/usr/bin/env python3
"""Scientific common-support correction of GloNAF vs HOSTS focal-host specificity.

The 201+ of 284 WGSRPD3 region prevalence stratum contains more documented
HOSTS hosts than available *non-documented* plants. No unrestricted prevalence
matching is possible there. We report all-host descriptive evidence, and run
a separate null on the complete frequency bins with actual nonhost support.
Do NOT relabel unmatched alien plants as confirmed nonhosts.
"""
from __future__ import annotations

from bisect import bisect_left
from collections import Counter,defaultdict
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from analyze_documented_host_specificity import (
    GLO_NAF,TDWG,HOSTS,PANEL,BINS,
    focal_hosts,glonaf_index,prevalence_bins,pair_jaccard,host_freq_matched_null,
)
from analyze_matched_glonaf_consumer_history import produce_matched_regions
from analyze_colonial_legacy_resource_beta import (
    load as load_colonial_butterfly,
    orthogonal_controls,residual,correlation,
)

ROOT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
LOCK=ROOT/"documented_host_specificity_lock_v2.json"
OUT=ROOT/"documented_host_specificity_v2.json"


def overlap_supported_mask(rawmask,frequency,buckets):
    """Select ALL members of jointly representable bins; never cherry-pick taxa."""
    raw_indices=[]
    work=rawmask
    while work:
        b=work & -work
        i=b.bit_length()-1
        if frequency[i]>0:
            raw_indices.append(i)
        work^=b
    if not raw_indices:
        raise RuntimeError("No source-active HOSTS taxa for prevalence matching")
    true_set=set(raw_indices)
    focal_bin=defaultdict(list)
    for i in raw_indices:
        focal_bin[bisect_left(BINS,int(frequency[i]))].append(i)
    kept=[]
    dropped=[]
    source={}
    for group,actual in sorted(focal_bin.items()):
        available=[i for i in buckets.get(group,[]) if i not in true_set]
        supported=len(available)>=len(actual)
        source[str(group)]={
          "focal_HOSTS_source_taxa":len(actual),
          "alternative_not_documented_source_taxa":len(available),
          "prevalence_upper_bound":BINS[group],
          "support":supported,
        }
        (kept if supported else dropped).extend(actual)
    if not kept:
        raise RuntimeError("No source-supported GloNAF frequency bins")
    allowed=0
    for i in kept:
        allowed |= 1<<i
    return allowed,{
       "supported_frequency_bins":sorted(int(k) for k,v in source.items() if v["support"]),
       "excluded_frequency_bins":sorted(int(k) for k,v in source.items() if not v["support"]),
       "eligible_original_taxa":len(raw_indices),
       "retained_comparable_taxa":len(kept),
       "excluded_taxa_due_to_no_pseudohost_support":len(dropped),
       "source_frequency_bins":source,
    }


def association(mask,full_masks,other_masks,data):
    """Same covariate sets and GloNAF source regions for each host definition."""
    h=[m&mask for m in full_masks]
    ia0,ib0=data["ia"],data["ib"]
    valid=np.fromiter((bool(h[int(i)]|h[int(j)]) for i,j in zip(ia0,ib0)),
                      bool,count=len(ia0))
    if int(valid.sum())<20000:
        raise RuntimeError("Known host subcohort has insufficient fixed region-pair support")
    ia,ib=ia0[valid],ib0[valid]
    x,no_host_pair=pair_jaccard(h,ia,ib)
    alien_background,missing_background=pair_jaccard(
        other_masks,ia,ib,undefined_zero=True)
    if no_host_pair or not np.all(np.isfinite(x)):
        raise RuntimeError("Unexpected undefined HOSTS host pair support")
    rich=np.asarray([m.bit_count() for m in full_masks])
    minrich=np.log1p(np.minimum(rich[ia],rich[ib]))
    maxrich=np.log1p(np.maximum(rich[ia],rich[ib]))
    covar=[q[valid] for q in data["control_sets"]["colonial_vs_resource_beta_after_flora"]]
    colonial=data["exposure"][valid]
    y=data["resource_delta"][valid]
    common=orthogonal_controls(covar+[colonial,minrich,maxrich,alien_background])
    yy=residual(y,common)
    xr=residual(x,common)
    rho=correlation(xr,yy)
    if rho is None:
        raise RuntimeError("Documented host prediction is constant after controls")
    other=orthogonal_controls(covar+[colonial,minrich,maxrich,x])
    other_rho=correlation(residual(alien_background,other),residual(y,other))
    stats={
      "tested_region_nodes":len(data["ids"]),
      "tested_region_pairs":len(ia),
      "omitted_pairs_host_both_empty":int((~valid).sum()),
      "mean_host_jaccard":float(np.mean(x)),
      "mean_other_alien_jaccard":float(np.mean(alien_background)),
      "mean_butterfly_resource_beta_delta":float(np.mean(y)),
      "observed_documented_host_rho_conditioning_on_other_aliens":rho,
      "other_alien_rho_conditioning_on_documented_hosts":other_rho,
      "other_alien_pairs_both_empty_treated_zero":int(missing_background),
    }
    return stats,{"basis":common,"response_residual":yy,
                  "ia":ia,"ib":ib,"observed":rho}


def analyse(draws=199,seed=31031):
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock["schema"]!="chocho_GloNAF_documented_HOSTS_prevalence_common_support_lock_v2":
        raise RuntimeError("Common-support design is not source-locked")
    data=load_colonial_butterfly()
    if len(data["ids"])!=284:
        raise RuntimeError("Source fixed panel not 284 TDWG3 regions")
    source_masks,_,flora_audit=produce_matched_regions(GLO_NAF,TDWG,(0.05,))
    full_masks=[source_masks[.05][code] for code in data["ids"]]
    host_names,host_audit=focal_hosts(HOSTS,PANEL)
    canonical,source_name_audit=glonaf_index(GLO_NAF)
    if source_name_audit["distinct_glonaf_names"]!=flora_audit["glonaf"]["distinct_species"]:
        raise RuntimeError("Original source bit indexing mismatch")
    matched=sorted(host_names.intersection(canonical))
    positions={i for n in matched for i in canonical[n]}
    true=0
    for j in positions:
        true |= 1<<j
    if len(matched)<20:
        raise RuntimeError("Named HOSTS–GloNAF plant intersection insufficient")
    other=[m&~true for m in full_masks]
    frequency,buckets=prevalence_bins(
         full_masks,source_name_audit["distinct_glonaf_names"])
    full_stats,_=association(true,full_masks,other,data)
    restricted,support=overlap_supported_mask(true,frequency,buckets)
    restricted_stats,internals=association(restricted,full_masks,other,data)
    fake,matched_profile,available,focal=host_freq_matched_null(
        restricted,full_masks,buckets,frequency,seed=seed,draws=draws)
    null=[]
    undef=[]
    for mask in fake:
        mm=[v&mask for v in full_masks]
        x,missing=pair_jaccard(
             mm,internals["ia"],internals["ib"],undefined_zero=True)
        rho=correlation(residual(x,internals["basis"]),
                        internals["response_residual"])
        if rho is None:
            raise RuntimeError("Geography-matched pseudohost Jaccard unidentifiable")
        null.append(rho)
        undef.append(missing)
    observed=restricted_stats["observed_documented_host_rho_conditioning_on_other_aliens"]
    return {
       "schema":"chocho_documented_HOSTS_focal_host_GloNAF_beta_common_support_v2",
       "status":"SOURCE_MATCHED_LOWER_BOUND_FOCAL_HOST_SPECIFICITY_SENSITIVITY",
       "source":{
          "host":"fixed HOSTS 808e0b869f9ec1adf8efff87cf6a395adda103e0",
          "glonaf":"Lenzner 2022, Zenodo 10.5281/zenodo.6640264",
          "method_lock_sha256":hashlib.sha256(LOCK.read_bytes()).hexdigest(),
          "butterfly_focal_species":239,"WGSRPD3_regions":284
       },
       "joined_names":{
          **host_audit,**source_name_audit,
          "focal_HOSTS_binomial_names_found_in_GloNAF":len(matched),
          "source_raw_bit_positions_for_known_HOSTS_candidates":len(positions),
          "source_region_active_known_taxa":sum(frequency[i]>0 for i in positions),
          "source_focal_HOSTS_names_not_found_in_GloNAF":len(host_names-set(matched))
       },
       "all_documented_host_case_without_strict_matched_null":full_stats,
       "overlap_only_prevalence_support":support,
       "matched_frequency_control":restricted_stats,
       "matched_null":{
          "draws":draws,"seed":seed,
          "retained_raw_source_taxa_by_prevalence_bin":matched_profile,
          "null_available_source_taxa":available,
          "focal_known_hosts_excluded_from_control_pool":focal,
          "median":float(np.median(null)),
          "q025":float(np.quantile(null,.025)),
          "q975":float(np.quantile(null,.975)),
          "observed":observed,
          "upper_tail":(1+sum(x>=observed for x in null))/(len(null)+1),
          "two_sided_abs_tail":(1+sum(abs(x)>=abs(observed) for x in null))/(len(null)+1),
          "mean_fake_both_empty_pairs_coded_zero":float(np.mean(undef)),
          "maximum_fake_both_empty_pairs_coded_zero":int(max(undef)),
       },
       "limits":[
          "Focal HOSTS aliases are exact normalized Genus epithet; synonym resolution remains incomplete.",
          "Recorded in HOSTS is not equivalent to suitable local host; noncatalogued does NOT imply true nonhost.",
          "Highly widespread known-host bin(s) lacked sufficient alternative plants and were removed as WHOLE strata from strict matched test.",
          "All-host and common-support tests have different plant predictors and potentially pair support; cannot claim full-taxon generality from overlap test.",
          "GloNAF source and WCVP source measure different plant universes and modern distribution, not pre-colonial time changes.",
          "Potential resource Jaccard based on same HOSTS input creates source-linked association; independent field butterfly occurrence is missing.",
          "Plant bins preserve only coarse 284-region occurrence frequency; neither plant family nor phylogeny nor spatial autocorrelation is conserved.",
          "Fake-host double-zero pairs are coded zero and reported, not strictly comparable to a fully positive-union observed host pair.",
          "All pair endpoints reused; tail fractions are exploratory not causal independent-sample p-values.",
       ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--draws",type=int,default=199)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    if not(19<=args.draws<=999):
        p.error("19..999 pseudohost draws")
    out=analyse(args.draws)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in out.items() if k!="limits"},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
