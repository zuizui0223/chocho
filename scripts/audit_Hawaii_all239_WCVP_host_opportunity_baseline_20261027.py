#!/usr/bin/env python3
"""Artificial host-richness matched resource-source baseline for 10 Hawaii aliens.

Original 239 butterfly globally documented HOSTS host portfolios.
Do 10 source-identified independently established Hawaii aliens have any
WCVP Hawaii INTRODUCED WILD accepted hosts beyond an exchangeability
reference with equal original documented host-list size strata?
This is NOT an ecological arrival/establishment p-value or negative census.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

from analyze_hosts_location_wcvp_origin_20261027 import (
    accepted_pairs,read_rows,sha256,wcvp_range
)
from audit_Hawaii_established_2024_butterfly_WCVP_host_coverage_20261027 import (
    classify_plants,HAWAII_WCVP,
)

SCHEMA="chocho_original239_Hawaii_intro_host_opportunity_background_20261027_v01"
STRATA=("1-5","6-20","21-50","51+")

def stratum(n):
    if n<=0:raise ValueError("Original species accepted HOSTS host list must be nonempty")
    if n<=5:return "1-5"
    if n<=20:return "6-20"
    if n<=50:return "21-50"
    return "51+"


def exact_hypergeom(total,support,draw):
    if not 0<=support<=total or not 0<=draw<=total:
        raise ValueError("Impossible source study host list matched taxon margins")
    return {
      k:math.comb(support,k)*math.comb(total-support,draw-k)/math.comb(total,draw)
      for k in range(max(0,draw-(total-support)),min(support,draw)+1)
    }


def exact_bin_preserving_reference(rows,positive):
    population=defaultdict(list)
    marked=Counter()
    for row in rows:
        population[stratum(row["source_accepted_host_count"])].append(row)
    for row in positive:
        marked[stratum(row["source_accepted_host_count"])]+=1
    mean=0.
    pmf={0:1.}
    source_bin=[]
    for label in STRATA:
        pop=population[label]
        n=len(pop)
        k=sum(row["hawaii_WCVP_introduced_wild_host_count"]>0 for row in pop)
        m=marked[label]
        if m>n:
            raise ValueError("Positive Hawaii species not in source richness stratum")
        dist=exact_hypergeom(n,k,m)
        mean += m*k/n if n else 0.
        source_bin.append({
            "original_source_host_breadth_bin":label,
            "all_original_239_source_species":n,
            "with_any_introduced_Hawaii_wild_host":k,
            "source_Zenodo_Hawaii_established_species_in_bin":m,
            "source_Zenodo_positive_with_introduced_Hawaii_host":sum(
                row["hawaii_WCVP_introduced_wild_host_count"]>0
                for row in positive
                if stratum(row["source_accepted_host_count"])==label),
            "source_mathematical_expected_with_introduced_host":m*k/n if n else 0.,
        })
        tmp=defaultdict(float)
        for i,pi in pmf.items():
            for j,pj in dist.items():
                tmp[i+j]+=pi*pj
        pmf=dict(tmp)
    mass=sum(pmf.values())
    if not math.isclose(mass,1,abs_tol=1e-9):
        raise ValueError("Artificial source matched host strata probabilities do not sum to one")
    observed=sum(x["hawaii_WCVP_introduced_wild_host_count"]>0 for x in positive)
    ordered=sorted(pmf)
    def quantile(q):
        c=0.
        for v in ordered:
            c+=pmf[v]
            if c>=q-1e-12:
                return v
        return ordered[-1]
    return {
        "original_source_bins":source_bin,
        "source_Zenodo_10_original_positive_aliens":len(positive),
        "source_observed_10_with_Hawaii_wild_introduced_host":observed,
        "source_all239_host_breadth_bin_matched_artificial_expected":
            float(mean),
        "source_artificial_95pct_interval":[quantile(.025),quantile(.975)],
        "source_artificial_right_tail_NOT_ecological_p":float(
            sum(p for k,p in pmf.items() if k>=observed)),
        "source_artificial_possible_support":[ordered[0],ordered[-1]],
        "reference_is_noncausal_missing_insect_arrival_exposure":True,
    }


def audit(root,positive):
    meta=json.loads((root/"SOURCE_MANIFEST.json").read_text())
    if meta.get("inputs_count")!=11:
        raise ValueError("Not original 11-member author-only source")
    for name,digest in meta["sha256"].items():
        if sha256(root/name)!=digest:
            raise ValueError("Original author WCVP constituent SHA changed")
    species=read_rows(root/"inputs/anthropogenic_species_metrics.csv")
    focal={r["species"] for r in species}
    if len(focal)!=239 or len(species)!=239:
        raise ValueError("Original 239 source cohort drift")
    accepted,qc,_=accepted_pairs(
        root/"inputs/original/insect_host_accepted.csv",focal)
    if qc["source_ambiguous_accepted_id_names"]:
        raise ValueError("Source accepted-host ambiguity")
    hosts=defaultdict(set)
    for (sp,raw),pid in accepted.items():
        hosts[sp].add(pid)
    if sum(len(v) for v in hosts.values())!=2596:
        raise ValueError("Original 2596 accepted butterfly-host links changed")
    native,current=wcvp_range(
        root/"inputs/original/native_extant_nondoubtful_wgsrpd3.csv",
        root/"inputs/original/extant_nondoubtful_wgsrpd3.csv")
    original10=positive["individual_source_positive_Hawaii_butterflies"]
    if len(original10)!=10:
        raise ValueError("Original 10 Hawaii alien record panel changed")
    ids={r["species"] for r in original10}
    if len(ids)!=10 or not ids <= focal:
        raise ValueError("Hawaii positive species not included in original 239 cohort")
    allrows=[]
    positives=[]
    compare={r["species"]:r for r in original10}
    for sp in sorted(focal):
        counts=classify_plants(hosts[sp],native,current)
        record={
            "species":sp,
            "source_accepted_host_count":len(hosts[sp]),
            "hawaii_WCVP_native_wild_host_count":counts["NATIVE_HAW_WILD"],
            "hawaii_WCVP_introduced_wild_host_count":counts["INTRODUCED_HAW_WILD"],
            "hawaii_source_unlisted_wild_host_count":counts["NOT_RECORDED_HAW_WILD"]
        }
        if sp in ids:
            prior=compare[sp]
            if (record["hawaii_WCVP_native_wild_host_count"],
                record["hawaii_WCVP_introduced_wild_host_count"],
                record["hawaii_source_unlisted_wild_host_count"])!=(
                prior["WCVP_Hawaii_wild_host_plants_native"],
                prior["WCVP_Hawaii_wild_host_plants_introduced"],
                prior["known_host_species_without_WCVP_wild_presence_Hawaii"]):
                raise ValueError("Previously verified ten Hawaii source host counts drifted")
            positives.append(record)
        allrows.append(record)
    status=Counter()
    for row in allrows:
        n=row["hawaii_WCVP_native_wild_host_count"]
        i=row["hawaii_WCVP_introduced_wild_host_count"]
        if n and i:status["BOTH_NATIVE_AND_INTRODUCED_WILD_HOSTS"]+=1
        elif n:status["ONLY_WCVP_NATIVE_WILD_HOSTS"]+=1
        elif i:status["ONLY_WCVP_INTRODUCED_WILD_HOSTS"]+=1
        else:status["NO_HOST_IN_WCVP_WILD_HAWAII_MAP"]+=1
    if sum(status.values())!=239:
        raise ValueError("Original 239 cohort host source coverage identity failed")
    ref=exact_bin_preserving_reference(allrows,positives)
    return {
        "schema":SCHEMA,
        "status":"SOURCE_ORIGINAL239_HAWAII_HOST_OPPORTUNITY_STRUCTURAL_COVERAGE_NOT_INVASION_FAILURE",
        "original_239_source_host_coverage":dict(status),
        "original_239_source_any_introduced_Hawaii_wild_host":sum(
            row["hawaii_WCVP_introduced_wild_host_count"]>0 for row in allrows),
        "original_239_source_any_native_Hawaii_wild_host":sum(
            row["hawaii_WCVP_native_wild_host_count"]>0 for row in allrows),
        "source_orig_239_host_count":len(allrows),
        "confirmed_10_positive_richness_stratified_reference":ref,
        "source_original_positive_total_accepted_host_plants":sum(
            x["source_accepted_host_count"] for x in positives),
        "source_original_all239_total_accepted_host_plants":sum(
            x["source_accepted_host_count"] for x in allrows),
        "scientific_gates":[
            "The 229 other source butterflies lack matching documented 2024 Hawaii ALIEN ESTABLISHMENT records, but are NOT evidenced immigration failures.",
            "If all 239 species already have introduced host plants in HAW, the original 10/10 documented established status has no structural discriminating power.",
            "Conditioning an artificial hypergeometric source null on original host richness does not control insect transport opportunity, first arrival, Hawaii insect nativity or field detection.",
            "The source mathematical right tail is NOT an ecological p value or an estimated invasion-risk coefficient.",
            "WCVP accepts wild native/naturalized plants, not every cultivated-only Hawaiian garden food resource.",
            "The prior source positive 10/10 does not show feeding on these host plants in Hawaiian wild populations.",
        ],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--Hawaii-positive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    positive=json.loads(args.__dict__["Hawaii_positive"].read_text())
    if positive.get("schema")!=\
        "chocho_zenodo_2024_10_aliens_Hawaii_original_WCVP_hosts_20261027_v01":
        raise ValueError("Original ten positive Hawaii source host matrix missing")
    r=audit(args.source_root,positive)
    r["source_hawaii_positive_summary_SHA256"]=sha256(args.__dict__["Hawaii_positive"])
    r["source_original_source_manifest_SHA256"]=sha256(args.source_root/"SOURCE_MANIFEST.json")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":r["status"],
        "all239_Hawaii_WCVP_origin_classes":r["original_239_source_host_coverage"],
        "all239_any_introduced_Hawaii_host":
            r["original_239_source_any_introduced_Hawaii_wild_host"],
        "all239_any_native_Hawaii_host":
            r["original_239_source_any_native_Hawaii_wild_host"],
        "source_bin_matched_artificial_reference":
            r["confirmed_10_positive_richness_stratified_reference"],
    },sort_keys=True))


if __name__=="__main__":
    main()
