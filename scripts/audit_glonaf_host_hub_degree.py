#!/usr/bin/env python3
"""Interpret 86 broadly catalogued focal larval-host source IDs, descriptively.

Do broadly GloNAF-indexed botanical hosts connect to more of the fixed239
HOSTS butterflies? Only source taxonomy/candidate resource network, no actual
local feeding, no bird/insect occupancy, and no p-values. Results are source
taxonomy audit and consumer-degree counts. No original third-party row data
are republished; a small top-10 public species name list may aid interpretation.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import defaultdict,Counter
import csv
import hashlib
import json
from pathlib import Path
import statistics

from analyze_documented_host_specificity import (
    GLO_NAF,TDWG,HOSTS,PANEL,BINS,binomial,glonaf_index,prevalence_bins
)
from analyze_matched_glonaf_consumer_history import produce_matched_regions
from analyze_colonial_legacy_resource_beta import load as load_colonial_butterfly

BASE=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
PREV=BASE/"documented_host_prevalence_leave_bin_out_v1.json"
OUT=BASE/"glonaf_documented_host_source_hubs_v1.json"


def consumer_links():
    with PANEL.open("r",encoding="utf8",newline="") as f:
        panel={r["species"].strip() for r in csv.DictReader(f)}
    if len(panel)!=239:
        raise RuntimeError("No fixed 239 GEB butterfly names")
    links=defaultdict(set)
    with HOSTS.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        for r in reader:
            species=(r["Insect Genus"].strip()+" "+r["Insect Species"].strip())
            if species not in panel:
                continue
            name=binomial(r["Hostplant Genus"],r["Hostplant Species"])
            if name:
                links[name].add(species)
    return links


def summarize(names,links,rawids,freq):
    degrees=[len(links.get(n,set())) for n in names]
    butterflies={sp for n in names for sp in links.get(n,set())}
    # Sort independently of butterfly resource outcome; names are from
    # independently curated HOSTS taxonomic assertions.
    top=sorted(((n,len(links.get(n,set()))) for n in names),
               key=lambda pair:(-pair[1],pair[0]))[:10]
    return {
        "raw_GloNAF_taxon_IDs":len(rawids),
        "unique_canonical_genus_species_aliases":len(names),
        "aliases_with_documented_focal239_consumers":sum(x>0 for x in degrees),
        "focal239_butterflies_with_at_least_one_documented_host_in_stratum":len(butterflies),
        "fraction_of_239_butterflies_with_a_documented_host_in_stratum":len(butterflies)/239,
        "mean_documented_focal_consumer_degree_per_alias":statistics.mean(degrees) if degrees else None,
        "median_documented_focal_consumer_degree_per_alias":statistics.median(degrees) if degrees else None,
        "max_documented_focal_consumer_degree_per_alias":max(degrees) if degrees else None,
        "min_source_region_prevalence":min((freq[j] for j in rawids),default=None),
        "max_source_region_prevalence":max((freq[j] for j in rawids),default=None),
        "top_10_binomial_aliases_by_documented_consumer_degree":top,
    }


def evaluate():
    frozen=json.loads(PREV.read_text(encoding="utf8"))
    if frozen["schema"]!="chocho_posthoc_documented_host_prevalence_influence_v1":
        raise RuntimeError("Prevalence influence source baseline not fixed")
    data=load_colonial_butterfly()
    if len(data["ids"])!=284:
        raise RuntimeError("Not fixed 284 region cohort")
    geographies,_,original=produce_matched_regions(GLO_NAF,TDWG,(.05,))
    full=[geographies[.05][code] for code in data["ids"]]
    aliases,proof=glonaf_index(GLO_NAF)
    if proof["distinct_glonaf_names"]!=original["glonaf"]["distinct_species"]:
        raise RuntimeError("GloNAF source IDs not preserved")
    consumer=consumer_links()
    matched=sorted(set(consumer)&set(aliases))
    active_ids={j for name in matched for j in aliases[name]}
    freq,_=prevalence_bins(full,proof["distinct_glonaf_names"])
    matched_active={j for j in active_ids if freq[j]>0}
    if len(matched_active)!=frozen["source_active_documented_host_IDs"]:
        raise RuntimeError("Matched known-host active taxa drift")
    raw_alias=defaultdict(set)
    for name in matched:
        for j in aliases[name]:
            if j in matched_active:
                raw_alias[j].add(name)
    if any(len(aa)!=1 for aa in raw_alias.values()):
        raise RuntimeError("GloNAF source species ID is ambiguously linked to host aliases")
    strata=defaultdict(list)
    for j in matched_active:
        strata[bisect_left(BINS,int(freq[j]))].append(j)
    result={}
    for band in sorted(strata):
        rawids=sorted(strata[band])
        names={next(iter(raw_alias[j])) for j in rawids}
        result[str(band)] = {
            "source_region_prevalence_range_inclusive":[
                1 if band==0 else BINS[band-1]+1,BINS[band]],
            **summarize(names,consumer,rawids,freq),
        }
    pivotal=result.get("7")
    if pivotal is None or pivotal["raw_GloNAF_taxon_IDs"]!=86:
        raise RuntimeError("Critical post-hoc bin7 source 86 taxon IDs drift")
    return {
       "schema":"chocho_fixed239_documented_HOSTS_host_hub_strata_v1",
       "status":"PUBLIC_SOURCE_BOTANICAL_AND_HOST_NETWORK_INTERPRETATION_NOT_OCCUPANCY",
       "source":{
           "original_GloNAF_Zenodo_md5":"16f5585d2500676951a4b664bbf0ffde",
           "HOSTS_fixed_source_commit":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
           "frozen_prevalence_influence_sha256":hashlib.sha256(PREV.read_bytes()).hexdigest(),
           "fixed_butterflies":239,"GloNAF_matched_regions":284,
       },
       "known_focal_host_strata":result,
       "high_prevalence_101_200_region_stratum":pivotal,
       "interpretation_limits":[
           "Documented degree is number of focal239 butterfly species named in HOSTS; sampling and taxon biases remain.",
           "Source genus-species aliasing does not resolve all botanical synonyms, so GloNAF source plant IDs may duplicate species names.",
           "Geographic prevalence in a modern GloNAF alien plant catalogue is not introduction age or historical maritime import.",
           "The butterfly potential-resource outcome is constructed from HOSTS and WCVP, not independent observed butterfly occurrence.",
           "A highly geographic host does not imply it is used by all documented butterflies in every region.",
           "The 86 source plant IDs were examined post-hoc as the strongest influence band; no causal or original confirmatory mechanism follows.",
       ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=OUT)
    args=ap.parse_args()
    obj=evaluate()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",
                           encoding="utf8")
    print(json.dumps({k:v for k,v in obj.items() if k!="known_focal_host_strata"},
                     indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
