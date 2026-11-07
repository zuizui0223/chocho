#!/usr/bin/env python3
"""Offline evidence integrity audit: original verified WPI OD and separate 27 UNVERIFIED what-if.

The full 645581-row independent raw ZIP reparse ran in scientific CI. This
additional gate checks permanently published source totals and all 1275 rows
without relying on third-party ZIP availability during review or CI reruns.
"""
from __future__ import annotations
import csv
import json
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/"island51_zenodo2020_portvoyage_audit_20261107"

def read(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def verify():
    old=read(PARENT/"original51_all1275_direct2020_source_pair_observation.csv")
    novel=read(HERE/"original51_verified_vs_UNVERIFIED_geonear_alias_2020_OD_pair_scenarios.csv")
    parent=json.loads((PARENT/"source_matched_voyage_receipt.json").read_text())
    source=json.loads((HERE/"source_only_geonear27_whatif_receipt.json").read_text())
    audit=json.loads((HERE/"independent_raw_source_audit.json").read_text())
    if (len(old)!=1275 or len(novel)!=1275 or
        parent["source_original_194_WPI_ports"]!=194 or
        parent["source_original_51_landmasses"]!=51 or
        parent["verified_external_endpoint_nodes"]!=46 or
        source["complete_original_2020_source_triplet_rows"]!=645581 or
        source["ADDITIONAL_UNVERIFIED_geography_only_endpoint_nodes"]!=27):
        raise ValueError("Prior verified or hypothetical source input cohort changed")
    if (source["status"]!="SOURCE_ONLY_HYPOTHETICAL_27_NEARBY_NODE_COVERAGE_CALCULATED_NOT_VERIFIED_TRAFFIC" or
        audit["status"]!="ORIGINAL_645581_SOURCE_ROWS_AND_1275_VERIFIED_VS_UNVERIFIED_PAIRS_RECONSTRUCTED"):
        raise ValueError("Published source-only or raw-source independent audit failed")
    seen=set()
    original_sum=0
    hypothetical_sum=0
    new_pairs=0
    source_pairs=0
    first_pair_extra=0
    for left,right in zip(old,novel):
        ka=(left["island_a"],left["island_b"])
        kb=(right["island_a"],right["island_b"])
        if ka!=kb or ka in seen:
            raise ValueError("Original source island pair order/identity altered or duplicated")
        seen.add(ka)
        v=int(left["direct_observed_2020_source_voyages_in_either_direction"])
        w=int(right["previous_VERIFIED_source_events"])
        h=int(right["hypothetical_27_geonear_alias_source_events_NOT_VERIFIED"])
        delta=int(right["additional_source_events_due_to_unverified_port_identity"])
        if v!=w or h<w or h-w!=delta or h<0:
            raise ValueError("Original source verified triplet counts were altered")
        if (right["previous_verified_direct_pair"].lower()=="true")!=(v>0) or (
            right["pair_exists_only_if_unverified_port_candidates_are_ASSUMED_same"].lower()=="true")!=(v==0 and h>0):
            raise ValueError("Source false-identity indicators inconsistent")
        if (right["zero_source_index_is_not_no_ship_traffic"].lower()!="true" or
            right["historical_host_or_butterfly_colonization_not_inferred"].lower()!="true"):
            raise ValueError("Unobserved source route erroneously treated as an ecological absence/mechanism")
        source_pairs+=int(v>0)
        new_pairs+=int(v==0 and h>0)
        original_sum+=v
        hypothetical_sum+=h
        if ka==("15","26"):
            first_pair_extra=delta
    hpair=sum(int(x["hypothetical_27_geonear_alias_source_events_NOT_VERIFIED"])>0 for x in novel)
    if (len(seen)!=1275 or source_pairs!=12 or
        source_pairs!=parent["distinct_original_landmass_unordered_pairs_with_observed_2020_voyage"] or
        original_sum!=133 or hpair!=28 or new_pairs!=16 or
        hypothetical_sum!=298 or hypothetical_sum-original_sum!=165 or first_pair_extra!=99):
        raise ValueError("Source scientific 12/28, 133/298 or additional 16/165 values do not reproduce")
    derived={
      "previous_verified_unordered_island_pairs":source_pairs,
      "hypothetical_unordered_island_pairs":hpair,
      "additional_island_pair_types_requiring_unverified_node_assumptions":new_pairs,
      "previous_verified_cross_island_source_events":original_sum,
      "hypothetical_cross_island_source_events":hypothetical_sum,
    }
    for key,value in derived.items():
        if source[key]!=value:
            raise ValueError("Source published receipt not equal to 1275 full-row table: "+key)
    if (source["ACTUAL_SOURCE_GATE_REMAINS_FAILED_WITHOUT_VERIFIED_IDENTITIES"] is not True or
        source["hypothetical_15islands_80pairs_threshold_reached_if_all_aliases_true"] is not False or
        source["unverified_aliases_are_not_promoted_to_original_WPI_matches"] is not True or
        audit["source_hypothetical_unverified_island_links"]!=28 or
        audit["publisher_source_rows_verified"]!=645581):
        raise ValueError("Source-only no-biological-fit/false-port-identity boundary was lost")
    result={
        "status":"ALL_1275_SOURCE_PAIRS_AND_SOURCE_RAW_AUDIT_RECEIPTS_RECONSTRUCTED_OFFLINE",
        "source_original_ports":194,
        "original_island_pair_rows":1275,
        "previous_verified_island_pairs":12,
        "hypothetical_UNVERIFIED_island_pairs":28,
        "source_extra_island_pair_types_only_if_false_names_assumed":16,
        "source_original_verified_cross_island_events":133,
        "source_hypothetical_cross_island_events":298,
        "prespecified_80_source_verified_island_pairs_UNMET":True,
        "all_unverified_geonear_candidates_still_UNVERIFIED":True,
    }
    print("ZENODO27_PUBLISHED_OFFLINE_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    verify()
