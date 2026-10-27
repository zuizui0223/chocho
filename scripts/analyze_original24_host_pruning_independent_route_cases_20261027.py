#!/usr/bin/env python3
"""Explore recipient-local host repertoire pruning in 24 already documented aliens.

Input source: independently audited original chocho #85 Zenodo 2024 alien
establishment positives × accepted HOSTS/WCVP *wild* plant source region.
The 24 source cases are NOT independent introduction attempts or direct
observed local larvae, and all 35 Zenodo records are positive-only.

Contrasting pathway evidence: Galápagos Calpodes ethlius (independently
classified plant-associated contaminant, not verified carrier plant/lot)
versus Okinawa Leptosia nina (a coarse Nansei botanical proxy and local
wild larva on Gyoboku, not known introduced-plant arrival pathway).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

SCHEMA="chocho_24_receiver_host_pruning_independent_pathway_cases_20261027_v01"
SOURCE_PROFILE_SHA="5ef806a12f8750572c5303bce40ea5ffa58ad56138ae8cdb75806ff7cbbb8efc"
SOURCE_PROFILE_SCHEMA="chocho_original239_zenodo35_recipient_host_opportunity_20261027_v01"
EXPECTED_COUNTS={
    "MAIN_COMPARABLE":24,
    "SENSITIVITY_BOTANICAL_PROXY":3,
    "UNSCORABLE_ORIGINAL_SOURCE_GEOGRAPHY":8,
}
EXTERNAL_COLS=(
    "consumer_species","recipient_zone","source_geography_grade",
    "first_source_butterfly_year","consumer_source_establishment_status",
    "source_arrival_intent","arrival_pathway_original_literal",
    "arrival_pathway_evidence_grade","arrival_pathway_source_url",
    "independent_host_plant_or_genus","independent_plant_first_record_year",
    "independent_plant_region_origin_status","independent_plant_intent",
    "independent_plant_introduction_pathway","independent_plant_source_url",
    "regional_larval_host_taxon_grade","larval_host_source_evidence_grade",
    "independent_local_host_evidence_url",
    "same_carrier_and_actual_larval_host_species_verified",
    "scientific_uncertainty",
)
FIXED_CASES={
    ("Calpodes ethlius","Galapagos"):"EXACT_WCVP_LEVEL3_SOURCE_REGION",
    ("Leptosia nina","Okinawa"):"BOTANICAL_NANSEI_PROXY_NOT_EXACT_OKINAWA",
}


def sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_profile(path, *, enforce_source=True):
    if enforce_source and sha_file(path)!=SOURCE_PROFILE_SHA:
        raise ValueError("Previously source-verified original 35 alien recipient profile SHA changed")
    obj=json.loads(Path(path).read_text())
    if obj.get("schema")!=SOURCE_PROFILE_SCHEMA:
        raise ValueError("Not original 239 accepted WCVP + Zenodo 2024 35 positive input")
    original=obj.get("source",{})
    if (original.get("original_2024_positive_species_recipient_records")!=35
        or original.get("original_2024_matched_aliens")!=18
        or original.get("original_239_butterflies")!=239
        or original.get("original_accepted_butterfly_host_pairs")!=2596):
        raise ValueError("Original 239/18/35/2596 source panel changed")
    src=obj["geographic_source_guard95_and_sensitivity90"]["0.95"]
    if src["source_positive35_mapped_geography_counts"]!=EXPECTED_COUNTS or \
       len(src["positive_original35_each_with_scope_and_host_counts"])!=35:
        raise ValueError("Original verified geographical source cohorts no longer match")
    return obj


def count_local_host_source(rows):
    """Original accepted wild HOSTS global link count by source region, no field use."""
    strict=[]
    unresolved=[]
    for row in rows:
        grade=row.get("geographic_status")
        if grade!="MAIN_COMPARABLE":
            unresolved.append({
                "species":row["species"],
                "original_2024_recipient_zone":row.get(
                    "independent_2024_recipient_zone",row.get("recipient_zone")),
                "source_geographic_evidence_grade":grade,
                "source_geographic_mapping_status":row.get("source_mapping_status"),
                "not_scored_as_absent":True,
            })
            continue
        sp=row["species"]
        region=row["independent_2024_recipient_zone"]
        intro=int(row["introduced_only_wild_source_hosts"])+int(row[
            "source_plant_native_and_introduced_in_country_atoms"])
        native=int(row["native_only_wild_source_hosts"])+int(row[
            "source_plant_native_and_introduced_in_country_atoms"])
        both=int(row["source_plant_native_and_introduced_in_country_atoms"])
        unlisted=int(row["unlisted_or_below_area_gate_known_hosts"])
        global_hosts=int(row["accepted_original_known_global_larval_hosts"])
        if not global_hosts>0 or not intro>0 or intro+native-both+unlisted!=global_hosts:
            raise ValueError("Previously source-verified botanical host decomposition changed")
        strict.append({
            "species":sp,
            "recipient_zone":region,
            "source_original_global_accepted_host_species":global_hosts,
            "source_current_recipient_introduced_wild_host_species":intro,
            "source_current_recipient_native_wild_host_species":native,
            "source_current_recipient_both_native_and_introduced_host_species":both,
            "source_current_recipient_no_WCVP_wild_record_or_area_below_gate":unlisted,
            "introduced_wild_host_fraction_of_original_globally_known_hosts":
                intro/global_hosts,
            "any_WCVP_native_wild_host_alternative":native>0,
            "exploratory_local_potential_redundancy_class":(
                "ONE_INTRODUCED_WILD_HOST" if intro==1 else
                "TWO_OR_THREE_INTRODUCED_WILD_HOSTS" if intro<=3 else
                "FOUR_OR_FIVE_INTRODUCED_WILD_HOSTS" if intro<=5 else
                "SIX_OR_MORE_INTRODUCED_WILD_HOSTS"
            ),
            "source_host_range_not_identical_to_historical_plant_distribution":True,
        })
    if len(strict)!=24 or len(unresolved)!=11:
        raise ValueError("Original strict 24 and separately coarse 11 mapping changed")
    positive_species=set(x["species"] for x in strict)
    bysp=defaultdict(list)
    for r in strict:
        bysp[r["species"]].append(r)
    groups=[]
    for species,grp in sorted(bysp.items()):
        if len(grp)<2:
            continue
        values=[x["source_current_recipient_introduced_wild_host_species"]
                for x in grp]
        groups.append({
            "species":species,
            "number_of_distinct_source_receiving_areas":len(grp),
            "original_accepted_global_larval_host_species_count":
                grp[0]["source_original_global_accepted_host_species"],
            "source_min_current_introduced_host_species":min(values),
            "source_max_current_introduced_host_species":max(values),
            "source_max_divided_by_min_current_introduced_hosts":
                max(values)/min(values),
            "source_receiving_areas_intro_host_species":[
                {"zone":x["recipient_zone"],
                 "accepted_introduced_wild_host_species":
                    x["source_current_recipient_introduced_wild_host_species"],
                 "any_native_wild_host":x["any_WCVP_native_wild_host_alternative"]}
                for x in sorted(grp,key=lambda q:q["recipient_zone"])
            ],
            "constant_original_species_HOSTS_list_across_receiving_areas":True,
        })
    ks=[r["source_current_recipient_introduced_wild_host_species"] for r in strict]
    fracs=[r["introduced_wild_host_fraction_of_original_globally_known_hosts"]
           for r in strict]
    total_global=sum(r["source_original_global_accepted_host_species"] for r in strict)
    status=Counter(r["exploratory_local_potential_redundancy_class"] for r in strict)
    counts={
        "n_original_strict_geographic_source_positive_cases":len(strict),
        "n_distinct_original_butterfly_taxa_in_24_cases":len(positive_species),
        "n_distinct_source_receiving_areas":len(
            set(x["recipient_zone"] for x in strict)),
        "n_with_only_1_introduced_wild_host":sum(k==1 for k in ks),
        "n_with_at_most_2_introduced_wild_hosts":sum(k<=2 for k in ks),
        "n_with_at_most_3_introduced_wild_hosts":sum(k<=3 for k in ks),
        "n_with_at_most_5_introduced_wild_hosts":sum(k<=5 for k in ks),
        "n_with_at_most_10_introduced_wild_hosts":sum(k<=10 for k in ks),
        "n_with_any_native_wild_host_alternative":sum(
            r["any_WCVP_native_wild_host_alternative"] for r in strict),
        "source_local_intro_species_median":statistics.median(ks),
        "source_local_intro_species_min":min(ks),
        "source_local_intro_species_max":max(ks),
        "source_local_intro_over_original_global_host_richness_median":
            statistics.median(fracs),
        "source_intro_plant_species_memberships_summed_over_24_positive_cases":
            sum(ks),
        "source_global_known_host_plant_species_memberships_over_24_positive_cases":
            total_global,
        "membership_weighted_source_local_intro_fraction":sum(ks)/total_global,
        "source_pruning_classes":dict(status),
        "n_original_taxa_with_two_or_more_strict_alien_recipient_geographies":
            len(groups),
        "source_yearly_first_establishment_not_cargo_interception":True,
        "no_actual_butterfly_host_consumption_inferred":True,
    }
    return counts,strict,unresolved,groups


def load_cases(path, rows):
    with Path(path).open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or [])!=EXTERNAL_COLS:
            raise ValueError("Independent published-pathway evidence source schema changed")
        cases=list(reader)
    if len(cases)!=2 or {
        (r["consumer_species"],r["recipient_zone"]):r["source_geography_grade"]
        for r in cases
    }!=FIXED_CASES:
        raise ValueError("Two independent original source case identities changed")
    bykey={(r["species"],r.get("recipient_zone",r.get("independent_2024_recipient_zone"))):r
           for r in rows}
    evidence=[]
    for row in cases:
        sp,region=row["consumer_species"],row["recipient_zone"]
        key=sp,region
        if key not in bykey:
            raise ValueError("Independent original arrival case absent from original positive 35")
        src=bykey[key]
        if row["source_geography_grade"]=="EXACT_WCVP_LEVEL3_SOURCE_REGION":
            if src.get("geographic_status")!="MAIN_COMPARABLE":
                raise ValueError("True source independent Galapagos beneficiary not mapped")
            if (int(src["introduced_only_wild_source_hosts"]),
                int(src["native_only_wild_source_hosts"]),
                int(src["accepted_original_known_global_larval_hosts"]))!=(3,0,5):
                raise ValueError("Calpodes accepted Galapagos original host source changed")
        else:
            if src.get("geographic_status")!="SENSITIVITY_BOTANICAL_PROXY":
                raise ValueError("Okinawa botanical proxy cannot enter primary 24 cases")
            if (int(src["introduced_only_wild_source_hosts"]),
                int(src["native_only_wild_source_hosts"]),
                int(src["accepted_original_known_global_larval_hosts"]))!=(0,2,10):
                raise ValueError("Nansei Leptosia 0 introduced/2 native accepted hosts source changed")
        for field in ("arrival_pathway_source_url","independent_plant_source_url",
                      "independent_local_host_evidence_url"):
            u=urlparse(row[field])
            if u.scheme!="https" or not u.netloc:
                raise ValueError("Original independent host/pathway source URL missing")
        if row["same_carrier_and_actual_larval_host_species_verified"]!="UNKNOWN":
            raise ValueError("No original source evidence links a verified carrier plant species to local larval development")
        b_year=row["first_source_butterfly_year"]
        plant_year=row["independent_plant_first_record_year"]
        if b_year and not (1700<=int(b_year)<=2026):
            raise ValueError("Invalid original first insect source year")
        if plant_year and not (1700<=int(plant_year)<=2026):
            raise ValueError("Invalid independent plant first record source year")
        case={
            "consumer_species":sp,
            "independently_established_recipient":region,
            "source_evidence_grade":row["source_geography_grade"],
            "source_first_record_of_consumer_year":int(b_year) if b_year else None,
            "independent_reported_arrival_pathway_class":row[
                "arrival_pathway_original_literal"],
            "arrival_pathway_evidence_grade":row["arrival_pathway_evidence_grade"],
            "independent_host_plant_taxon":row["independent_host_plant_or_genus"],
            "source_first_record_of_named_plant_year":int(plant_year)
                if plant_year else None,
            "first_record_temporal_order_does_not_prove_absence":{
                "butterfly_report_precedes_this_botanical_first_report":
                    bool(plant_year and b_year and int(b_year)<int(plant_year)),
                "first_report_is_not_actual_introduction_date":True,
            },
            "source_larval_evidence_grade":row["larval_host_source_evidence_grade"],
            "same_carrier_host_species_documented_to_have_transported_this_consumer":False,
            "carrier_host_identification_status":"UNVERIFIED_NOT_ZERO",
            "sources":{
                "independent_arrival_category":row["arrival_pathway_source_url"],
                "plant_status_and_first_record":row["independent_plant_source_url"],
                "larval_use_or_regional_taxon":row["independent_local_host_evidence_url"],
            },
            "interpretive_constraints":row["scientific_uncertainty"],
        }
        evidence.append(case)
    return evidence


def run(args):
    doc=load_profile(args.original_profile)
    src=doc["geographic_source_guard95_and_sensitivity90"]["0.95"]
    raw=src["positive_original35_each_with_scope_and_host_counts"]
    counts,strict,unscorable,repeated=count_local_host_source(raw)
    evidence=load_cases(args.independent_cases,raw)
    if (counts["n_with_at_most_3_introduced_wild_hosts"]!=8
        or counts["source_local_intro_species_median"]!=8
        or counts["n_with_any_native_wild_host_alternative"]!=9):
        raise ValueError("Original 24 host-support case margins changed from independently validated source")
    if len(repeated)!=5 or len(strict)!=24:
        raise ValueError("Original within-species 5 repeat positive taxon evidence changed")
    report={
        "schema":SCHEMA,
        "status":"EXPLORATORY_24_LOCAL_WILD_HOST_PRUNING_PLUS_INDEPENDENT_NONCAUSAL_PATHWAY_CASES",
        "source":{
            "original_input_verified_35_recipient_profile_SHA256":
                sha_file(args.original_profile),
            "original_source_18_species_35_positive_receiving_areas":True,
            "original_global_butterfly_source_cohort":239,
            "original_accepted_butterfly_host_plant_pairs":2596,
            "original_strict_geographic_botanical_zones":7,
            "original_2024_species_alien_establishment_positive_only":True,
            "source_independent_route_case_ledger_SHA256":
                sha_file(args.independent_cases),
        },
        "aggregate_source_original24":counts,
        "source_recorded_24_accepted_host_geography_profiles":strict,
        "source_original35_other_11_geographically_proxy_or_unscorable":unscorable,
        "within_butterfly_across_multiple_receiving_areas":repeated,
        "independent_two_invasion_route_context_cases":evidence,
        "interpretation_limits":[
            "All accepted globally documented HOSTS plant links are potential geographically available resource taxa only, not wild stage-specific feeding locally.",
            "Source original 24 species×recipient cases include repeated butterflies and are all alien establishment POSITIVE; no failed-arrival cohort or interpretable causal p values.",
            "Source global-HOSTS breadth is literature-dependent and the locally reduced WCVP wild list omits planted-only ornamentals, ephemeral cultivation and unknown taxonomic synonyms.",
            "An original WCVP region-wide introduced host is not confirmed present in the actual exact site of insect reproduction or at first insect arrival date.",
            "Source single wild-host matches do not prove biological obligate dependence or extinction sensitivity, and plant taxa are not weighted by standing leaf biomass.",
            "Galapagos introduction-pathway CATEGORY contaminant on plants is not observed named live nursery shipment nor actual insect stage.",
            "A 1991 first botanical collection label is not proof plant Canna indica absent in 1987 when Calpodes ethlius was first recorded.",
            "Okinawa-to-Nansei is an incomplete botanical proxy, and the photographed local Crateva host is not independently voucher-confirmed as modern accepted species.",
            "Okinawa arrival via humans versus natural adult dispersal has no source-verified transport record.",
            "The categorical contrast of two selected island cases cannot estimate arrival-route frequencies or a general orphan-resource percentage.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":report["status"],
        "source_24_pruning_summary":counts,
        "source_same_species_multiple_receiving_regions":repeated,
        "two_qualitative_pathway_cases":evidence,
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-profile",type=Path,required=True)
    p.add_argument("--independent-cases",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
