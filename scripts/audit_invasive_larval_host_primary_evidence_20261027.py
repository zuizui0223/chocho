#!/usr/bin/env python3
"""Source-verified larval-stage evidence screen for 13 invasive butterfly events.

The evidence file is curated from primary studies, official study abstracts,
WCVP/POWO and explicit scientific records. It is NOT a second independent
statistical sample of butterflies or a rate of evolutionary host adoption.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

SCHEMA = "chocho_13_invasion_event_larval_host_independent_evidence_20261027_v1"
FROZEN_ORIGINAL_EVENTS = frozenset({
    ("Pieris rapae", "CA"), ("Pieris rapae", "NZ"), ("Pieris rapae", "AU"),
    ("Thymelicus lineola", "CA"), ("Pieris brassicae", "CL"),
    ("Danaus plexippus", "US-HI"), ("Danaus plexippus", "AU"),
    ("Danaus plexippus", "NZ"), ("Papilio demoleus", "DO"),
    ("Papilio demoleus", "PR"), ("Papilio xuthus", "US-HI"),
    ("Cacyreus marshalli", "ES"), ("Chilades pandava", "GU"),
})
EVENT_COLS = ("species", "in_original_239", "native_region_statement",
              "recipient_region", "recipient_code", "first_year_lower",
              "first_year_upper", "establishment_by_year", "date_interpretation",
              "arrival_mechanism_status", "source_confidence",
              "candidate_host_context", "primary_evidence_url", "interpretive_limits")
EVIDENCE_COLS = (
    "record_id", "butterfly", "recipient_code", "evidence_site",
    "observation_year", "host_name", "host_taxonomic_rank",
    "recipient_host_origin", "host_origin_evidence_url", "evidence_tier",
    "developmental_endpoint", "field_derived_local_host_record",
    "geographic_scope", "study_citation", "evidence_source_url", "evidence_note",
)
DIRECT = {"FIELD_LARVA_HOST_REARED", "FIELD_EGG_AND_LARVA_HOST",
          "FIELD_LARVAL_DAMAGE_OR_OVIPOSITION"}
NO_FIELD = {"CONTROLLED_ONLY_LARVAL_COMPLETION",
            "CONTROLLED_ONLY_OVIPOSITION",
            "REGIONAL_HOST_LIFE_CYCLE_REVIEW",
            "ADULT_OR_ESTABLISHMENT_ONLY"}
ORIGINS = {"NATIVE_AT_RECIPIENT", "ALIEN_NATURALIZED_RECORDED",
           "ALIEN_CULTIVATED_NOT_ESTABLISHED", "PLANT_STATUS_UNRESOLVED"}
STRICT_LOCAL = {"RECIPIENT_SITE", "SAME_RECEIVING_COUNTRY_DIFFERENT_SITE",
                "SAME_RECEIVING_TERRITORY_DIFFERENT_ISLAND"}
NONSTRICT = {"NATIVE_RANGE_COMPARATOR",
             "RECIPIENT_SITE_OUTSIDE_ORIGINAL_13_EVENTS",
             "EXPERIMENT_OUTSIDE_ORIGINAL_13_EVENTS",
             "MULTICOUNTRY_REVIEW_NOT_MALLORCA_EVENT",
             "RECEIVING_COUNTRY_REVIEW_NO_SITE_PAIR",
             "RECIPIENT_SITE_NO_HOST_LARVA",
             "RECEIVING_COUNTRY_NO_HOST_LARVA"}


def read_file(path, cols):
    with Path(path).open(encoding="utf-8", newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or [])!=cols:
            raise ValueError("Source evidence schema or first-record CSV columns changed")
        return list(reader)


def ensure_url(url, *, missing_allowed=False):
    if missing_allowed and not url:
        return
    parsed=urlparse(url)
    if parsed.scheme!="https" or not parsed.netloc or not parsed.path:
        raise ValueError("Unverified study URL or plant source URL")


def score(event_rows, evidence_rows):
    keys={(r["species"], r["recipient_code"]):r for r in event_rows}
    if set(keys) != FROZEN_ORIGINAL_EVENTS or len(event_rows)!=13 or \
       len({r["species"] for r in event_rows})!=8:
        raise ValueError("Original frozen 13 events / eight species altered")
    core={r["species"] for r in event_rows if r["in_original_239"]=="true"}
    if core!={"Pieris rapae","Thymelicus lineola","Pieris brassicae",
             "Danaus plexippus","Papilio demoleus"}:
        raise ValueError("Five original-239 invasive candidate identities changed")
    if len(evidence_rows)!=15 or len({r["record_id"] for r in evidence_rows})!=15:
        raise ValueError("Curated publication-stage evidence roster changed")
    grouped=defaultdict(list)
    alltiers=Counter()
    known_species={r["species"] for r in event_rows}
    source_ids=set()
    for r in evidence_rows:
        taxon, loc=r["butterfly"],r["recipient_code"]
        if taxon not in known_species:
            raise ValueError("New butterfly species selected by favorable result")
        if r["host_taxonomic_rank"] not in ("species","genus"):
            raise ValueError("Unsupported local host taxonomic resolution")
        if r["evidence_tier"] not in DIRECT|NO_FIELD or \
           r["recipient_host_origin"] not in ORIGINS:
            raise ValueError("Unknown evidence tier or plant-native origin")
        if r["geographic_scope"] not in STRICT_LOCAL|NONSTRICT:
            raise ValueError("Unexpected geographical scope")
        if (r["field_derived_local_host_record"]=="true") != \
                (r["evidence_tier"] in DIRECT):
            raise ValueError("Field stage requirement disagrees with source evidence label")
        if r["evidence_tier"] in DIRECT and (
                not r["observation_year"] or
                r["developmental_endpoint"].startswith("no_")):
            raise ValueError("Direct local field stage missing source evidence year or stage")
        if loc=="GU" and taxon=="Chilades pandava" and \
                r["host_name"]=="Cycas micronesica" and \
                r["recipient_host_origin"]!="NATIVE_AT_RECIPIENT":
            raise ValueError("Guam native cycad incorrectly treated as introduced")
        if taxon=="Cacyreus marshalli" and r["host_name"].startswith("Geranium ") \
                and r["evidence_tier"] in DIRECT:
            raise ValueError("Experimental native European Geranium is not confirmed wild larval use")
        if r["host_taxonomic_rank"]=="genus" and r["recipient_host_origin"] in (
                "NATIVE_AT_RECIPIENT","ALIEN_NATURALIZED_RECORDED"):
            raise ValueError("Genus-only host source cannot prove exact WCVP plant origin")
        ensure_url(r["evidence_source_url"])
        ensure_url(r["host_origin_evidence_url"],missing_allowed=True)
        if r["observation_year"]:
            year=int(r["observation_year"])
            if not 1800<=year<=2026:
                raise ValueError("Implausible source observation year")
            if (taxon,loc) in keys and year<int(keys[taxon,loc]["first_year_lower"]):
                raise ValueError("Source larval study predates original reported butterfly arrival")
        if r["geographic_scope"] in STRICT_LOCAL and (taxon,loc) not in keys:
            raise ValueError("Strict local evidence cannot silently add a new invasion event")
        alltiers[r["evidence_tier"]]+=1
        grouped[(taxon,loc)].append(r)
        source_ids.add(r["study_citation"])
    entries=[]
    core_direct=set()
    any_direct=set()
    field_host_routes=Counter()
    for k,event in sorted(keys.items()):
        sources=grouped.get(k,[])
        strict=[r for r in sources if r["geographic_scope"] in STRICT_LOCAL
                and r["evidence_tier"] in DIRECT]
        direct_origins=set()
        site_tiers=sorted({r["evidence_tier"] for r in strict})
        for x in strict:
            direct_origins.add(x["recipient_host_origin"])
            field_host_routes[x["recipient_host_origin"]]+=1
        is_direct=bool(strict)
        if is_direct:
            any_direct.add(k[0])
            if event["in_original_239"]=="true":
                core_direct.add(k[0])
        entries.append({
            "butterfly":k[0],
            "recipient_region_code":k[1],
            "original_239_study_species":event["in_original_239"]=="true",
            "literature_evidence_rows_total":len(sources),
            "event_compatible_field_host_evidence":is_direct,
            "original_field_evidence_tiers":site_tiers,
            "strict_field_evidence_host_origin_classes":sorted(direct_origins),
            "event_field_host_evidence_record_count":len(strict),
            "independently_sourced_larval_field_studies_may_overlap_HOSTS":True,
            "no_field_host_evidence_is_NOT_ecological_absence":not is_direct,
        })
    context=[r for r in evidence_rows if
             (r["butterfly"],r["recipient_code"]) not in keys]
    controlled_native=[
        r for r in evidence_rows
        if r["evidence_tier"]=="CONTROLLED_ONLY_LARVAL_COMPLETION"
        and r["recipient_host_origin"]=="NATIVE_AT_RECIPIENT"
    ]
    count=Counter(x["event_compatible_field_host_evidence"] for x in entries)
    return {
        "schema":SCHEMA,
        "status":"STAGE_RESOLVED_PRIMARY_SOURCE_CASE_SCREEN_NOT_GLOBAL_ADOPTION_ESTIMATE",
        "original_event_species":8,
        "original_event_count":13,
        "original_event_239_butterfly_species":5,
        "curated_independent_publication_evidence_rows":len(evidence_rows),
        "distinct_study_citations":len(source_ids),
        "all_source_evidence_tiers":dict(sorted(alltiers.items())),
        "original_events_with_strict_field_larval_host_evidence":count[True],
        "original_events_without_source_qualified_field_larval_host_evidence":count[False],
        "original_239_butterfly_species_with_event_qualified_field_evidence":len(core_direct),
        "all_species_with_event_qualified_field_evidence":len(any_direct),
        "event_qualified_field_source_host_origin_by_row":dict(field_host_routes),
        "native_European_Geranium_lab_development_studies_without_wild_host_confirmation":
            len(controlled_native),
        "separate_native_and_new_site_source_rows_outside_original_13_events":
            len(context),
        "original_events":entries,
        "case_evidence_source_independence":"PRIMARY_RESEARCH_CITATION_RETRIEVED_SEPARATELY_BUT_ORIGINAL_HOSTS_CITATION_OVERLAP_NOT_IDENTIFIABLE",
        "critical_stage_interpretation":{
            "field_larva_reared":"Field-collected larvae supported on source-listed host under captive rearing; not evidence free-living adults emerged.",
            "field_oviposition":"Observed eggs/larvae in receiving habitat; no demonstrated development to adult without separate rearing.",
            "controlled_lab_only":"Potential host use demonstrated by captive experiments; never promoted to observed wild consumer-resource niche.",
            "genus_host":"A genus-only host report cannot be linked to one introduced/nonnative WCVP plant species.",
            "unknown_literature":"Failure to find a site-specific host record is missing evidence, not an absence of host use.",
            "source_dependence":"A separately retrieved study may also underpin a historical HOSTS catalogue entry; the HOSTS CSV lacks per-row source citations, so formal statistical independence cannot be asserted."
        },
        "chronology_discrepancies":[
            {"species":"Papilio demoleus","receiving":"DO",
             "original_2004_peer_reviewed_specimen":"2004-03-29",
             "conflicting_EPP0_global_distribution_first_year":2001,
             "original_specimen_source":"https://bvearmb.do/handle/123456789/2202",
             "conflicting_EPPO_source":"https://gd.eppo.int/taxon/PAPIDE/distribution/DO",
             "resolution":"FIRST_DOCUMENTED_PEER_REVIEWED_CAPTURE_2004_BUT_EARLIER_RETROSPECTIVE_EPPO_DATE_UNRESOLVED"},
            {"species":"Pieris rapae","receiving":"NZ",
             "first_early_source_reported":1929,
             "establishment_source_reported":1930,
             "first_record_source":"https://www.nzbutterflies.org.nz/species-info/cabbage-white/",
             "1930_encyclopedia_source":"https://teara.govt.nz/en/1966/butterfly-white",
             "resolution":"FIRST_OBSERVATION_AND_ESTABLISHMENT_DATES_NOT_INTERCHANGEABLE"}
        ],
        "limitations":[
            "This is an intentionally selected invasion-case evidence panel, not a representative 239-species survey.",
            "No recipient field study is assumed independent of original compiled HOSTS citation content.",
            "Native Guam host use does not by itself prove a genetically evolved host switch; ecological fitting is an hypothesis.",
            "Strong laboratory Geranium performance does not mean European native Geranium hosts are colonized in the field.",
            "A field-collected larva reared to pupation is not demonstrated free-ranging full life-cycle completion.",
            "Missing local records are never converted into zero biological host use.",
            "No true climate-niche shifts, direct competitive suppression or larval population abundance changes were measured."
        ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events",required=True,type=Path)
    ap.add_argument("--evidence",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    events=read_file(args.events,EVENT_COLS)
    evidence=read_file(args.evidence,EVIDENCE_COLS)
    x=score(events,evidence)
    x["source_sha256"]={
        "event_register":hashlib.sha256(args.events.read_bytes()).hexdigest(),
        "curated_stage_evidence":hashlib.sha256(args.evidence.read_bytes()).hexdigest()
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(x,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "status":x["status"],
        "original_13_events":x["original_event_count"],
        "strict_field_host_supported_events":x["original_events_with_strict_field_larval_host_evidence"],
        "original_in_panel_species_with_local_field_evidence":x["original_239_butterfly_species_with_event_qualified_field_evidence"],
        "regional_non_field_or_external_comparator_records":x["separate_native_and_new_site_source_rows_outside_original_13_events"],
        "controlled_native_Geranium_success_without_field_evidence":x["native_European_Geranium_lab_development_studies_without_wild_host_confirmation"],
        "host_origin_breakdown_field_rows":x["event_qualified_field_source_host_origin_by_row"],
    },sort_keys=True))


if __name__=="__main__":
    main()
