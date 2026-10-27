#!/usr/bin/env python3
"""Posthoc new-source discovery update WITHOUT rewriting first 13-event audit.

Keep exact first-pass 3/13 source-cohort result and add explicit evidence
discovered after first-pass completion: two site-level recipient records,
two recipient-country-without-specific-site reports and two orthogonal
ecological interpretation controls (native co-consumer, native-range use).

Neither study count nor reporting fraction estimates global host switching.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
from urllib.parse import urlparse

from audit_invasive_larval_host_primary_evidence_20261027 import (
    score, read_file, EVENT_COLS, EVIDENCE_COLS,
    FROZEN_ORIGINAL_EVENTS,
)

SCHEMA="chocho_invasive_13_event_posthoc_larval_evidence_addendum_20261027_v1"
COLS=("addendum_id","consumer_taxon","recipient_code","event_relation",
      "evidence_location","field_observation_year_or_blank","host_name",
      "host_taxonomic_certainty","host_origin_status",
      "independent_juvenile_evidence_class","site_source_quality",
      "original_study_source","evidence_url","host_origin_kew_url",
      "scientific_qualification")
EXACT_SITE={
    "FIELD_LARVAL_FEEDING_RETROSPECTIVE",
    "FIELD_LARVAE_ON_HOST_CROPS_DIRECT_STUDY",
}
COUNTRY_ONLY={
    "COUNTRY_FIELD_EGGS_AND_LARVAE_FROM_HOST",
    "COUNTRY_FIELD_CATERPILLAR_FOODPLANT_OBSERVED",
}
CONTEXT={
    "NATIVE_AND_INVASIVE_CATERPILLAR_HOST_CO_OCCURRENCE",
    "NATIVE_CONSUMER_USE_ON_ALIEN_PLANT",
}
LEVELS={
    "INVASION_RECIPIENT_SITE",
    "INVASION_RECIPIENT_COUNTRY_UNLOCATED",
    "OTHER_CONSUMER_SHARED_HOST_EVIDENCE",
    "NATIVE_BUTTERFLY_REGION_COMPARATOR",
}
ORIGINS={"INTRODUCED_PLANT_WCVP_Kew","INTRODUCED_CULTIGEN_Kew"}
EXACT_ORIGINAL_ADDENDUM={
    "A01":("Pieris brassicae","CL","INVASION_RECIPIENT_SITE"),
    "A02":("Pieris rapae","NZ","INVASION_RECIPIENT_SITE"),
    "A03":("Papilio demoleus","PR","INVASION_RECIPIENT_COUNTRY_UNLOCATED"),
    "A04":("Papilio demoleus","DO","INVASION_RECIPIENT_COUNTRY_UNLOCATED"),
    "A05":("Tatochila autodice blanchardi","CL",
           "OTHER_CONSUMER_SHARED_HOST_EVIDENCE"),
    "A06":("Pieris brassicae nepalensis","IN",
           "NATIVE_BUTTERFLY_REGION_COMPARATOR"),
}


def parse_addendum(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        r=csv.DictReader(f)
        if tuple(r.fieldnames or ())!=COLS:
            raise ValueError("Posthoc larval evidence addendum schema has changed")
        rows=list(r)
    if len(rows)!=6 or set(r["addendum_id"] for r in rows)!=set(EXACT_ORIGINAL_ADDENDUM):
        raise ValueError("Not all six previously located addendum sources retained")
    for r in rows:
        key=r["addendum_id"]
        identity=(r["consumer_taxon"],r["recipient_code"],r["event_relation"])
        if identity!=EXACT_ORIGINAL_ADDENDUM[key] or r["event_relation"] not in LEVELS:
            raise ValueError("Original posthoc evidence source case identity changed")
        if r["host_origin_status"] not in ORIGINS:
            raise ValueError("Original accepted host origin Kew evidence absent")
        if r["independent_juvenile_evidence_class"] not in EXACT_SITE|COUNTRY_ONLY|CONTEXT:
            raise ValueError("Unverified larval evidence class")
        if r["event_relation"]=="INVASION_RECIPIENT_SITE" and \
           r["independent_juvenile_evidence_class"] not in EXACT_SITE:
            raise ValueError("Direct field locality original type mismatch")
        if r["event_relation"]=="INVASION_RECIPIENT_COUNTRY_UNLOCATED" and \
           r["independent_juvenile_evidence_class"] not in COUNTRY_ONLY:
            raise ValueError("Unverified site should not be represented as located")
        if r["event_relation"] not in (
            "INVASION_RECIPIENT_SITE","INVASION_RECIPIENT_COUNTRY_UNLOCATED"
        ) and r["independent_juvenile_evidence_class"] not in CONTEXT:
            raise ValueError("Native consumer context misrepresented as invaded event")
        if r["event_relation"].startswith("INVASION") and \
           (r["consumer_taxon"],r["recipient_code"]) not in FROZEN_ORIGINAL_EVENTS:
            raise ValueError("New opportunistic butterfly invasion event added")
        for col in ("evidence_url","host_origin_kew_url"):
            u=urlparse(r[col])
            if u.scheme!="https" or not u.netloc:
                raise ValueError("Missing original independent evidence URL")
        year=r["field_observation_year_or_blank"]
        if year and not 1900<=int(year)<=2026:
            raise ValueError("Impossible field observation year")
        if key in ("A03","A04") and year:
            raise ValueError("Specific field year for country-level sources not verified")
        if key=="A01" and (
            r["host_name"]!="Tropaeolum majus" or
            r["field_observation_year_or_blank"]!="1971"):
            raise ValueError("Chilean early resource observation altered")
        if key=="A05" and r["consumer_taxon"]=="Pieris brassicae":
            raise ValueError("Native Chilean competitor incorrectly relabelled invader")
        if key=="A06" and r["field_observation_year_or_blank"]!="2018":
            raise ValueError("Native range comparator cannot predate 1971 host invasion")
    return rows


def synthesize(first_pass, addendum):
    if (first_pass.get("original_events_with_strict_field_larval_host_evidence")!=3
        or first_pass.get("original_event_count")!=13
        or first_pass.get("original_239_butterfly_species_with_event_qualified_field_evidence")!=1):
        raise ValueError("Prior first-pass original 13 event finding unexpectedly changed")
    events={(v["butterfly"],v["recipient_region_code"]):v
            for v in first_pass["original_events"]}
    if set(events)!=FROZEN_ORIGINAL_EVENTS:
        raise ValueError("Original first-pass frozen receiving-area key roster changed")
    strict={k for k,r in events.items() if r["event_compatible_field_host_evidence"]}
    additions_exact={ (r["consumer_taxon"],r["recipient_code"]) for r in addendum
                      if r["independent_juvenile_evidence_class"] in EXACT_SITE}
    additions_country={ (r["consumer_taxon"],r["recipient_code"]) for r in addendum
                        if r["independent_juvenile_evidence_class"] in COUNTRY_ONLY}
    if additions_exact & strict or additions_country & strict or \
       additions_exact & additions_country or len(additions_exact)!=2 or \
       len(additions_country)!=2:
        raise ValueError("First-pass and posthoc new locality/field evidence overlap")
    updated_site=strict|additions_exact
    updated_region=updated_site|additions_country
    context=[r for r in addendum if r["independent_juvenile_evidence_class"] in CONTEXT]
    if len(context)!=2:
        raise ValueError("Both distinct native and invaded butterfly sources must be retained")
    route=defaultdict(set)
    for r in addendum:
        if (r["consumer_taxon"],r["recipient_code"]) in updated_region:
            route[r["recipient_code"]].add(r["host_origin_status"])
    native_original_239={
        "Pieris rapae","Thymelicus lineola","Pieris brassicae",
        "Danaus plexippus","Papilio demoleus"
    }
    per_event=[]
    for key in sorted(events):
        tier=("EXACT_FIELD_STUDY_OR_FIRST_PASS_SITE" if key in updated_site else
              "COUNTRY_FIELD_EVIDENCE_UNLOCATED" if key in additions_country else
              "NO_CONFIRMED_SOURCE_LOCAL_FIELD_STAGE_IN_THIS_CURATED_SET")
        per_event.append({
            "butterfly":key[0],"recipient_code":key[1],
            "original_239":key[0] in native_original_239,
            "source_local_field_evidence_precision":tier,
            "source_field_host_use_NOT_establishment_proof":True,
        })
    return {
        "schema":SCHEMA,
        "status":"POSTHOC_LITERATURE_DISCOVERY_SOURCE_EVIDENCE_LEVELS_NOT_SAMPLING_PREVALENCE",
        "original_event_species":8,
        "original_13_invasion_recipient_events":13,
        "original_239_butterflies_in_candidate_cohort":5,
        "original_first_pass_literature_rows":15,
        "original_first_pass_field_supported_site_events":3,
        "posthoc_addendum_cited_source_rows":6,
        "posthoc_exact_site_additional_events":len(additions_exact),
        "posthoc_country_only_immatures_additional_events":len(additions_country),
        "site_resolved_any_field_larval_host_evidence_events":len(updated_site),
        "country_or_site_resolved_any_field_larval_host_evidence_events":len(updated_region),
        "events_remaining_without_recipient_field_juvenile_host_evidence":
             len(events)-len(updated_region),
        "original_239_invasion_candidate_species_with_site_field_evidence":
            len({k[0] for k in updated_site if k[0] in native_original_239}),
        "original_239_invasion_candidate_species_with_country_field_evidence":
            len({k[0] for k in updated_region if k[0] in native_original_239}),
        "context_other_native_competitor_and_native_consumer_host_reports":len(context),
        "case_status":per_event,
        "new_source_evidence_summary":[
            {"id":r["addendum_id"],"consumer":r["consumer_taxon"],
             "region":r["recipient_code"],"host":r["host_name"],
             "host_origin_Kew":r["host_origin_status"],
             "evidence_class":r["independent_juvenile_evidence_class"],
             "event_relation":r["event_relation"],
             "field_observation_year_or_blank":r["field_observation_year_or_blank"],
             "primary_source_link":r["evidence_url"],
             "plant_source_link":r["host_origin_kew_url"],
             "limitation":r["scientific_qualification"]}
            for r in addendum
        ],
        "negative_inference_guards":[
            "Counts describe how many of the hand-curated original events have any located source evidence; NOT a population prevalence estimate.",
            "Two country-level immature records cannot be promoted to event-site records without their original field locality.",
            "1971 Chile larval use of introduced T. majus is historical source evidence, not newly evolved host adaptation.",
            "1998 indigenous Tatochila vs introduced Pieris cohabitation on introduced Tropaeolum is potential resource sharing; not proven competitive suppression.",
            "2018 India P. brassicae nepalensis also using T. majus is later than Chile 1971 and cannot demonstrate the plant was utilized by ancestors before invasion.",
            "Kew present-day introduced geography does not date the arrival of that host plant.",
            "Source-host and consumer native ranges differ and some original studies may also underpin HOSTS intermediate records."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--events",type=Path,required=True)
    p.add_argument("--firstpass-evidence",type=Path,required=True)
    p.add_argument("--addendum",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    first=score(read_file(args.events,EVENT_COLS),
                read_file(args.firstpass_evidence,EVIDENCE_COLS))
    rows=parse_addendum(args.addendum)
    result=synthesize(first,rows)
    result["source_sha256"]={
        "original_13_event_manifest":hashlib.sha256(args.events.read_bytes()).hexdigest(),
        "original_precommit_15_stage_rows":hashlib.sha256(args.firstpass_evidence.read_bytes()).hexdigest(),
        "posthoc_six_literature_rows":hashlib.sha256(args.addendum.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "status":result["status"],
        "site_resolved_events":result["site_resolved_any_field_larval_host_evidence_events"],
        "site_or_country_resolved_events":result["country_or_site_resolved_any_field_larval_host_evidence_events"],
        "site_resolved_239species":result["original_239_invasion_candidate_species_with_site_field_evidence"],
        "country_resolved_239species":result["original_239_invasion_candidate_species_with_country_field_evidence"],
        "unresolved_original_events":result["events_remaining_without_recipient_field_juvenile_host_evidence"]
    },sort_keys=True))


if __name__=="__main__":
    main()
