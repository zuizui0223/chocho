"""Non-network validity and counterexample tests for independent larval-host stage panel."""
from __future__ import annotations
import csv
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_invasive_larval_host_primary_evidence_20261027 import (
    score, read_file, EVENT_COLS, EVIDENCE_COLS,
)

EVENTS=ROOT/"data/exploratory/invasive_butterfly_verified_entry_events_20261027.csv"
STUDIES=ROOT/"data/exploratory/invasive_butterfly_independent_larval_stages_20261027.csv"


def loaded():
    return read_file(EVENTS,EVENT_COLS),read_file(STUDIES,EVIDENCE_COLS)


def test_exact_original_event_cohort_and_epistemic_hierarchy():
    events,studies=loaded()
    assert len(events)==13 and len(studies)==15
    x=score(events,studies)
    assert x["original_event_count"]==13
    assert x["original_event_239_butterfly_species"]==5
    assert x["curated_independent_publication_evidence_rows"]==15
    assert len(x["original_events"])==13
    assert (x["original_events_with_strict_field_larval_host_evidence"]+
            x["original_events_without_source_qualified_field_larval_host_evidence"])==13
    assert "FIRST_DOCUMENTED_PEER_REVIEWED_CAPTURE_2004" in \
        x["chronology_discrepancies"][0]["resolution"]
    assert x["case_evidence_source_independence"].startswith("PRIMARY_RESEARCH_CITATION")


def test_native_island_cycad_is_not_relabelled_as_moved_plant_resource():
    events,studies=loaded()
    result=score(events,studies)
    row=next(r for r in result["original_events"]
             if r["butterfly"]=="Chilades pandava" and
                r["recipient_region_code"]=="GU")
    assert row["event_compatible_field_host_evidence"]
    assert row["strict_field_evidence_host_origin_classes"]==[
        "NATIVE_AT_RECIPIENT"
    ]
    assert "CONTROLLED_ONLY_OVIPOSITION" not in row["original_field_evidence_tiers"]


def test_controlled_only_geranium_adult_emergence_is_not_wild_host_use():
    events,studies=loaded()
    x=score(events,studies)
    assert x["native_European_Geranium_lab_development_studies_without_wild_host_confirmation"]==2
    sp=[r for r in x["original_events"]
        if r["butterfly"]=="Cacyreus marshalli"]
    assert sp and not any(r["event_compatible_field_host_evidence"] for r in sp)
    assert x["separate_native_and_new_site_source_rows_outside_original_13_events"]>=3


def test_regional_review_or_adult_citrus_never_promoted_to_larval_field():
    e,s=loaded()
    x=score(e,s)
    for region in ("DO","PR"):
        p=next(r for r in x["original_events"]
               if r["butterfly"]=="Papilio demoleus" and
                  r["recipient_region_code"]==region)
        assert not p["event_compatible_field_host_evidence"]
        assert p["no_field_host_evidence_is_NOT_ecological_absence"]
    for r in x["original_events"]:
        assert "ADULT_OR_ESTABLISHMENT_ONLY" not in r["original_field_evidence_tiers"]


def test_wrong_geographic_or_botanical_promotions_fail():
    events,studies=loaded()
    corrupt=[dict(r) for r in studies]
    for row in corrupt:
        if row["record_id"]=="E01":
            row["recipient_host_origin"]="ALIEN_NATURALIZED_RECORDED"
    with pytest.raises(ValueError,match="Guam native cycad"):
        score(events,corrupt)
    corrupt=[dict(r) for r in studies]
    for row in corrupt:
        if row["record_id"]=="E05":
            row["evidence_tier"]="FIELD_LARVA_HOST_REARED"
            row["field_derived_local_host_record"]="true"
    with pytest.raises(ValueError,match="not confirmed wild"):
        score(events,corrupt)
    corrupt=[dict(r) for r in studies]
    for row in corrupt:
        if row["record_id"]=="E04":
            row["recipient_host_origin"]="ALIEN_NATURALIZED_RECORDED"
    with pytest.raises(ValueError,match="Genus-only host"):
        score(events,corrupt)


def test_an_unverified_location_cannot_create_a_new_original_event():
    e,s=loaded()
    corrupt=[dict(r) for r in s]
    for row in corrupt:
        if row["record_id"]=="E04":
            row["geographic_scope"]="RECIPIENT_SITE"
    with pytest.raises(ValueError,match="cannot silently add"):
        score(e,corrupt)


def test_original_event_identifiers_are_not_mutated_by_new_literature():
    events,studies=loaded()
    corrupt=[dict(e) for e in events]
    corrupt[0]["recipient_code"]="XX"
    with pytest.raises(ValueError,match="Original frozen"):
        score(corrupt,studies)
