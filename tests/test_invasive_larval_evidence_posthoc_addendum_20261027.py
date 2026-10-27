"""No-network repeatable stage-source check for posthoc original 13 invasion case extension."""
from __future__ import annotations
import csv,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_invasive_larval_evidence_posthoc_addendum_20261027 import (
    parse_addendum,synthesize
)
from audit_invasive_larval_host_primary_evidence_20261027 import (
    read_file,score,EVENT_COLS,EVIDENCE_COLS
)
EVENTS=ROOT/"data/exploratory/invasive_butterfly_verified_entry_events_20261027.csv"
BASE=ROOT/"data/exploratory/invasive_butterfly_independent_larval_stages_20261027.csv"
ADD=ROOT/"data/exploratory/invasive_butterfly_larval_evidence_posthoc_addendum_20261027.csv"


def prepared():
    src=score(read_file(EVENTS,EVENT_COLS),read_file(BASE,EVIDENCE_COLS))
    add=parse_addendum(ADD)
    return src,add


def test_original_first_pass_remains_fixed_and_new_evidence_is_labeled_posthoc():
    base,add=prepared()
    assert base["original_events_with_strict_field_larval_host_evidence"]==3
    assert len(add)==6
    out=synthesize(base,add)
    assert out["original_first_pass_field_supported_site_events"]==3
    assert out["posthoc_exact_site_additional_events"]==2
    assert out["posthoc_country_only_immatures_additional_events"]==2
    assert out["site_resolved_any_field_larval_host_evidence_events"]==5
    assert out["country_or_site_resolved_any_field_larval_host_evidence_events"]==7
    assert out["events_remaining_without_recipient_field_juvenile_host_evidence"]==6
    assert out["original_239_invasion_candidate_species_with_site_field_evidence"]==3
    assert out["original_239_invasion_candidate_species_with_country_field_evidence"]==4
    assert len(out["case_status"])==13
    assert len(out["new_source_evidence_summary"])==6
    assert "POSTHOC" in out["status"]


def test_chile_introduced_tropaeolum_is_host_ecology_not_evolved_new_ability():
    _,add=prepared()
    a=next(x for x in add if x["addendum_id"]=="A01")
    assert a["consumer_taxon"]=="Pieris brassicae"
    assert a["recipient_code"]=="CL"
    assert a["field_observation_year_or_blank"]=="1971"
    assert a["host_name"]=="Tropaeolum majus"
    assert a["host_origin_status"]=="INTRODUCED_PLANT_WCVP_Kew"
    b=next(x for x in add if x["addendum_id"]=="A06")
    assert b["consumer_taxon"]=="Pieris brassicae nepalensis"
    assert b["field_observation_year_or_blank"]=="2018"
    assert b["event_relation"]=="NATIVE_BUTTERFLY_REGION_COMPARATOR"


def test_puerto_rico_is_not_automatically_source_verified_at_specimen_site():
    base,add=prepared()
    out=synthesize(base,add)
    x=next(q for q in out["case_status"] if
           q["butterfly"]=="Papilio demoleus" and q["recipient_code"]=="PR")
    assert x["source_local_field_evidence_precision"]== \
        "COUNTRY_FIELD_EVIDENCE_UNLOCATED"
    assert "SOURCE_SITE" not in x["source_local_field_evidence_precision"]
    y=next(q for q in out["case_status"] if
           q["butterfly"]=="Papilio demoleus" and q["recipient_code"]=="DO")
    assert y["source_local_field_evidence_precision"]== \
        "COUNTRY_FIELD_EVIDENCE_UNLOCATED"


def test_shared_host_indigenous_consumer_does_not_change_original_event_count():
    base,add=prepared()
    out=synthesize(base,add)
    assert len([q for q in out["new_source_evidence_summary"] if
                q["event_relation"]=="OTHER_CONSUMER_SHARED_HOST_EVIDENCE"])==1
    assert out["context_other_native_competitor_and_native_consumer_host_reports"]==2
    assert out["original_13_invasion_recipient_events"]==13


def test_false_pre_1971_native_host_repertoire_chronology_rejected(monkeypatch):
    _,src=prepared()
    rows=[dict(x) for x in src]
    rows[5]["field_observation_year_or_blank"]="1965"
    # parse_addendum validates source identities, chronology, and source URLs;
    # write mutated evidence into scratch clone rather than silently editing source.
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as d:
        p=Path(d)/"mutated.csv"
        with p.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        with pytest.raises(ValueError,match="cannot predate"):
            parse_addendum(p)


def test_wrong_firstpass_results_cannot_be_hidden_in_new_study():
    first,add=prepared()
    first=dict(first)
    first["original_events_with_strict_field_larval_host_evidence"]=4
    with pytest.raises(ValueError,match="Prior first-pass"):
        synthesize(first,add)


def test_unverified_plants_do_not_become_introduced_without_Kew_source():
    _,rows=prepared()
    rows=[dict(r) for r in rows]
    rows[0]["host_origin_status"]="PLANT_STATUS_UNRESOLVED"
    with pytest.raises(ValueError,match="accepted host origin"):
        # parse stored CSV, not the in-memory result
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"wrong.csv"
            with p.open("w",newline="",encoding="utf-8") as out:
                w=csv.DictWriter(out,fieldnames=list(rows[0]))
                w.writeheader();w.writerows(rows)
            parse_addendum(p)
