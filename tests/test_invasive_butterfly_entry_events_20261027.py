"""Source-independent controls for invasion event chronology and original 239 membership."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from validate_invasive_butterfly_entry_events_20261027 import evaluate, read_register

REGISTER = ROOT / "data/exploratory/invasive_butterfly_verified_entry_events_20261027.csv"

def fake_original_239():
    events = read_register(REGISTER)
    core = {r["species"] for r in events if r["in_original_239"] == "true"}
    assert len(core) == 5
    return sorted(core) + [f"Calibration species {n}" for n in range(234)]


def test_current_curated_register_unique_and_has_explicit_statuses():
    events = read_register(REGISTER)
    assert len(events) == 13
    assert len({r["species"] for r in events}) == 8
    assert len({(r["species"], r["recipient_code"]) for r in events}) == 13
    r = evaluate(REGISTER, fake_original_239())
    assert r["candidate_in_original_239"] == 5
    assert r["candidate_outside_original_239"] == 3
    assert r["distinct_species_recipient_event_count"] == 13
    assert r["source_original_239_species_verified"] == 239


def test_wrong_239_panel_membership_fails_before_publication():
    source = fake_original_239()
    source[source.index("Papilio demoleus")] = "Wrong replacement species"
    with pytest.raises(ValueError, match="eligibility drift"):
        evaluate(REGISTER, source)


def test_recipient_year_range_distinguishes_first_records_and_establishment():
    events = read_register(REGISTER)
    x = next(r for r in events
             if r["species"] == "Pieris rapae" and r["recipient_code"] == "AU")
    assert (x["first_year_lower"], x["first_year_upper"]) == ("1929", "1929")
    assert x["establishment_by_year"] == "1937"
    x = next(r for r in events
             if r["species"] == "Cacyreus marshalli")
    assert x["first_year_lower"] == "1987"
    assert x["establishment_by_year"] == "1989"


def test_butterfly_own_arrival_is_not_assumed_from_introduced_host_status():
    events = read_register(REGISTER)
    x = [r for r in events if r["species"] == "Danaus plexippus"]
    assert len(x) == 3
    assert all("uncertain" in row["arrival_mechanism_status"]
               or "possible" in row["arrival_mechanism_status"] for row in x)


def test_no_false_original_239_inclusion_for_added_case_studies():
    events = read_register(REGISTER)
    off = {r["species"] for r in events
           if r["in_original_239"] == "false"}
    assert off == {"Cacyreus marshalli", "Chilades pandava", "Papilio xuthus"}
