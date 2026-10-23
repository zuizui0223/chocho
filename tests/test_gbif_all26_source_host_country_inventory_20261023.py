"""Source-only 26-positive HOSTS/WCVP catalogue; no biological inference."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import pytest

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/derive_gbif_all26_source_host_country_inventory_20261023.py"
spec=importlib.util.spec_from_file_location("source_inventory_all26",SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_source_name_gate_never_guesses_infraspecific_or_synonym():
    assert m.name_gate(["Citrus maxima"])=="EXACT_SOURCE_BINOMIAL_NEEDS_GBIF_TAXON_RESOLUTION"
    assert m.name_gate(["Citrus maxima (Burm.) Merr."])=="NOT_EXACT_SPECIES_BINOMIAL"
    assert m.name_gate(["Citrus maxima","Citrus grandis"])=="AMBIGUOUS_ACCEPTED_ID_NAME"
    assert m.name_gate([])=="AMBIGUOUS_ACCEPTED_ID_NAME"

def test_source_names_not_juvenile_outcome_selection():
    assert len(m.FOCAL)==5
    assert m.SCHEMA.endswith("_v1")
    assert m.EXPECTED_ZIP_SHA=="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"

def test_no_silent_acceptance_of_wrong_source():
    with pytest.raises(ValueError):
        m.census({"species":["A a"]},[],{},{"status":"NOT_COMPLETE"},
                 {},{},locked=True)
