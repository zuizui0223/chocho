"""Offline tests for source-locked GBIF dataset registry identities."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("dataset_registry_followup",ROOT/"scripts/audit_gbif_dataset_registry_20261020.py")
mod=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(mod)

def test_six_unique_provenance_dataset_ids():
    keys=mod.candidate_keys(ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_provenance_facet_qc_20261020.json")
    assert len(keys)==6
    assert "50c9509d-22c7-4a22-a47d-8c48425ef4a7" in keys
    assert "6c8d87be-829d-4364-b2bc-578258761cc2" in keys

def test_registry_key_and_title_must_validate():
    key="50c9509d-22c7-4a22-a47d-8c48425ef4a7"
    row=mod.parse_metadata(json.dumps({"key":key,"title":"iNaturalist Research-grade Observations","type":"OCCURRENCE"}).encode(),key)
    assert row["dataset_key"]==key and row["title"].startswith("iNaturalist")
    with pytest.raises(ValueError,match="mismatch"):
        mod.parse_metadata(json.dumps({"key":"mismatched","title":"Wrong"}).encode(),key)
    with pytest.raises(ValueError,match="title"):
        mod.parse_metadata(json.dumps({"key":key,"title":""}).encode(),key)
