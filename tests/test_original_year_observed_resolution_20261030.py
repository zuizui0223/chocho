"""Replay observed local raw-year resolution against the existing source-locked five-link screen.

The private original Parquet is NOT bundled. CI validates the restricted
selected-year derivative, source identity receipt and original resolver.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030"
SCRIPT = BASE / "resolve_original_island_chain_years.py"
SCREEN = BASE / "dated_host_island_chain_screen.json"
RESULT = BASE / "full_original_island_chain_years_observed_20261030.json"
RECEIPT = BASE / "original_year_private_source_audit_20261030.json"

def git_blob(path):
    b = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()

def load():
    spec = importlib.util.spec_from_file_location("original_year_resolution_real_snapshot", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def test_original_source_identity_and_replay_without_private_raw():
    # Anonymous GEB bundle intentionally substitutes exploratory source provenance.
    # Retain full strict SHA/source checks in normal paper CI only.
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory original-island source intentionally redacted in anonymous GEB bundle")
    assert git_blob(SCRIPT) == "f60e9d3c7b75b032727dbd1b8c5c8055ae7f4b82"
    assert git_blob(SCREEN) == "1455d93681cf99c4e7d891873a08e5109ee0cc44"
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    snapshot = json.loads(RESULT.read_text(encoding="utf-8"))
    assert snapshot["local_original_year_parquet_sha256"] == receipt["local_original_year_parquet_sha256"]
    assert len(snapshot["local_original_year_parquet_sha256"]) == 64
    assert (receipt["rows"],receipt["unique_species"],receipt["unique_island_ids"]) == (83573,224,3384)
    assert receipt["source_kind"].endswith("NOT_PUBLISHED")
    mod = load()
    source_mapping = mod.species_keys_from_followup(ROOT)
    source_screen = json.loads(SCREEN.read_text(encoding="utf-8"))
    rows = []
    for link in snapshot["rows"]:
        key = source_mapping[link["butterfly"]]
        for endpoint in link["endpoints"]:
            rows.extend({"speciesKey": key, "landmass_id": endpoint["island"], "year": int(y)}
                        for y in endpoint["all_indexed_years_from_original_year_table"])
    actual = mod.full_year_resolution(source_screen, source_mapping, rows)
    expect = {k:v for k,v in snapshot.items() if k != "local_original_year_parquet_sha256"}
    assert actual == expect

def test_near_island_pairs_are_chronologically_supported_but_have_no_nonfocal_control():
    snapshot = json.loads(RESULT.read_text(encoding="utf-8"))
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    assert snapshot["fully_resolved_pair_count"] == 5
    assert snapshot["distinct_order_and_prior_re_record_count"] == 5
    links = snapshot["rows"]
    near = [r for r in links if r["island_port_representative_distance_km"] < 500]
    far = [r for r in links if r["island_port_representative_distance_km"] > 10000]
    assert (len(near),len(far)) == (3,2)
    assert all(r["earlier_indexed_island"] == "276" for r in near)
    assert {r["butterfly"]: r["other_island_recorded_after_years"] for r in near} == {
        "Colias croceus":4,"Colias hyale":4,"Nymphalis antiopa":6}
    rerec_downstream = 0
    for r in near:
        downstream = next(e for e in r["endpoints"] if e["island"] != r["earlier_indexed_island"])
        t0 = downstream["first_indexed_year_from_original_year_table"]
        rerec_downstream += any(t0 < y <= t0+5
                                for y in downstream["all_indexed_years_from_original_year_table"])
    assert rerec_downstream == 2
    baselines = receipt["near_island_pair_baselines"]
    assert [(b["both_islands_first_indexed_2011_2022"],b["other_nonfocal_new_butterflies_same_pair"])
            for b in baselines] == [(2,0),(1,0)]
    assert receipt["cannot_fit_host_causal_effect_or_confirm_ship_route_from_this_snapshot"] is True
