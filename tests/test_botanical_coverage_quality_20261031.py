"""Public CI validates source-locked GBIF botanical index summary; private zip is optional."""
from __future__ import annotations
from pathlib import Path
import json,math
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_host_controls_20261031"
def test_botanical_source_quality_reconciles():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Private botanical index data excluded from GEB blinded review material")
    q=json.loads((D/"botanical_coverage_quality_20261031.json").read_text())
    s=json.loads((D/"daily_botanical_port_pilot_20261031.json").read_text())
    assert q["source_artifact_gzip_sha256"]==s["taxon_year_archive"]["sha256"]
    assert q["source_species_year_facet_rows"]==s["taxon_year_archive"]["rows"]==4843
    assert q["unique_species_keys_across_nine_ports"]==1447
    assert q["source_api_counts"]["total_indexed_botanical_records"]==sum(
        x["plant_records_indexed"] for x in s["yearly_facets_recording_support"])==14742
    assert q["source_api_counts"]["records_without_species_key"]==sum(
        x["records_without_species_key"] for x in s["yearly_facets_recording_support"])==661
    assert q["source_api_counts"]["records_assigned_species_key"]==14081
    by_site={p["island_id"]:p for p in q["sites"]}
    assert len(by_site)==s["pilot_islands"]==9
    assert sum(x["records_with_species_key"] for x in q["sites"])==14081
    x=by_site["81"]
    assert x["unique_species_keys"]==4
    assert x["top_species_records"]==164
    assert math.isclose(x["top_species_fraction"],164/171,abs_tol=1e-6)
    assert x["distinct_recorded_years"]==3
    assert x["inverse_simpson_effective_species"]<2
    assert len(q["pairs"])==7 and q["observed_host_absence_proven"] is False
    assert q["original_port_counts_exactly_matched"]==2
    assert q["original_port_counts_differed"]==7
