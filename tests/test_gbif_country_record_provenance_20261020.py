"""No-network tests for exact source-locked GBIF dataset/year/basis facets."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts/audit_gbif_country_record_provenance_20261020.py"
SPEC = importlib.util.spec_from_file_location("country_facet_quality", SRC)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)

def sample():
    return json.dumps({"limit": 0, "results": [], "count": 9,
        "facets": [
            {"field": "DATASET_KEY", "counts": [
                {"name": "data-1", "count": 7}, {"name": "data-2", "count": 2}]},
            {"field": "YEAR", "counts": [
                {"name": "2020", "count": 4}, {"name": "2024", "count": 5}]},
            {"field": "BASIS_OF_RECORD", "counts": [
                {"name": "HUMAN_OBSERVATION", "count": 9}]}]}).encode()

def test_locked_seven_candidates_and_api_filters():
    cohort = mod.load_candidates(
        ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_record_turnover_summary_20261020.json",
        ROOT / "data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json",
        ROOT / "data/exploratory/gbif_country_20261020/gbif_country_counts.csv")
    assert len(cohort) == 7
    assert sum(p["added_host_country"] for p in cohort) == 3
    assert {(p["species"], p["country"], p["taxon_key"], p["frozen_late_records"])
            for p in cohort if p["frozen_late_records"] >= 5} == {
        ("Aguna asander", "DO", 1944574, 9),
        ("Papilio demoleus", "IL", 1938069, 86)}
    for pair in cohort:
        u = urlparse(mod.query_url(pair))
        x = parse_qs(u.query)
        assert x["taxonKey"] == [str(pair["taxon_key"])]
        assert x["year"] == ["2018,2025"]
        assert x["country"] == [pair["country"]]
        assert x["limit"] == ["0"]
        assert set(x["facet"]) == {"datasetKey","year","basisOfRecord"}
        assert set(x["basisOfRecord"]) == set(mod.BASIS)

def test_three_facets_and_concentration_summary():
    x = mod.validate_response(sample())
    assert x["record_count"] == 9
    assert x["dataset_count"] == 2
    assert x["years_recorded"] == ["2020","2024"]
    assert x["largest_dataset_fraction"] == pytest.approx(7/9)
    assert x["largest_year_fraction"] == pytest.approx(5/9)

def test_reject_missing_or_mismatched_facets_and_truncated_response():
    v = json.loads(sample())
    v["facets"].pop()
    with pytest.raises(ValueError, match="Missing"):
        mod.validate_response(json.dumps(v).encode())
    v = json.loads(sample())
    v["facets"][0]["counts"][0]["count"] = 8
    with pytest.raises(ValueError, match="Incomplete"):
        mod.validate_response(json.dumps(v).encode())
    with pytest.raises(ValueError):
        mod.validate_response(b'{"incomplete":')

def test_saved_provenance_facet_audit_complete_and_not_host_use():
    result = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_provenance_facet_qc_20261020.json"
    report = json.loads(result.read_text(encoding="utf-8"))
    assert report["schema"] == "gbif_country_record_provenance_facet_qc_v1"
    assert report["frozen_candidates"] == report["live_facets_complete"] == 7
    assert len(report["results"]) == 7
    assert sum(x["added_host_country"] for x in report["results"]) == 3
    assert all(x["status"] == "COMPLETE_LIVE_FACET_QC" and
               x["source_drift_records"] == 0 and
               sum(x["facet_counts"]["BASISOFRECORD"].values()) == x["frozen_late_records"]
               for x in report["results"])
    assert all(set(x["facet_counts"]["BASISOFRECORD"]) == {"HUMAN_OBSERVATION"}
               for x in report["results"])
    bykey = {(x["species"],x["country"]):x for x in report["results"]}
    aguna = bykey[("Aguna asander","DO")]
    papilio = bykey[("Papilio demoleus","IL")]
    assert aguna["dataset_count"] == 1 and len(aguna["years_recorded"]) == 6
    assert papilio["dataset_count"] == 3 and papilio["facet_counts"]["YEAR"]["2025"] == 61
    assert aguna["added_host_country"] == 1 and papilio["added_host_country"] == 0
    assert report["quality_audit_only_not_establishment"] is True
