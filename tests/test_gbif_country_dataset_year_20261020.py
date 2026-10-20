"""No-network contract tests for exact 7-cell × 11-dataset year source audit."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts/audit_gbif_country_dataset_year_20261020.py"
SPEC = importlib.util.spec_from_file_location("dataset_year_qc", SRC)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)
QC = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_provenance_facet_qc_20261020.json"
REG = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_dataset_registry_qc_20261020.json"


def mock(years):
    return json.dumps({"count": sum(years.values()), "limit": 0, "results": [],
                       "facets": [{"field": "YEAR", "counts": [
                           {"name": str(y), "count": n} for y, n in years.items()]}]}).encode()


def test_source_seven_cells_eleven_dataset_combinations_and_filter_lock():
    cohort, meta = mod.load_fixed_pairs(QC, REG)
    assert meta["seven_candidate_cells"] == 7
    assert meta["dataset_country_cells"] == 11
    assert meta["dataset_registry_keys"] == 6
    assert len(cohort) == 11
    assert {p["source_country_total"] for p in cohort
            if (p["species"], p["country"]) == ("Papilio demoleus", "IL")} == {86}
    for pair in cohort:
        args = parse_qs(urlparse(mod.query_url(pair)).query)
        assert args["year"] == ["2018,2025"]
        assert args["country"] == [pair["country"]]
        assert args["datasetKey"] == [pair["dataset_key"]]
        assert args["taxonKey"] == [str(pair["taxon_key"])]
        assert set(args["basisOfRecord"]) == set(mod.BASIS)
        assert args["hasCoordinate"] == ["true"]
        assert args["hasGeospatialIssue"] == ["false"]
        assert args["occurrenceStatus"] == ["PRESENT"]
        assert args["limit"] == ["0"]
        assert args["facet"] == ["year"]


def test_validate_year_facet_and_reject_partial_duplicates():
    assert mod.validate_year_facet(mock({"2023": 2, "2025": 3})) == {
        "record_count": 5, "year_counts": {"2023": 2, "2025": 3}}
    base = json.loads(mock({"2024": 2, "2025": 3}))
    base["facets"][0]["counts"][1]["count"] = 2
    with pytest.raises(ValueError, match="sum"):
        mod.validate_year_facet(json.dumps(base).encode())
    base = json.loads(mock({"2024": 2, "2025": 3}))
    base["facets"][0]["counts"].append({"name": "2024", "count": 1})
    with pytest.raises(ValueError, match="Duplicated"):
        mod.validate_year_facet(json.dumps(base).encode())
    with pytest.raises(ValueError):
        mod.validate_year_facet(b'{"broken":')
    with pytest.raises(ValueError, match="out-of-range"):
        mod.validate_year_facet(mock({"2026": 1}))


def test_joint_country_year_reconciliation_does_not_mask_dataset_drift():
    pairs, meta = mod.load_fixed_pairs(QC, REG)
    fetched = []
    for p in pairs:
        # Toy deterministic conditional split: put *all* each dataset's
        # counts into one year => no source country year marginal match.
        r = {**p, "status": "COMPLETE_MATCHED_SOURCE_COUNT", "record_count": p["source_dataset_count"],
             "year_counts": {"2025": p["source_dataset_count"]}}
        fetched.append(r)
    doc = mod.assemble(pairs, fetched, meta)
    assert doc["status"] == "INCOMPLETE_NO_INFERENCE"
    assert doc["complete_dataset_queries"] == 11
    assert doc["complete_country_cells"] < 7
    fetched[0]["status"] = "INDEX_DRIFT_NOT_USED_FOR_SOURCE_INFERENCE"
    doc = mod.assemble(pairs, fetched, meta)
    assert doc["status"] == "INCOMPLETE_NO_INFERENCE"
    assert any(x["status"] == "SOURCE_DRIFT_OR_UNAVAILABLE_NO_INFERENCE" for x in doc["country_cells"])


def test_saved_eleven_dataset_year_query_matrix_reconciles_original():
    p = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_dataset_year_matrix_20261020.json"
    assert p.is_file()
    output = json.loads(p.read_text(encoding="utf-8"))
    assert output["schema"] == "gbif_country_seven_candidate_dataset_by_year_qc_v1"
    assert output["status"] == "ALL_SEVEN_SOURCE_CELLS_RECONCILED"
    assert output["complete_country_cells"] == 7
    assert output["complete_dataset_queries"] == 11
    source = json.loads(QC.read_text(encoding="utf-8"))
    reference = {(r["species"],r["country"]):r for r in source["results"]}
    assert len(output["country_cells"]) == 7
    for item in output["country_cells"]:
        row = reference[(item["species"],item["country"])]
        assert item["country_years"] == row["facet_counts"]["YEAR"]
        assert item["record_count"] == row["frozen_late_records"]
        assert sum(d["record_count"] for d in item["datasets"]) == item["record_count"]
        assert all(d["status"] == "COMPLETE_MATCHED_SOURCE_COUNT" for d in item["datasets"])
    israel = output["primary_two_case_readout"]["israel_papilio"]
    assert israel["record_count"] == 86 and israel["2025_total"] == 61
    assert israel["2025_by_dataset_title"] == {
        "iNaturalist Research-grade Observations":33,
        "INPA - Observations":27,"eButterfly Surveys":1}
    assert israel["inpa_2025_contribution"] == 27
    assert israel["non_inpa_year_marginals"]["2024"] == 9
    assert israel["non_inpa_year_marginals"]["2025"] == 34
    assert output["outcomes_cannot_establish_butterfly_colonization"] is True

