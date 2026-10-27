"""Offline tests for the 73-case exact-source GBIF point acquisition gate."""
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_butterfly_73_point_source_20261027 import (
    query_url, valid_metadata, summarize_one, retrieve_case, make_report,
)


def case():
    return {
        "species": "Test butterfly",
        "country": "CU",
        "resource_class": "novel",
        "gbif_accepted_species_key": 1944574,
        "original_gbif_record_count": 2,
    }


def record(idnum=1, year=2023, stage="Adult", uncertainty=500):
    return {
        "key": idnum, "decimalLatitude": 23.0 + idnum / 10000.,
        "decimalLongitude": -81.0, "year": year,
        "speciesKey": 1944574, "countryCode": "CU",
        "coordinateUncertaintyInMeters": uncertainty,
        "basisOfRecord": "HUMAN_OBSERVATION",
        "lifeStage": stage,
        "datasetKey": "TEST_DATA",
    }


class Fetcher:
    def __init__(self, documents):
        self.docs = documents
        self.calls = 0
        self.retries = 0

    def get(self, url):
        self.calls += 1
        if not self.docs:
            raise AssertionError("Unexpected external query")
        return self.docs.pop(0), "3" * 64


def test_query_identifies_one_source_taxon_and_original_bases():
    url = query_url(case(), 0)
    q = parse_qs(urlparse(url).query)
    assert q["taxonKey"] == ["1944574"]
    assert q["country"] == ["CU"]
    assert q["year"] == ["2010,2025"]
    assert q["hasGeospatialIssue"] == ["false"]
    assert len(q["basisOfRecord"]) == 5
    assert q["limit"] == ["300"]
    with pytest.raises(ValueError):
        query_url(case(), 1)
    with pytest.raises(ValueError):
        query_url(case(), -300)
    with pytest.raises(ValueError):
        query_url({**case(), "resource_class": "unknown"}, 0)


def test_source_point_metadata_is_not_imputed_when_missing():
    good, status = valid_metadata(record(), 1944574, "CU")
    assert status == "valid" and good["strict"]
    assert good["stage"] == "adult"
    assert good["reported_establishment_annotation"] is False
    missing, status = valid_metadata({**record(), "lifeStage": None,
                                      "coordinateUncertaintyInMeters": None},
                                     1944574, "CU")
    assert status == "valid" and not missing["strict"]
    assert missing["stage"] == "missing_stage"
    wrong, status = valid_metadata({**record(), "speciesKey": 9999},
                                   1944574, "CU")
    assert wrong is None and status == "accepted_species_key_mismatch"


def test_full_index_not_wrongly_treated_as_independent_populations():
    r = summarize_one(case(), [record(1), record(2, stage="Larva")], [
        {"offset": 0, "records": 2, "source_body_sha256": "h"}], 2)
    assert r["full_source_status"] == "COMPLETE_NEW_LIVE_GBIF_POINTS"
    assert r["validated_point_indices"] == 2
    assert r["point_source_reporting_year_counts"] == {"2018_2025": 2}
    assert r["point_source_stage_labels"] == {"adult": 1, "juvenile_source_label": 1}
    assert r["point_0p01degree_reporting_cells"] == 1
    assert r["climate_annotation_status"].startswith("NOT_MEASURED")
    with pytest.raises(ValueError, match="duplicate"):
        summarize_one(case(), [record(1), record(1)], [], 2)


def test_exhaustive_page_and_overbudget_gates():
    valid, ok = retrieve_case(case(), Fetcher([
        {"count": 2, "results": [record(1), record(2)]}
    ]), remaining=100)
    assert valid["full_source_status"] == "COMPLETE_NEW_LIVE_GBIF_POINTS"
    assert valid["new_live_indexed_count"] == 2
    cap = retrieve_case(case(), Fetcher([
        {"count": 5001, "results": []}
    ]), remaining=40000)
    assert cap["full_source_status"] == "NOT_EVALUABLE_PREDECLARED_RECORD_CAP"
    bad = retrieve_case(case(), Fetcher([
        {"count": 3, "results": [record(1), record(2)]}
    ]), remaining=40000)
    assert bad["full_source_status"] == "INCOMPLETE_GBIF_SOURCE_QUERY"
    partial = retrieve_case(case(), Fetcher([
        {"count": 2, "results": [record(1), record(1)]}
    ]), remaining=40000)
    assert partial["full_source_status"] == "INCOMPLETE_GBIF_SOURCE_QUERY"


def test_no_inference_from_partial_source_snapshot():
    source = {"schema": "test"}
    complete = {
        "full_source_status": "COMPLETE_NEW_LIVE_GBIF_POINTS",
        "resource_country_class": "analogue", "butterfly": "A", "country": "CU",
        "validated_point_indices": 2, "new_live_indexed_count": 2,
        "reported_butterfly_coord_uncertainty_le1000m": 1,
        "point_source_reporting_year_counts": {"2018_2025": 2},
        "point_source_stage_labels": {"adult": 2},
    }
    reader = Fetcher([])
    r = make_report(source, [complete], reader)
    assert r["status"] == "INCOMPLETE_73_CASE_QUERY_SNAPSHOT"
    assert r["cases_complete"] == 1
    assert r["point_climate_not_measured"] is True
