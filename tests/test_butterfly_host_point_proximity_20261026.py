"""No-network regression tests for fixed butterfly/host point validation."""

import json
import math
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_butterfly_host_point_proximity_20261026 import (
    distance_km, evaluate, fetch_pair, records_url, source_inputs, valid_point,
)

BASE = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020"
BUTTERFLIES = BASE / "gbif_positive48_gadm_lifestage_facets_20261020.json"
PLANTS = BASE / "gbif_all26_host_gadm_source_20261023.json"
INVENTORY = BASE / "gbif_all26_source_host_country_inventory_20261023.json"


def test_exact_original_source_roster_and_names():
    cases, plants, meta = source_inputs(INVENTORY, BUTTERFLIES, PLANTS)
    assert len(cases) == 26 and len(plants) == 44
    assert meta["frozen_case_count"] == 26
    assert meta["frozen_unique_host_species"] == 24
    assert meta["frozen_butterfly_source_count"] == 467
    assert meta["frozen_plant_source_count"] == 3595
    assert sum(len(r["hosts"]) for r in cases.values()) == 46


def test_query_preserves_exact_gbif_country_accepted_taxon_and_basis():
    source = json.loads(BUTTERFLIES.read_text())
    r = source["observations"][0]
    url = records_url(r["query_url"], r["taxon_key"], r["country"])
    x = parse_qs(urlparse(url).query)
    assert x["taxonKey"] == [str(r["taxon_key"])]
    assert x["country"] == [r["country"]]
    assert x["year"] == ["2010,2025"]
    assert len(x["basisOfRecord"]) == 5
    assert "facet" not in x
    assert x["limit"] == ["300"] and x["offset"] == ["0"]
    with pytest.raises(ValueError):
        records_url(r["query_url"], r["taxon_key"] + 1, r["country"])
    with pytest.raises(ValueError):
        records_url(r["query_url"], r["taxon_key"], "XX")


def make_point(lat, lon, year, uncertainty=None, stage=None, text=None):
    return {"lat": lat, "lon": lon, "year": year, "uncertainty": uncertainty,
            "strict": uncertainty is not None and uncertainty <= 1000,
            "stage_juvenile": stage == "Larva",
            "annotation_text": (text or "").lower() if stage == "Larva" else "",
            "dataset": "synthetic"}


def test_haversine_distance_and_temporal_matching():
    b = make_point(30, -70, 2021, 200, "Larva",
                   "Larva observed feeding on Bauhinia variegata.")
    near = make_point(30.01, -70, 2023, 100)
    old = make_point(30, -70, 2016, 100)
    assert 1 < distance_km(b, near) < 1.2
    result = evaluate([b], [near, old], ["Bauhinia variegata"])
    assert result["three_years"]["1km"] == 0
    assert result["three_years"]["5km"] == 1
    assert result["same_year"]["5km"] == 0
    assert result["all_years"]["1km"] == 1
    assert result["strict_known_uncertainty_1km"]["5km_3years"] == 1
    assert result["literal_juvenile_host_action_annotations"] == 1
    assert result["butterfly_with_any_year_matched_plant"] == 1


def test_unknown_coordinate_uncertainty_not_treated_as_known_precise():
    raw = {"key": 23, "decimalLatitude": 10., "decimalLongitude": 20.,
           "year": 2020, "countryCode": "CU", "speciesKey": 3333,
           "lifeStage": "Adult"}
    p, label = valid_point(raw, 3333, "CU")
    assert label == "valid" and not p["strict"]
    assert p["uncertainty"] is None
    raw["coordinateUncertaintyInMeters"] = 500
    p, label = valid_point(raw, 3333, "CU")
    assert label == "valid" and p["strict"]
    raw["speciesKey"] = 3334
    assert valid_point(raw, 3333, "CU")[1] == "accepted_species_key_mismatch"
    raw["speciesKey"] = 3333
    raw["decimalLatitude"] = 999
    assert valid_point(raw, 3333, "CU")[1] == "invalid_coordinate_or_year"


class DummyThrottle:
    def __init__(self, batches):
        self.batches = batches
        self.i = 0

    def fetch(self, url):
        if self.i >= len(self.batches):
            raise AssertionError("Unexpected additional GBIF source page")
        data = self.batches[self.i]
        self.i += 1
        return data, "0" * 64


def make_query():
    b = json.loads(BUTTERFLIES.read_text())["observations"][0]
    return {"gbif_species_key": b["taxon_key"], "country": b["country"],
            "original_faceted_index_count": b["record_count"],
            "source_url": b["query_url"]}


def test_complete_and_bounded_zero_index_responses():
    q = make_query()
    x, points = fetch_pair(q, DummyThrottle([{"count": 0, "results": []}]), 100)
    assert x["status"] == "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT"
    assert x["count"] == 0 and not points
    x, points = fetch_pair(q, DummyThrottle([{"count": 5001, "results": []}]), 12000)
    assert x["status"] == "OVER_BOUNDED_ACQUISITION_CAP"
    assert not points


def test_page_shortfall_or_duplicate_results_fail_closed():
    q = make_query()
    x, p = fetch_pair(q, DummyThrottle([{"count": 2, "results": [{"key": 1}]}]), 12000)
    assert x["status"] == "INCOMPLETE_GBIF_SOURCE"
    assert not p
    x, p = fetch_pair(q, DummyThrottle([{"count": 2, "results": [
        {"key": 1}, {"key": 1}
    ]}]), 12000)
    assert x["status"] == "INCOMPLETE_GBIF_SOURCE" and not p


def test_valid_point_needs_year_and_matching_country():
    p, reason = valid_point({"decimalLatitude": 32., "decimalLongitude": -67.,
                             "year": 2020, "speciesKey": 5,
                             "countryCode": "CU"}, 5, "CU")
    assert reason == "valid"
    _, reason = valid_point({"decimalLatitude": 32., "decimalLongitude": -67.,
                             "year": 2020, "speciesKey": 5,
                             "countryCode": "GB"}, 5, "CU")
    assert reason == "country_mismatch"
