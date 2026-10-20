"""Nontrivial scientific invariants for the independent country-grain analysis."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_gbif_country_host_opportunity_20261020 import (
    classify_country, collision_species, summarize_strata, build_strata,
)


def country():
    return {"country": "AA", "area": 100.0, "level1": "1",
            "material": {0: 90.0, 1: 5.0},
            "guard": {0: 90.0, 1: 5.0, 2: 1.0},
            "atoms": [{"members": {0}, "area": 90.}, {"members": {1}, "area": 5.},
                      {"members": {2}, "area": 1.}]}


def test_minor_native_guard_prevents_false_outside_native_country():
    result = classify_country(country(), native=1 << 2, added=1, climate=3)
    assert result["native_guard_present"] is True
    assert result["added_country"] is False
    assert result["native_fraction"] == pytest.approx(.01)


def test_unmapped_area_is_not_renormalized_out_of_climate_fraction():
    c = country()
    c["material"] = {0: 89.0, 1: 6.0}
    c["atoms"] = [{"members": {0}, "area": 89.}, {"members": {1}, "area": 6.},
                  {"members": {2}, "area": 1.}]
    result = classify_country(c, native=0, added=1, climate=1)
    assert result["climate_fraction"] == pytest.approx(.89)
    assert result["added_climate_fraction"] == pytest.approx(.89)
    assert result["added_country"] is True


def test_overlapping_regions_count_geographic_area_once():
    c = {"area": 100., "material": {0: 50., 1: 50.}, "guard": {0: 50., 1: 50.},
         "atoms": [{"members": {0}, "area": 30.}, {"members": {1}, "area": 30.},
                   {"members": {0, 1}, "area": 20.}]}
    x = classify_country(c, native=0, added=3, climate=3)
    assert x["added_fraction"] == pytest.approx(.8)
    assert x["climate_fraction"] == pytest.approx(.8)
    assert x["added_climate_fraction"] == pytest.approx(.8)


def test_shared_accepted_key_excludes_every_alias_instead_of_arbitrary_choice():
    records = [{"species": "Genus a", "taxon_key": 3},
               {"species": "Genus b", "taxon_key": 3},
               {"species": "Other c", "taxon_key": 4},
               {"species": "Unknown d", "taxon_key": None}]
    assert collision_species(records) == {"Genus a", "Genus b"}


def test_known_conditional_comparison_and_all_positive_no_information():
    strata = [{"species": "Genus a", "genus": "Genus", "stratum": "1",
               "N": 4, "A": 2, "K": 2, "O": 2, "expected": 1.0}]
    result = summarize_strata(strata, label="synthetic", null_draws=199, bootstrap_draws=99)
    assert result["observed_minus_expected_pp"] == pytest.approx(50.0)
    assert result["information"]["nonzero_conditional_variance_strata"] == 1
    assert result["status"] == "BELOW_PROTOCOL_INFORMATION_THRESHOLD"
    all_positive = [{**strata[0], "K": 4, "O": 2, "expected": 2.0}]
    assert summarize_strata(all_positive, label="saturated", null_draws=9,
                            bootstrap_draws=9)["information"]["nonzero_conditional_variance_strata"] == 0


def test_continuous_score_with_no_area_variation_has_no_information():
    cells = [{"country": code, "area_score": 0., "positive": int(i == 0)}
             for i, code in enumerate(("AA", "BB"))]
    strata = [{"species": "Genus a", "genus": "Genus", "stratum": "1",
               "N": 2, "A": 1, "K": 1, "O": 1, "expected": .5, "cells": cells}]
    result = summarize_strata(strata, label="constant_area", null_draws=9,
                              bootstrap_draws=9, continuous=True)
    assert result["observed_minus_expected_pp"] == 0
    assert result["continuous_score_informative_strata"] == 0
    assert result["status"] == "BELOW_PROTOCOL_INFORMATION_THRESHOLD"


def test_focal_observation_cannot_satisfy_its_own_effort_gate():
    countries = {"AA": {"country": "AA", "area": 100., "level1": "1", "material": {0: 100.}},
                 "BB": {"country": "BB", "area": 100., "level1": "1", "material": {1: 100.}}}
    exposures = {("Genus a", "AA", .05): {"native_guard_present": False,
                   "added_country": True, "added_fraction": 1., "climate_fraction": 1.,
                   "all_material_climate": True, "added_climate_fraction": 1.},
                 ("Genus a", "BB", .05): {"native_guard_present": False,
                   "added_country": False, "added_fraction": 0., "climate_fraction": 1.,
                   "all_material_climate": True, "added_climate_fraction": 0.}}
    counts = {("Genus a", "AA"): 1, ("Other b", "BB"): 1}
    rows, audit = build_strata(countries, exposures, counts, targets={"Genus a"},
        effort_species={"Genus a", "Other b"}, gate=1, q=.05, record_minimum=1)
    assert rows == []
    assert audit["candidate_country_species_cells"] == 1


def test_strict_added_exposure_drops_intermediate_country_instead_of_calling_negative():
    countries = {code: {"country": code, "area": 100., "level1": "1", "material": {0: 100.}}
                 for code in ("AA", "BB", "CC")}
    exposures = {("Genus a", code, .05): {"native_guard_present": False,
                  "added_country": frac > 0, "added_fraction": frac,
                  "climate_fraction": 1., "all_material_climate": True,
                  "added_climate_fraction": frac}
                 for code, frac in (("AA", 1.), ("BB", 0.), ("CC", .2))}
    counts = {("Other b", code): 1 for code in countries}
    rows, audit = build_strata(countries, exposures, counts, targets={"Genus a"},
        effort_species={"Other b"}, gate=1, q=.05, record_minimum=1,
        geography_mode="added90")
    assert len(rows) == 1
    assert rows[0]["N"] == 2 and rows[0]["A"] == 1
    assert audit["excluded_intermediate_added_cells"] == 1
