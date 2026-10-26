"""Synthetic, source-independent regression tests for niche opportunity vs index."""

import hashlib
import math
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_butterfly_niche_entry_country_20261026 import (
    blob_digest, opportunity_class, geography_summary, exact_convolution,
    compare, PRIMARY,
)


def test_git_blob_hash_uses_real_nul_not_printable_escape():
    b = b"test"
    expected = hashlib.sha1(b"blob 4" + bytes([0]) + b).hexdigest()
    assert blob_digest(b) == expected


def test_native_host_climate_purity_and_strict_area_cutoffs():
    assert opportunity_class(.049, 0.) == "low_or_zero_added"
    assert opportunity_class(.1, .08) == "analogue"
    assert opportunity_class(.1, .02) == "novel"
    assert opportunity_class(.1, .05) == "mixed"
    assert opportunity_class(.1, .085, purity=.90) == "mixed"
    assert opportunity_class(.1, .015, purity=.90) == "mixed"
    assert opportunity_class(.1, .09, purity=.90) == "analogue"
    with pytest.raises(ValueError, match="Impossible"):
        opportunity_class(.1, .2)
    with pytest.raises(ValueError):
        opportunity_class(float("nan"), .01)
    with pytest.raises(ValueError):
        opportunity_class(0, 0, purity=.5)


def make_row(species, country, kind, *, native=False, level1="1"):
    return {
        "species": species, "country": country, "level1": level1,
        "country_km2": 100.0, "native": native,
        "added": .20 if kind != "none" else 0.,
        "analog": .20 if kind == "analogue" else 0.,
        "novel": .20 if kind == "novel" else 0.,
    }


def source_fixture(include_effort=True):
    rows = [
        make_row("Test butterfly", "A", "novel"),
        make_row("Test butterfly", "B", "analogue"),
        make_row("Test butterfly", "C", "novel", native=True),
        make_row("Test butterfly", "D", "none"),
    ]
    count = Counter({("Test butterfly", "A"): 2, ("Test butterfly", "C"): 3})
    if include_effort:
        for country in ("A", "B"):
            for other in ("Effort one", "Effort two", "Effort three"):
                count[(other, country)] = 1
    return rows, {("all_accepted_bases", "pooled"): count}, {
        "all_accepted_bases": {"Test butterfly", "Effort one", "Effort two", "Effort three"}
    }


def test_native_countries_and_low_area_are_explicitly_excluded():
    rows, counts, complete = source_fixture()
    result = compare(rows, counts, complete, **PRIMARY)
    assert result["status"] == "INSUFFICIENT_NICHE_ENTRY_INFORMATION"
    assert result["strata_with_both_resource_classes"] == 1
    assert result["candidate_country_species_cells"] == 2
    assert result["novel_candidate_cells"] == 1
    assert result["positive_novel_cells"] == 1
    assert result["positive_analog_cells"] == 0
    assert math.isclose(result["expected_positive_novel_cells"], 0.5)
    assert math.isclose(result["exact_conditional_reference"]["right_tail_not_biological_p"], 0.5)
    assert result["source_diagnostics"]["native_host_country_excluded"] == 1
    assert result["source_diagnostics"]["source_low_or_zero_added"] == 1


def test_no_own_species_effort_leakage():
    rows, counts, complete = source_fixture(include_effort=False)
    result = compare(rows, counts, complete, **PRIMARY)
    assert result["candidate_country_species_cells"] == 0
    assert result["strata_with_both_resource_classes"] == 0
    assert result["source_diagnostics"]["excluded_insufficient_other_butterfly_records"] == 2


def test_exact_stratified_conditional_null():
    strata = [
        {"N": 4, "A": 2, "K": 1},
        {"N": 3, "A": 1, "K": 0},
    ]
    val = exact_convolution(strata, 1)
    assert math.isclose(val["right_tail_not_biological_p"], 0.5)
    assert math.isclose(val["left_tail_not_biological_p"], 1)
    assert val["conditional_support"] == [0, 1]
    assert exact_convolution([], 0)["right_tail_not_biological_p"] == 1.


def test_area_accounting_is_weighted_by_country_area_not_equal_cells():
    rows = [
        make_row("A b", "X", "novel"),
        make_row("A b", "Y", "analogue"),
        make_row("A b", "Z", "novel", native=True),
    ]
    rows[0]["country_km2"] = 300.
    d = geography_summary(rows)
    assert math.isclose(d["novel_share_area_weighted"], .75)
    assert d["source_cells_native_guard_excluded"] == 1
    assert d["source_cells_outside_native"] == 2
