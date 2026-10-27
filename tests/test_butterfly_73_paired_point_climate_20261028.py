"""Frozen aggregate provenance / pair-thinning invariants. No internet or raw points."""
from __future__ import annotations

import sys
from copy import deepcopy
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_butterfly_73_paired_point_climate_20261028 import (
    SOURCE_CSV_SHA256, DEFAULT_SOURCE, read_cases, compare, run,
)


def test_frozen_case_checksum_and_source_accounting():
    rows, digest = read_cases(DEFAULT_SOURCE)
    assert digest == SOURCE_CSV_SHA256
    assert len(rows) == 73
    assert len({r["species"] for r in rows}) == 20
    assert sum(r["n"] for r in rows) == 6076
    out = run(rows, digest)
    assert out["frozen_cohort"]["single_record_cases"] == 18
    assert out["frozen_cohort"]["cases_ge5_records"] == 39
    assert out["frozen_cohort"]["top_two_record_counts"] == 4418
    assert out["frozen_cohort"]["top_two_record_share"] == pytest.approx(4418 / 6076)


def test_paired_species_does_not_count_records_as_individual_trials():
    rows, _ = read_cases(DEFAULT_SOURCE)
    p = compare(rows, "area95", 1, boot_iterations=100)
    assert p["paired_species"] == 10
    assert p["paired_species_countries_each_side_ge2"] == 3
    assert (p["positive_species"], p["negative_species"], p["tied_species"]) == (7, 1, 2)
    assert p["mean_equal_species_difference"] == pytest.approx(0.4242, abs=0.0001)
    assert p["by_class"]["novel"]["cases"] == 32
    assert p["by_class"]["analogue"]["cases"] == 41
    thinned = compare(rows, "area95", 5, boot_iterations=100)
    assert thinned["paired_species"] == 2
    assert thinned["cases_retained"] == 39
    assert thinned["mean_equal_species_difference"] == pytest.approx(0.3381, abs=0.0001)
    assert {s["species"] for s in thinned["paired_species_details"]} == {"Mechanitis lysimnia", "Zizula hylax"}


def test_strict_coordinates_and_climate_thresholds_are_not_interchanged():
    rows, _ = read_cases(DEFAULT_SOURCE)
    strict = compare(rows, "area95", 1, strict=True, boot_iterations=100)
    assert strict["paired_species"] == 9
    assert strict["by_class"]["novel"]["points"] == 260
    assert strict["by_class"]["analogue"]["points"] == 2918
    area90 = compare(rows, "area90", 1, boot_iterations=100)
    assert area90["paired_species"] == 10
    assert area90["mean_equal_species_difference"] > .6
    equal = compare(rows, "equal95", 1, boot_iterations=100)
    assert equal["mean_equal_species_difference"] > .4


def test_reject_corrupt_aggregates_and_no_unsafe_zero_fill():
    rows, _ = read_cases(DEFAULT_SOURCE)
    bad = deepcopy(rows)
    bad[0]["strict_n"] = bad[0]["n"] + 1
    # Reproduce read_cases validations using a temporary malformed frozen CSV.
    import csv
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "bad.csv"
        keys = list(rows[0].keys())
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(bad)
        with pytest.raises(ValueError, match="Invalid denominator"):
            read_cases(path, verify_frozen=False)
        with pytest.raises(ValueError, match="SHA256 mismatch"):
            read_cases(path)
    with pytest.raises(ValueError):
        compare(rows, "wrong-climate-product")
    with pytest.raises(ValueError):
        compare(rows, "area95", 0)
