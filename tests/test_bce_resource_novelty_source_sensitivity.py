"""Unit checks for source-inventory dependence of host-range novelty accounting."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/analyze_bce_resource_novelty_source_sensitivity.py"
spec = importlib.util.spec_from_file_location("bce_host_novelty", PATH)
assert spec and spec.loader
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_native_absorption_and_introduced_redundancy_are_distinct() -> None:
    native = {"a": {"N"}, "b": {"M"}}
    contemporary = {"a": {"N", "Q"}, "b": {"M", "N", "Q"}}
    data = audit.count_species_opportunity({"a", "b"}, native, contemporary)
    assert data["native"] == 2
    assert data["contemporary"] == 3
    assert data["H"] == 3
    assert data["A"] == 1  # new host b in a's native region N
    assert data["R"] == 1  # two introduced host plants support Q
    assert data["added"] == 1
    assert data["H"] == data["A"] + data["R"] + data["added"]


def test_extra_native_host_can_reclassify_old_added_region() -> None:
    base_n = {"a": {"N"}}
    base_c = {"a": {"N", "Q"}}
    augmented_n = {"a": {"N"}, "b": {"Q"}}
    augmented_c = {"a": {"N", "Q"}, "b": {"Q", "R"}}
    initial = audit.count_species_opportunity({"a"}, base_n, base_c)
    later = audit.count_species_opportunity({"a", "b"}, augmented_n, augmented_c)
    assert initial["added_regions"] == {"Q"}
    assert later["added_regions"] == {"R"}
    assert initial["added_regions"] - later["added_regions"] == {"Q"}
    assert later["native_regions"] - initial["native_regions"] == {"Q"}
    assert later["H"] == 2
    assert later["A"] == 1
    assert later["R"] == 0
    assert later["added"] == 1


def test_zero_added_and_denominator_honesty() -> None:
    native = {"a": {"A"}}
    contemporary = {"a": {"A"}}
    x = audit.count_species_opportunity({"a"}, native, contemporary)
    y = audit.summation([x])
    assert y["native_cells"] == 1
    assert y["candidate_incidences"] == 0
    assert y["native_absorption_fraction"] is None
    assert y["distinct_new_region_per_candidate_fraction"] is None


def test_conflicting_wcvp_plant_ids_must_not_be_silently_merged() -> None:
    with pytest.raises(ValueError, match="WCVP range mismatch"):
        audit.combine_host_ranges(
            {"a": {"N"}}, {"a": {"Q"}}, {"a"}
        )
    out = audit.combine_host_ranges(
        {"a": {"N"}}, {"a": {"N"}, "b": {"Q"}}, {"a", "b"}
    )
    assert out["a"] == {"N"}
    assert out["b"] == {"Q"}


def test_frozen_expected_inventory_receipts_guard_against_cherry_picking() -> None:
    assert audit.SOURCE_SHA256 == (
        "f46d18703fe5cc8d97699af6f084f652a704e46b0f9e8c68e32dd70decc32c53"
    )
    assert audit.EXPECTED_BASE["candidate_incidences"] == 75346
    assert audit.EXPECTED_BASE["native_absorbed_incidences"] == 46387
    assert audit.EXPECTED_AUGMENTED["distinct_added_cells"] == 14967
    assert audit.EXPECTED_AUGMENTED["native_cells"] == 29852
