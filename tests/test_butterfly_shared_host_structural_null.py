"""Synthetic structural-null tests; never pretend toy outputs are ecology."""
from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analyze_butterfly_shared_host_structural_null.py"
spec = importlib.util.spec_from_file_location("structural_null", SCRIPT)
assert spec and spec.loader
null = importlib.util.module_from_spec(spec)
spec.loader.exec_module(null)


def test_exact_shared_host_region_and_single_support_are_distinct() -> None:
    native_host = {"a": 1 << 0, "b": 1 << 1}
    introduced = {"a": 1 << 2, "b": 1 << 2}
    host_lists = [["a", "b"], ["a"]]
    native_species = [(1 << 0) | (1 << 1), 1 << 0]
    common = [(0, 1, ["a"], 1 << 0)]
    x = null.site_contribution_stats(
        host_lists, native_host, introduced, common, native_species
    )
    assert x["added_units"] == 2
    assert x["single_units"] == 1
    assert x["single_fraction"] == 0.5
    assert x["introduced_host_butterfly_region_incidences"] == 3
    assert x["introduced_incidences_in_native_supported_cells"] == 0
    assert x["introduced_incidences_in_newly_supported_cells"] == 3
    assert x["redundant_introduced_host_incidences_in_new_cells"] == 1
    assert x["mean_introduced_hosts_per_new_cell"] == 1.5
    assert x["novel_exact_shared_pair_regions"] == 1
    assert x["shared_pair_units_per_added_unit"] == 0.5
    assert x["added_cell_total_other_shared_butterflies"] == 2
    assert x["added_cell_mean_other_shared_butterflies"] == 1.0
    assert x["added_cell_fraction_with_any_shared_butterfly"] == 1.0
    assert x["added_cell_fraction_without_shared_butterfly"] == 0.0



def test_no_co_user_when_butterflies_gain_distinct_hosts_in_same_region() -> None:
    # Both butterflies acquire a plant resource in region 2, but no EXACT
    # plant is shared. The opportunity envelopes overlap, host use does not.
    native_host = {"a": 1 << 0, "b": 1 << 1}
    introduced = {"a": 1 << 2, "b": 1 << 2}
    result = null.site_contribution_stats(
        [["a"], ["b"]],
        native_host,
        introduced,
        [],  # no exact shared host pair
        [1 << 0, 1 << 1],
    )
    assert result["added_units"] == 2
    assert result["single_units"] == 2
    assert result["novel_exact_shared_pair_regions"] == 0
    assert result["added_cell_mean_other_shared_butterflies"] == 0
    assert result["added_cell_fraction_without_shared_butterfly"] == 1


def test_added_cell_exposures_can_differ_from_novel_shared_pair_regions() -> None:
    # Shared host was already accessible to both butterflies in region 0,
    # and an introduced host from the same pair reaches region 2 where
    # butterfly 0 has no native resources but butterfly 1 has native
    # resources through plant b. The pair's new shared-host region counts
    # once, the newly gained-cell peer exposure counts only butterfly 0.
    n = {"a": 1 << 0, "b": 1 << 2}
    c = {"a": 1 << 2, "b": 0}
    result = null.site_contribution_stats(
        [["a"], ["a", "b"]],
        n, c,
        [(0, 1, ["a"], 1 << 0)],
        [1 << 0, (1 << 0) | (1 << 2)],
    )
    assert result["added_units"] == 1
    assert result["novel_exact_shared_pair_regions"] == 1
    assert result["added_cell_total_other_shared_butterflies"] == 1
    assert result["added_cell_mean_other_shared_butterflies"] == 1
    assert result["added_cell_fraction_without_shared_butterfly"] == 0
    assert result["introduced_host_butterfly_region_incidences"] == 2
    assert result["introduced_incidences_in_native_supported_cells"] == 1
    assert result["introduced_incidences_in_newly_supported_cells"] == 1
    assert result["redundant_introduced_host_incidences_in_new_cells"] == 0

def test_swap_preserves_both_geographic_margins() -> None:
    native = {"a": {0}, "b": {1}, "c": {4}}
    intro_original = {"a": {2}, "b": {3}, "c": {2}}
    family = {"a": "F", "b": "F", "c": "F"}
    level1 = {0: "A", 1: "A", 2: "A", 3: "A", 4: "A"}
    intro = {h: set(cols) for h, cols in intro_original.items()}
    groups, cumulative, total = null.make_groups(native, intro, family, level1)
    assert total == 3
    assert cumulative[-1] == 3
    accepted = null.attempt_range_swaps(
        groups, cumulative, native, intro, 300, random.Random(73)
    )
    assert accepted > 0
    assert {h: len(v) for h, v in intro.items()} == {
        h: len(v) for h, v in intro_original.items()
    }
    for region in (2, 3):
        assert sum(region in u for u in intro.values()) == sum(
            region in u for u in intro_original.values()
        )
    assert all(not intro[h] & native[h] for h in intro)
    assert all((h, u) in [(x, y) for g in groups for x, y in g]
               for h, vals in intro.items() for u in vals)


def test_family_and_continent_blocks_are_not_crossed() -> None:
    native = {"a": set(), "b": set(), "c": set()}
    intro = {"a": {1}, "b": {2}, "c": {3}}
    groups, cumulative, total = null.make_groups(
        native, intro, {"a": "F", "b": "G", "c": "F"},
        {1: "A", 2: "A", 3: "B"}
    )
    assert total == 3
    assert groups == []
    assert cumulative == []
    assert null.attempt_range_swaps(groups, cumulative, native, intro, 10,
                                    random.Random(8)) == 0


def test_exclusion_rejects_spurious_introduced_native_overlap() -> None:
    with pytest.raises(ValueError, match="overlap"):
        null.make_groups({"a": {1}}, {"a": {1}}, {"a": "F"}, {1: "A"})


def test_percentile_interpolation_and_evaluability_flag() -> None:
    assert null.percentiles([1., 2., 3., 4., 5.]) == pytest.approx(
        (1.1, 3., 4.9)
    )
    base = {
        "species_hosts": [["a"], ["a"]],
        "native_bits": {"a": 1},
        "introduced_sets": {"a": {1}},
        "native_sets": {"a": {0}},
        "shared_pairs": [(0, 1, ["a"], 1)],
        "native_species": [1, 1],
        "groups": [],
        "cumulative": [],
        "total_introduced_edges": 1,
    }
    with pytest.raises(ValueError, match="empirical reconstruction drift"):
        null.run_null(base, reps=19, seed=1, burn_factor=1, step_factor=1)
