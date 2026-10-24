"""Source-agnostic numerical and bitset tests for the gain-support plant null."""

import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_butterfly_gain_support_plant_null_20261024 import (
    individual_resource_states, family_region_margins, margin_assert,
    ess_estimate, regression_states,
)
from analyze_butterfly_shared_host_structural_null import (
    make_groups, attempt_range_swaps,
)


def test_single_support_duplicate_bits():
    # Species one: two newly opened cells, one has two distinct hosts.
    model = {
        "species": ["A a", "B b"],
        "species_hosts": [["a", "b"], ["a"]],
        "native_species": [1 << 0, 1 << 0],
    }
    introduced = {"a": {1}, "b": {1, 2}}
    out = individual_resource_states(model, introduced)
    assert out == [
        ("A a", 2, 1, 1),
        ("B b", 1, 1, 1),
    ]


def test_family_continent_swap_preserves_exact_margins():
    native = {"a": {2}, "b": {3}, "c": set()}
    introduced = {
        "a": {0, 1},
        "b": {0, 1},
        "c": {0, 1},
    }
    fam = {"a": "Fam", "b": "Fam", "c": "Fam"}
    level1 = {u: "1" for u in range(4)}
    groups, cumulative, E = make_groups(native, introduced, fam, level1)
    model = {"native_sets": native, "host_families": fam, "level1": level1}
    ref = family_region_margins(model, introduced)
    assert E == 6
    swaps = attempt_range_swaps(
        groups, cumulative, native, introduced, 300, random.Random(42)
    )
    assert swaps >= 0
    margin_assert(model, introduced, ref)
    for h, regions in introduced.items():
        assert regions.isdisjoint(native[h])


def test_ess_and_regression_source_guards():
    e = ess_estimate([math.sin(i*2.7) for i in range(199)])
    assert 0 < e["ess"] <= 199
    assert e["median"] is not None
    constant = ess_estimate([0.15]*199)
    assert constant["ess"] == 0

    # Model rank and cohort gate must never invent a small-sample beta.
    model = {
        "frozen_metrics": {
            "A a": {"native": 1, "host_count": 1, "host_family_count": 1, "adequate": True},
            "B b": {"native": 1, "host_count": 1, "host_family_count": 1, "adequate": True},
        }
    }
    try:
        regression_states(model, [
            ("A a",1,1,1),("B b",1,0,1)
        ])
    except ValueError as e:
        assert "full-rank" in str(e)
    else:
        raise AssertionError("Small singular panels may not return zero-valued effects")
