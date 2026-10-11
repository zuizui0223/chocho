"""Safety and invariance checks for the continent x plant-family geographic null."""
from __future__ import annotations

import random
import runpy
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NS = runpy.run_path(str(REPO / "scripts/analyze_butterfly_resource_redundancy_geonull.py"))


def test_double_swap_preserves_family_region_and_host_continent_margins():
    edges = {
        ("Fabaceae", "1"): [
            ("h1", "EUR1"), ("h2", "EUR2"),
            ("h3", "EUR1"), ("h4", "EUR2"),
        ],
        ("Poaceae", "1"): [("h5", "EUR1")],
    }
    introduced = {
        "h1": {"EUR1"}, "h2": {"EUR2"}, "h3": {"EUR1"}, "h4": {"EUR2"},
        "h5": {"EUR1"},
    }
    native = {"h1": {"EUR3"}, "h2": {"EUR3"}, "h3": {"EUR3"}}
    before = NS["degree_signature"](edges)
    # Each block in the same dictionary is backed by the original mutable list.
    swappable = NS["build_swappable_blocks"](edges)
    attempted, accepted = NS["swaps"](swappable, native, introduced, random.Random(42), 200)
    assert attempted > 0
    assert accepted > 0
    assert before == NS["degree_signature"](edges)
    assert all(not (introduced[h] & native.get(h, set())) for h in introduced)
    assert sum(map(len, introduced.values())) == 5
    assert len({(h, region) for h, regs in introduced.items() for region in regs}) == 5


def test_resource_singletons_only_count_cells_outside_native_union():
    species = ["butterfly_A", "butterfly_B"]
    hosts = {"butterfly_A": {"h1", "h2"}, "butterfly_B": {"h3"}}
    native = {"h1": {"NATIVE"}, "h2": set(), "h3": set()}
    introductions = {"h1": {"I1"}, "h2": {"I1", "I2"}, "h3": {"I3"}}
    group = {"both": species, "single": ["butterfly_B"]}
    summary, species_metrics = NS["geographic_resources"](
        species, hosts, native, introductions, group
    )
    assert species_metrics["butterfly_A"] == (1, 2, 1, 3, 0, 1)
    assert species_metrics["butterfly_B"] == (0, 1, 1, 1, 0, 0)
    assert summary["both"]["introduced_host_region_incidences"] == 4
    assert summary["both"]["introduced_incidences_already_covered_by_native_hosts"] == 0
    assert summary["both"]["duplicated_introduced_incidences_outside_native_resource"] == 1
    assert summary["both"]["added_units"] == 3
    assert summary["both"]["single_host_added_units"] == 2
    assert summary["both"]["single_host_share"] == 2 / 3
    assert summary["single"]["single_host_share"] == 1


def test_permutation_p_value_tail_matches_definition():
    observed = {"single_host_share": 0.8}
    draws = [{"single_host_share": v} for v in (0.2, 0.3, 0.4, 0.85)]
    result = NS["summarize_metric"](observed, draws, "single_host_share")
    # Monte Carlo +1 correction: one random value >= 0.8
    assert result["p_observed_ge_null"] == 2 / 5
    assert result["p_observed_le_null"] == 4 / 5
