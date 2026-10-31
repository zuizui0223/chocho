"""Validate published species-level wing / host / port gap diagnostics.

Tests use existing public source aggregates; no restricted HOSTS edges or GBIF
raw coordinates are needed. These tests do not validate the causal hypothesis.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031"
SPEC = importlib.util.spec_from_file_location("butterfly_maritime_trait_gap", D / "analyze_trait_gap.py")
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_rank_ties_and_partial_rank_well_defined():
    assert MODULE.rank([2.0, 1.0, 2.0]) == [2.5, 1.0, 2.5]
    assert math.isclose(MODULE.spearman([1, 2, 3], [3, 2, 1]), -1)
    assert MODULE.spearman([1, 1, 1], [1, 2, 3]) is None


def test_all225_and_both_period_scope_are_retained():
    matched, all225 = MODULE.make_screen()
    assert len(all225) == 225
    assert len({r["taxon_key"] for r in all225}) == 225
    assert len(matched) == 193
    assert len([r for r in all225 if r["disposition"] == "MATCHED"]) == 193
    assert sum(r["n"] for r in matched) == 3136
    assert sum(r["far500"] for r in matched) == 164
    assert sum(r["far500"] > 0 for r in matched) == 86
    assert sum(r["port"] > 0 for r in matched) == 149
    assert all(0 <= r["far500"] <= r["far100"] <= r["n"] for r in matched)


def test_public_snapshot_recomputed_without_statistical_refitting():
    matched, _ = MODULE.make_screen()
    v = MODULE.report(matched, n_boot=0)
    snapshot = json.loads((D / "trait_gap_results.json").read_text(encoding="utf-8"))
    assert v["source_git_blob_pins"] == snapshot["source_git_blob_pins"]
    for key in ("n_candidate_species", "n_matched_species", "n_index_events",
                "n_record_gaps_500km", "n_species_with_gap_500km",
                "n_species_with_port_candidate_host_index"):
        assert v[key] == snapshot[key]
    for k, d in snapshot["correlations"].items():
        assert math.isclose(v["correlations"][k]["point"], d["point"],
                            rel_tol=0, abs_tol=1e-9)
    for threshold, sens in snapshot["minimum_event_sensitivity"].items():
        assert v["minimum_event_sensitivity"][threshold]["species"] == sens["species"]
        for k, a in sens["correlations"].items():
            b = v["minimum_event_sensitivity"][threshold]["correlations"][k]
            if a is None:
                assert b is None
            else:
                assert math.isclose(a, b, rel_tol=0, abs_tol=1e-9)


def test_deterministic_800_species_resamples_match_frozen_intervals():
    matched, _ = MODULE.make_screen()
    recomputed = MODULE.report(matched, n_boot=800)
    archived = json.loads((D / "trait_gap_results.json").read_text(encoding="utf-8"))
    assert recomputed["bootstrap_draws"] == archived["bootstrap_draws"] == 800
    for k, published in archived["correlations"].items():
        recalculated = recomputed["correlations"][k]
        assert recalculated["bootstrap_n"] == published["bootstrap_n"] == 800
        for expected, actual in zip(published["species_resampled_95"],
                                    recalculated["species_resampled_95"]):
            assert math.isclose(expected, actual, rel_tol=0, abs_tol=1e-8), k
