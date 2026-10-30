"""Consequential mask, geometry, fixed-margin kernel and Monte Carlo checks."""
import importlib.util
import itertools
import math
from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("scipy")
PATH = Path(__file__).resolve().parents[1] / "scripts/analyze_island_recorded_sources_20261029.py"
sys.path.insert(0, str(PATH.parent))


def module():
    assert PATH.exists(), "Recorded-source diagnostic is not implemented"
    spec = importlib.util.spec_from_file_location("recorded_sources", PATH)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def enumerate_fiber(mask, row_margins, column_margins):
    """Independent exhaustive enumeration, used only on tiny test fixtures."""
    cells = list(zip(*np.where(mask)))
    result = []
    for chosen in itertools.combinations(cells, sum(row_margins)):
        state = np.zeros(mask.shape, dtype=bool)
        for row, column in chosen:
            state[row, column] = True
        if np.array_equal(state.sum(1), row_margins) and np.array_equal(state.sum(0), column_margins):
            result.append(state)
    return result


def exact_kernel(m, mask, states):
    """Enumerate ALL ordered positive-edge draws, including invalid proposals."""
    keys = {s.tobytes(): i for i, s in enumerate(states)}
    transition = np.zeros((len(states), len(states)))
    for i, state in enumerate(states):
        edges = list(zip(*np.where(state)))
        max_length = min(np.count_nonzero(state.sum(0)), np.count_nonzero(state.sum(1)), len(edges))
        lengths = m.cycle_length_probabilities(max_length)
        transition[i, i] += .5
        if not lengths:
            transition[i, i] += .5
        for length, probability in lengths.items():
            ordered_draws = list(itertools.permutations(edges, length))
            for chosen in ordered_draws:
                proposal = state.copy()
                m.apply_cycle(proposal, mask, chosen)
                transition[i, keys[proposal.tobytes()]] += .5 * probability / len(ordered_draws)
    return transition


def test_structural_zero_fiber_needs_three_cycle_and_has_symmetric_lazy_kernel():
    m = module()
    mask = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=bool)
    states = enumerate_fiber(mask, [1, 1, 1], [1, 1, 1])
    assert len(states) == 2
    transition = exact_kernel(m, mask, states)
    # k=2 is rejected. Of six ordered triples, three have the allowed orientation.
    # lazy .5 * k3 weight .15/(.80+.15) * allowed orientation .5 = 3/76.
    assert transition[0, 1] == pytest.approx(3 / 76)
    assert np.allclose(transition, transition.T)
    assert np.allclose(transition.sum(1), 1)
    assert np.allclose(np.ones(2) / 2 @ transition, np.ones(2) / 2)
    assert np.all(np.diag(transition) >= .5)
    assert all(m.fiber_has_move(s, mask) for s in states)


@pytest.mark.parametrize("row_margins,column_margins", [([1, 1, 1], [1, 1, 1]), ([2, 1, 1], [2, 1, 1])])
def test_all_ordered_draws_preserve_uniform_stationarity_on_irregular_mask(row_margins, column_margins):
    m = module()
    mask = np.array([[1, 1, 1], [1, 1, 1], [1, 0, 1]], dtype=bool)
    states = enumerate_fiber(mask, row_margins, column_margins)
    assert len(states) >= 4
    transition = exact_kernel(m, mask, states)
    assert np.allclose(transition, transition.T, atol=1e-13)
    assert np.allclose(transition.sum(1), 1)
    reachable = {0}
    while True:
        expanded = reachable | {j for i in reachable for j in range(len(states)) if transition[i, j] > 0}
        if expanded == reachable:
            break
        reachable = expanded
    assert len(reachable) == len(states)


def test_cycle_rejects_repeated_rows_without_changing_matrix_or_resampling():
    m = module()
    state = np.array([[1, 1, 0], [0, 0, 1]], dtype=bool)
    before = state.copy()
    assert m.apply_cycle(state, np.ones_like(state), [(0, 0), (0, 1)]) == "repeated_row_or_column"
    assert np.array_equal(state, before)
    assert set(m.cycle_length_probabilities(7)) == {2, 3, 4, 5, 6, 7}
    assert sum(m.cycle_length_probabilities(7).values()) == pytest.approx(1)


def test_attempted_time_chain_includes_invalid_cycles_and_lazy_holds():
    m = module()
    mask = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=bool)
    state = np.eye(3, dtype=bool)
    values = np.array([[1, 3, 0], [0, 1, 3], [3, 0, 1]], dtype=float)
    samples, report = m.sample_chain(state, mask, values, seed=9017, burn_in=100, samples=2000, thin=10)
    assert report["attempted_steps"] == 20100
    assert sum(report["outcomes"].values()) == 20100
    assert .03 < report["accepted_moves"] / 20100 < .05
    assert set(samples) == {1.0, 3.0}
    assert abs(float(np.mean(samples)) - 2) < .12
    assert report["stored_margin_checks"] == 2000


def test_irregular_mask_chain_agrees_with_exhaustive_uniform_statistic():
    m = module()
    mask = np.array([[1, 1, 1], [1, 1, 1], [1, 0, 1]], dtype=bool)
    states = enumerate_fiber(mask, [1, 1, 1], [1, 1, 1])
    values = np.array([[0, 2, 5], [7, 3, 1], [6, 4, 8]], dtype=float)
    expected = np.mean([np.mean(values[s]) for s in states])
    samples, _ = m.sample_chain(states[0], mask, values, seed=3721, burn_in=1000, samples=6000, thin=10)
    assert np.mean(samples) == pytest.approx(expected, abs=.08)


def risk_fixture():
    roster = [{"island_key": key, "representative_lon": lon, "representative_lat": 0}
              for key, lon in [("a", -20), ("b", -10), ("c", -5), ("d", 120)]]
    counts = {(r["island_key"], w): {"x": 1, "y": 1, "z": 1}
              for r in roster for w in ("2010_2017", "2018_2025")}
    counts["b", "2010_2017"]["s"] = 2
    counts["a", "2018_2025"]["s"] = 1
    counts["b", "2018_2025"]["s"] = 2
    counts["d", "2018_2025"]["w"] = 1
    distances = {tuple(sorted((a, b))): (10., 11.) for a, b in itertools.combinations("abcd", 2)}
    return roster, counts, {"s": "Original species"}, distances


def test_risk_mask_uses_fixed_late_four_taxa_and_excludes_baseline_presence():
    m = module()
    risk = m.build_risk_set(*risk_fixture())
    rows = risk["cells"]
    assert {r["target_island"] for r in rows} == {"a", "d"}
    assert sum(r["late_positive"] for r in rows) == 1
    assert {r["longitude_block"] for r in rows} == {"longitude_1", "longitude_3"}
    assert all(r["source_island_count"] == 1 for r in rows)
    assert all(r["sampled_nearest_source_island"] == "b" for r in rows)
    assert all(r["nearest_source_lower_km"] == 10 for r in rows)
    assert all(r["nearest_source_upper_km"] == 11 for r in rows)
    for row in rows:
        ri, ci = row["row_index"], row["column_index"]
        assert risk["mask"][ri, ci]
        # Every other species/block row is structurally forbidden at this target.
        other_rows = [i for i, label in enumerate(risk["row_labels"]) if label[1] != row["longitude_block"]]
        assert not risk["mask"][other_rows, ci].any()


def test_declared_source_and_late_thresholds_have_separate_risk_and_outcome_roles():
    m = module()
    fixture = risk_fixture()
    source = m.build_risk_set(*fixture, source_min_records=3)
    assert source["cells"] == []
    late = m.build_risk_set(*fixture, late_min_records=2)
    assert len(late["cells"]) == 2
    assert not late["state"].any()


def test_geometry_pair_uses_coasts_not_representative_centers():
    pytest.importorskip("geopandas")
    m = module()
    from measure_island_port_geography_20261029 import unit_xyz
    samples = {"a": unit_xyz([[0, 0], [10, 0]]), "b": unit_xyz([[10.1, 0]])}
    measured = m.measure_pair_distances(samples, step_km=1)
    lo, hi = measured["a", "b"]
    assert hi == pytest.approx(11.119508, rel=1e-6)
    assert hi - lo == pytest.approx(1)


def test_deterministic_fiber_and_constant_statistic_never_claim_sampling_support():
    m = module()
    state = np.eye(3, dtype=bool)
    assert not m.fiber_has_move(state, state.copy())
    diagnostic = m.assess_samples(np.ones((4, 2000)), observed=1.)
    assert not diagnostic["qualified_tail"]
    assert diagnostic["qualified_tail_probability"] is None
    assert "degenerate_statistic" in diagnostic["gate_failures"]
    zero = m.run_null(np.zeros((2, 2), dtype=bool), np.ones((2, 2), dtype=bool), np.ones((2, 2)))
    assert zero["summary"]["status"] == "NO_POSITIVE_OUTCOMES"
    assert zero["samples"] == []


def test_chain_diagnostics_detect_nonmixing_and_gate_tiny_monte_carlo_tails():
    m = module()
    rng = np.random.default_rng(19)
    independent = rng.normal(size=(4, 4000))
    ok = m.assess_samples(independent, observed=0.)
    assert ok["split_rhat"] < 1.05
    assert ok["pooled_ess"] >= 1000
    assert ok["tail_indicator_ess"] >= 1000
    assert ok["qualified_tail"]
    tiny = m.assess_samples(independent, observed=-100.)
    assert not tiny["qualified_tail"]
    assert "fewer_than_20_samples_in_one_tail" in tiny["gate_failures"]
    drifted = independent + np.arange(4)[:, None] * 3
    bad = m.assess_samples(drifted, observed=0.)
    assert bad["split_rhat"] >= 1.05
    assert not bad["qualified_tail"]


def test_complete_zero_snapshot_is_distinct_from_missing_or_duplicate_queries():
    m = module()
    roster = [{"island_key": "a"}]
    zero = {"status": "COMPLETE", "count": 0, "first_page_recheck_passed": True,
            "facets": {f: {} for f in ("SPECIES_KEY", "DATASET_KEY", "YEAR", "BASIS_OF_RECORD")},
            "unassigned_facet_records": {f: 0 for f in ("SPECIES_KEY", "DATASET_KEY", "YEAR", "BASIS_OF_RECORD")}}
    queries = [{**zero, "island_key": "a", "window": w} for w in ("2010_2017", "2018_2025")]
    assert m.validate_occurrence_snapshot(roster, {"results": queries}) == {
        ("a", "2010_2017"): {}, ("a", "2018_2025"): {}}
    for changed in ([queries[0]], queries + [queries[0]], [{**queries[0], "status": "MISSING_QUERY_FAILED"}, queries[1]]):
        with pytest.raises(ValueError):
            m.validate_occurrence_snapshot(roster, {"results": changed})
