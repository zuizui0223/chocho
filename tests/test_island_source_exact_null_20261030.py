"""Exact labeled-table law, independent unranking and numerical-bound checks."""
import csv
from functools import lru_cache
import gzip
import importlib.util
import itertools
import math
from pathlib import Path
import sys

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/analyze_island_source_exact_null_20261030.py"


@lru_cache(None)
def module():
    assert PATH.exists(), "Exact conditional-table analysis is not implemented"
    spec = importlib.util.spec_from_file_location("island_source_exact_null", PATH)
    result = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = result
    spec.loader.exec_module(result)
    return result


def tiny_tables(shape):
    for bits in itertools.product((False, True), repeat=math.prod(shape)):
        yield np.array(bits, dtype=bool).reshape(shape)


def enumerate_fiber(mask, row_margins, column_margins):
    return [s for s in tiny_tables(mask.shape)
            if not np.any(s & ~mask)
            and np.array_equal(s.sum(1), row_margins)
            and np.array_equal(s.sum(0), column_margins)]


def test_every_three_by_three_mask_and_feasible_margin_fiber_has_exact_rank_bijection():
    m = module()
    states = list(tiny_tables((3, 3)))
    margins = [(tuple(s.sum(1)), tuple(s.sum(0))) for s in states]
    checked_tables = 0
    for mask in states:
        groups = {}
        for state, key in zip(states, margins):
            if not np.any(state & ~mask):
                groups.setdefault(key, []).append(state)
        for (rows, columns), expected in groups.items():
            fiber = m.ExactBlockFiber(mask, rows, columns)
            assert fiber.count == len(expected)
            actual = [fiber.unrank(k) for k in range(fiber.count)]
            assert {s.tobytes() for s in actual} == {s.tobytes() for s in expected}
            assert len({s.tobytes() for s in actual}) == fiber.count
            assert np.array_equal(fiber.occupancy_counts(), np.sum(expected, axis=0))
            checked_tables += len(expected)
    assert checked_tables == 3 ** 9


def test_structural_zero_six_cycle_and_labeled_unequal_margins():
    m = module()
    mask = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=bool)
    fiber = m.ExactBlockFiber(mask, [1, 1, 1], [1, 1, 1])
    assert fiber.count == 2
    assert np.array_equal(fiber.occupancy_counts(), mask.astype(int))
    unequal = np.array([[1, 1, 1], [1, 1, 0], [0, 1, 1]], dtype=bool)
    expected = enumerate_fiber(unequal, [2, 1, 1], [2, 1, 1])
    actual = m.ExactBlockFiber(unequal, [2, 1, 1], [2, 1, 1])
    assert actual.count == len(expected) == 2
    assert {actual.unrank(k).tobytes() for k in range(actual.count)} == {s.tobytes() for s in expected}


@pytest.mark.parametrize("shape", [(0, 0), (0, 3), (3, 0), (2, 3)])
def test_empty_zero_margin_fiber_contains_one_empty_table(shape):
    m = module()
    fiber = m.ExactBlockFiber(np.ones(shape, dtype=bool), [0] * shape[0], [0] * shape[1])
    assert fiber.count == 1
    assert fiber.unrank(0).shape == shape
    assert not fiber.unrank(0).any()
    assert not fiber.occupancy_counts().any()
    with pytest.raises(ValueError, match="rank"):
        fiber.unrank(1)


def test_infeasible_hall_constraint_is_zero_and_cannot_be_sampled():
    m = module()
    mask = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]], dtype=bool)
    fiber = m.ExactBlockFiber(mask, [1, 1, 1, 1], [2, 1, 1, 0])
    assert fiber.count == 0
    with pytest.raises(ValueError, match="rank"):
        fiber.unrank(0)


@pytest.mark.parametrize("rows,columns", [([1, 0], [0, 0]), ([-1, 1], [0, 0]), ([.5, .5], [1, 0])])
def test_invalid_margins_are_rejected_without_coercing_or_rebalancing(rows, columns):
    with pytest.raises(ValueError):
        module().ExactBlockFiber(np.ones((2, 2), dtype=bool), rows, columns)


def test_count_caps_fail_instead_of_returning_truncated_completion_counts():
    m = module()
    mask = np.ones((3, 3), dtype=bool)
    with pytest.raises(m.DPResourceLimitError, match="states"):
        m.ExactBlockFiber(mask, [1] * 3, [1] * 3, max_states=1)
    with pytest.raises(m.DPResourceLimitError, match="transitions"):
        m.ExactBlockFiber(mask, [1] * 3, [1] * 3, max_transitions=1)


def test_centered_moments_and_occupancies_match_independent_exhaustive_scores():
    m = module()
    mask = np.array([[1, 1, 1], [1, 1, 0], [0, 1, 1]], dtype=bool)
    fiber = m.ExactBlockFiber(mask, [1, 1, 1], [1, 1, 1])
    values = np.array([[1, 5, 8], [2, 6, 0], [0, 3, 9]], dtype=float)
    expected = [math.fsum(values[s]) for s in enumerate_fiber(mask, [1] * 3, [1] * 3)]
    mean = math.fsum(expected) / len(expected)
    variance = math.fsum((s - mean) ** 2 for s in expected) / len(expected)
    observed = fiber.moments(values)
    assert observed[0] == pytest.approx(mean, abs=1e-12)
    assert observed[1] == pytest.approx(variance, abs=1e-12)
    assert math.fsum((fiber.occupancy_counts() / fiber.count * values).ravel()) == pytest.approx(mean)
    before = [fiber.unrank(k).tobytes() for k in range(fiber.count)]
    changed = fiber.moments(values ** 2)
    assert changed[0] != observed[0]
    assert [fiber.unrank(k).tobytes() for k in range(fiber.count)] == before
    shifted = fiber.moments(values + 1e10)
    assert shifted[1] == pytest.approx(variance, abs=1e-5)


def test_product_rank_bijection_and_integer_rank_sampling_with_replacement():
    m = module()
    first = m.ExactBlockFiber(np.ones((2, 2), dtype=bool), [1, 1], [1, 1])
    second = m.ExactBlockFiber(np.ones((2, 3), dtype=bool), [1, 1], [1, 1, 0])
    expected = set(itertools.product(
        [first.unrank(k).tobytes() for k in range(first.count)],
        [second.unrank(k).tobytes() for k in range(second.count)]))
    actual = {tuple(s.tobytes() for s in m.unrank_product([first, second], k)) for k in range(4)}
    assert actual == expected
    draws = list(m.draw_global_ranks([2, 2], seed=7723, samples=20))
    assert draws == list(m.draw_global_ranks([2, 2], seed=7723, samples=20))
    assert len(draws) == 20 and len(set(draws)) < 20
    large = list(m.draw_global_ranks([2 ** 60 + 3, 17], seed=9723, samples=20))
    assert all(type(k) is int and 0 <= k < (2 ** 60 + 3) * 17 for k in large)
    assert any(k > 2 ** 53 for k in large)


def test_per_table_validation_catches_structural_zeros_and_margin_changes():
    m = module()
    state = np.eye(3, dtype=bool)
    m.validate_table(state, state, [1] * 3, [1] * 3)
    with pytest.raises(ValueError, match="mask"):
        m.validate_table(np.roll(state, 1, axis=0), state, [1] * 3, [1] * 3)
    with pytest.raises(ValueError, match="margin"):
        m.validate_table(np.zeros_like(state), state, [1] * 3, [1] * 3)


def test_inclusive_tail_and_exact_binomial_upper_bound_do_not_treat_zero_hits_as_zero_probability():
    m = module()
    zero = m.tail_summary(0, 100_000)
    assert zero["one_sided_95_upper"] == pytest.approx(2.9956874019427947e-5, rel=1e-12)
    assert zero["simultaneous_three_scenario_95_upper"] > zero["one_sided_95_upper"]
    assert zero["qualified_bound"]
    assert m.tail_summary(10, 10)["one_sided_95_upper"] == 1
    from scipy.stats import binom
    upper = m.tail_summary(3, 10)["one_sided_95_upper"]
    assert binom.cdf(3, 10, upper) == pytest.approx(.05, abs=1e-12)
    assert m.in_closer_tail(1 + .5e-12, 1)
    assert not m.in_closer_tail(1 + 2e-12, 1)
    for count, total in [(-1, 10), (11, 10), (0, 0)]:
        with pytest.raises(ValueError):
            m.tail_summary(count, total)


def test_deterministic_gzip_preserves_every_rank_as_decimal_integer(tmp_path):
    m = module()
    rows = [{"rank": 2 ** 91 + 7, "T": "7.123456789012345"}, {"rank": 0, "T": "8.5"}]
    first, second = tmp_path / "first.csv.gz", tmp_path / "second.csv.gz"
    m.write_gzip_csv(first, rows, ["rank", "T"])
    m.write_gzip_csv(second, rows, ["rank", "T"])
    assert first.read_bytes() == second.read_bytes()
    with gzip.open(first, "rt", newline="") as handle:
        saved = list(csv.DictReader(handle))
    assert [int(r["rank"]) for r in saved] == [2 ** 91 + 7, 0]
    assert saved[0]["T"] == rows[0]["T"]


def test_input_hash_gate_rejects_tampering_before_loading_data(tmp_path):
    m = module()
    source = tmp_path / "pinned.csv"
    source.write_text("original\n")
    pins = {"pinned.csv": m.sha(source.read_bytes())}
    m.verify_pinned_files(tmp_path, pins)
    source.write_text("changed\n")
    with pytest.raises(ValueError, match="Pinned input changed"):
        m.verify_pinned_files(tmp_path, pins)


def test_pinned_loader_and_block_counts_reproduce_frozen_labeled_inputs():
    m = module()
    _, risks = m.load_inputs()
    expected = {
        "primary": (2385, 66, [1, 6870383478, 91827, 5541566504]),
        "source_records_ge3": (2197, 62, [1, 1403180478, 18727, 2678695685]),
        "late_records_ge2": (2385, 35, [1, 4692, 1, 952792]),
    }
    for name, (cells, positives, counts) in expected.items():
        risk = risks[name]
        assert int(risk["mask"].sum()) == cells
        assert int(risk["state"].sum()) == positives
        blocks = m.build_blocks(risk)
        assert [b["fiber"].count for b in blocks] == counts
        for rank in [0, math.prod(counts) - 1]:
            state = m.table_from_global_rank(blocks, rank, risk["mask"].shape)
            m.validate_table(state, risk["mask"], risk["row_margins"], risk["column_margins"])


def test_sample_runner_preserves_forced_contribution_denominator_every_rank_and_every_margin(tmp_path):
    m = module()
    mask = np.zeros((4, 4), dtype=bool)
    mask[:3, :3], mask[3, 3] = True, True
    values = np.array([[1, 2, 4, 0], [8, 16, 32, 0], [64, 128, 256, 0], [0, 0, 0, 1000]], dtype=float)
    state = np.eye(4, dtype=bool)
    labels = [(f"s{i}", "longitude_1" if i < 3 else "longitude_0") for i in range(4)]
    cells = [{"row_index": i, "column_index": j, "longitude_block": labels[i][1],
              "species_key": f"s{i}", "original_name": f"Species {i}", "target_island": f"island{j}",
              "late_positive": int(state[i, j])} for i, j in zip(*np.where(mask))]
    risk = {"analysis": "toy", "mask": mask, "state": state, "values": values,
            "row_margins": state.sum(1), "column_margins": state.sum(0),
            "row_labels": labels, "column_labels": [f"island{j}" for j in range(4)], "cells": cells,
            "observed_T": 1273 / 4}
    blocks = m.build_blocks(risk)
    summary, occupancy_rows = m.summarize_fiber(risk, blocks)
    assert summary["fiber_cardinality_decimal"] == "6"
    assert summary["null_mean_T"] == pytest.approx((1000 + 511 / 3) / 4)
    assert summary["support"]["forced_positive_cells"] == 1
    assert summary["support"]["variable_positive_cells"] == 3
    assert len(occupancy_rows) == 10
    saved = tmp_path / "toy_samples.csv.gz"
    result = m.sample_scenario(risk, blocks, seed=1729, samples=40, path=saved)
    assert result["stored_samples"] == result["stored_margin_checks"] == 40
    with gzip.open(saved, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 40
    assert len({r["global_rank"] for r in rows}) <= 6
    tail_count = 0
    for row in rows:
        table = m.table_from_global_rank(blocks, int(row["global_rank"]), mask.shape)
        m.validate_table(table, mask, state.sum(1), state.sum(0))
        expected = math.fsum(values[table]) / 4
        assert float(row["mean_log1p_source_gap_km"]) == expected
        assert table[3, 3]
        tail_count += expected <= risk["observed_T"] + 1e-12
    assert result["lower_tail_count_including_ties"] == tail_count
