"""Geometric fixtures for exclusion, sampling bounds and pruning safety."""
import importlib.util
from pathlib import Path
import pytest

# Dedicated analysis runs execute these fixtures; the lightweight paper test
# environment does not install the exploratory GIS stack.
pytest.importorskip("geopandas")


def module():
    path = Path(__file__).parents[1] / "scripts/measure_named_island_neighbor_gaps_20261029.py"
    assert path.exists(), "Named-island neighbor measurement must exist"
    spec = importlib.util.spec_from_file_location("neighbor_gaps", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_self_is_excluded_and_nearest_other_unit_is_retained():
    m = module()
    samples = {0: m.unit_xyz([[0, 0]]), 1: m.unit_xyz([[1, 0]]),
               2: m.unit_xyz([[10, 0]])}
    result = m.nearest_other(0, samples, 1.0)
    assert result["neighbor_index"] == 1
    assert 111 < result["upper_km"] < 112
    assert 110 < result["lower_km"] < 111
    assert result["candidates_evaluated"] == 1


def test_large_sample_cap_cannot_prune_a_close_coast_by_its_center():
    m = module()
    # Unit 2's first sample is far away, but its other coast sample is closest.
    samples = {0: m.unit_xyz([[0, 0]]), 1: m.unit_xyz([[2, 0]]),
               2: m.unit_xyz([[40, 0], [0.1, 0]])}
    result = m.nearest_other(0, samples, 1.0)
    assert result["neighbor_index"] == 2
    assert 11 < result["upper_km"] < 12
    assert result["lower_km"] <= result["upper_km"]
