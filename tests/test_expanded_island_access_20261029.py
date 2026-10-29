import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("geopandas")
PATH = Path(__file__).resolve().parents[1] / "scripts/measure_expanded_island_access_20261029.py"
sys.path.insert(0, str(PATH.parent))


def module():
    assert PATH.exists(), "Expanded coast-access diagnostic is not implemented"
    spec = importlib.util.spec_from_file_location("expanded_access", PATH)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_candidate_area_rule_excludes_focal_and_small_nearer_landmass():
    m = module()
    areas = {0: 5000, 1: 9, 2: 1000, 3: 900}
    assert m.candidate_indices(areas, 0, 5) == [1, 2, 3]
    assert m.candidate_indices(areas, 0, 1000) == [2]
    points = {i: m.unit_xyz([[lon, 0]]) for i, lon in [(0, 0), (1, .01), (2, 5), (3, 2)]}
    caps = {i: m.spherical_cap(p, 1) for i, p in points.items()}
    result = m.measure_candidates(0, m.candidate_indices(areas, 0, 1000), points, caps, 1, {})
    assert result["neighbor_index"] == 2
    assert 554 < result["lower_km"] < result["upper_km"] < 557
    assert result["candidates_evaluated"] + result["candidates_pruned"] == 1


def test_candidate_pruning_keeps_far_center_with_close_coast():
    m = module()
    points = {0: m.unit_xyz([[0, 0]]), 1: m.unit_xyz([[10, 0], [.1, 0]]),
              2: m.unit_xyz([[.5, 0]])}
    caps = {i: m.spherical_cap(p, 1) for i, p in points.items()}
    result = m.measure_candidates(0, [1, 2], points, caps, 1, {})
    assert result["neighbor_index"] == 1
    assert result["upper_km"] == pytest.approx(11.119508, rel=1e-6)
