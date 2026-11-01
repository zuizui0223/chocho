"""Prevent optimistic island-pair pseudo-replication and retrospective time claims."""
import importlib.util
import json
import math
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "provenance/reviewer_defenses/results/island_observational_robustness_20261031"
SCRIPT = DIR / "loio_pair_baseline.py"
RECEIPT = DIR / "loio_pair_baseline_lambda1.json"


def test_held_out_island_baseline_replays():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory source not part of anonymous GEB reviewer bundle")
    spec = importlib.util.spec_from_file_location("loio_butterflies", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    actual = m.evaluate(lam=1.0)
    frozen = json.loads(RECEIPT.read_text(encoding="utf8"))
    for key in ("status", "islands", "pairs", "folds",
                "training_pairs_per_fold", "evaluation_pairs_per_fold",
                "feature_groups", "limitations"):
        assert actual[key] == frozen[key], key
    assert "new" not in " ".join(sum(actual["feature_groups"].values(), []))
    assert "2001-2024" in actual["data"]
    for m_name, stats in frozen["model_stats"].items():
        for k, value in stats.items():
            assert math.isclose(actual["model_stats"][m_name][k], value,
                                rel_tol=0, abs_tol=1e-10), (m_name, k)
    for name, stats in frozen["comparisons"].items():
        for k, value in stats.items():
            v = actual["comparisons"][name][k]
            if isinstance(value, (int, float)):
                assert math.isclose(v, value, rel_tol=0, abs_tol=1e-10), (name, k)
            else:
                assert len(v) == len(value)
                for a, b in zip(v, value):
                    assert math.isclose(a, b, rel_tol=0, abs_tol=1e-10)
    assert actual["comparisons"]["M1 vs M0"]["islands_improved"] == 30
    assert actual["comparisons"]["E1 vs E0"]["islands_improved"] == 29


def test_held_out_island_regularization_direction_stable():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory source not part of anonymous GEB reviewer bundle")
    spec = importlib.util.spec_from_file_location("loio_butterflies2", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for lam in (0, 10):
        o = m.evaluate(lam=lam)
        assert 1.0 < o["comparisons"]["M1 vs M0"]["mse_reduction_pct"] < 3.0
        assert 1.0 < o["comparisons"]["E1 vs E0"]["mse_reduction_pct"] < 2.0
        assert o["folds"] == o["islands"] == 51
