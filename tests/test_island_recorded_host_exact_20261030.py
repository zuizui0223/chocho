"""Focused synthetic and frozen-data checks of exact recorded-host co-recording."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/analyze_island_recorded_host_exact_20261030.py"
spec = importlib.util.spec_from_file_location("recorded_host_exact", PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_two_by_two_manual_enumeration():
    # Each row has one positive; each target island has one positive.
    # Exactly two binary tables: perfect matching (H=2) or cross (H=0).
    rows = [
        (1, ((0, 1, 3.0), (1, 0, 5.0))),
        (1, ((0, 0, 7.0), (1, 1, 11.0))),
    ]
    assert module.block_distribution(rows, (1, 1)) == {0: 1, 2: 1}
    moments = module.block_moments(rows, (1, 1))
    assert moments["count"] == 2
    assert moments["expected_H"] == 1.0
    assert moments["variance_H"] == 1.0
    score = module.summary_for_distribution({0: 1, 2: 1}, 2)
    assert score["inclusive_upper_tail_count"] == 1
    assert score["inclusive_upper_tail_probability"] == 0.5


def test_forced_zero_and_forced_positive():
    assert module.block_distribution([], ()) == {0: 1}
    rows = [(1, ((0, 1, 8.0),))]
    assert module.block_distribution(rows, (1,)) == {1: 1}
    assert module.block_moments(rows, (1,))["variance_H"] == 0.0
    assert module.convolve({0: 2, 1: 3}, {1: 2, 2: 1}) == {
        1: 4, 2: 8, 3: 3
    }


def test_distance_tilt_synthetic():
    # When lower-distance tables also have H=2, a proximity tilt raises P(H>=2).
    scores = [("tiny", [(2, 1.0), (0, 3.0)])]
    at_zero = module.weighted_profile(scores, 0.0, 2)
    stronger = module.weighted_profile(scores, 0.5, 2)
    assert math.isclose(at_zero["probability_H_at_least_observed"], 0.5, abs_tol=1e-12)
    assert stronger["probability_H_at_least_observed"] > 0.5
    assert stronger["null_expected_distance_score_D"] < at_zero["null_expected_distance_score_D"]
    assert math.isclose(stronger["normalized_score_probability_sum"], 1.0, abs_tol=1e-12)


def test_frozen_reproducible_fiber_and_no_drift():
    # Anonymous-stage intentionally excludes exploratory provenance/data.
    # Only that explicitly documented export stage can skip this data test.
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Original exploratory input is intentionally absent from anonymous GEB snapshot")
    result = module.run(ROOT)
    s = result["host_recording_statistic"]
    assert result["original_candidate_cells"] == 2385
    assert result["original_later_positive_cells"] == 66
    assert s["observed_recorded_host_cells"] == 20
    assert s["inclusive_upper_tail_count"] == 6487980
    assert s["exact_fiber_count"] == 147117600
    assert s["score_support_min"] == 12
    assert s["score_support_max"] == 22
    assert math.isclose(s["null_expected_recorded_host_cells"], 17.06018423356553, abs_tol=1e-10)
    assert math.isclose(s["inclusive_upper_tail_probability"], 0.04410063785706129, abs_tol=1e-12)
    assert len(result["blocks"]) == 12
    assert result["informative_host_blocks"] == 3
    assert sum(result["score_distribution"].values()) == 147117600
    assert math.isclose(result["distance_score_companion"]["null_H_D_correlation"],
                        -0.1488477943734855, abs_tol=1e-9)
    loo = result["leave_one_informative_block_out"]
    assert math.isclose(loo["longitude_1:tropical"]["inclusive_upper_tail_probability"],
                        0.06428571428571428, abs_tol=1e-12)
    assert math.isclose(loo["longitude_3:north_extratropical"]["inclusive_upper_tail_probability"],
                        0.21338738532983137, abs_tol=1e-12)
    assert math.isclose(loo["longitude_3:tropical"]["inclusive_upper_tail_probability"],
                        0.06644276415602213, abs_tol=1e-12)
    tilt = result["distance_tilt_sensitivity"]
    assert tilt["status"] == "distance_moment_matched"
    assert tilt["not_a_composite_null_p_value"] is True
    assert math.isclose(tilt["fitted_lambda"], 0.7611080093702449, abs_tol=1e-7)
    assert math.isclose(tilt["fitted_lambda_profile"]["null_expected_recorded_host_score_H"],
                        18.049327866145973, abs_tol=1e-7)
    assert math.isclose(tilt["fitted_lambda_profile"]["probability_H_at_least_observed"],
                        0.1305459571904408, abs_tol=1e-7)
    assert math.isclose(tilt["uniform"]["probability_H_at_least_observed"],
                        s["inclusive_upper_tail_probability"], abs_tol=1e-12)
    assert result["not_biological_host_absence"] is True
    assert result["not_confirmatory"] is True
