"""Only assumed variance and requested precision; NEVER experimental evidence."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

P = Path(__file__).resolve().parents[1] / "scripts/plan_kyoto_six_arm_precision_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("kyoto_prec", P)
assert SPEC and SPEC.loader
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def test_hypothetical_0point2_sd_and_0point1_halfwidth():
    r = m.scenario(0.2, 0.1)
    assert r["independent_complete_blocks_illustrative"] == 21
    assert r["minimum_whole_cages_for_primary_FOUR_arm"] == 84
    assert r["minimum_whole_cages_for_extended_SIX_arm"] == 126
    assert r["additional_cages_to_estimate_identity_arm_comparisons"] == 42
    assert r["calculated_conservative_halfwidth"] <= 0.1
    assert r["not_power_or_empirical_variance"] is True


def test_illustrative_precision_bound_monotone_and_extends_to_six_arms():
    for sd in (0.10, 0.20, 0.30, 0.40):
        for h in (0.1, 0.15, 0.2):
            v = m.scenario(sd, h)
            assert v["independent_complete_blocks_illustrative"] >= 10
            assert v["minimum_whole_cages_for_extended_SIX_arm"] == (
                v["independent_complete_blocks_illustrative"] * 6
            )
            assert v["calculated_conservative_halfwidth"] <= h + 1e-12
        assert m.scenario(sd, 0.1)["independent_complete_blocks_illustrative"] >= (
            m.scenario(sd, 0.2)["independent_complete_blocks_illustrative"]
        )


def test_no_invented_prior_sd_or_pilot_inference_claim():
    v = m.grid([0.2], [0.1, 0.2])
    assert v["status"] == "HYPOTHETICAL_PILOT_PLANNING_ONLY_ZERO_OBSERVED_VARIANCE"
    assert v["primary_estimand"] == "S_competitor_x_food"
    assert not v["authorizations_or_animal_work_completed"]
    assert len(v["scenarios"]) == 2
    assert all(c["not_power_or_empirical_variance"] for c in v["scenarios"])


@pytest.mark.parametrize("sd,h", [
    (0, 0.1), (-1, 0.1), (0.2, 0), (0.2, -0.1),
    (float("nan"), 0.1), (0.2, float("inf")), (100000, 1e-7),
])
def test_invalid_and_extreme_requests_fail_closed(sd, h):
    with pytest.raises(ValueError):
        m.scenario(sd, h)


def test_sample_size_not_power_and_10_is_only_illustrative_floor():
    r = m.scenario(0.05, 0.25)
    assert r["independent_complete_blocks_illustrative"] == 10
    assert math.isclose(r["conservative_critical_value_t_df9_for_B_gte10"], 2.2621571627409915)
