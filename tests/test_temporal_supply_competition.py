from __future__ import annotations

import math

import pytest

from butterfly_specialization_ecology.temporal_supply_competition import (
    evaluate, fraction_unmet, realized_food, reversal_interval, supply_regimes,
)

TOY = {"A": (0.75, 0.25), "B": (0.25, 0.75)}


def test_timing_treatments_preserve_total_supply_and_solo_access():
    arms = supply_regimes(TOY, 0.2)
    assert arms["early_shortage"] == pytest.approx((0.8, 1.0))
    assert arms["distributed_shortage"] == pytest.approx((0.9, 0.9))
    assert arms["late_shortage"] == pytest.approx((1.0, 0.8))
    for supply in arms.values():
        assert sum(supply) == pytest.approx(1.8)
        for demand in TOY.values():
            assert all(d <= r for d, r in zip(demand, supply))


def test_allocation_conserves_food_and_exposes_assumed_shortage():
    received = realized_food(TOY, (0.8, 1.0))
    assert received["A"] == pytest.approx((0.60, 0.25))
    assert received["B"] == pytest.approx((0.20, 0.75))
    assert sum(received[s][0] for s in received) == pytest.approx(0.8)
    assert fraction_unmet(TOY["A"], received["A"]) == pytest.approx(0.15)
    assert fraction_unmet(TOY["B"], received["B"]) == pytest.approx(0.05)


def test_assumed_calendar_demand_asynchrony_can_reverse_relative_harm():
    r = evaluate(TOY, 0.2)
    assert r["relative_disadvantage_reverses"] is True
    early = r["regimes"]["early_shortage"]
    even = r["regimes"]["distributed_shortage"]
    late = r["regimes"]["late_shortage"]
    assert early["relative_disadvantage_A_minus_B"] == pytest.approx(+0.30)
    assert even["relative_disadvantage_A_minus_B"] == pytest.approx(0.0)
    assert late["relative_disadvantage_A_minus_B"] == pytest.approx(-0.30)
    assert early["illustrative_mixed_emergence_probability"]["A"] == pytest.approx(
        0.8 * math.exp(-0.45)
    )
    assert late["illustrative_mixed_emergence_probability"]["B"] == pytest.approx(
        0.8 * math.exp(-0.45)
    )
    assert all(
        all(value == pytest.approx(0) for value in x["solo_unmet_fraction"].values())
        for x in r["regimes"].values()
    )


def test_no_shortage_and_same_relative_demand_time_preclude_reversal():
    assert evaluate(TOY, 0.0)["relative_disadvantage_reverses"] is False
    aligned = {"A": (0.6, 0.4), "B": (1.2, 0.8)}
    assert reversal_interval(aligned) is None
    assert evaluate(aligned, 0.2)["relative_disadvantage_reverses"] is False
    for arm in evaluate(aligned, 0.2)["regimes"].values():
        assert arm["relative_disadvantage_A_minus_B"] == pytest.approx(0.0)


def test_sensitivity_interval_and_counterexamples():
    bounds = reversal_interval(TOY)
    assert bounds is not None
    assert bounds["lower_exclusive"] == pytest.approx(1 / 3)
    assert bounds["upper_exclusive"] == pytest.approx(3)
    assert not evaluate(TOY, 0.2, shortage_sensitivity={"A": 0.2, "B": 3})[
        "relative_disadvantage_reverses"
    ]
    assert not evaluate(TOY, 0.2, shortage_sensitivity={"A": 12, "B": 3})[
        "relative_disadvantage_reverses"
    ]


def test_unequal_total_food_demand_does_not_by_itself_prevent_reversal():
    hypothetical = {"A": (1.95, 0.65), "B": (0.25, 0.75)}
    r = evaluate(hypothetical, 0.2)
    assert r["seasonal_supply_total"] == pytest.approx(3.4)
    assert r["relative_disadvantage_reverses"] is True


@pytest.mark.parametrize("demands,shortage", [
    (TOY, -0.1),
    (TOY, float("nan")),
    (TOY, 0.25),
    ({"A": (1.0, 0.0), "B": (1.0, 1.0)}, 0.1),
    ({"A": (1.0, 1.0), "C": (1.0, 1.0)}, 0.1),
])
def test_invalid_inputs_fail_closed(demands, shortage):
    with pytest.raises(ValueError):
        evaluate(demands, shortage)


def test_invalid_fitness_assumptions_fail_closed():
    with pytest.raises(ValueError):
        evaluate(TOY, 0.2, baseline_emergence={"A": 1.2, "B": 0.8})
    with pytest.raises(ValueError):
        evaluate(TOY, 0.2, shortage_sensitivity={"A": -1, "B": 3})
