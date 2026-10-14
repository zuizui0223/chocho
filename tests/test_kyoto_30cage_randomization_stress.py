"""Toy fixed-potential outcome stress tests; NO observed butterfly source data."""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
P = ROOT / "scripts/calibrate_kyoto_30cage_randomization.py"
SPEC = importlib.util.spec_from_file_location("kyoto_randomization_stress", P)
assert SPEC and SPEC.loader
model = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(model)

from randomize_kyoto_30cage_nearbalanced import assign


@pytest.mark.parametrize("scenario,expected", [
    ("global_no_effects", 0),
    ("nonzero_main_effects_without_interaction", 0),
    ("greenhouse_varying_main_effects_without_interaction", 0),
    ("opposing_greenhouse_interactions_zero_average", 0),
    ("hypothetical_positive_resource_rescue", .15),
])
def test_hypothetical_finite_population_truth_is_not_field_result(scenario, expected):
    p = model.make_fixed_potentials(scenario)
    assert len(p) == 30
    assert all(len(x) == 4 and min(x.values()) >= 0 and max(x.values()) <= 1
               for x in p.values())
    assert model.true_average_interaction(p) == pytest.approx(expected, abs=1e-10)


def test_nonzero_main_effects_cannot_be_mistaken_for_interaction_effect():
    inventory = model._inventory()
    a = assign(inventory, seed=20261014)
    null = model.evaluate_one_allocation(
        a, model.make_fixed_potentials("global_no_effects")
    )
    nuisance = model.evaluate_one_allocation(
        a, model.make_fixed_potentials("nonzero_main_effects_without_interaction")
    )
    assert nuisance["greenhouse_adjusted_estimate"] == pytest.approx(
        null["greenhouse_adjusted_estimate"], abs=1e-11
    )
    assert nuisance["estimated_cage_model_SE"] == pytest.approx(
        null["estimated_cage_model_SE"], abs=1e-11
    )
    assert nuisance["raw_unblocked_estimate"] == pytest.approx(
        null["raw_unblocked_estimate"], abs=1e-11
    )


def test_randomization_calibration_does_not_emit_pvalues_or_empirical_power_claims():
    x = model.calibrate(
        "greenhouse_varying_main_effects_without_interaction",
        n_randomizations=64, seed=20261014
    )
    assert x["truth_known_only_because_potential_outcomes_are_invented"] == pytest.approx(0)
    assert x["hypothetical_original_whole_cages"] == 30
    assert x["Monte_Carlo_allocation_draws"] == 64
    assert x["not_butterfly_sampling_power_or_validity_estimate"]
    assert "p_value" not in x
    assert not any("significant" in k.lower() for k in x)
    assert 0 <= x["fraction_nominal_model_t_intervals_cover_hypothetical_fixed_population_D"] <= 1
    assert x["empirical_SD_of_adjusted_randomization_distribution"] >= 0
    assert math.isfinite(x["mean_model_based_cage_SE"])


def test_positive_effect_hypothetical_parameter_is_not_empirically_estimated():
    x = model.calibrate("hypothetical_positive_resource_rescue", n_randomizations=64)
    assert x["truth_known_only_because_potential_outcomes_are_invented"] == pytest.approx(.15)
    assert abs(x["mean_greenhouse_adjusted_primary"] - .15) < .05
    assert x["adjusted_root_mean_squared_error"] >= 0
    assert x["raw_root_mean_squared_error"] >= 0


def test_small_set_of_randomized_schedules_requires_an_explicit_minimum():
    with pytest.raises(ValueError, match="16"):
        model.calibrate("global_no_effects", n_randomizations=1)
    with pytest.raises(ValueError, match="16"):
        model.calibrate("global_no_effects", n_randomizations=10001)
    with pytest.raises(ValueError, match="Unknown"):
        model.calibrate("unobserved species", n_randomizations=64)
    with pytest.raises(ValueError, match="30"):
        model.true_average_interaction({})
    with pytest.raises(ValueError, match="Four complete"):
        model.true_average_interaction({
            cage: {"NATURAL_NONE":.5} for cage in range(30)
        })


def test_deterministic_seed_replays_full_omission_randomization():
    one = model.calibrate("global_no_effects", n_randomizations=64, seed=12)
    two = model.calibrate("global_no_effects", n_randomizations=64, seed=12)
    assert one == two
    assert one["truth_known_only_because_potential_outcomes_are_invented"] == pytest.approx(0)


def test_design_outcomes_are_distinct_from_known_zero_average_heterogeneity():
    x = model.calibrate(
        "opposing_greenhouse_interactions_zero_average",
        n_randomizations=64, seed=11
    )
    assert x["truth_known_only_because_potential_outcomes_are_invented"] == pytest.approx(0)
    # Site-specific interactions cancel *on average*, but are nonzero locally.
    p = model.make_fixed_potentials("opposing_greenhouse_interactions_zero_average")
    local = [
        r["SUPPORTED_ADD_S"] - r["SUPPORTED_NONE"] -
        r["NATURAL_ADD_S"] + r["NATURAL_NONE"]
        for r in p.values()
    ]
    assert min(local) < 0 < max(local)
    assert x["not_butterfly_sampling_power_or_validity_estimate"]
