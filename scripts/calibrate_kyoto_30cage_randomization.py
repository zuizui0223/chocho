#!/usr/bin/env python3
"""Design-based CALIBRATION with fabricated fixed potential adult-risk curves.

This is a deterministic *method stress-test*, NOT reanalysis of real butterflies.
Potential values are arbitrary bounded expected focal-adult success proportions
for 30 cages × four treatment assignments, not collected adult counts.

The 10-greenhouse near-balanced randomizer is REPLAYED in each draw, including
a fresh random choice of the two eligible omission patterns, greenhouse
assignment, and within-greenhouse assignment. The primary D estimate is
refitted using greenhouse fixed effects. We compare observed estimates with
the finite-population average interaction, known ONLY in simulated scenarios.

Crucial distinction: a Fisher permutation under the GLOBAL SHARP NULL
(all four potential outcomes identical) cannot directly test the WEAK null
of zero average factorial interaction when competitor and food main effects
are nonzero. This tool does NOT issue randomization-test p-values or confuse
the two hypotheses. Its model-t coverage frequencies are synthetic diagnostic
frequencies, NOT biological power or actual statistical validity guarantees.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

from analyze_kyoto_four_arm_greenhouse_bibd import regression
from randomize_kyoto_30cage_nearbalanced import ARMS, assign, design_counts

T_DF17 = 2.1098155778331806
SCENARIOS = (
    "global_no_effects",
    "nonzero_main_effects_without_interaction",
    "greenhouse_varying_main_effects_without_interaction",
    "opposing_greenhouse_interactions_zero_average",
    "hypothetical_positive_resource_rescue",
)


def _inventory() -> list[dict[str, str]]:
    return [
        {"cage_id": f"SIM:G{g:02d}:C{c}", "greenhouse_id": f"SIM:G{g:02d}"}
        for g in range(10) for c in range(3)
    ]


def make_fixed_potentials(scenario: str) -> dict[str, dict[str, float]]:
    if scenario not in SCENARIOS:
        raise ValueError("Unknown explicitly hypothetical simulation scenario")
    outcomes = {}
    for gh in range(10):
        angle = 2 * math.pi * gh / 10
        for cage in range(3):
            cid = f"SIM:G{gh:02d}:C{cage}"
            # An arbitrary deterministic intercept with greenhouse variation
            # plus within-greenhouse cage heterogeneity, NOT a prior estimate.
            baseline = .51 + .08 * math.sin(angle) + .033 * (cage - 1)
            competitor, support, interaction = 0.0, 0.0, 0.0
            if scenario != "global_no_effects":
                competitor, support = -.11, .12
            if scenario in (
                "greenhouse_varying_main_effects_without_interaction",
                "opposing_greenhouse_interactions_zero_average",
            ):
                competitor += .15 * math.sin(angle)
                support += .09 * math.cos(angle)
            if scenario == "opposing_greenhouse_interactions_zero_average":
                interaction = .20 * math.sin(angle)
            if scenario == "hypothetical_positive_resource_rescue":
                interaction = .15
            by_arm = {
                name: baseline + competitor * comp + support * food
                + interaction * comp * food
                for name, (comp, food) in ARMS.items()
            }
            if min(by_arm.values()) < 0 or max(by_arm.values()) > 1:
                raise ValueError("Toy success probabilities fall outside 0..1")
            outcomes[cid] = by_arm
    return outcomes


def true_average_interaction(potentials: dict[str, dict[str, float]]) -> float:
    if len(potentials) != 30:
        raise ValueError("30 fixed original hypothetical cages required")
    for row in potentials.values():
        if set(row) != set(ARMS):
            raise ValueError("Four complete potential responses are essential for this audit")
    return statistics.mean(
        row["SUPPORTED_ADD_S"] - row["SUPPORTED_NONE"]
        - row["NATURAL_ADD_S"] + row["NATURAL_NONE"]
        for row in potentials.values()
    )


def evaluate_one_allocation(allocation: list[dict], potentials: dict) -> dict:
    if len(allocation) != 30 or not design_counts(allocation)["design_valid"]:
        raise ValueError("Expected the complete original contrast-optimal allocation")
    assigned = {
        row["cage_id"]: potentials[row["cage_id"]][row["treatment"]]
        for row in allocation
    }
    pseudo_outcomes = {
        cage: {"focal_adult_low": value, "focal_eventual_adult_upper": value}
        for cage, value in assigned.items()
    }
    fit = regression(allocation, pseudo_outcomes)
    if fit["n_residual_df"] != 17:
        raise RuntimeError("Fixed greenhouse treatment design lost 17 residual df")
    if not math.isclose(fit["analytical_contrast_variance_multiplier"], .60, abs_tol=1e-10):
        raise RuntimeError("Random allocation no longer has exact 0.6 information")
    means = {
        arm: statistics.mean(
            assigned[row["cage_id"]] for row in allocation if row["treatment"] == arm
        )
        for arm in ARMS
    }
    raw = (means["SUPPORTED_ADD_S"] - means["SUPPORTED_NONE"]
           - means["NATURAL_ADD_S"] + means["NATURAL_NONE"])
    adjusted = fit["greenhouse_adjusted_observed_adult_contrast"]
    se = fit["model_based_standard_error_given_homoscedastic_cage_residuals"]
    return {
        "raw_unblocked_estimate": raw,
        "greenhouse_adjusted_estimate": adjusted,
        "nominal_independent_cage_model_95pct_interval": [
            adjusted - T_DF17 * se,
            adjusted + T_DF17 * se,
        ],
        "estimated_cage_model_SE": se,
    }


def calibrate(
    scenario: str,
    *, n_randomizations: int = 512,
    seed: int = 20261014,
) -> dict:
    if scenario not in SCENARIOS:
        raise ValueError("Unknown synthetic scenario name")
    if type(n_randomizations) is not int or not 16 <= n_randomizations <= 10000:
        raise ValueError("Require 16..10000 Monte Carlo allocations; no pseudo-precision")
    if type(seed) is not int:
        raise ValueError("integer seed required")
    toy = make_fixed_potentials(scenario)
    true_d = true_average_interaction(toy)
    invent = _inventory()
    rng = random.Random(seed)
    estimates = []
    for _ in range(n_randomizations):
        # All randomization layers are rerun, not merely cage labels permuted
        # within one fixed omission pattern.
        allocation = assign(invent, seed=rng.getrandbits(64))
        estimates.append(evaluate_one_allocation(allocation, toy))

    adjusted = [r["greenhouse_adjusted_estimate"] for r in estimates]
    raw = [r["raw_unblocked_estimate"] for r in estimates]
    se = [r["estimated_cage_model_SE"] for r in estimates]
    interval_contains_true = [
        r["nominal_independent_cage_model_95pct_interval"][0] - 1e-12 <= true_d <=
        r["nominal_independent_cage_model_95pct_interval"][1] + 1e-12
        for r in estimates
    ]
    if any(not math.isfinite(x) for x in adjusted + raw + se):
        raise RuntimeError("Nonfinite model calibration")
    return {
        "scenario": scenario,
        "truth_known_only_because_potential_outcomes_are_invented": true_d,
        "hypothetical_original_whole_cages": len(toy),
        "Monte_Carlo_allocation_draws": n_randomizations,
        "randomization_seed": seed,
        "mean_greenhouse_adjusted_primary": statistics.mean(adjusted),
        "adjusted_bias_relative_to_hypothetical_finite_population_D": (
            statistics.mean(adjusted) - true_d
        ),
        "adjusted_root_mean_squared_error": math.sqrt(
            sum((x - true_d)**2 for x in adjusted)/n_randomizations
        ),
        "empirical_SD_of_adjusted_randomization_distribution": statistics.stdev(adjusted),
        "mean_unblocked_raw_primary": statistics.mean(raw),
        "raw_bias_relative_to_hypothetical_finite_population_D": statistics.mean(raw) - true_d,
        "raw_root_mean_squared_error": math.sqrt(
            sum((x - true_d)**2 for x in raw)/n_randomizations
        ),
        "mean_model_based_cage_SE": statistics.mean(se),
        "fraction_nominal_model_t_intervals_cover_hypothetical_fixed_population_D": (
            sum(interval_contains_true)/n_randomizations
        ),
        "not_butterfly_sampling_power_or_validity_estimate": True,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--draws", type=int, default=512)
    ap.add_argument("--seed", type=int, default=20261014)
    ap.add_argument("--receipt", required=True, type=Path)
    args = ap.parse_args()
    out = {
        "schema": "chocho_kyoto_nearbalanced_randomization_stress_test_v01",
        "status": "HYPOTHETICAL_FIXED_POTENTIAL_OUTCOMES_NOT_ANIMAL_RESULTS",
        "design": "4-arm, 10 source greenhouses x 3 whole cages, same two eligible omission classes",
        "null_hypothesis_warning": (
            "No interaction is a WEAK composite null with possible nonzero "
            "competitor and foliage main effects. A naive Fisher randomization "
            "test of NO EFFECT of any arm uses a different GLOBAL SHARP null. "
            "This workflow intentionally emits NO p-value."
        ),
        "methodological_prior_art": [
            "Wu & Ding (2021) JASA DOI 10.1080/01621459.2020.1750415",
            "Muralidharan et al. (2020) NBER DOI 10.3386/w26562",
        ],
        "scenarios": [
            calibrate(scenario, n_randomizations=args.draws, seed=args.seed)
            for scenario in SCENARIOS
        ],
        "no_new_butterfly_cages_larvae_or_adults_observed": True,
        "no_real_current_greenhouse_inventory_or_family_effect_variance": True,
        "GEB_main_not_touched": True,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
