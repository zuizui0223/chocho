#!/usr/bin/env python3
"""Future independent 10-greenhouse 30-cage Kyoto butterfly ITT analysis.

No field data are bundled; tests use synthetic fixtures. The randomized four
treatments remain S-added yes/no × accessible host-foliage support yes/no.
The original 30 greenhouse trial is PRIOR ART, not new adult-emergence data.

The new near-balanced allocation is source-inventory/seed checked. Per-cage
initial focal A larval cohorts and added S larvae are separately tracked.
Food-shortage observations are post-treatment diagnostics only.

An OLS fixed-greenhouse estimate is reported with MODEL-BASED uncertainty,
not a guaranteed interval under heterogeneous cage or family effects.
Raw randomized-arm risks are also reported without greenhouse conditioning.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math

import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

from analyze_kyoto_four_arm_greenhouse_bibd import (
    ASSIGNMENT_COLUMNS, FATE_COLUMNS, RESOURCE_COLUMNS,
    read_csv, nonnegative, fates_by_cage, foliage_by_cage, regression,
)
from randomize_kyoto_four_arm_greenhouse_bibd import ARMS, ALLOCATION_COLUMNS
from randomize_kyoto_30cage_nearbalanced import (
    B, K, DESIGN, VARIANCE_FACTOR_IF_HOMOSCEDASTIC, assign, design_counts,
)

T_CRITICAL_DF17 = 2.1098155778331806
PRIMARY = "(SUPPORTED_ADD_S-SUPPORTED_NONE)-(NATURAL_ADD_S-NATURAL_NONE)"


def check_assignment(allocation: list[dict], info: list[dict]) -> dict:
    if len(allocation) != B * K:
        raise ValueError("Exactly 30 original whole-cage assignments required")
    original = {}
    for row in allocation:
        cage = str(row["cage_id"]).strip()
        if not cage or cage in original:
            raise ValueError("Original cage IDs missing or duplicated")
        if row["design_version"] != DESIGN or row["treatment"] not in ARMS:
            raise ValueError("Not this prospective near-balanced design")
        if (nonnegative(row["competitor_present"], "original S assignment", integer=True),
                nonnegative(row["food_support"], "original support assignment", integer=True)) != ARMS[row["treatment"]]:
            raise ValueError("Original factors do not match allocated arm")
        if not row["greenhouse_id"].strip():
            raise ValueError("Blank original greenhouse ID")
        nonnegative(row["random_seed"], "random seed", integer=True)
        original[cage] = row
    if len({r["random_seed"] for r in original.values()}) != 1:
        raise ValueError("Allocation seeds disagree")
    structure = design_counts(list(original.values()))
    if not structure["design_valid"]:
        raise ValueError("Original 30-cage near-balanced incidence was modified")
    seed = int(allocation[0]["random_seed"])
    replay = assign(
        [{"cage_id": row["cage_id"], "greenhouse_id": row["greenhouse_id"]}
         for row in original.values()],
        seed=seed,
    )
    for expected in replay:
        recorded = original[expected["cage_id"]]
        if any(str(expected[col]) != str(recorded[col]) for col in ALLOCATION_COLUMNS):
            raise ValueError("Original random seed replay mismatch: no posthoc cage relabeling")

    if len(info) != len(allocation):
        raise ValueError("Original randomized whole cages cannot be omitted from analysis")
    cages: dict[str, dict] = {}
    for row in info:
        cage = row["cage_id"]
        if not cage or cage not in original or cage in cages:
            raise ValueError("Missing, unrandomized or duplicate experimental cage")
        if any(row[col] != original[cage][col]
               for col in ("greenhouse_id", "treatment", "competitor_present", "food_support")):
            raise ValueError("Source greenhouse, treatment or randomized factors changed")
        native = nonnegative(row["focal_a_initial_n"], "focal A original count", integer=True, positive=True)
        exotic = nonnegative(row["added_s_initial_n"], "added S original count", integer=True)
        assigned_s = ARMS[row["treatment"]][0]
        if (exotic > 0) != bool(assigned_s):
            raise ValueError("Original S density violates assigned S absence/presence")
        n_visits = nonnegative(row["expected_n_food_visits"], "planned foliage visits", integer=True, positive=True)
        if n_visits < 2:
            raise ValueError("Predeclare repeated foliage observation opportunities")
        floor = nonnegative(row["food_floor_cm2"], "food floor", positive=True)
        if not row["donor_batch_id"].strip():
            raise ValueError("Independent provenance of donor foliage required")
        cages[cage] = {**row, "focal_a_initial_n": native,
                       "added_s_initial_n": exotic,
                       "expected_n_food_visits": n_visits,
                       "food_floor_cm2": floor}
    if len({r["focal_a_initial_n"] for r in cages.values()}) != 1:
        raise ValueError("Keep initially assigned focal cohort size fixed")
    assigned_s_values = {r["added_s_initial_n"] for r in cages.values() if r["added_s_initial_n"] > 0}
    if len(assigned_s_values) != 1:
        raise ValueError("Added S cohort size must be predeclared equal in all positive arms")
    by_greenhouse: dict[str, list[dict]] = defaultdict(list)
    for cage in cages.values():
        by_greenhouse[cage["greenhouse_id"]].append(cage)
    for greenhouse, rr in by_greenhouse.items():
        if len(rr) != K:
            raise ValueError("Incomplete original greenhouse, no outcome-dependent dropped cages")
        for key in ("expected_n_food_visits", "food_floor_cm2", "donor_batch_id"):
            if len({str(r[key]) for r in rr}) != 1:
                raise ValueError(f"Greenhouse {greenhouse}: inconsistent initial {key}")
    return cages


def heteroscedastic_cage_sensitivity(cages: dict, outcomes: dict, adjusted: dict) -> dict:
    """HC3 cage-level sandwich sensitivity, NOT a small-sample exact CI.

    With fixed greenhouse intercepts, a correctly coded covariance calculation
    must use the same exact unit-level contrast weights as the primary fit.
    This does not handle correlated treatment responses among greenhouses or
    prove significance with only 10 original site blocks.
    """
    original=sorted(cages.values(),key=lambda r:(r["greenhouse_id"],r["cage_id"]))
    gh=sorted({row["greenhouse_id"] for row in original})
    other=[arm for arm in ARMS if arm!="NATURAL_NONE"]
    X=np.array([
        [int(r["greenhouse_id"]==block) for block in gh] +
        [int(r["treatment"]==arm) for arm in other]
        for r in original
    ],dtype=float)
    y=np.array([outcomes[r["cage_id"]]["focal_adult_low"] for r in original],dtype=float)
    xtxi=np.linalg.inv(X.T@X)
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    errors=y-X@beta
    leverages=np.einsum("ij,jk,ik->i",X,xtxi,X)
    if np.any(leverages>=1-1e-10) or np.any(leverages<0):
        raise ValueError("Unsupported leverage for HC3 in small incomplete block model")
    exact_weights=np.array([
        adjusted["whole_cage_linear_contrast_weights"][r["cage_id"]]
        for r in original
    ])
    variance=float(np.sum((exact_weights*errors/(1-leverages))**2))
    if variance < 0 or not math.isfinite(variance):
        raise ValueError("Numerical problem in HC3 variance calculation")
    return {
        "hc3_cage_level_primary_SE":math.sqrt(variance),
        "max_hat_leverage":float(max(leverages)),
        "small_sample_cage_heteroscedasticity_sensitivity_only":True,
        "NOT_greenhouse_cluster_robust_or_validated_biological_inference":True,
    }


def analyze(allocation: list[dict], info: list[dict],
            fates: list[dict], foliage: list[dict]) -> dict:
    cages = check_assignment(allocation, info)
    outcomes = fates_by_cage(cages, fates)
    food = foliage_by_cage(cages, foliage)
    adjusted = regression(list(cages.values()), outcomes)
    if adjusted["n_residual_df"] != 17:
        raise RuntimeError("Original 10x3 design must have 17 residual degrees of freedom")
    multiplier = adjusted["analytical_contrast_variance_multiplier"]
    if not math.isclose(multiplier, VARIANCE_FACTOR_IF_HOMOSCEDASTIC, abs_tol=1e-11):
        raise RuntimeError("Near-balanced primary contrast information no longer equals 0.6")
    estimate = adjusted["greenhouse_adjusted_observed_adult_contrast"]
    se = adjusted["model_based_standard_error_given_homoscedastic_cage_residuals"]
    adjusted["nominal_t17_95pct_interval_under_cage_error_assumptions"] = [
        estimate - T_CRITICAL_DF17 * se,
        estimate + T_CRITICAL_DF17 * se,
    ]
    adjusted["model_is_not_design_based_precision_or_empirical_population_generalization"] = True
    by_arm = {}
    risks = {}
    for arm in ARMS:
        ids = [c for c, row in cages.items() if row["treatment"] == arm]
        original_focal = sum(outcomes[c]["original_focal_A"] for c in ids)
        adult = sum(outcomes[c]["observed_flight_capable_focal_adults"] for c in ids)
        pending_diapause = sum(outcomes[c]["living_diapause_focal"] for c in ids)
        pending_other = sum(
            outcomes[c]["living_pupa_pending_focal"] +
            outcomes[c]["other_focal_unknown"] for c in ids
        )
        risks[arm] = adult / original_focal
        by_arm[arm] = {
            "original_randomized_whole_cages": len(ids),
            "original_assigned_focal_A_larvae": original_focal,
            "confirmed_flight_capable_focal_A_adults": adult,
            "living_diapausing_focal_A_pupae": pending_diapause,
            "other_living_or_unknown_focal_A": pending_other,
            "raw_observed_flight_capable_adult_risk": risks[arm],
            "worst_case_possible_eventual_adult_upper": (
                adult + pending_diapause + pending_other
            ) / original_focal,
            "food_access_visits_measured_below_floor": sum(
                food[c]["measured_edible_area_below_original_floor"] for c in ids
            ),
            "missing_food_access_measurements": sum(
                food[c]["missing_edible_area_visits"] for c in ids
            ),
            "food_access_measurement_not_a_survival_selection_rule": True,
        }
    raw_primary = (risks["SUPPORTED_ADD_S"] - risks["SUPPORTED_NONE"]
                   - risks["NATURAL_ADD_S"] + risks["NATURAL_NONE"])
    omitted = []
    for greenhouse in sorted({c["greenhouse_id"] for c in cages.values()}):
        reduced = [c for c in cages.values() if c["greenhouse_id"] != greenhouse]
        fit = regression(reduced, outcomes)
        omitted.append({
            "greenhouse_left_out": greenhouse,
            "other_greenhouses_remaining": 9,
            "adjusted_primary_for_descriptive_sensitivity": fit[
                "greenhouse_adjusted_observed_adult_contrast"
            ],
        })
    reduced_estimates=[
        row["adjusted_primary_for_descriptive_sensitivity"] for row in omitted
    ]
    delete_mean=sum(reduced_estimates)/len(reduced_estimates)
    gh_jackknife_se=math.sqrt(
        (len(reduced_estimates)-1)/len(reduced_estimates) *
        sum((value-delete_mean)**2 for value in reduced_estimates)
    )
    hc3=heteroscedastic_cage_sensitivity(cages,outcomes,adjusted)
    return {
        "schema": "chocho_kyoto_30cage_nearbalanced_original_focal_adults_v01",
        "status": "FUTURE_AUTHENTIC_FIELD_DATA_REQUIRED_NO_BUTTERFLY_RESULTS",
        "randomization_scope": "10 original distinct greenhouses x 3 whole cages each",
        "design": design_counts(allocation),
        "primary_estimand": PRIMARY,
        "raw_intention_to_treat_arm_mean_difference_in_differences": raw_primary,
        "primary_greenhouse_fixed_effect_adjusted": adjusted,
        "cage_heteroscedasticity_HC3_sensitivity_only":hc3,
        "ten_greenhouse_delete_one_jackknife_SE_descriptive_only":gh_jackknife_se,
        "jackknife_SE_is_NOT_randomization_exact_or_power_evidence":True,
        "per_arm_original_cohort_outcomes": by_arm,
        "leave_one_greenhouse_out_description_only": omitted,
        "original_focal_A_count_present_in_all_randomized_cages": True,
        "original_added_S_count_separately_followed": True,
        "no_postrandomization_food_or_survivor_selection": True,
        "greenhouse_specific_treatment_effects_not_identifiable_with_one_cage_per_cell": True,
        "source_maternal_families_independent_not_auto_verified_by_code": True,
        "two_stage_randomization_applied_prior_to_outcomes": True,
        "source_host_is_native_plant_not_invasive_host_distribution_evidence": True,
        "GEB_main_untouched": True,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ("allocation", "assignments", "fates", "resources", "receipt"):
        ap.add_argument("--" + arg, required=True, type=Path)
    args = ap.parse_args()
    inputs = {
        "allocation": args.allocation, "assignments": args.assignments,
        "fates": args.fates, "resources": args.resources,
    }
    out = analyze(
        read_csv(args.allocation, ALLOCATION_COLUMNS),
        read_csv(args.assignments, ASSIGNMENT_COLUMNS),
        read_csv(args.fates, FATE_COLUMNS),
        read_csv(args.resources, RESOURCE_COLUMNS),
    )
    out["authentic_source_sha256"] = {
        key: hashlib.sha256(path.read_bytes()).hexdigest()
        for key, path in inputs.items()
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
