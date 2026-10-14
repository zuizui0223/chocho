#!/usr/bin/env python3
"""POST-HOC source-only butterfly food-intake shape comparison at matched ages.

Excludes missing intervals rather than assuming no consumption. Original
figshare feed records are analyzed via reconstruct_kyoto_late_instar_intake.
NO competitive survival, randomized food pulses, or generalized species-level
phenology effects are estimated.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import reconstruct_kyoto_late_instar_intake as original

PRIMARY_AGES = (8, 10, 12, 14, 16, 18)
ALT_WINDOWS = ((8, 10, 12, 14, 16), (10, 12, 14, 16, 18))
EARLY_AGES = (8, 10, 12)
LATE_AGES = (14, 16, 18)


def _mean(values):
    if not values:
        raise ValueError("empty vector is not an estimable response")
    return sum(values) / len(values)


def _profile_for_species(events: list[dict], *, ages: tuple[int, ...], expected_ids: int) -> dict:
    if len(set(ages)) != len(ages) or not ages or tuple(sorted(ages)) != ages:
        raise ValueError("unique sorted age bins required")
    age_set = set(ages)
    by_id: dict[str, dict[int, float]] = defaultdict(dict)
    observed_ids = set()
    for record in events:
        identifier = str(record["original_larva_id"])
        day = record["end_day_since_founder_hatch"]
        mass = record["consumed_leaf_dry_mg"]
        if not identifier or not isinstance(day, int) or not math.isfinite(mass) or mass < 0:
            raise ValueError("invalid source age/mass/ID")
        observed_ids.add(identifier)
        if day not in age_set:
            continue
        if day in by_id[identifier]:
            raise ValueError("duplicate original larva x age observation")
        by_id[identifier][day] = mass

    if not (1 < expected_ids) or len(observed_ids) > expected_ids:
        raise ValueError("original individual IDs not reconcilable with source count")
    included = []
    exclusions = []
    for identifier in sorted(observed_ids):
        present = by_id.get(identifier, {})
        missing = sorted(age_set - set(present))
        mass_sum = sum(present.values())
        if missing or mass_sum <= 0:
            exclusions.append({
                "original_id": identifier,
                "missing_ages": missing,
                "positive_matched_mass": mass_sum > 0,
            })
            continue
        included.append({
            "original_id": identifier,
            "total_dry_mg_in_window": mass_sum,
            "age_shares": {str(day): present[day] / mass_sum for day in ages},
            "age_center_of_mass": sum(day * present[day] / mass_sum for day in ages),
        })
    if len(included) < 2:
        raise ValueError("fewer than two individual profiles with complete matched observation")

    per_age = {str(t): _mean([x["age_shares"][str(t)] for x in included]) for t in ages}
    if abs(sum(per_age.values()) - 1) > 1e-10:
        raise ValueError("normalized demand shares fail conservation")

    weighted_mass_by_age = {str(day): sum(x["age_shares"][str(day)] * x["total_dry_mg_in_window"]
                                        for x in included) for day in ages}
    weighted_total = sum(weighted_mass_by_age.values())
    pooled = {day: mass / weighted_total for day, mass in weighted_mass_by_age.items()}

    return {
        "total_tracked_ids_in_original": expected_ids,
        "unique_original_ids_with_any_2day_matched_intake": len(observed_ids),
        "included_complete_matched_age_ids": len(included),
        "excluded_ids_and_reason": exclusions,
        "each_larva_equal_weight_normalized_age_profile": per_age,
        "mass_weighted_normalized_age_profile": pooled,
        "individual_age_centers": [x["age_center_of_mass"] for x in included],
        "mean_age_center_of_mass": _mean([x["age_center_of_mass"] for x in included]),
        "median_age_center_of_mass": statistics.median([x["age_center_of_mass"] for x in included]),
        "individual_equal_weight_age_shares": included,
    }


def profile_comparison(source_by_species: dict[str, dict], *, ages: tuple[int, ...]) -> dict:
    profiles = {
        sp: _profile_for_species(
            source_by_species[sp]["interval_events"],
            ages=ages,
            expected_ids=source_by_species[sp]["n_single_larva_ids"],
        ) for sp in ("A", "S")
    }
    av = profiles["A"]["each_larva_equal_weight_normalized_age_profile"]
    sv = profiles["S"]["each_larva_equal_weight_normalized_age_profile"]
    overlap = sum(min(av[str(day)], sv[str(day)]) for day in ages)
    tv = 0.5 * sum(abs(av[str(day)] - sv[str(day)]) for day in ages)
    if not math.isclose(1 - overlap, tv, abs_tol=1e-10):
        raise ValueError("total-variation/overlap identity failed")

    pa = profiles["A"]["mass_weighted_normalized_age_profile"]
    ps = profiles["S"]["mass_weighted_normalized_age_profile"]

    return {
        "age_bin_end_days_since_source_species_founder_start": list(ages),
        "by_species": profiles,
        "temporal_demand_profile_overlap": overlap,
        "total_variation_between_normalized_profiles": tv,
        "mass_weighted_temporal_overlap_sensitivity": sum(min(pa[str(t)], ps[str(t)]) for t in ages),
        "A_minus_S_mean_age_center_of_mass_days": (
            profiles["A"]["mean_age_center_of_mass"] - profiles["S"]["mean_age_center_of_mass"]
        ),
    }


def primary_with_leave_one_out(source_by_species: dict[str, dict]) -> dict:
    primary = profile_comparison(source_by_species, ages=PRIMARY_AGES)
    a_rows = primary["by_species"]["A"]["individual_equal_weight_age_shares"]
    s_rows = primary["by_species"]["S"]["individual_equal_weight_age_shares"]
    if not set(EARLY_AGES).isdisjoint(LATE_AGES) or set(EARLY_AGES + LATE_AGES) != set(PRIMARY_AGES):
        raise ValueError("the two prespecified windows must partition the support")
    a_early = sum(primary["by_species"]["A"]["each_larva_equal_weight_normalized_age_profile"][str(t)]
                  for t in EARLY_AGES)
    s_early = sum(primary["by_species"]["S"]["each_larva_equal_weight_normalized_age_profile"][str(t)]
                  for t in EARLY_AGES)
    difference = a_early - s_early
    alternatives = []
    for drop, rows in (("A", a_rows), ("S", s_rows)):
        for idx in range(len(rows)):
            survivors = [row for j, row in enumerate(rows) if j != idx]
            if len(survivors) < 1:
                raise ValueError("leave-one-out not supported")
            early = _mean([sum(x["age_shares"][str(t)] for t in EARLY_AGES) for x in survivors])
            value = early - s_early if drop == "A" else a_early - early
            alternatives.append({
                "dropped_species": drop,
                "original_id": rows[idx]["original_id"],
                "A_minus_S_early_share": value,
            })
    profile_result = {
        **primary,
        "early_support_end_days": list(EARLY_AGES),
        "late_support_end_days": list(LATE_AGES),
        "A_early_share": a_early,
        "S_early_share": s_early,
        "A_minus_S_early_share": difference,
        "leave_one_individual_out": {
            "min_A_minus_S_early_share": min(x["A_minus_S_early_share"] for x in alternatives),
            "max_A_minus_S_early_share": max(x["A_minus_S_early_share"] for x in alternatives),
            "sign_consistent_with_full_sample": all(
                (x["A_minus_S_early_share"] > 0) == (difference > 0)
                and (x["A_minus_S_early_share"] < 0) == (difference < 0)
                for x in alternatives
            ),
            "iterations": alternatives,
        },
    }
    if 0 < a_early < 1 and 0 < s_early < 1 and a_early != s_early:
        bounds = sorted((s_early / a_early, (1 - s_early) / (1 - a_early)))
        profile_result["hypothetical_reversal_condition"] = {
            "possible_under_illustrative_scarcity_allocation": True,
            "direction_if_sensitivities_equal": (
                "early-shortage relatively harms A; late-shortage relatively harms S"
                if a_early > s_early
                else "early-shortage relatively harms S; late-shortage relatively harms A"
            ),
            "k_A_over_k_S_open_interval": bounds,
            "warning": "This is only a derived model corridor, not observed competition or emergence.",
        }
    else:
        profile_result["hypothetical_reversal_condition"] = {
            "possible_under_illustrative_scarcity_allocation": False,
            "warning": "No admissible nonsingular corridor in these aggregate profiles.",
        }
    return profile_result


def analyze(reconstruction: dict) -> dict:
    if reconstruction.get("status") != "SOURCE_RECONSTRUCTED_DESCRIPTIVE_ONLY":
        raise ValueError("no verified original reconstruction")
    results = reconstruction["per_species"]
    if any(not results[s]["table1_dry_intake_within_2pct"] for s in ("A", "S")):
        raise ValueError("independent Table 1 mass calibration reproduction failed")
    primary = primary_with_leave_one_out(results)
    secondary = {}
    for ages in ALT_WINDOWS:
        cmp = profile_comparison(results, ages=ages)
        secondary["-".join(map(str, ages))] = {
            "ages": list(ages),
            "complete_A": cmp["by_species"]["A"]["included_complete_matched_age_ids"],
            "complete_S": cmp["by_species"]["S"]["included_complete_matched_age_ids"],
            "overlap": cmp["temporal_demand_profile_overlap"],
            "total_variation": cmp["total_variation_between_normalized_profiles"],
            "A_minus_S_age_center_days": cmp["A_minus_S_mean_age_center_of_mass_days"],
        }
    return {
        "schema": "chocho_kyoto_matched_biological_age_intake_shape_v01",
        "status": "POSTHOC_WITHIN_SOURCE_DESCRIPTOR_NOT_INDEPENDENT_ECOLOGICAL_EFFECT",
        "original_MD5s": reconstruction["source_identity"],
        "primary_six_bin_complete_case": primary,
        "alternate_age_window_sensitivity": secondary,
        "study_limitations": [
            "one separately reared original 2016 cohort per species",
            "2-day resolution, incomplete later ages and complete-case attrition",
            "individual larvae are not independent maternal/population replicates",
            "demand amounts conditional on measured third-to-fifth instars",
            "age since separate hatch groups is not wild plant phenology",
            "no experimental change in host-resource supply timing",
            "no contemporaneous true edible resource availability in mixed cages",
            "no adult eclosion, survival, diapause or field demographic estimate",
        ],
        "main_GEB_Paper_unchanged": True,
    }


def run():
    raw = original.audit()
    raw["status"] = "SOURCE_RECONSTRUCTED_DESCRIPTIVE_ONLY"
    return analyze(raw)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    try:
        outcome = run()
    except Exception as exc:
        outcome = {
            "schema": "chocho_kyoto_matched_biological_age_intake_shape_v01",
            "status": "SOURCE_OR_MATCHED_SUPPORT_BLOCKED",
            "reason": f"{type(exc).__name__}: {str(exc)[:260]}",
        }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(outcome, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(outcome, indent=2, ensure_ascii=False))
    if outcome["status"] == "SOURCE_OR_MATCHED_SUPPORT_BLOCKED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
