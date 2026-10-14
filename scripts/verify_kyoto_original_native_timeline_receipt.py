#!/usr/bin/env python3
"""Fail closed if the pinned ORIGINAL 2012 Kyoto source no longer reproduces the receipt.

The expected file is a scientific SOURCE provenance claim, not pilot outcomes.
No estimates are refitted, no MD5 mismatch tolerated, no raw source corrections.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def verify(original: dict, frozen: dict) -> dict:
    if original.get("status") != "AUTHOR_ORIGINAL_SOURCE_DESCRIPTIVE_POSTHOC_NOT_NEW_CAUSAL_MECHANISM":
        raise ValueError("Original-source audit withheld the primary descriptive outcome")
    if frozen.get("status") != "MD5_VERIFIED_POSTHOC_SOURCE_DESCRIPTIVE_NO_NEW_ANIMAL_EXPERIMENT":
        raise ValueError("Frozen numeric receipt has unexpected biological status")
    p=frozen["provenance"]
    if original["original_source_md5"]!=p["original_file_md5"]:
        raise ValueError("Wrong original Figshare CSV identity or checksum")
    if (original["whole_cages"],original["cage_x_visit_original_rows"]) != (
        p["independent_original_cages"], p["cage_repeated_censuses"]
    ):
        raise ValueError("Original source independent cage/visit unit mismatch")
    a=frozen["stage_accounting"]
    s=original["source_stage_total_and_cumulative_reconciliation"]["counts"]
    required={
        "native_total_a_equals_stages_plus_cumul_ap":
            "native_total_equals_current_all_stages_plus_cumulative_pupae",
        "introduced_total_s_equals_stages_plus_cumul_sp":
            "introduced_total_equals_current_all_stages_plus_cumulative_pupae",
        "source_rows_where_native_total_plus_cumul_exceeds_original_cohort":
            "native_total_plus_cumulative_pupae_exceeds_original_assigned",
    }
    for frozen_key, source_key in required.items():
        if a[frozen_key]!=s[source_key]:
            raise ValueError(f"Frozen original stage reconciliation failed: {frozen_key}")
    if original["matched_exact_native_density_early_cage_count"] != a["exact_A_matched_density_cages"]:
        raise ValueError("Initial native matched group size changed")
    early=original["matched_density_day7_15_21_source_accounting"]
    for label,key in (
        ("day_7_15_21_matched_cage_by_census_rows","whole_cage_x_early_census_rows"),
        ("all_48_early_rows_total_equals_stages_plus_cumulative","source_total_equals_stage_count_plus_pupation"),
        ("all_48_early_rows_cumulative_pupae_zero","no_A_pupa_already_recorded"),
    ):
        expected=a[label]
        expected=early["whole_cage_x_early_census_rows"] if expected is True and label.startswith("all_") else expected
        if expected!=early[key]:
            raise ValueError(f"Early original stage/accounting gate failed: {label}")
    blocks=original["native_density_matched_strata"]
    treatment=frozen["originally_published_density_treatment_structure"]
    if [int(k) for k in blocks] != treatment["exact_matched_native_A_initial_counts"]:
        raise ValueError("Source native-density strata do not reproduce original treatment structure")
    if sum(2*len((x["lower_S_initial_n"],x["higher_S_initial_n"])) for x in blocks.values()) !=16:
        raise ValueError("Pair count is not the expected original independent cages")
    for key,pair in blocks.items():
        if [pair["lower_S_initial_n"],pair["higher_S_initial_n"]] != treatment["lower_and_higher_S_by_A"][key]:
            raise ValueError("Original S high/low source density contrast shifted")
        if pair["cages_per_density_cell"] != treatment["cages_per_combination"]:
            raise ValueError("Incorrect original cage sampling unit")
    for day in frozen["day_comparisons"]:
        label=str(day["day"])
        observed=original["equal_native_stratum_weight_summary_by_day"].get(label)
        if observed is None:
            raise ValueError("Missing previously frozen original census day")
        pairs_to_compare=(
            ("raw_total_a_fraction_high_S_minus_low_S","mean_higher_minus_lower_S_raw_total_A_fraction_across_strata"),
            ("pupation_fraction_high_S_minus_low_S","mean_higher_minus_lower_S_pupation_fraction_across_A_density_strata"),
        )
        for expected,actual in pairs_to_compare:
            if not math.isclose(day[expected],observed[actual],abs_tol=1e-12,rel_tol=0):
                raise ValueError(f"Original fixed-native effect changed on day {label}: {expected}")
        observed_strata=observed["source_total_A_fraction_stratum_differences"]
        if len(observed_strata)!=len(day["raw_by_original_A"]) or any(
            not math.isclose(a,b,abs_tol=1e-12,rel_tol=0)
            for a,b in zip(observed_strata,day["raw_by_original_A"])
        ):
            raise ValueError(f"Original per-A stratum timing changed on day {label}")
    observed_anomalies={
        (str(row["source_cage_id"]),row["species"],row["day"],
         row["previous_cumulative_pupation"],row["current_source_cumulative_pupation"])
        for row in original["source_data_quality_anomalies"]
    }
    expected_anomalies={
        (x["cage"],"S" if x["species"]=="Sericinus montela" else "A",
         x["day"],x["previous"],x["current"])
        for x in frozen["source_cumulative_pupation_decreases_at_last_census"]
    }
    if observed_anomalies != expected_anomalies:
        raise ValueError("Original source day-48 cumulative anomalies changed")
    if not all(original["terminal_fraction_withheld_for_species_with_cumulative_anomaly"].values()):
        raise ValueError("Original invalid terminal cumulative pupation was reinterpreted")
    if frozen["new_butterfly_experimental_results"] or not original["GEB_PR38_untouched"]:
        raise ValueError("Unwarranted promotion to independent biological results or main GEB changes")
    return {
        "schema":"chocho_kyoto_source_timeline_frozen_integrity_v01",
        "status":"FULL_SOURCE_MD5_AND_FROZEN_TEMPORAL_ORIGINAL_VALUES_REPRODUCED",
        "source_independent_cages":30,
        "matched_original_native_density_cages":16,
        "original_7_15_21_source_rows_checked":48,
        "anomalous_last_day_original_entries_kept":len(observed_anomalies),
        "comparative_days_checked":[x["day"] for x in frozen["day_comparisons"]],
        "source_only_not_adult_emergence_or_resource_causality":True,
        "main_GEB_untouched":True,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original",type=Path,required=True)
    ap.add_argument("--frozen",type=Path,required=True)
    args=ap.parse_args()
    original=json.loads(args.original.read_text(encoding="utf-8"))
    frozen=json.loads(args.frozen.read_text(encoding="utf-8"))
    print(json.dumps(verify(original,frozen),indent=2))


if __name__=="__main__":
    main()
