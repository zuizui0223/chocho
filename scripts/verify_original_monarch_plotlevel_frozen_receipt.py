#!/usr/bin/env python3
"""Independent guarded replay of already-published 2021 plot-level outcome counts."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

LABELS = {
    "missing_assigned_no_observed_adult": "missing_as_no_recorded_adult",
    "missing_possible_one_adult": "missing_as_one_possible_adult",
}
COMPARISONS = {
    "ambient_host_difference": "observed_plot_paired_host_difference_ambient",
    "elevated_host_difference": "observed_plot_paired_host_difference_elevated",
    "interaction": "observed_temperature_by_species_adult_emergence_risk_interaction",
}
INTERVALS = {
    "interaction_bootstrap_95pct": "interaction",
    "ambient_host_difference_bootstrap_95pct": "ambient_host_difference",
    "elevated_host_difference_bootstrap_95pct": "elevated_host_difference",
}


def close(x, y):
    return math.isclose(x, y, abs_tol=1e-12, rel_tol=0)


def verify(original: dict, frozen: dict) -> dict:
    if original.get("status") != "PUBLISHED_PRIOR_ART_POSTHOC_DESCRIPTIVE_NO_NEW_BUTTERFLY_FINDING":
        raise ValueError("Original result incorrectly promoted to new empirical finding")
    if original["source"]["original_git_blob_sha1"] != frozen["source"]["git_blob_sha1"]:
        raise ValueError("Scientific author Git blob has drifted")
    if original["source"]["original_csv_sha256"] != frozen["source"]["sha256"]:
        raise ValueError("Original raw source SHA256 has drifted")
    f = frozen["data_units"]
    u = original["source_units"]
    for key,source in (
        ("original_larvae_assigned","n_original_larvae"),
        ("distinct_host_plants","n_distinct_independent_host_plants"),
        ("experimental_temperature_plot_units","n_distinct_experimental_plots"),
        ("larvae_per_plant","larvae_per_host_plant"),
    ):
        if f[key] != u[source]:
            raise ValueError("Original plant/plot unit of analysis mismatch")
    if (u["original_temperature_plot_randomization_documented"] or
            frozen["data_units"]["original_plot_temperature_assignment_randomized_documented"] or
            not frozen["no_extrapolation"]["temperature_plot_randomization_not_explicitly_documented_in_authors_methods"]):
        raise ValueError("Undocumented thermal plot allocation cannot be called randomized")
    if set(u["plots_by_temperature"]) != {"ambient", "elevated"} or set(u["plots_by_temperature"].values()) != {15}:
        raise ValueError("Temperature treatment no longer represented by 15 plots")
    for key,source in (
        ("known_adult_fates","both_observed"),
        ("adult_eclosions","adult_1"),
        ("nonadult_known","adult_0"),
        ("unknown","ambiguous_stage_fate"),
        ("pupated_in_known_fates","pupa_1"),
    ):
        if f[key] != u["original_status"][source]:
            raise ValueError(f"Original life history fate count changed: {key}")
    if original["pupae_vs_adult_recorded_count_difference"] != f["postpupation_loss_among_known_pupae"]:
        raise ValueError("Pupal-to-adult original fate gap changed")
    if len(u["uncertain_original_larvae"]) != 1:
        raise ValueError("Not exactly one original uncertain adult")
    unknown = u["uncertain_original_larvae"][0]
    if unknown["original_larva_ID"] != f["unknown_individual_id"] or (
        unknown["Temp"] + "|" + unknown["Milkweed"] + "|" + unknown["OE_treatment"]
    ) != f["unknown_condition"]:
        raise ValueError("Original missing adult identity or randomized arm changed")
    eight=original["source_original_eight_factorial_cells"]
    if len(eight) != 8 or len(frozen["eight_cells"]) != 8:
        raise ValueError("Missing original factorial cell")
    for row in frozen["eight_cells"]:
        label="|".join([row["temp"],row["host"],row["inoc"]])
        if label not in eight:
            raise ValueError("Missing original host-temperature-OE treatment")
        for frozen_key, original_key in (
            ("n","larvae_originally_assigned"),("known","known_adult_outcomes"),
            ("adults","adult_eclosions"),("pupae","source_pupation_successes"),
        ):
            if row[frozen_key] != eight[label][original_key]:
                raise ValueError(f"Original 8-cell count mismatch: {label}, {frozen_key}")
    if sum(row["adults"] for row in frozen["eight_cells"]) != f["adult_eclosions"]:
        raise ValueError("Frozen pupation/adult count not conserved")
    for frozen_label,src_label in LABELS.items():
        desired=frozen["plot_paired_host_risk_difference"][frozen_label]
        result=original["primary_original_cohort_plot_stratified_missing_fate_bounds"][src_label]
        for key,source in COMPARISONS.items():
            if not close(result[source], desired[key]):
                raise ValueError(f"Plot-level result mismatch: {frozen_label}, {key}")
        ci=result["stratified_whole_plot_bootstrap_95_intervals"]
        for key,source in INTERVALS.items():
            if len(ci[source]) != 2 or any(
                not close(a,b) for a,b in zip(ci[source], desired[key])
            ):
                raise ValueError(f"Source plot bootstrap mismatch: {frozen_label}, {key}")
        if result["draws"] != frozen["plot_paired_host_risk_difference"]["draws"] or (
            result["seed"] != frozen["plot_paired_host_risk_difference"]["bootstrap_seed"]
        ):
            raise ValueError("Bootstrap randomization or number of plots changed")
    if not original["GEB_main_and_PR38_unmodified"] or not frozen["no_extrapolation"]["GEB_main_and_pr38_unchanged"]:
        raise ValueError("Unapproved promotion to GEB main")
    if not frozen["no_extrapolation"]["not_a_new_confirmatory_fitness_interaction_result"]:
        raise ValueError("Wrong original published result claim")
    return {
        "schema":"chocho_monarch_plotlevel_original_frozen_consistency_v01",
        "status":"AUTHOR_BLOB_AND_ALL_EIGHT_CELL_ADULT_FATES_AND_PLOT_BOOTSTRAP_REPRODUCED",
        "experiment_original_year":2021,
        "original_individual_assigned":240,
        "original_temperature_treatment_plot_units":30,
        "original_distinct_plants":120,
        "unknown_fates_preserved":1,
        "source_published_priors_not_new_biological_effect":True,
        "main_GEB_untouched":True,
    }



def verify_assigned_OE_plot_reproduction(original:dict, manifest:dict) -> dict:
    if manifest.get("status") != (
        "AUTHOR_ORIGINAL_PLOT_SOURCE_OE_POSTHOC_PLUS_INDEPENDENT_2019_PUBLISHED_PHENOTYPE_NOT_NOVEL"
    ):
        raise ValueError("Invalid source status for OE and reproductive state comparison")
    source=manifest["original_Ragonese"]
    if original["source"]["original_git_blob_sha1"]!=source["source_git_blob"]:
        raise ValueError("OE manifest original author repository Git source mismatch")
    if original["source_units"]["n_original_larvae"]!=source["n_larvae"] or (
        original["source_units"]["n_distinct_experimental_plots"]!=source["n_original_plots"]
    ):
        raise ValueError("OE manifest is not on original 30-plot and 240-larva source")
    if not source["temperature_plot_location_random_assignment_not_explicitly_documented"]:
        raise ValueError("Unverified temperature randomization cannot be asserted")
    if not manifest["frozen_OE_contrasts"]["no_OE_diagnosis_based_posttreatment_selection"]:
        raise ValueError("Post-randomization infection diagnoses cannot define arms")
    for case, original_label in LABELS.items():
        # Unlike 2019 adult reproductive status, assigned OE branches share
        # original Ragonese plot IDs, so their 95% intervals must be
        # resampled together within temperature and checked in the source.
        stored=manifest["frozen_OE_contrasts"][original_label]
        audited=(original["primary_original_cohort_plot_stratified_missing_fate_bounds"]
                 [original_label]["OE_inoculation_stratified_posthoc_sensitivity"])
        if audited["plot_resampling_units_per_temperature"]!=15 or (
            audited["bootstrap_draws"]!=manifest["frozen_OE_contrasts"]["n_resamples"]
            or audited["seed"]!=manifest["frozen_OE_contrasts"]["seed"]
        ):
            raise ValueError("Original OE plot/seed matching failed")
        for key in ("infected","control","difference_control_minus_infected"):
            audited_value=(
                audited["observed_host_by_warming_interaction_by_assigned_OE"][key]
                if key!="difference_control_minus_infected"
                else audited["difference_in_host_warming_interactions_control_minus_infected"]
            )
            if not close(audited_value,stored[key]["point"]):
                raise ValueError(f"Wrong source OE-stratified effect: {case} {key}")
            bootstrap=audited["paired_whole_plot_bootstrap_95pct"][key]
            if len(bootstrap)!=2 or any(
                not close(actual,frozen)
                for actual,frozen in zip(bootstrap,stored[key]["interval_95"])
            ):
                raise ValueError(f"Wrong source OE plot cluster interval: {case} {key}")
            if not (bootstrap[0]<=0<=bootstrap[1]):
                raise ValueError("OE interaction interval was promoted to significance")
    prior=manifest["published_reproductive_state_prior_art"]
    if prior["article_doi"]!="10.3390/insects10080253" or (
        prior["field_diet_experiment"]["adult_reproductive_status_analyzed"]
        ==prior["field_diet_experiment"]["adult_eclosions"]
    ):
        raise ValueError("2019 263 analyzed adult subset cannot be merged with 295 emerged")
    if not prior["cannot_transfer_reproductive_fraction_to_Ragonese2021_adult_outcomes"]:
        raise ValueError("Cross-source nonlinked reproductive status falsely multiplied")
    if not manifest["main_GEB_PR38_unchanged"]:
        raise ValueError("Scientific audit cannot silently revise GEB submission")
    return {
        "schema":"chocho_monarch_original_OE_strata_reproduction_gate_verified_v01",
        "status":"AUTHOR_SOURCE_OE_PLOT_CONTRAST_EXACT_WITH_CROSS_STUDY_ENDPOINT_STOP",
        "original_BLOB_verified":True,
        "both_2019_and_2021_sources_cannot_be_multiplied":True,
        "original_assigned_plot_OE_triple_difference_not_novel":True,
        "GEB_submission_unmodified":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original",required=True,type=Path)
    p.add_argument("--frozen",required=True,type=Path)
    p.add_argument("--oe-reproductive-manifest",type=Path)
    a=p.parse_args()
    orig=json.loads(a.original.read_text(encoding="utf-8"))
    frozen=json.loads(a.frozen.read_text(encoding="utf-8"))
    results={"original_adult_plot":verify(orig,frozen)}
    if a.oe_reproductive_manifest is not None:
        manifest=json.loads(a.oe_reproductive_manifest.read_text(encoding="utf-8"))
        results["assigned_OE_and_2019_endpoint_stop"]=verify_assigned_OE_plot_reproduction(orig,manifest)
    print(json.dumps(results,indent=2))


if __name__=="__main__":
    main()
