#!/usr/bin/env python3
"""Source-blind checks for a future Kyoto leaf-stage pilot.

No butterfly, leaf, or facility measurements are included in this repository.
Only user-provided authentic measurements are accepted. This is a feasibility
and data-integrity audit, NOT proof of physiological accessibility or fitness.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

HOST = "Aristolochia debilis"
SPECIES = {"Atrophaneura alcinous", "Sericinus montela"}
STAGES = {"young", "mature"}
LEAVES = [
    "donor_plant_id", "leaf_batch_id", "host_species", "leaf_stage",
    "sample_date", "shoot_node", "whole_plant_age_years",
    "leaf_area_cm2", "fresh_mass_g", "dry_mass_g",
]
ASSAYS = [
    "cup_id", "source_block_id", "maternal_family_id", "species",
    "instar_class", "assigned_leaf_stage", "source_donor_plant_id",
    "assay_leaf_batch_id", "initial_larvae_n", "alive_end_n",
    "observation_minutes", "feeding_confirmed",
]
FACILITIES = ["greenhouse_id", "eligible_whole_cages"]


def read_rows(path: Path, expected: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected:
            raise ValueError(f"{path}: require exact columns {expected}")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError(f"{path}: empty/invalid source rows")
    return rows


def number(value, field, *, integer=False, zero=False, blank=False):
    if str(value or "").strip() == "" and blank:
        return None
    try:
        n = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{field}: missing or invalid measurement") from exc
    if not math.isfinite(n) or n < 0 or (not zero and n == 0):
        raise ValueError(f"{field}: invalid nonpositive/nonfinite measurement")
    if integer and n != int(n):
        raise ValueError(f"{field}: requires a whole count")
    return int(n) if integer else n


def leaf_calibration(rows: list[dict[str, str]]) -> dict:
    batch_ids = set()
    donor_stage: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    dates = set()
    for row in rows:
        donor, batch = row["donor_plant_id"].strip(), row["leaf_batch_id"].strip()
        stage = row["leaf_stage"].strip().lower()
        if not donor or not batch or batch in batch_ids:
            raise ValueError("blank donor, duplicate or blank calibration leaf-batch ID")
        batch_ids.add(batch)
        if row["host_species"] != HOST or stage not in STAGES:
            raise ValueError("calibration requires original Aristolochia debilis young/mature leaves")
        try:
            sampled = date.fromisoformat(row["sample_date"])
        except ValueError as exc:
            raise ValueError("original leaf date must be ISO YYYY-MM-DD") from exc
        dates.add(sampled.isoformat())
        if str(row["shoot_node"]).strip():
            number(row["shoot_node"], "shoot_node", integer=True)
        number(row["whole_plant_age_years"], "whole_plant_age_years", blank=True)
        area = number(row["leaf_area_cm2"], "leaf_area_cm2")
        wet = number(row["fresh_mass_g"], "fresh_mass_g")
        dry = number(row["dry_mass_g"], "dry_mass_g")
        if dry > wet + 1e-9:
            raise ValueError("calibrated dry tissue mass exceeds original fresh mass")
        donor_stage[donor][stage].append({
            "dry_fraction": dry / wet,
            "dry_mg_per_cm2": dry * 1000 / area,
            "wet_mg_per_cm2": wet * 1000 / area,
        })
    matched = [d for d, by in donor_stage.items() if set(by) == STAGES]
    def mean(records, key):
        return statistics.mean([x[key] for x in records])
    donor_summaries = {
        donor: {
            stage: {
                "n_original_leaf_batches": len(recs),
                "mean_dry_fraction": mean(recs, "dry_fraction"),
                "mean_dry_mg_per_cm2": mean(recs, "dry_mg_per_cm2"),
            }
            for stage, recs in sorted(by.items())
        }
        for donor, by in sorted(donor_stage.items())
    }
    stage_summaries = {}
    for stage in sorted(STAGES):
        donors_here = [donor_summaries[d][stage] for d in donor_summaries if stage in donor_summaries[d]]
        stage_summaries[stage] = {
            "independent_donor_plants": len(donors_here),
            "donor_equal_mean_dry_fraction": (
                statistics.mean(x["mean_dry_fraction"] for x in donors_here) if donors_here else None
            ),
            "donor_equal_mean_dry_mg_per_cm2": (
                statistics.mean(x["mean_dry_mg_per_cm2"] for x in donors_here) if donors_here else None
            ),
        }
    paired_diff = [
        donor_summaries[d]["young"]["mean_dry_fraction"] -
        donor_summaries[d]["mature"]["mean_dry_fraction"]
        for d in matched
    ]
    return {
        "n_original_leaf_batches": len(rows),
        "independent_original_donor_plants": len(donor_stage),
        "donors_paired_young_and_mature": len(matched),
        "all_matched_original_donor_plant_ids": sorted(matched),
        "n_original_sampling_dates": len(dates),
        "by_stage": stage_summaries,
        "paired_within_donor_young_minus_mature_dry_fraction_mean": (
            statistics.mean(paired_diff) if paired_diff else None
        ),
        "donor_stage_records": donor_summaries,
        "structural_matched_donor_gate": len(matched) >= 2,
        "not_a_claim_that_two_donors_support_statistical_generality": True,
        "not_a_dry_mass_measurement_from_2016_csv": True,
    }


def feeding_assay(rows: list[dict[str, str]], calibration: dict) -> dict:
    donor_ids = set(calibration["donor_stage_records"])
    seen_cups, seen_foliage = set(), set()
    arms: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        cup, foliage = row["cup_id"].strip(), row["assay_leaf_batch_id"].strip()
        family = row["maternal_family_id"].strip()
        block = row["source_block_id"].strip()
        donor = row["source_donor_plant_id"].strip()
        stage = row["assigned_leaf_stage"].strip().lower()
        species = row["species"].strip()
        instar = row["instar_class"].strip()
        if not cup or cup in seen_cups or not foliage or foliage in seen_foliage:
            raise ValueError("cup and fresh assay leaf-batch IDs must be unique")
        seen_cups.add(cup)
        seen_foliage.add(foliage)
        if not block or not family or donor not in donor_ids:
            raise ValueError("missing source block/maternal family/unknown plant donor")
        if stage not in STAGES or stage not in calibration["donor_stage_records"][donor]:
            raise ValueError("assay foliage stage is uncalibrated for its original donor")
        if species not in SPECIES or instar not in {"L1", "L3"}:
            raise ValueError("species/larval instar not in the protocol")
        n = number(row["initial_larvae_n"], "initial_larvae_n", integer=True)
        alive = number(row["alive_end_n"], "alive_end_n", integer=True, zero=True)
        if alive > n:
            raise ValueError("more alive larvae than originally assigned")
        mins = number(row["observation_minutes"], "observation_minutes")
        observed = row["feeding_confirmed"].strip().lower()
        if observed not in {"yes", "no", "unknown"}:
            raise ValueError("feeding_confirmed must be yes/no/unknown, not inferred from mortality")
        arms[f"{species}|{instar}|{stage}"].append({
            "source_block_id": block,
            "maternal_family_id": family,
            "source_donor_plant_id": donor,
            "started": n, "alive": alive, "minutes": mins, "observed": observed,
        })
    out = {}
    for name, records in sorted(arms.items()):
        eligible = [r for r in records if r["observed"] != "unknown"]
        out[name] = {
            "whole_assay_cups": len(records),
            "independent_maternal_family_ids": len({r["maternal_family_id"] for r in records}),
            "independent_donor_plant_ids": len({r["source_donor_plant_id"] for r in records}),
            "original_larvae_started": sum(r["started"] for r in records),
            "alive_at_declared_observation": sum(r["alive"] for r in records),
            "feeding_yes": sum(r["observed"] == "yes" for r in eligible),
            "feeding_no_after_prescribed_observation": sum(r["observed"] == "no" for r in eligible),
            "feeding_unknown": sum(r["observed"] == "unknown" for r in records),
            "cup_feeding_confirmed_fraction_among_observed": (
                sum(r["observed"] == "yes" for r in eligible) / len(eligible) if eligible else None
            ),
            "not_an_individual_level_or_population_level_survival_effect": True,
        }
    mandatory = {f"{species}|L1|{stage}" for species in SPECIES for stage in STAGES}
    milestone = mandatory.issubset(out) and all(
        out[name]["independent_maternal_family_ids"] >= 2
        and out[name]["independent_donor_plant_ids"] >= 2
        for name in mandatory
    )
    return {
        "by_species_instar_leaf_stage": out,
        "required_neonate_2species_by_2leafstage_cells": sorted(mandatory),
        "four_core_neonate_cells_present_with_minimum_source_repetition": milestone,
        "not_pilot_power_or_ecological_confirmatory_inference": True,
        "unknown_feeding_is_not_reclassified_as_no": True,
        "same_source_cups_not_independent_maternal_families": True,
    }


def facility_audit(rows: list[dict[str, str]]) -> dict:
    seen = set()
    cap = {}
    for row in rows:
        greenhouse = row["greenhouse_id"].strip()
        if not greenhouse or greenhouse in seen:
            raise ValueError("duplicate/blank greenhouse inventory")
        seen.add(greenhouse)
        cap[greenhouse] = number(
            row["eligible_whole_cages"], "eligible_whole_cages", integer=True, zero=True
        )
    return {
        "independent_greenhouse_ids": len(cap),
        "original_whole_cage_capacity_by_greenhouse": dict(sorted(cap.items())),
        "current_eligible_total_whole_cages": sum(cap.values()),
        "greenhouses_supporting_at_least_four_cages": sum(v >= 4 for v in cap.values()),
        "greenhouses_supporting_at_least_six_cages": sum(v >= 6 for v in cap.values()),
        "greenhouse_local_complete_four_arm_blocks_capacity": sum(v // 4 for v in cap.values()),
        "greenhouse_local_complete_six_arm_blocks_capacity": sum(v // 6 for v in cap.values()),
        "greenhouse_incomplete_block_design_not_analyzed_here": True,
        "not_a_claim_about_facility_currently_available_without_real_inventory": True,
    }


def audit(leaf, assay, facilities):
    c = leaf_calibration(leaf)
    a = feeding_assay(assay, c)
    f = facility_audit(facilities)
    return {
        "schema": "chocho_kyoto_leaf_stage_pilot_feasibility_v01",
        "status": "REAL_INPUT_STRUCTURAL_AUDIT_ONLY_NO_FITNESS_CAUSAL_ESTIMATE",
        "donor_stage_mass_calibration": c,
        "no_choice_neonate_feeding_pilot": a,
        "verified_current_facility_capacity": f,
        "operational_gate": (
            "STRUCTURAL_PILOT_SUPPORTED_NOT_STUDY_APPROVED"
            if c["structural_matched_donor_gate"] and
               a["four_core_neonate_cells_present_with_minimum_source_repetition"]
            else "PILOT_INPUT_INCOMPLETE_OR_INSUFFICIENT"
        ),
        "interpretation": (
            "Leaf stage and species-by-instar feeding are experimental exposure "
            "screening variables, NOT independently observed adult reproduction. "
            "Maternal families/donor plants and future cage/greenhouse structure "
            "must be replicated beyond this structural gate."
        ),
        "true_demographic_fitness_response_estimated": False,
        "preliminary_2023_neonate_study_reproduced": False,
        "butterflies_reared_by_script": False,
        "do_not_release_exotic_butterflies": True,
        "GEB_main_unchanged": True,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--leaf-calibration", type=Path, required=True)
    p.add_argument("--feeding-assay", type=Path, required=True)
    p.add_argument("--facility-inventory", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    args = p.parse_args()
    paths = {
        "leaf_calibration": args.leaf_calibration,
        "feeding_assay": args.feeding_assay,
        "facility_inventory": args.facility_inventory,
    }
    output = audit(
        read_rows(args.leaf_calibration, LEAVES),
        read_rows(args.feeding_assay, ASSAYS),
        read_rows(args.facility_inventory, FACILITIES),
    )
    output["source_input_sha256"] = {
        name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
