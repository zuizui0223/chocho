#!/usr/bin/env python3
"""Analyze a FUTURE six-arm randomized whole-cage butterfly experiment.

No field observations are bundled. ITT uses original recipients only, while
originally assigned additional A or S larvae have separate mandatory fates.
Unknown and living-diapause pupae stay in denominators and are neither coded as
confirmed death nor automatically as eventual flight-capable adults.

Primary planned contrast (original 2023 biological puzzle):
  (SUPPORTED_ADD_S - SUPPORTED_NONE) -
  (NATURAL_ADD_S   - NATURAL_NONE)

Key secondary contrast (is any rescue specific to the added species?):
  (SUPPORTED_ADD_S - SUPPORTED_ADD_A) -
  (NATURAL_ADD_S   - NATURAL_ADD_A)

Both estimate effects of assigned TREATMENT STRATEGIES, not pure leaf biomass
or direct heterospecific interference. Added A eats more per larva than
added S in the original 2023 trials: equal headcount is not equal feeding.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from randomize_kyoto_six_arm_identity_food import ARMS, ALLOCATION_COLUMNS, DESIGN

ALLOCS = list(ALLOCATION_COLUMNS)
INFO = [
    "cage_id", "block_id", "treatment", "food_support", "added_group",
    "recipient_a_initial_n", "added_larvae_initial_n", "food_floor_cm2",
    "expected_n_resource_visits", "donor_source_batch_id", "study_season",
]
FATES = ["cage_id", "larva_id", "cohort_role", "species", "fate"]
VISITS = [
    "cage_id", "observation_datetime", "recipient_a_alive_n",
    "recipient_a_late_instar_alive_n", "added_a_alive_n", "added_s_alive_n",
    "accessible_leaf_area_cm2", "edible_fresh_mass_g",
    "accessible_leaf_area_added_cm2", "leaf_stage", "added_leaf_stage",
    "leaf_source_batch_id", "sham_handling_done",
]
STAGES = {"young", "mature", "mixed", "none", "unknown"}
SUCCESS = "flight_capable_adult"
FAILURE = {"dead_immature", "pupa_confirmed_dead", "adult_nonflight"}
UNRESOLVED = {"lost_to_followup", "other_documented",
              "diapause_alive", "pupa_alive_pending_emergence"}
FATE_VALUES = {SUCCESS} | FAILURE | UNRESOLVED
NAME_A = "Atrophaneura alcinous"
NAME_S = "Sericinus montela"

# Within-block contrast coefficients.
COEFFS = {
    "identity_x_food_secondary": {
        "SUPPORTED_ADD_S": +1, "SUPPORTED_ADD_A": -1,
        "NATURAL_ADD_S": -1, "NATURAL_ADD_A": +1,
    },
    "identity_natural_S_minus_A": {"NATURAL_ADD_S": +1, "NATURAL_ADD_A": -1},
    "identity_supported_S_minus_A": {"SUPPORTED_ADD_S": +1, "SUPPORTED_ADD_A": -1},
    "S_competitor_natural": {"NATURAL_ADD_S": +1, "NATURAL_NONE": -1},
    "S_competitor_supported": {"SUPPORTED_ADD_S": +1, "SUPPORTED_NONE": -1},
    "A_competitor_natural": {"NATURAL_ADD_A": +1, "NATURAL_NONE": -1},
    "A_competitor_supported": {"SUPPORTED_ADD_A": +1, "SUPPORTED_NONE": -1},
    "S_competitor_x_food": {
        "SUPPORTED_ADD_S": +1, "SUPPORTED_NONE": -1,
        "NATURAL_ADD_S": -1, "NATURAL_NONE": +1,
    },
    "A_competitor_x_food": {
        "SUPPORTED_ADD_A": +1, "SUPPORTED_NONE": -1,
        "NATURAL_ADD_A": -1, "NATURAL_NONE": +1,
    },
}
PRIMARY = "S_competitor_x_food"


def csv_rows(path: Path, columns: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != columns:
            raise ValueError(f"{path}: exact columns required: {columns}")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError(f"{path}: empty or malformed source rows")
    return rows


def numeric(v: str, field: str, *, missing: bool = False) -> float | None:
    if str(v or "").strip() == "" and missing:
        return None
    try:
        num = float(v)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field}: missing or invalid numeric value") from exc
    if not math.isfinite(num) or num < 0:
        raise ValueError(f"{field}: cannot be nonfinite or negative")
    return num


def int_value(v: str, field: str, *, positive: bool = False) -> int:
    x = numeric(v, field)
    if x != int(x) or (positive and x == 0):
        raise ValueError(f"{field}: expected nonnegative integer" + (" >0" if positive else ""))
    return int(x)


def assignment(allocation: list[dict], info: list[dict]) -> tuple[dict, dict]:
    alloc = {}
    for row in allocation:
        cage = row["cage_id"].strip()
        if not cage or cage in alloc:
            raise ValueError("Duplicate or blank original randomized cage_id")
        if row["treatment"] not in ARMS:
            raise ValueError("Unknown six-arm treatment")
        expected = ARMS[row["treatment"]]
        if (int_value(row["food_support"], "randomized food_support"),
                row["added_group"]) != expected:
            raise ValueError("Randomized treatment factor inconsistent with arm")
        if row["design_version"] != DESIGN:
            raise ValueError("Wrong original six-arm design version")
        seed = row["random_seed"].strip()
        try:
            int(seed)
        except ValueError as exc:
            raise ValueError("No auditable allocation seed") from exc
        alloc[cage] = row
    if len({r["random_seed"] for r in allocation}) != 1:
        raise ValueError("Different random seeds inside one allocation receipt")
    cages: dict[str, dict] = {}
    for row in info:
        cage = row["cage_id"].strip()
        if not cage or cage in cages:
            raise ValueError("Duplicate or blank study cage_id")
        if cage not in alloc:
            raise ValueError("No original assignment for experimental cage")
        original = alloc[cage]
        for field in ("block_id", "treatment", "food_support", "added_group"):
            if row[field] != original[field]:
                raise ValueError("Post-randomization treatment or block change")
        focal = int_value(row["recipient_a_initial_n"], "recipient_a_initial_n", positive=True)
        added = int_value(row["added_larvae_initial_n"], "added_larvae_initial_n")
        group = row["added_group"]
        if (added == 0) != (group == "none"):
            raise ValueError("Assigned competitor identity and added headcount conflict")
        floor = numeric(row["food_floor_cm2"], "food_floor_cm2")
        if floor <= 0:
            raise ValueError("Predeclared edible food floor must be >0")
        nvisit = int_value(row["expected_n_resource_visits"], "expected_n_resource_visits", positive=True)
        if nvisit < 2:
            raise ValueError("Require repeated food-access observations")
        if not row["donor_source_batch_id"].strip() or not row["study_season"].strip():
            raise ValueError("Missing source-donor provenance or experimental season")
        cages[cage] = {
            **row, "recipient_a_initial_n": focal,
            "added_larvae_initial_n": added, "food_floor_cm2": floor,
            "expected_n_resource_visits": nvisit,
        }
    if set(cages) != set(alloc):
        raise ValueError("Missing randomized whole cage cannot be dropped from analysis")
    blocks = defaultdict(list)
    for cage in cages.values():
        blocks[cage["block_id"]].append(cage)
    for block, rows in blocks.items():
        count = Counter(r["treatment"] for r in rows)
        if set(count) != set(ARMS) or len(set(count.values())) != 1:
            raise ValueError(f"{block}: all six arms must be balanced and retained")
        for name in ("recipient_a_initial_n", "expected_n_resource_visits",
                     "food_floor_cm2", "study_season", "donor_source_batch_id"):
            if len(set(r[name] for r in rows)) != 1:
                raise ValueError(f"{block}: {name} must be assigned equally within source block")
        added_n = {r["added_larvae_initial_n"] for r in rows if r["added_group"] != "none"}
        if len(added_n) != 1 or next(iter(added_n)) < 1:
            raise ValueError(f"{block}: added A and S must have same prespecified HEADCOUNT")
    return cages, blocks


def fates_for_cages(cages: dict, fates: list[dict]) -> dict:
    by_cage = defaultdict(list)
    used_ids = set()
    for row in fates:
        cage, larva = row["cage_id"], row["larva_id"]
        if cage not in cages:
            raise ValueError("Fate from unknown randomized cage")
        if not larva or larva in used_ids:
            raise ValueError("Original larva IDs must be nonblank and globally unique")
        used_ids.add(larva)
        if row["cohort_role"] not in ("focal", "added"):
            raise ValueError("No reliable original focal/added cohort designation")
        if row["fate"] == "pupa_no_adult":
            raise ValueError("Ambiguous pupa_no_adult: distinguish living diapause from death")
        if row["fate"] not in FATE_VALUES:
            raise ValueError("Unknown adult/pupal fate, cannot complete original cohort")
        spec = NAME_A if (row["cohort_role"] == "focal"
                          or cages[cage]["added_group"] == "A") else NAME_S
        if row["species"] != spec:
            raise ValueError("Added species or focal larva identity does not match assignment")
        by_cage[cage].append(row)
    result = {}
    for cage, contract in cages.items():
        by_role = defaultdict(list)
        for row in by_cage.get(cage, []):
            by_role[row["cohort_role"]].append(row["fate"])
        if len(by_role["focal"]) != contract["recipient_a_initial_n"]:
            raise ValueError(f"{cage}: original assigned focal A cohort not completely followed")
        if len(by_role["added"]) != contract["added_larvae_initial_n"]:
            raise ValueError(f"{cage}: original additional larvae missing or reassigned")
        focal_fates = Counter(by_role["focal"])
        added_fates = Counter(by_role["added"])
        original_n = contract["recipient_a_initial_n"]
        adult = focal_fates[SUCCESS]
        pending = sum(focal_fates[x] for x in UNRESOLVED)
        result[cage] = {
            "original_focal_A": original_n,
            "original_added": contract["added_larvae_initial_n"],
            "observed_flight_capable_focal_adults": adult,
            "focal_lower_adult_fraction": adult / original_n,
            "focal_upper_eventual_adult_fraction": (adult + pending) / original_n,
            "focal_known_dead_or_nonflight": sum(focal_fates[x] for x in FAILURE),
            "focal_living_diapausing": focal_fates["diapause_alive"],
            "focal_living_nonemerged_pupae": focal_fates["pupa_alive_pending_emergence"],
            "focal_unresolved_other": sum(
                focal_fates[x] for x in ("lost_to_followup", "other_documented")),
            "focal_fate_unknown_for_eventual_emergence": pending,
            "added_fate_counts": dict(added_fates),
        }
    return result


def visits_for_cages(cages: dict, visits: list[dict]) -> dict:
    rows_by_cage = defaultdict(list)
    seen = set()
    for row in visits:
        cage = row["cage_id"]
        if cage not in cages:
            raise ValueError("Food observation outside original randomized cage")
        stamp = row["observation_datetime"]
        try:
            dt = datetime.fromisoformat(stamp)
        except ValueError as exc:
            raise ValueError("Food observations need ISO 8601 timestamp") from exc
        if dt.tzinfo is None:
            raise ValueError("Food observation timestamp missing timezone")
        key = (cage, dt.isoformat())
        if key in seen:
            raise ValueError("Repeated observation at same cage/date")
        seen.add(key)
        c = cages[cage]
        alive = int_value(row["recipient_a_alive_n"], "recipient_a_alive_n")
        late = int_value(row["recipient_a_late_instar_alive_n"], "recipient_a_late_instar_alive_n")
        extra_a = int_value(row["added_a_alive_n"], "added_a_alive_n")
        extra_s = int_value(row["added_s_alive_n"], "added_s_alive_n")
        if not (0 <= late <= alive <= c["recipient_a_initial_n"]):
            raise ValueError("Late-instar or total focal A count exceeds original randomized cohort")
        if extra_a + extra_s > c["added_larvae_initial_n"]:
            raise ValueError("Extra larva count exceeds originally allocated additions")
        if (c["added_group"] != "A" and extra_a) or (c["added_group"] != "S" and extra_s):
            raise ValueError("Added larva species does not correspond to original assigned treatment")
        leaf_area = numeric(row["accessible_leaf_area_cm2"], "accessible_leaf_area_cm2", missing=True)
        mass = numeric(row["edible_fresh_mass_g"], "edible_fresh_mass_g", missing=True)
        added_area = numeric(row["accessible_leaf_area_added_cm2"], "accessible_leaf_area_added_cm2", missing=True)
        leaf_stage = row["leaf_stage"].strip().lower()
        added_stage = row["added_leaf_stage"].strip().lower()
        if leaf_stage not in STAGES or (added_stage and added_stage not in STAGES):
            raise ValueError("Unrecognized young/mature/mixed leaf developmental state")
        if added_area is not None and added_area > 0 and not added_stage:
            raise ValueError("Food supplementation needs contemporaneous added_leaf_stage")
        if added_area is not None and added_area > 0 and c["food_support"] == "0":
            raise ValueError("Natural-access arm must not be given accessible added foliage")
        if row["sham_handling_done"] not in ("0", "1"):
            raise ValueError("Record whether matched sham handling was performed")
        if not row["leaf_source_batch_id"].strip():
            raise ValueError("Source of actual foliage not recorded")
        rows_by_cage[cage].append({
            "stamp": dt.isoformat(), "alive": alive, "late": late,
            "extra_a": extra_a, "extra_s": extra_s,
            "area": leaf_area, "fresh_mass": mass, "added_area": added_area,
            "leaf_stage": leaf_stage, "added_stage": added_stage,
            "donor_batch": row["leaf_source_batch_id"],
            "sham": row["sham_handling_done"] == "1",
        })
    out = {}
    for cage, c in cages.items():
        series = rows_by_cage.get(cage, [])
        if len(series) != c["expected_n_resource_visits"]:
            raise ValueError(f"{cage}: missing required scheduled resource observation")
        observed = [s for s in series if s["area"] is not None]
        late = [s for s in series if s["late"] > 0]
        late_observed = [s for s in late if s["area"] is not None]
        floor = c["food_floor_cm2"]
        out[cage] = {
            "scheduled_visits": len(series),
            "known_accessible_area_visits": len(observed),
            "area_missing_visits": len(series) - len(observed),
            "food_floor_breach_visits": sum(s["area"] < floor for s in observed),
            "late_living_focal_stage_visits": len(late),
            "late_living_focal_with_area": len(late_observed),
            "late_living_focal_food_floor_breach_visits": sum(
                s["area"] < floor for s in late_observed),
            "late_stage_food_area_missing": len(late) - len(late_observed),
            "visits_fresh_mass_observed": sum(s["fresh_mass"] is not None for s in series),
            "total_accessible_area_added_cm2": sum(
                s["added_area"] or 0 for s in series),
            "added_tissue_age_counts": dict(Counter(
                s["added_stage"] for s in series
                if s["added_area"] is not None and s["added_area"] > 0)),
            "sham_handling_reported_visits": sum(s["sham"] for s in series),
            "donor_batch_ids_used": sorted({s["donor_batch"] for s in series}),
            "interpretation": (
                "Postrandomization diagnostics ONLY. Food floor attained is NOT an "
                "eligibility criterion; low survivor numbers can themselves cause "
                "high leaf availability; late-instar observations condition on "
                "survival and cannot prove previous absence of food shortage."
            ),
        }
    return out


def _risk_by_arm(block_rows: list[dict], fates: dict, field: str) -> dict[str, float]:
    return {
        arm: sum(fates[r["cage_id"]][field] for r in block_rows if r["treatment"] == arm)
             / sum(r["treatment"] == arm for r in block_rows)
        for arm in ARMS
    }


def _contrasts(risks: dict[str, float]) -> dict[str, float]:
    return {name: sum(multiplier * risks[arm] for arm, multiplier in weights.items())
            for name, weights in COEFFS.items()}


def _missing_bounds(low: dict, high: dict) -> dict:
    # For every signed arm coefficient, choose lower for + and upper for -
    # for the lower bound; reverse for upper. Never blindly count all missing
    # cohorts as successes for a difference-in-differences upper bound.
    return {
        name: [
            sum(v * (low[arm] if v > 0 else high[arm]) for arm, v in coefs.items()),
            sum(v * (high[arm] if v > 0 else low[arm]) for arm, v in coefs.items()),
        ]
        for name, coefs in COEFFS.items()
    }


def _mean(vals: list[float]) -> float:
    if not vals:
        raise ValueError("No source block contrast")
    return sum(vals) / len(vals)


def _block_bootstrap(rows: list[dict], *, seed: int = 20261014, draws: int = 3000):
    # A reproducible descriptive cluster/bootstrap interval, NOT a design-
    # based power guarantee. A weakly clustered pilot must withhold this.
    if len(rows) < 10:
        return None
    rng = random.Random(seed)
    samples = []
    for _ in range(draws):
        samples.append(_mean([
            rows[rng.randrange(len(rows))]["contrasts"][PRIMARY]
            for __ in range(len(rows))
        ]))
    samples.sort()
    return [samples[int(draws * .025)], samples[int(draws * .975)]]


def analyze(allocation: list[dict], info: list[dict], fates: list[dict],
            visits: list[dict]) -> dict:
    cages, blocks = assignment(allocation, info)
    ends = fates_for_cages(cages, fates)
    access = visits_for_cages(cages, visits)
    by_block = []
    for block, members in sorted(blocks.items()):
        risks_low = _risk_by_arm(members, ends, "focal_lower_adult_fraction")
        risks_high = _risk_by_arm(members, ends, "focal_upper_eventual_adult_fraction")
        b = _missing_bounds(risks_low, risks_high)
        by_block.append({
            "block_id": block,
            "n_independent_cages": len(members),
            "risks_known_adult_by_arm": risks_low,
            "risks_eventual_adult_upper_by_arm": risks_high,
            "contrasts": _contrasts(risks_low),
            "eventual_adult_attrition_bounds": b,
        })
    contrasts = {
        name: _mean([r["contrasts"][name] for r in by_block]) for name in COEFFS
    }
    worst = {
        name: [
            _mean([r["eventual_adult_attrition_bounds"][name][0] for r in by_block]),
            _mean([r["eventual_adult_attrition_bounds"][name][1] for r in by_block]),
        ]
        for name in COEFFS
    }
    if not all(worst[k][0] - 1e-12 <= contrasts[k] <= worst[k][1] + 1e-12
               for k in COEFFS):
        raise RuntimeError("Signed missing-fate bounds do not enclose observed contrasts")
    per_arm = {}
    for arm in ARMS:
        ids = [r["cage_id"] for r in cages.values() if r["treatment"] == arm]
        n = sum(ends[i]["original_focal_A"] for i in ids)
        k = sum(ends[i]["observed_flight_capable_focal_adults"] for i in ids)
        pending = sum(ends[i]["focal_fate_unknown_for_eventual_emergence"] for i in ids)
        late = sum(access[i]["late_living_focal_with_area"] for i in ids)
        floor_breach = sum(access[i]["late_living_focal_food_floor_breach_visits"] for i in ids)
        per_arm[arm] = {
            "randomized_whole_cages": len(ids),
            "original_focal_larvae": n,
            "verified_flight_capable_adults": k,
            "living_diapausing_focal_pupae": sum(ends[i]["focal_living_diapausing"] for i in ids),
            "other_living_focal_pupae": sum(ends[i]["focal_living_nonemerged_pupae"] for i in ids),
            "unresolved_focal_fates": pending,
            "observed_eventual_adult_probability_lower": k / n,
            "eventual_adult_probability_upper": (k + pending) / n,
            "food_access_occasion_count": sum(access[i]["scheduled_visits"] for i in ids),
            "missing_leaf_area_visits": sum(access[i]["area_missing_visits"] for i in ids),
            "late_focal_area_observed_visits": late,
            "late_focal_food_floor_breach_visits": floor_breach,
            "late_focal_floor_breach_fraction_observed": floor_breach / late if late else None,
            "sham_handling_reported": sum(access[i]["sham_handling_reported_visits"] for i in ids),
            "accessible_supplemented_leaf_area_cm2": sum(
                access[i]["total_accessible_area_added_cm2"] for i in ids),
        }
    return {
        "schema": "chocho_kyoto_six_arm_identity_x_food_response_v01",
        "data_state": "USER_SUPPLIED_FUTURE_ANIMAL_OUTCOMES_ONLY_NO_RESULT_BUNDLED",
        "randomization_unit": "one whole cage/enclosure, not plants or larvae within it",
        "initial_recipient_cohort": NAME_A,
        "n_randomized_cages": len(cages),
        "n_source_blocks": len(blocks),
        "per_arm": per_arm,
        "block_level": by_block,
        "primary_planned_estimand": PRIMARY,
        "primary_observed_flight_capable_adult_risk_difference": contrasts[PRIMARY],
        "other_planned_contrasts": {k: v for k, v in contrasts.items() if k != PRIMARY},
        "worst_case_eventual_emergence_attrition_intervals": worst,
        "independent_block_bootstrap_primary_interval_if_10plus_blocks": _block_bootstrap(by_block),
        "all_original_focal_fates_resolved": all(
            ends[c]["focal_fate_unknown_for_eventual_emergence"] == 0 for c in cages),
        "all_original_added_larvae_accounted_for": True,
        "extra_larvae_have_same_HEADCOUNT_not_same_consumption": True,
        "primary_treatments_are_intention_to_treat": True,
        "primary_DOES_NOT_condition_on_achieved_food_availability": True,
        "primary_scientific_question": (
            "Does assigned food support alleviate the TOTAL intervention effect "
            "of adding Sericinus to the original Atrophaneura recipient cohort?"
        ),
        "secondary_scientific_question": (
            "Is any support modification different when the same number "
            "of added larvae are conspecific Atrophaneura instead?"
        ),
        "conspecific_vs_heterospecific_identity_not_consumption_matched": True,
        "resource_shortage_vs_quality_vs_direct_effect_disentangled": False,
        "living_diapausing_pupae_are_not_dead": True,
        "consumer_species_introduced_host_plant_is_native": True,
        "external_GEB_main_not_modified": True,
        "interpretive_stop": (
            "This analysis cannot label the support strategy a pure leaf-mass "
            "intervention, infer physiological tolerance from host consumption, "
            "or claim an introduced PLANT effect. Bootstrap, if returned, is "
            "only a block-level descriptive sensitivity and NOT power evidence."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for field in ("allocation", "assignments", "fates", "resource", "receipt"):
        p.add_argument("--" + field, type=Path, required=True)
    args = p.parse_args()
    paths = {k: getattr(args, k) for k in ("allocation", "assignments", "fates", "resource")}
    input_data = [
        csv_rows(paths["allocation"], ALLOCS),
        csv_rows(paths["assignments"], INFO),
        csv_rows(paths["fates"], FATES),
        csv_rows(paths["resource"], VISITS),
    ]
    receipt = analyze(*input_data)
    receipt["input_file_sha256"] = {
        label: hashlib.sha256(path.read_bytes()).hexdigest()
        for label, path in paths.items()
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
