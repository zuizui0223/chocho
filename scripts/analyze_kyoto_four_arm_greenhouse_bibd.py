#!/usr/bin/env python3
"""Future Kyoto *four-arm, greenhouse-fixed-effects* BIBD causal ITT analysis.

Eight independent greenhouses, 3 whole cages per greenhouse, 4 treatments,
each treatment seen 6 times, every pair co-occurring in 4 greenhouses.
Unlike the existing 4-arm COMPLETE-block reader this supports the original
published 3-cage greenhouse footprint without pseudoreplication.

Each initially assigned focal Atrophaneura individual is retained through
observed adult flight, confirmed death, diapause, pending pupa or unresolved.
Original added Sericinus cohorts are tracked separately. Post-treatment leaf
availability is a diagnostic, never a reason to select "successful" cages.

The primary is (supported S - supported none) - (natural S - natural none).
The OLS fixed-greenhouse descriptive interval uses an illustrative
homoscedastic independent cage residual assumption; it is NOT a measured
fitness effect without actual source data, nor a guaranteed precision level.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

from randomize_kyoto_four_arm_greenhouse_bibd import (
    ALLOCATION_COLUMNS, ARMS, DESIGN, B, R, LAMBDA, counts,
)

ASSIGNMENT_COLUMNS = [
    "cage_id", "greenhouse_id", "treatment", "competitor_present",
    "food_support", "focal_a_initial_n", "added_s_initial_n",
    "expected_n_food_visits", "food_floor_cm2", "donor_batch_id",
]
FATE_COLUMNS = ["cage_id", "larva_id", "species", "cohort_role", "fate"]
RESOURCE_COLUMNS = [
    "cage_id", "observation_datetime", "focal_a_alive_n",
    "focal_a_late_instar_alive_n", "added_s_alive_n",
    "accessible_edible_leaf_area_cm2", "accessible_leaf_added_area_cm2",
    "leaf_development_stage", "added_leaf_stage", "leaf_source_batch_id",
    "sham_handling_done",
]
FOCAL = "Atrophaneura alcinous"
ADDED = "Sericinus montela"
SUCCESS = "flight_capable_adult"
FAIL = {"dead_immature", "pupa_confirmed_dead", "adult_nonflight"}
PENDING = {"diapause_alive", "pupa_alive_pending_emergence",
           "lost_to_followup", "other_documented"}
VALID_FATES = {SUCCESS} | FAIL | PENDING
VALID_LEAF = {"young", "mature", "mixed", "none", "unknown"}
T_CRITICAL_DF13 = 2.1603686564610127


def read_csv(path: Path, schema: list[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != schema:
            raise ValueError(f"Original data must have exact columns {schema}")
        rows = list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError(f"Empty/invalid CSV {path}")
    return rows


def nonnegative(v, name: str, *, integer=False, positive=False, missing=False):
    if str(v or "").strip() == "" and missing:
        return None
    try:
        x = float(v)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}: missing/non-numeric original data") from exc
    if not math.isfinite(x) or x < 0 or (positive and x <= 0):
        raise ValueError(f"{name}: negative/nonfinite/nonpositive")
    if integer and int(x) != x:
        raise ValueError(f"{name}: non-integer count")
    return int(x) if integer else x


def validate_assignment(allocation: list[dict], info: list[dict]) -> dict:
    if len(allocation) != B * 3:
        raise ValueError("Required 24 original whole-cage assignments")
    if len({r["cage_id"] for r in allocation}) != len(allocation):
        raise ValueError("Duplicate original cage allocation")
    if len({r["random_seed"] for r in allocation}) != 1:
        raise ValueError("Allocation lacks one preserved random seed")
    original = {}
    for r in allocation:
        if r["design_version"] != DESIGN or r["treatment"] not in ARMS:
            raise ValueError("BIBD allocation belongs to another experiment or unknown arm")
        if (nonnegative(r["competitor_present"], "competitor", integer=True),
                nonnegative(r["food_support"], "food support", integer=True)) != ARMS[r["treatment"]]:
            raise ValueError("Original random assignment and factor labels disagree")
        if not r["cage_id"].strip() or not r["greenhouse_id"].strip():
            raise ValueError("Missing original whole-cage/greenhouse ID")
        nonnegative(r["random_seed"], "random seed", integer=True)
        original[r["cage_id"]] = r
    if not counts(original.values())["complete_v4_b8_k3_r6_lambda4"]:
        raise ValueError("Original assignment is not a connected balanced BIBD")
    if len(info) != len(original):
        raise ValueError("Missing randomized whole cage cannot be excluded")
    cages = {}
    for r in info:
        cage = r["cage_id"]
        if not cage or cage in cages or cage not in original:
            raise ValueError("Missing, duplicate or unrandomized cage information")
        for k in ("greenhouse_id", "treatment", "competitor_present", "food_support"):
            if r[k] != original[cage][k]:
                raise ValueError("Postrandomization treatment or greenhouse change")
        n = nonnegative(r["focal_a_initial_n"], "focal count", integer=True, positive=True)
        added = nonnegative(r["added_s_initial_n"], "added S count", integer=True)
        expected, support = ARMS[r["treatment"]]
        if bool(added) != bool(expected):
            raise ValueError("Competitor originally assigned does not match added S count")
        visits = nonnegative(r["expected_n_food_visits"], "scheduled visits", integer=True, positive=True)
        if visits < 2:
            raise ValueError("At least two independent scheduled foliage occasions")
        floor = nonnegative(r["food_floor_cm2"], "food floor", positive=True)
        if not r["donor_batch_id"].strip():
            raise ValueError("Original source donor plant batch required")
        cages[cage] = {
            **r, "focal_a_initial_n": n, "added_s_initial_n": added,
            "expected_n_food_visits": visits, "food_floor_cm2": floor,
        }
    if len({v["focal_a_initial_n"] for v in cages.values()}) != 1:
        raise ValueError("Focal A initial cohort size must be fixed across cages")
    added_positive = {v["added_s_initial_n"] for v in cages.values()
                      if v["added_s_initial_n"] > 0}
    if len(added_positive) != 1:
        raise ValueError("Added S density not fixed before treatments")
    by_gh = defaultdict(list)
    for v in cages.values():
        by_gh[v["greenhouse_id"]].append(v)
    for gh, members in by_gh.items():
        for k in ("expected_n_food_visits", "food_floor_cm2", "donor_batch_id"):
            if len({r[k] for r in members}) != 1:
                raise ValueError(f"{gh}: preassigned resource observation/plant batch mismatched")
    return cages


def fates_by_cage(cages: dict, rows: list[dict]) -> dict:
    observed = defaultdict(lambda: defaultdict(list))
    ids = set()
    for r in rows:
        cage, ident = r["cage_id"], r["larva_id"]
        if cage not in cages or not ident or ident in ids:
            raise ValueError("Unassigned cage or duplicate/missing original larva")
        ids.add(ident)
        role = r["cohort_role"]
        if role not in {"focal", "added"}:
            raise ValueError("Role must be original focal A or added S")
        correct = FOCAL if role == "focal" else ADDED
        if r["species"] != correct:
            raise ValueError("Wrong original focal versus added butterfly species")
        if r["fate"] not in VALID_FATES:
            raise ValueError("Unknown or ambiguous pupa_no_adult outcome")
        observed[cage][role].append(r["fate"])
    result = {}
    for cage, contract in cages.items():
        x = observed[cage]["focal"]
        y = observed[cage]["added"]
        if len(x) != contract["focal_a_initial_n"] or len(y) != contract["added_s_initial_n"]:
            raise ValueError("Every originally randomized larva must have fate row")
        ct = Counter(x)
        low = ct[SUCCESS] / len(x)
        possible = sum(ct[k] for k in PENDING)
        result[cage] = {
            "original_focal_A": len(x),
            "original_added_S": len(y),
            "observed_flight_capable_focal_adults": ct[SUCCESS],
            "focal_adult_low": low,
            "focal_eventual_adult_upper": (ct[SUCCESS] + possible)/len(x),
            "living_diapause_focal": ct["diapause_alive"],
            "living_pupa_pending_focal": ct["pupa_alive_pending_emergence"],
            "other_focal_unknown": ct["lost_to_followup"] + ct["other_documented"],
            "confirmed_death_or_nonflight": sum(ct[k] for k in FAIL),
            "added_S_fate_counts": dict(Counter(y)),
        }
    return result


def foliage_by_cage(cages: dict, rows: list[dict]) -> dict:
    out = defaultdict(list)
    seen = set()
    for r in rows:
        cid = r["cage_id"]
        if cid not in cages:
            raise ValueError("Foliage observation belongs to unassigned cage")
        try:
            dt = datetime.fromisoformat(r["observation_datetime"])
        except ValueError as exc:
            raise ValueError("Resource measurement needs ISO date with offset") from exc
        if dt.tzinfo is None:
            raise ValueError("Timezone needed for original resource observation")
        key=(cid,dt.isoformat())
        if key in seen:
            raise ValueError("Duplicate cage resource observation")
        seen.add(key)
        c = cages[cid]
        alive = nonnegative(r["focal_a_alive_n"], "living focal", integer=True)
        late = nonnegative(r["focal_a_late_instar_alive_n"], "late focal", integer=True)
        added = nonnegative(r["added_s_alive_n"], "living introduced S", integer=True)
        if not 0 <= late <= alive <= c["focal_a_initial_n"] or added > c["added_s_initial_n"]:
            raise ValueError("Living native/introduced counts exceed original allocation")
        area = nonnegative(r["accessible_edible_leaf_area_cm2"], "accessible leaf area", missing=True)
        supplemented = nonnegative(r["accessible_leaf_added_area_cm2"], "additional edible area", missing=True)
        if r["leaf_development_stage"] not in VALID_LEAF or (
            r["added_leaf_stage"] and r["added_leaf_stage"] not in VALID_LEAF
        ):
            raise ValueError("Unknown stage of actually edible host leaf")
        if supplemented is not None and supplemented > 0:
            if c["food_support"] != "1" or not r["added_leaf_stage"]:
                raise ValueError("Accessible supplementation mismatches original intervention")
        if r["sham_handling_done"] not in ("0","1") or not r["leaf_source_batch_id"].strip():
            raise ValueError("Missing sham handling or original leaf donor identification")
        out[cid].append({
            "area": area, "added_area": supplemented, "late": late, "alive": alive,
            "donor_id": r["leaf_source_batch_id"],
            "sham": r["sham_handling_done"] == "1",
        })
    summary = {}
    for cage, c in cages.items():
        v = out[cage]
        if len(v) != c["expected_n_food_visits"]:
            raise ValueError("A randomized cage has missing scheduled food observation")
        observed = [r for r in v if r["area"] is not None]
        late_obs = [r for r in observed if r["late"] > 0]
        summary[cage] = {
            "scheduled_visits": len(v),
            "missing_edible_area_visits": len(v) - len(observed),
            "measured_edible_area_below_original_floor": sum(
                r["area"] < c["food_floor_cm2"] for r in observed),
            "late_living_focal_measured_visits": len(late_obs),
            "late_living_focal_below_floor_visits": sum(
                r["area"] < c["food_floor_cm2"] for r in late_obs),
            "sham_handling_true_visits": sum(r["sham"] for r in v),
            "post_treatment_diagnostic_only": True,
        }
    return summary


def regression(cages: list[dict], adult_fates: dict) -> dict:
    """Fixed-greenhouse adjusted 2x2 factorial contrast on cage-level risks.

    Allowed to inspect 7-greenhouse leave-one-greenhouse-out sensitivity, but
    validated original study still requires all eight balanced greenhouses.
    """
    order = sorted(cages, key=lambda r:(r["greenhouse_id"],r["cage_id"]))
    gh = sorted({r["greenhouse_id"] for r in order})
    base = "NATURAL_NONE"
    other = [a for a in ARMS if a != base]
    if len(order) < len(gh)+len(other)+1:
        raise ValueError("Not enough cages to estimate block-adjusted response")
    X = np.array([
        [int(r["greenhouse_id"] == b) for b in gh] +
        [int(r["treatment"] == a) for a in other]
        for r in order
    ], dtype=float)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError("Incomplete block design is disconnected; contrast not identified")
    y_low = np.array([adult_fates[r["cage_id"]]["focal_adult_low"] for r in order])
    y_up = np.array([adult_fates[r["cage_id"]]["focal_eventual_adult_upper"] for r in order])
    contrast=np.array([0.0]*len(gh)+[-1.0,-1.0,+1.0])
    xtxi=np.linalg.inv(X.T@X)
    betahat=np.linalg.lstsq(X,y_low,rcond=None)[0]
    observed=float(contrast@betahat)
    w=X@xtxi@contrast
    if not math.isclose(observed,float(w@y_low),abs_tol=1e-9):
        raise RuntimeError("Adjusted whole-cage weights fail estimator check")
    low=float(np.where(w >= 0, w*y_low, w*y_up).sum())
    high=float(np.where(w >= 0, w*y_up, w*y_low).sum())
    residuals=y_low-X@betahat
    df=len(order)-X.shape[1]
    if df <= 0:
        raise ValueError("Model has no independent residual degrees of freedom")
    residual_sd=float(np.sqrt(sum(residuals**2)/df))
    multiplier=float(contrast@xtxi@contrast)
    se=math.sqrt(multiplier)*residual_sd
    # t(13) is valid only for the 24-cage original design.
    model_interval=([observed-T_CRITICAL_DF13*se,observed+T_CRITICAL_DF13*se]
                    if df==13 else None)
    return {
        "n_greenhouses":len(gh), "n_independent_whole_cages":len(order),
        "n_residual_df":df,
        "greenhouse_adjusted_observed_adult_contrast":observed,
        "possible_eventual_emergence_contrast_lower":low,
        "possible_eventual_emergence_contrast_upper":high,
        "whole_cage_linear_contrast_weights": {
            r["cage_id"]:float(weight) for r,weight in zip(order,w)
        },
        "analytical_contrast_variance_multiplier":multiplier,
        "residual_cage_risk_SD_under_fixed_greenhouse_additivity":residual_sd,
        "model_based_standard_error_given_homoscedastic_cage_residuals":se,
        "nominal_t_interval_given_homoscedastic_independent_cage_errors":model_interval,
        "model_based_interval_is_NOT_design_based_and_NOT_evidence_of_power":True,
    }


def analyze(allocation,info,fates,resources):
    cages=validate_assignment(allocation,info)
    ends=fates_by_cage(cages,fates)
    foliage=foliage_by_cage(cages,resources)
    fit=regression(list(cages.values()),ends)
    if fit["n_residual_df"]!=13 or abs(fit["analytical_contrast_variance_multiplier"]-.75)>1e-9:
        raise RuntimeError("Frozen BIBD model must have 13 residual df and variance multiplier 0.75")
    gh=sorted({r["greenhouse_id"] for r in cages.values()})
    loo=[]
    for omitted in gh:
        trial=[r for r in cages.values() if r["greenhouse_id"] != omitted]
        value=regression(trial,ends)
        loo.append({"greenhouse_left_out":omitted,
                    "adjusted_primary":value["greenhouse_adjusted_observed_adult_contrast"]})
    per_arm={}
    for arm in ARMS:
        original=[cid for cid,c in cages.items() if c["treatment"]==arm]
        total=sum(ends[cid]["original_focal_A"] for cid in original)
        success=sum(ends[cid]["observed_flight_capable_focal_adults"] for cid in original)
        unresolved=sum(
            ends[cid]["living_diapause_focal"]+
            ends[cid]["living_pupa_pending_focal"]+
            ends[cid]["other_focal_unknown"] for cid in original
        )
        per_arm[arm]={
            "independent_whole_cages":len(original),
            "original_focal_A_larvae":total,
            "verified_flight_capable_adults":success,
            "living_diapause_A":sum(ends[cid]["living_diapause_focal"] for cid in original),
            "other_unknown_eventual_adult_fates":unresolved-
                sum(ends[cid]["living_diapause_focal"] for cid in original),
            "raw_observed_adult_fraction":success/total,
            "raw_eventual_adult_fraction_upper":(success+unresolved)/total,
            "food_visits_missing_leaf_area":sum(
                foliage[cid]["missing_edible_area_visits"] for cid in original),
            "measured_food_floor_breach_visits":sum(
                foliage[cid]["measured_edible_area_below_original_floor"] for cid in original),
        }
    return {
        "schema":"chocho_kyoto_four_arm_BIBD_original_cohort_ITT_v01",
        "state":"FUTURE_NEW_REAL_FIELD_DATA_REQUIRED__SYNTHETIC_TESTS_ARE_NOT_RESULTS",
        "design_properties":counts(allocation),
        "original_greenhouse_unit_not_greenhouse_identical_to_entire_cage":True,
        "primary_estimand":"(SUPPORTED_ADD_S-SUPPORTED_NONE)-(NATURAL_ADD_S-NATURAL_NONE)",
        "primary_block_adjusted":fit,
        "per_treatment_raw_arm_fates":per_arm,
        "leave_one_original_greenhouse_out":loo,
        "leave_one_greenhouse_effect_min":min(x["adjusted_primary"] for x in loo),
        "leave_one_greenhouse_effect_max":max(x["adjusted_primary"] for x in loo),
        "no_pseudo_replication_of_nested_pots_or_larvae":True,
        "no_dropping_cages_for_postrandomization_food_observations":True,
        "no_invasive_host_plant_effect":True,
        "source_population_generalization_not_established":True,
        "cannot_identify_leaf_biomass_only_vs_plant_quality_vs_contact":True,
        "main_GEB_unchanged":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("allocation","assignments","fates","resources","receipt"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    sourcefiles={
        "allocation":a.allocation, "assignments":a.assignments,
        "fates":a.fates, "resources":a.resources,
    }
    out=analyze(
        read_csv(a.allocation,ALLOCATION_COLUMNS),
        read_csv(a.assignments,ASSIGNMENT_COLUMNS),
        read_csv(a.fates,FATE_COLUMNS),
        read_csv(a.resources,RESOURCE_COLUMNS),
    )
    out["original_input_sha256"]={
        key:hashlib.sha256(path.read_bytes()).hexdigest()
        for key,path in sourcefiles.items()
    }
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    a.receipt.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
