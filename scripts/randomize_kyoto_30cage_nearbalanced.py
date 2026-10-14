#!/usr/bin/env python3
"""Prospective 10-greenhouse × 3-whole-cage near-balanced incomplete-block design.

NO recorded plant, butterfly or cage inventory is bundled. Only an actual,
pre-eligible 30-cage inventory is accepted. The *four-arm invader×accessible
foliage* contrast is prioritized; source-greenhouse effects remain blocked.

Because 30/4 is not an integer, an exact BIBD is impossible. We choose
omission totals (3,2,2,3) or (2,3,3,2) in the arm order
[NATURAL_NONE, NATURAL_ADD_S, SUPPORTED_NONE, SUPPORTED_ADD_S].
Each greenhouse has three distinct arms, each arm 7 or 8 cages; both
factor levels each have 15 total cages. The primary contrast has model-based
variance multiplier 0.6 sigma², verified independently by matrix tests.

Randomly choose the equally precise complementary omission pattern; randomly
permute omissions among the 10 original greenhouse IDs; randomly permute the
three arms among the three independently eligible cages within each greenhouse.
No post-hoc selection by plant or butterfly response is allowed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from randomize_kyoto_four_arm_greenhouse_bibd import (
    ARMS, ALLOCATION_COLUMNS, INVENTORY_COLUMNS,
)

DESIGN = "kyoto_10greenhouse_30cage_nearbalanced_v01"
B, K = 10, 3
POSITIVE_CONTRAST_ARMS = ("NATURAL_NONE", "SUPPORTED_ADD_S")
NEGATIVE_CONTRAST_ARMS = ("NATURAL_ADD_S", "SUPPORTED_NONE")
VARIANCE_FACTOR_IF_HOMOSCEDASTIC = 0.60
OMISSION_PATTERN_CANDIDATES = (POSITIVE_CONTRAST_ARMS, NEGATIVE_CONTRAST_ARMS)


def validate_inventory(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    if len(rows) != B * K:
        raise ValueError("Require exactly 30 independently eligible whole cages, 10 greenhouses x 3")
    all_cages = set()
    by_greenhouse: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if set(row) != set(INVENTORY_COLUMNS):
            raise ValueError("Only cage_id,greenhouse_id columns allowed; no invented outcomes")
        cage, greenhouse = row["cage_id"].strip(), row["greenhouse_id"].strip()
        if not cage or not greenhouse or cage in all_cages:
            raise ValueError("Missing/duplicate independently eligible whole-cage or greenhouse ID")
        all_cages.add(cage)
        by_greenhouse[greenhouse].append(cage)
    if len(by_greenhouse) != B or any(len(cages) != K for cages in by_greenhouse.values()):
        raise ValueError("Exactly 10 original greenhouses with three eligible cages each required")
    return by_greenhouse


def design_counts(rows: list[dict]) -> dict:
    by_greenhouse: dict[str, set[str]] = defaultdict(set)
    per_arm = Counter()
    seen = set()
    for row in rows:
        greenhouse, arm = row["greenhouse_id"], row["treatment"]
        cage = row["cage_id"]
        if arm not in ARMS or not greenhouse or not cage or cage in seen:
            raise ValueError("Invalid or duplicate near-balanced arm/cage/greenhouse")
        if arm in by_greenhouse[greenhouse]:
            raise ValueError("Within-greenhouse duplicate treatment")
        seen.add(cage)
        by_greenhouse[greenhouse].add(arm)
        per_arm[arm] += 1
    omissions = Counter()
    for greenhouse, selected in by_greenhouse.items():
        if len(selected) != K:
            raise ValueError(f"Greenhouse {greenhouse} must host three DISTINCT arms")
        for arm in set(ARMS) - selected:
            omissions[arm] += 1
    cooccurrences = {
        f"{a}|{b}": sum({a, b}.issubset(arms) for arms in by_greenhouse.values())
        for a, b in itertools.combinations(ARMS, 2)
    }
    totals = {
        "competitor_absent": sum(per_arm[a] for a, (s, _) in ARMS.items() if s == 0),
        "competitor_present": sum(per_arm[a] for a, (s, _) in ARMS.items() if s == 1),
        "natural_foliage": sum(per_arm[a] for a, (_, f) in ARMS.items() if f == 0),
        "supported_foliage": sum(per_arm[a] for a, (_, f) in ARMS.items() if f == 1),
    }
    omissions_correct = any(
        all(
            omissions[a] == (3 if a in three_omissions else 2)
            for a in ARMS
        )
        for three_omissions in OMISSION_PATTERN_CANDIDATES
    )
    return {
        "n_greenhouses": len(by_greenhouse),
        "n_whole_cages": len(rows),
        "replicates_by_arm": {a: per_arm[a] for a in ARMS},
        "omissions_by_arm": {a: omissions[a] for a in ARMS},
        "pair_cooccurrences": cooccurrences,
        "main_factor_cage_counts": totals,
        "design_valid": (
            len(by_greenhouse) == B
            and len(rows) == B * K
            and omissions_correct
            and set(per_arm.values()) == {7, 8}
            and set(totals.values()) == {15}
            and sorted(cooccurrences.values()) == [4, 5, 5, 5, 5, 6]
        ),
    }


def assign(rows: list[dict[str, str]], *, seed: int) -> list[dict]:
    if not isinstance(seed, int):
        raise ValueError("Seed must be an integer frozen before outcomes")
    blocks = validate_inventory(rows)
    rng = random.Random(seed)
    three_omission_arms = OMISSION_PATTERN_CANDIDATES[rng.randrange(2)]
    omissions = [a for a in ARMS for _ in range(3 if a in three_omission_arms else 2)]
    rng.shuffle(omissions)
    output = []
    for greenhouse, omitted in zip(sorted(blocks), omissions):
        included = [a for a in ARMS if a != omitted]
        rng.shuffle(included)
        for cage, arm in zip(sorted(blocks[greenhouse]), included):
            comp, support = ARMS[arm]
            output.append({
                "cage_id": cage, "greenhouse_id": greenhouse,
                "treatment": arm, "competitor_present": comp,
                "food_support": support, "random_seed": seed,
                "design_version": DESIGN,
            })
    if not design_counts(output)["design_valid"]:
        raise RuntimeError("Randomized 30-cage near-balanced incidence does not satisfy contract")
    return output


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", required=True, type=Path)
    ap.add_argument("--allocation", required=True, type=Path)
    ap.add_argument("--receipt", required=True, type=Path)
    ap.add_argument("--seed", required=True, type=int)
    args = ap.parse_args()
    with args.inventory.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != INVENTORY_COLUMNS:
            raise ValueError("Original inventory must use cage_id,greenhouse_id exactly")
        inventory = list(reader)
    assigned = assign(inventory, seed=args.seed)
    args.allocation.parent.mkdir(parents=True, exist_ok=True)
    with args.allocation.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALLOCATION_COLUMNS)
        writer.writeheader()
        writer.writerows(assigned)
    receipt = {
        "schema": "chocho_kyoto_30cage_nearbalanced_assignment_v01",
        "design": DESIGN,
        "status": "ASSIGNED_ONLY_IF_USER_SUPPLIED_VALID_SOURCE_INVENTORY",
        "random_seed": args.seed,
        "source_inventory_sha256": hashlib.sha256(args.inventory.read_bytes()).hexdigest(),
        "assigned_original_whole_cages_sha256": hashlib.sha256(args.allocation.read_bytes()).hexdigest(),
        "design_counts": design_counts(assigned),
        "conditional_primary_variance_multiplier": VARIANCE_FACTOR_IF_HOMOSCEDASTIC,
        "no_pilot_source_CV_or_empirical_fitness_observed": True,
        "present_facility_access_not_confirmed_by_historical_paper": True,
        "no_permits_or_containment_claimed": True,
        "GEB_submission_main_untouched": True,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
