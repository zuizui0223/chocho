#!/usr/bin/env python3
"""Randomize a future 4-treatment balanced incomplete-block Kyoto experiment.

Eight *preselected, independently eligible* greenhouses with exactly three
WHOLE cages each. Each treatment is omitted from two greenhouses, appears
in six, and co-occurs with each other treatment in four. This solves the
physical 3-cage/greenhouse block SIZE constraint, not power or access to
real greenhouses. The 10x3 original 2012 setup is a historical example only.

The randomized greenhouse omission pattern is generated before outcomes, then
the three allowed labels are shuffled among original cage IDs. No seedlings,
caterpillars, observations or field permit access are created by this script.
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

DESIGN = "kyoto_four_arm_greenhouse_bibd_v01"
ARMS = {
    "NATURAL_NONE": (0, 0),
    "NATURAL_ADD_S": (1, 0),
    "SUPPORTED_NONE": (0, 1),
    "SUPPORTED_ADD_S": (1, 1),
}
INVENTORY_COLUMNS = ["cage_id", "greenhouse_id"]
ALLOCATION_COLUMNS = [
    "cage_id", "greenhouse_id", "treatment", "competitor_present",
    "food_support", "random_seed", "design_version",
]
B = 8
K = 3
V = 4
R = 6
LAMBDA = 4


def validate(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    if len(rows) != B * K:
        raise ValueError(
            "BIBD requires exactly 24 PRESELECTED eligible whole cages "
            "from 8 greenhouses x 3; no automatic greenhouse selection"
        )
    ids = set()
    blocks: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if set(row) != set(INVENTORY_COLUMNS):
            raise ValueError("Exact cage_id,greenhouse_id inventory only: no outcomes or inferred capacity")
        cid, greenhouse = row["cage_id"].strip(), row["greenhouse_id"].strip()
        if not cid or not greenhouse or cid in ids:
            raise ValueError("Blank or duplicate original whole-cage/greenhouse identifier")
        ids.add(cid)
        blocks[greenhouse].append(cid)
    if len(blocks) != B or any(len(v) != K for v in blocks.values()):
        raise ValueError("Each of exactly 8 original eligible greenhouses requires 3 cages")
    return blocks


def counts(rows: list[dict]) -> dict:
    block_arms: dict[str, set[str]] = defaultdict(set)
    treatment_counts = Counter()
    for row in rows:
        arm = row["treatment"]
        if arm not in ARMS:
            raise ValueError("Unknown BIBD treatment")
        if arm in block_arms[row["greenhouse_id"]]:
            raise ValueError("Duplicate treatment within same greenhouse")
        block_arms[row["greenhouse_id"]].add(arm)
        treatment_counts[arm] += 1
    pair = {
        f"{a}|{b}": sum({a,b}.issubset(arms) for arms in block_arms.values())
        for a,b in itertools.combinations(ARMS, 2)
    }
    omitted = Counter()
    for greenhouse, arms in block_arms.items():
        if len(arms) != K:
            raise ValueError(f"Greenhouse {greenhouse} not assigned exactly 3 different arms")
        for absent in set(ARMS) - arms:
            omitted[absent] += 1
    return {
        "n_greenhouses": len(block_arms),
        "n_whole_cages": len(rows),
        "replicates_by_treatment": {a:treatment_counts[a] for a in ARMS},
        "omissions_by_treatment": {a:omitted[a] for a in ARMS},
        "pair_cooccurrences_by_treatment_pair": pair,
        "complete_v4_b8_k3_r6_lambda4": (
            len(rows) == 24 and len(block_arms) == 8
            and all(treatment_counts[a] == R and omitted[a] == 2 for a in ARMS)
            and all(n == LAMBDA for n in pair.values())
        ),
    }


def assign(inventory: list[dict[str, str]], *, seed: int) -> list[dict]:
    if not isinstance(seed, int):
        raise ValueError("Provide one integer seed fixed before observation")
    blocks = validate(inventory)
    rng = random.Random(seed)
    omission_labels = [arm for arm in ARMS for _ in range(2)]
    rng.shuffle(omission_labels)
    result: list[dict] = []
    for greenhouse, omit in zip(sorted(blocks), omission_labels):
        selected = [a for a in ARMS if a != omit]
        rng.shuffle(selected)
        for cage, arm in zip(sorted(blocks[greenhouse]), selected):
            competitor, food_support = ARMS[arm]
            result.append({
                "cage_id": cage, "greenhouse_id": greenhouse, "treatment": arm,
                "competitor_present": competitor, "food_support": food_support,
                "random_seed": seed, "design_version": DESIGN,
            })
    c = counts(result)
    if not c["complete_v4_b8_k3_r6_lambda4"]:
        raise RuntimeError("Randomized incomplete-block incidence is not balanced")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory", required=True, type=Path)
    p.add_argument("--allocation", required=True, type=Path)
    p.add_argument("--receipt", required=True, type=Path)
    p.add_argument("--seed", required=True, type=int)
    args = p.parse_args()
    with args.inventory.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != INVENTORY_COLUMNS:
            raise ValueError("Eligible inventory headings must be cage_id,greenhouse_id")
        inventory = list(reader)
    result = assign(inventory, seed=args.seed)
    args.allocation.parent.mkdir(parents=True, exist_ok=True)
    with args.allocation.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALLOCATION_COLUMNS)
        writer.writeheader()
        writer.writerows(result)
    summary = {
        "schema": "chocho_kyoto_four_arm_greenhouse_BIBD_allocation_v01",
        "status": "ONLY_AN_ALLOCATION_CONTRACT_UNLESS_REAL_CAGES_SUPPLIED",
        "design": DESIGN,
        "seed": args.seed,
        "source_inventory_sha256": hashlib.sha256(args.inventory.read_bytes()).hexdigest(),
        "allocation_sha256": hashlib.sha256(args.allocation.read_bytes()).hexdigest(),
        "design_counts": counts(result),
        "baseline_facility_from_published_2012_cages_not_verified_now": True,
        "no_animal_data_or_causal_effect_estimated": True,
        "not_a_power_calculation_or_sample_size_recommendation": True,
        "GEB_submission_main_unchanged": True,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
