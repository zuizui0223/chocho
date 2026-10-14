#!/usr/bin/env python3
"""Balanced, blocked SIX-ARM allocation of WHOLE CAGES, never nested host plants.

Only an already eligible inventory is accepted; no animal or outcome rows exist.
No sample size, field eligibility, source population, or biological effects are
invented. Keeping this separate from the published two-by-two experiment avoids
retrospective treatment reassignment.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

DESIGN = "kyoto_identity_food_six_arm_v01"
ARMS = {
    "NATURAL_NONE": (0, "none"),
    "NATURAL_ADD_A": (0, "A"),
    "NATURAL_ADD_S": (0, "S"),
    "SUPPORTED_NONE": (1, "none"),
    "SUPPORTED_ADD_A": (1, "A"),
    "SUPPORTED_ADD_S": (1, "S"),
}
ALLOCATION_COLUMNS = (
    "cage_id", "block_id", "treatment", "food_support", "added_group",
    "random_seed", "design_version",
)


def assign(rows: list[dict[str, str]], seed: int) -> list[dict]:
    if not isinstance(seed, int):
        raise ValueError("random seed must be an integer, recorded before assignment")
    if not rows:
        raise ValueError("No real, eligible cage inventory was supplied")
    blocks: dict[str, list[str]] = defaultdict(list)
    seen: set[str] = set()
    for row in rows:
        if set(row) != {"cage_id", "block_id"}:
            raise ValueError("Inventory must have only cage_id and block_id; no outcomes")
        cage, block = row["cage_id"].strip(), row["block_id"].strip()
        if not cage or not block:
            raise ValueError("blank cage ID or source block")
        if cage in seen:
            raise ValueError(f"duplicate original whole-cage identifier {cage}")
        seen.add(cage)
        blocks[block].append(cage)
    for block, cages in blocks.items():
        if len(cages) < 6 or len(cages) % 6:
            raise ValueError(
                f"block {block}: expected >=6 and divisible by 6 independently eligible whole cages"
            )
    rng = random.Random(seed)
    output: list[dict] = []
    for block in sorted(blocks):
        cages = sorted(blocks[block])
        per_arm = len(cages) // 6
        labels = [arm for arm in ARMS for _ in range(per_arm)]
        rng.shuffle(labels)
        for cage, label in zip(cages, labels):
            food, added = ARMS[label]
            output.append({
                "cage_id": cage,
                "block_id": block,
                "treatment": label,
                "food_support": food,
                "added_group": added,
                "random_seed": seed,
                "design_version": DESIGN,
            })
    counts = Counter((x["block_id"], x["treatment"]) for x in output)
    if any(counts[(block, arm)] != len(cages) // 6
           for block, cages in blocks.items() for arm in ARMS):
        raise RuntimeError("internal source-block balance check failed")
    return output


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory", type=Path, required=True)
    p.add_argument("--allocation", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    p.add_argument("--seed", type=int, required=True)
    args = p.parse_args()
    with args.inventory.open(encoding="utf-8-sig", newline="") as handle:
        data = csv.DictReader(handle)
        if data.fieldnames != ["cage_id", "block_id"]:
            raise ValueError("Inventory columns must be cage_id,block_id in this exact order")
        inventory = list(data)
    allocated = assign(inventory, args.seed)
    args.allocation.parent.mkdir(parents=True, exist_ok=True)
    with args.allocation.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ALLOCATION_COLUMNS))
        writer.writeheader()
        writer.writerows(allocated)
    receipt = {
        "schema": "kyoto_six_arm_blocked_randomization_receipt_v01",
        "design_version": DESIGN,
        "source_inventory_sha256": hashlib.sha256(args.inventory.read_bytes()).hexdigest(),
        "output_allocation_sha256": hashlib.sha256(args.allocation.read_bytes()).hexdigest(),
        "seed": args.seed,
        "independent_source_blocks": len(set(x["block_id"] for x in allocated)),
        "independently_randomized_whole_cages": len(allocated),
        "treatment_counts": dict(Counter(x["treatment"] for x in allocated)),
        "no_animal_observations_or_outcomes_generated": True,
        "no_sample_size_or_power_claim": True,
        "permissions_and_containment_required_before_trial": True,
        "GEB_submission_main_unchanged": True,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
