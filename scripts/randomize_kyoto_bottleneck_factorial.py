#!/usr/bin/env python3
"""Auditable blocked assignment for a future Kyoto butterfly pilot.

No biological data simulated; does not set sample size or assert an effect.
Input: original inventory CSV with unique real plant_id,block_id for plants
that have met experimental and legal eligibility requirements.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

DESIGN="kyoto_hidden_resource_bottleneck_v01"
ARMS={
    "NATURAL_NO_COMPETITOR":(0,0),
    "NATURAL_COMPETITOR":(1,0),
    "CLAMP_NO_COMPETITOR":(0,1),
    "CLAMP_COMPETITOR":(1,1),
}

def assign(rows:list[dict[str,str]],seed:int)->list[dict[str,str|int]]:
    if not rows:
        raise ValueError("no eligible experimental plants submitted")
    allowed={"plant_id","block_id"}
    if any(set(row)!=allowed for row in rows):
        raise ValueError("inventory must contain only plant_id,block_id; no outcome columns")
    seen=set()
    blocks=defaultdict(list)
    for row in rows:
        plant=row["plant_id"].strip()
        block=row["block_id"].strip()
        if not plant or not block:
            raise ValueError("blank plant or source block")
        if plant in seen:
            raise ValueError(f"duplicate plant identity: {plant}")
        seen.add(plant)
        blocks[block].append(plant)

    # Randomization performed after eligibility list is finalized; each complete
    # block contains equal numbers of each 2x2 treatment.
    rng=random.Random(seed)
    assignments=[]
    for block in sorted(blocks):
        plants=sorted(blocks[block])
        if len(plants)<4 or len(plants)%4:
            raise ValueError(f"block {block}: expected >=4 and multiple of 4 plants")
        replicate=len(plants)//4
        candidates=[name for name in ARMS for _ in range(replicate)]
        rng.shuffle(candidates)
        for plant,arm in zip(plants,candidates):
            competitor,clamp=ARMS[arm]
            assignments.append({
                "plant_id":plant,"block_id":block,
                "treatment":arm,"competitor_present":competitor,
                "food_clamp":clamp,"random_seed":seed,
                "design_version":DESIGN,
            })
    if len(assignments)!=len(rows):
        raise RuntimeError("assignments not exactly all eligible plants")
    assert all(set(Counter(x["treatment"] for x in assignments if x["block_id"]==b).values())=={len(blocks[b])//4} for b in blocks)
    return assignments

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory",type=Path,required=True)
    p.add_argument("--allocation",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--seed",type=int,required=True)
    args=p.parse_args()
    with args.inventory.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=["plant_id","block_id"]:
            raise ValueError("source list must have columns plant_id,block_id in that order")
        allocation=assign(list(reader),args.seed)
    args.allocation.parent.mkdir(parents=True,exist_ok=True)
    with args.allocation.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(allocation[0]))
        writer.writeheader()
        writer.writerows(allocation)
    receipt={
        "schema":"kyoto_factorial_randomization_receipt_v01",
        "design_version":DESIGN,
        "random_seed":args.seed,
        "n_source_blocks":len(set(x["block_id"] for x in allocation)),
        "n_independently_assigned_plants":len(allocation),
        "per_treatment":dict(Counter(x["treatment"] for x in allocation)),
        "biological_pilot_conducted":False,
        "biological_outcome_estimated":False,
        "sampling_assumption":"All plants in inventory eligible BEFORE treatment assignment. Block count is user-supplied, not a statistically justified target.",
        "GEB_PR38_untouched":True,
    }
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
