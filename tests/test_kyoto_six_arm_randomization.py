"""Entirely synthetic ledger fixtures; these are NOT butterfly observations."""
from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
path = ROOT / "scripts/randomize_kyoto_six_arm_identity_food.py"
spec = importlib.util.spec_from_file_location("kyoto6_random", path)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def eligible(blocks=2, per_block=6):
    return [{"cage_id": f"B{b}:C{c}", "block_id": f"B{b}"}
            for b in range(blocks) for c in range(per_block)]


def test_deterministic_six_arm_randomization_with_block_balance():
    a = m.assign(eligible(3, 12), seed=20261014)
    b = m.assign(list(reversed(eligible(3, 12))), seed=20261014)
    assert a == b
    assert len(a) == 36
    for group in range(3):
        counts = Counter(row["treatment"] for row in a if row["block_id"] == f"B{group}")
        assert dict(counts) == {name: 2 for name in m.ARMS}
    assert all(r["design_version"] == m.DESIGN for r in a)
    assert all("survival" not in r for r in a)
    assert any(row["added_group"] == "none" and row["food_support"] == 0 for row in a)


def test_distinct_seed_can_change_assignment_but_not_balance():
    a = m.assign(eligible(2, 6), seed=12)
    b = m.assign(eligible(2, 6), seed=17)
    assert a != b
    for rows in (a, b):
        assert Counter(x["treatment"] for x in rows) == Counter({k: 2 for k in m.ARMS})


@pytest.mark.parametrize("inventory,pattern", [
    ([], "No real"),
    ([{"cage_id": "X", "block_id": "B", "outcome": "1"}], "no outcomes"),
    ([{"cage_id": " ", "block_id": "B"}], "blank cage"),
    (eligible(1, 5), "divisible by 6"),
    (eligible(1, 7), "divisible by 6"),
])
def test_invalid_inventory_fails_closed(inventory, pattern):
    with pytest.raises(ValueError, match=pattern):
        m.assign(inventory, seed=9)


def test_duplicate_cage_id_rejected():
    inventory = eligible(1, 6)
    inventory[-1]["cage_id"] = inventory[0]["cage_id"]
    with pytest.raises(ValueError, match="duplicate original"):
        m.assign(inventory, seed=9)
