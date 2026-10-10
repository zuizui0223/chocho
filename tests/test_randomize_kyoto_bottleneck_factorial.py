"""Tests use fictitious PLANT IDENTITIES only; no simulated butterfly survival."""
import importlib.util
from collections import Counter,defaultdict
from pathlib import Path

script=Path(__file__).resolve().parents[1]/"scripts/randomize_kyoto_bottleneck_factorial.py"
s=importlib.util.spec_from_file_location("assignment",script)
assert s and s.loader
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

def inventory(block_count=4,plant_per_block=8):
    return [{"plant_id":f"TEST_PLANT_{b:02d}_{p:02d}","block_id":f"B{b:02d}"}
            for b in range(block_count) for p in range(plant_per_block)]

def test_block_balanced_2x2_without_outcome_data():
    x=m.assign(inventory(),seed=20261010)
    assert len(x)==32
    for block in {r["block_id"] for r in x}:
        assert Counter(z["treatment"] for z in x if z["block_id"]==block)=={a:2 for a in m.ARMS}
    assert len(set(r["plant_id"] for r in x))==len(x)
    assert {r["competitor_present"] for r in x}=={0,1}
    assert {r["food_clamp"] for r in x}=={0,1}

def test_stable_randomization_and_seed_matters():
    rows=inventory()
    a=m.assign(rows,20261010)
    b=m.assign(list(reversed(rows)),20261010)
    c=m.assign(rows,20261011)
    assert a==b
    assert [(r["plant_id"],r["treatment"]) for r in a]!=[(r["plant_id"],r["treatment"]) for r in c]

def test_no_uneven_or_duplicate_blocks():
    bads=[
        inventory(1,5),
        inventory(1,3),
        inventory(1,4)+[inventory(1,4)[0]],
        inventory(1,4)+[{"plant_id":"","block_id":"B2"}],
    ]
    for rows in bads:
        try:
            m.assign(rows,42)
        except ValueError:
            pass
        else:
            raise AssertionError("accepted an invalid original plant inventory")

def test_no_outcome_defined_at_randomization():
    rows=inventory(1,4)
    rows[0]["survived"]=1
    try:
        m.assign(rows,42)
    except ValueError as exc:
        assert "no outcome" in str(exc)
    else:
        raise AssertionError("outcome already known before randomization")
