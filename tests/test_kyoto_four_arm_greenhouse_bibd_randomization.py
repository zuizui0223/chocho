"""Mathematical and data-integrity checks ONLY; all adult fates are invented.

The 2012 original paper reported 10 greenhouses x 3 cages, but none of
these synthetic original IDs or individuals is an actual historical record.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
P=ROOT/"scripts/randomize_kyoto_four_arm_greenhouse_bibd.py"
SPEC=importlib.util.spec_from_file_location("kyoto_bibd_random",P)
assert SPEC and SPEC.loader
mod=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def inventory():
    return [
        {"cage_id":f"G{g}:C{c}","greenhouse_id":f"G{g}"}
        for g in range(1,9) for c in range(1,4)
    ]


def test_exact_four_arm_8x3_incomplete_design_and_pair_balance():
    a=mod.assign(inventory(),seed=20261014)
    b=mod.assign(list(reversed(inventory())),seed=20261014)
    assert a==b
    x=mod.counts(a)
    assert x["complete_v4_b8_k3_r6_lambda4"]
    assert x["n_greenhouses"]==8
    assert x["n_whole_cages"]==24
    assert set(x["replicates_by_treatment"].values())=={6}
    assert set(x["omissions_by_treatment"].values())=={2}
    assert len(x["pair_cooccurrences_by_treatment_pair"])==6
    assert set(x["pair_cooccurrences_by_treatment_pair"].values())=={4}
    assert len({r["treatment"] for r in a})==4
    assert all("fate" not in r for r in a)
    assert mod.ALLOCATION_COLUMNS==list(a[0])


def test_randomized_omissions_and_within_greenhouse_slots_are_repeatable():
    a=mod.assign(inventory(),seed=20261014)
    b=mod.assign(inventory(),seed=20261016)
    assert a != b
    for x in (a,b):
        observed=Counter(r["greenhouse_id"] for r in x)
        assert set(observed.values())=={3}
        assert all(len([v for v in x if v["greenhouse_id"]==g])==3 for g in observed)


@pytest.mark.parametrize("inventory_rows",[
    [],
    [{"cage_id":"x","greenhouse_id":"G1"}],
    [{"cage_id":f"X{i}","greenhouse_id":f"G{i//3}"} for i in range(30)],
    [{"cage_id":f"X{i}","greenhouse_id":f"G{i//4}"} for i in range(24)],
])
def test_no_silent_selection_or_false_historical_facility_eligibility(inventory_rows):
    with pytest.raises(ValueError,match="exactly 24|exactly 8"):
        mod.assign(inventory_rows,seed=11)


def test_duplicate_cage_and_outcome_column_fail_closed():
    rows=inventory()
    rows[-1]["cage_id"]=rows[0]["cage_id"]
    with pytest.raises(ValueError,match="duplicate"):
        mod.assign(rows,seed=1)
    rows=inventory()
    rows[0]["success"]=1
    with pytest.raises(ValueError,match="no outcomes"):
        mod.assign(rows,seed=1)


def test_randomizer_does_not_admit_unbalanced_guessed_allocation():
    a=mod.assign(inventory(),seed=1)
    a[0]["treatment"]=a[1]["treatment"]
    with pytest.raises(ValueError,match="Duplicate treatment"):
        mod.counts(a)


def test_design_incidence_exact_information_multiplier():
    import numpy as np
    a=mod.assign(inventory(),seed=20261014)
    incidence=np.array([
        [int(any(x["greenhouse_id"]==g and x["treatment"]==arm for x in a))
         for arm in mod.ARMS]
        for g in sorted({x["greenhouse_id"] for x in a})
    ],dtype=float)
    expected=np.full((4,4),4.0)
    np.fill_diagonal(expected,6)
    assert incidence.T@incidence==pytest.approx(expected)
    information=6*np.eye(4)-incidence.T@incidence/3
    c=np.array([+1,-1,-1,+1],dtype=float)
    assert c@np.linalg.pinv(information)@c==pytest.approx(.75)
    assert sum(c)==0


def test_source_based_design_receipt_math_remains_independently_reproducible():
    receipt = ROOT / "provenance/reviewer_defenses/results/kyoto_four_arm_greenhouse_bibd_information_20261014.json"
    x = json.loads(receipt.read_text(encoding="utf-8"))
    design = x["bibd"]
    assert x["status"] == "MATHEMATICAL_DESIGN_INFORMATION_NOT_BUTTERFLY_DATA"
    assert design["independent_greenhouses"] * design["whole_cages_each"] == 24
    assert design["times_each_treatment_assigned"] * design["treatments"] == 24
    assert design["greenhouses_per_treatment_pair"] * (design["treatments"] - 1) == (
        design["times_each_treatment_assigned"] * (design["whole_cages_each"] - 1)
    )
    assert x["primary_contrast_variance_multiplier_per_sigma_squared"] == pytest.approx(.75)
    assert x["extra_variance_fraction_vs_24cage_complete_blocks"] == pytest.approx(.125)
    assert x["relative_efficiency_vs_24cage_complete_blocks"] == pytest.approx(8/9)
    assert x["hypothetical_precision"]["model_based_nominal_95pct_halfwidth"] == pytest.approx(
        2.1603686564610127*.2*(.75**.5)
    )
    assert not x["newly_observed_butterfly_effect"]
    assert not x["current_greenhouse_inventory_provided"]
    assert not x["source_families_independently_replicated"]
    assert x["original_GEB_main_unchanged"]
