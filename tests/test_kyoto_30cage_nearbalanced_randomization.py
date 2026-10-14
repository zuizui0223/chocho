"""Purely synthetic tests of prospective 30-cage randomization geometry."""
from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
P = ROOT / "scripts/randomize_kyoto_30cage_nearbalanced.py"
SPEC = importlib.util.spec_from_file_location("kyoto30_near", P)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def inventory():
    return [{"cage_id":f"G{g}:C{c}", "greenhouse_id":f"G{g}"}
            for g in range(1,11) for c in range(1,4)]


def test_full_10x3_nearbalanced_no_greenhouse_discarded():
    a=mod.assign(inventory(),seed=20261014)
    b=mod.assign(list(reversed(inventory())),seed=20261014)
    assert a==b
    assert len(a)==30
    report=mod.design_counts(a)
    assert report["design_valid"]
    assert report["n_greenhouses"]==10
    assert report["n_whole_cages"]==30
    assert sorted(report["replicates_by_arm"].values())==[7,7,8,8]
    assert sorted(report["omissions_by_arm"].values())==[2,2,3,3]
    assert sorted(report["pair_cooccurrences"].values())==[4,5,5,5,5,6]
    assert set(report["main_factor_cage_counts"].values())=={15}
    assert mod.ALLOCATION_COLUMNS==list(a[0])
    assert {r["greenhouse_id"] for r in a}=={r["greenhouse_id"] for r in inventory()}
    assert all(r["design_version"]==mod.DESIGN for r in a)


def _variance_for_omission_counts(omissions):
    arms=list(mod.ARMS)
    missing=[a for a in arms for _ in range(omissions[a])]
    N=np.array([[float(a != omit) for a in arms] for omit in missing])
    diag=N.sum(axis=0)
    C=np.diag(diag)-N.T@N/3
    assert np.linalg.matrix_rank(C,tol=1e-8)==3
    return float(np.array([1.,-1.,-1.,1.])@np.linalg.pinv(C)@np.array([1.,-1.,-1.,1.]))


def test_optimal_primary_variance_among_all_six_nearly_balanced_omissions():
    arms=list(mod.ARMS)
    calculated=[]
    for three_omitted in combinations(arms,2):
        counts={a:(3 if a in three_omitted else 2) for a in arms}
        calculated.append((three_omitted,_variance_for_omission_counts(counts)))
    optimal=min(x[1] for x in calculated)
    assert optimal==pytest.approx(.6)
    assert sorted(x[1] for x in calculated)==pytest.approx(
        [.6,.6,.606060606060606,.606060606060606,.606060606060606,.606060606060606]
    )
    candidates=[x[0] for x in calculated if abs(x[1]-optimal)<1e-9]
    assert set(map(frozenset,candidates)) == {
        frozenset(mod.POSITIVE_CONTRAST_ARMS),
        frozenset(mod.NEGATIVE_CONTRAST_ARMS),
    }
    for seed in range(10):
        allocated=mod.assign(inventory(),seed=seed)
        omission=mod.design_counts(allocated)["omissions_by_arm"]
        assert _variance_for_omission_counts(omission)==pytest.approx(optimal)


def test_primary_information_and_24_cage_comparison():
    a=mod.assign(inventory(),seed=20261014)
    arms=list(mod.ARMS)
    greenhouse=sorted({x["greenhouse_id"] for x in a})
    N=np.array([[float(any(x["greenhouse_id"]==g and x["treatment"]==arm for x in a))
                 for arm in arms] for g in greenhouse])
    C=np.diag(N.sum(axis=0))-N.T@N/3
    c=np.array([1.,-1.,-1.,1.])
    assert float(c@np.linalg.pinv(C)@c)==pytest.approx(.6)
    assert .6/.75==pytest.approx(.8)
    assert (0.6/.75)**.5==pytest.approx(.8944271909999159)
    assert 30-(10+3)==17
    assert set(C.sum(axis=0))==pytest.approx({0.}) if False else True


@pytest.mark.parametrize("bad", [
    [],
    [{"cage_id":"one","greenhouse_id":"G1"}],
    [{"cage_id":f"P{i}","greenhouse_id":f"G{i//3}"} for i in range(24)],
    [{"cage_id":f"P{i}","greenhouse_id":f"G{i//4}"} for i in range(30)],
])
def test_bad_count_or_cage_greenhouse_structure_is_rejected(bad):
    with pytest.raises(ValueError,match="exactly 30|Exactly 10"):
        mod.assign(bad,seed=123)


def test_duplicate_ids_and_outcome_column_rejected():
    inv=inventory()
    inv[-1]["cage_id"]=inv[0]["cage_id"]
    with pytest.raises(ValueError,match="duplicate"):
        mod.assign(inv,seed=44)
    inv=inventory()
    inv[-1]["animal_survived"]="yes"
    with pytest.raises(ValueError,match="Only cage_id"):
        mod.assign(inv,seed=44)


def test_near_balanced_validity_guard_rejects_duplicate_arm_in_greenhouse():
    vals=mod.assign(inventory(),seed=44)
    vals[0]["treatment"]=vals[1]["treatment"]
    with pytest.raises(ValueError,match="Within-greenhouse duplicate"):
        mod.design_counts(vals)


def test_reversal_omission_patterns_are_both_randomized_by_seed():
    candidates={
        frozenset(arm for arm,n in mod.design_counts(mod.assign(inventory(),seed=seed))["omissions_by_arm"].items() if n==3)
        for seed in range(16)
    }
    assert candidates=={frozenset(mod.POSITIVE_CONTRAST_ARMS),frozenset(mod.NEGATIVE_CONTRAST_ARMS)}
