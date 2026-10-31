"""Exact gain-richness marginal preservation for 239 consumer resource geography."""
import importlib.util
from pathlib import Path
import random

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_global_resource_gain_fixed_margins.py"


def load():
    spec=importlib.util.spec_from_file_location("fixed_gain_margins",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_added_resources_swap_preserves_species_and_region_margins():
    m=load()
    # Every region has the same native species 0. Added resources are
    # species 1/2/3 and do not overwrite the native presence bit.
    native=[1,1,1,1]
    additions=[2,4,8,2]
    original=additions.copy()
    before=m.marginal_counts(native,additions,4)
    edges=m.additions_list(additions)
    swaps=m.swap_new_resources(native,additions,edges,random.Random(13),800)
    assert swaps>0
    assert m.marginal_counts(native,additions,4)==before
    assert sorted(edges)==sorted(m.additions_list(additions))
    assert all((a&n)==0 for a,n in zip(additions,native))
    assert sum(x.bit_count() for x in original)==sum(x.bit_count() for x in additions)
    assert sum(x.bit_count() for x in native)==4


def test_pair_similarity_has_fixed_native_reference():
    m=load()
    native=[1,1,1,1]
    added=[2,4,8,2]
    selected=[0,1,2,3]
    pairs,base=m.pair_baseline(native,selected)
    assert len(pairs)==6
    assert base==pytest.approx(1)
    observed=m.current_mean(native,added,pairs)-base
    assert observed<0
    _,counts=m.marginal_counts(native,added,4)
    assert counts==[0,2,1,1]


def test_immutable_resource_matrix_shares_paper_source():
    m=load()
    assert m.OUTPUT.name=="resource_gain_fixed_margins_v1.json"
    assert m.MATRIX.name.endswith(".matrix.csv.gz")
    assert m.RECEIPT.name=="beta_baseline_v3.json"



def test_within_world_division_swaps_preserve_species_gains_inside_each_division():
    m=load()
    native=[1]*6
    additions=[2,4,8,2,4,8]
    strata=["A","A","A","B","B","B"]
    edges=m.additions_list(additions)
    original=additions.copy()
    accepted=m.swap_new_resources(
        native,additions,edges,random.Random(17),2500,
        region_strata=strata)
    assert accepted>0
    for label in set(strata):
        former=[v if strata[i]==label else 0 for i,v in enumerate(original)]
        latter=[v if strata[i]==label else 0 for i,v in enumerate(additions)]
        assert m.marginal_counts(native,former,4)==m.marginal_counts(native,latter,4)
    assert sorted(edges)==sorted(m.additions_list(additions))
