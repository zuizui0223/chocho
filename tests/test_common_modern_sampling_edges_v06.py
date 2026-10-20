"""Unit tests for exact same-edge overlap between 0BP and 21ka modern anchors."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_common_modern_sampling_edges_v06 import paired_transition_counts

GOOD="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE"
BLOCKED="PHYSICAL_ANCHOR_BLOCKED"

def samples(labels,eligible=None):
    n=len(labels)
    if eligible is None: eligible=[True]*n
    return (np.array([0. if a else np.inf for a in eligible]),
            np.array(labels,dtype=int),
            np.array([GOOD if a else BLOCKED for a in eligible],dtype=object))

def test_paired_edge_transitions_are_on_identical_pairs_not_marginal_denominators():
    # Same->different for 0-1, different->same 1-2; others remain different.
    a=samples([1,1,2,3])
    b=samples([1,2,2,3])
    q=paired_transition_counts([0,1,2,3],
            [(0,1),(1,2),(2,3),(0,3)],a,b,100)
    assert q["edge_contactable_in_both_eras_n"]==4
    assert (q["same_to_same"],q["same_to_different"],
            q["different_to_same"],q["different_to_different"])==(0,1,1,2)
    assert q["same_component_paired_difference"]==pytest.approx(0)
    assert q["same_host_component_fraction_0BP_on_exact_common_pairs"]==pytest.approx(.25)
    assert q["same_host_component_fraction_21ka_on_exact_common_pairs"]==pytest.approx(.25)

def test_ice_blocked_later_pair_is_excluded_not_different():
    a=samples([1,1,2,3])
    b=samples([1,2,2,3],[True,False,True,True])
    q=paired_transition_counts([0,1,2,3],[(0,1),(0,2),(2,3)],a,b,100)
    assert q["edge_contactable_in_both_eras_n"]==2
    assert q["edge_contactable_only_0BP_n"]==1
    assert q["edge_contactable_only_21ka_n"]==0
    assert q["same_to_different"]==0
    assert q["different_to_different"]==2

def test_never_convert_unmatchable_pair_into_zero_difference():
    a=samples([1,1],[True,True]);b=samples([1,1],[False,True])
    q=paired_transition_counts([0,1],[(0,1)],a,b,100)
    assert q["edge_contactable_in_both_eras_n"]==0
    assert q["same_component_paired_difference"] is None
    assert q["same_host_component_fraction_21ka_on_exact_common_pairs"] is None

def test_radius_sensitivity_applied_to_both_eras_and_both_endpoints():
    good=np.array([GOOD,GOOD,GOOD],dtype=object)
    a=(np.array([20.,70.,125.]),np.array([1,1,2]),good)
    b=(np.array([25.,140.,40.]),np.array([1,2,1]),good)
    e=[(0,1),(0,2),(1,2)]
    q50=paired_transition_counts([0,1,2],e,a,b,50)
    q150=paired_transition_counts([0,1,2],e,a,b,150)
    assert q50["edge_contactable_in_both_eras_n"]==0
    assert q150["edge_contactable_in_both_eras_n"]==3

def test_duplicate_node_index_fails_closed():
    s=samples([1,1])
    with pytest.raises(ValueError):
        paired_transition_counts([0,0],[(0,0)],s,s,100)
