"""Invariant tests for comparing *identical* 0BP/21ka source original edges under host vs other-plant controls."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_0bp_matched_otherplant_patch_contact_v08 import full_quad_comparison
GOOD="PHYSICAL_LAND_WITH_NEAREST_HOST_CLIMATE"
BLOCKED="PHYSICAL_ANCHOR_BLOCKED"

def cloud(labels,ok=None):
    if ok is None:ok=[True]*len(labels)
    return (np.array([0. if x else np.inf for x in ok]),
            np.array(labels,dtype=int),
            np.array([GOOD if x else BLOCKED for x in ok],dtype=object))

def test_identical_4way_source_pairs_used_for_difference_in_difference():
    nodes=[0,1,2,3]
    edges=[(0,1),(1,2),(0,2),(2,3)]
    actual0=cloud([1,1,2,3])
    actual21=cloud([1,2,2,3])
    control0=cloud([1,1,1,3])
    control21=cloud([1,2,1,3])
    a=full_quad_comparison(nodes,edges,actual0,actual21,control0,control21,100)
    assert a['fourway_exact_common_original_edges_n']==4
    assert a['actual_0BP_same_fraction']==pytest.approx(.25)
    assert a['actual_21ka_same_fraction']==pytest.approx(.25)
    assert a['alternate_0BP_same_fraction']==pytest.approx(.75)
    assert a['alternate_21ka_same_fraction']==pytest.approx(.25)
    assert a['actual_paired_delta']==pytest.approx(0)
    assert a['alternate_paired_delta']==pytest.approx(-.5)
    assert a['actual_minus_alternate_paired_delta']==pytest.approx(.5)
    assert a['modern_0BP_baseline_same_patch_gap']==pytest.approx(-.5)

def test_unavailable_otherplant_anchor_must_not_be_counted_as_host_split():
    nodes=[0,1,2,3]
    edge=[(0,1),(0,2),(1,2),(2,3)]
    real0=cloud([1,1,2,3]); real21=cloud([1,2,2,3])
    alt0=cloud([1,1,2,3])
    alt21=cloud([1,1,2,3],[True,False,True,True])
    a=full_quad_comparison(nodes,edge,real0,real21,alt0,alt21,100)
    assert a['fourway_exact_common_original_edges_n']==2
    assert a['all_original_edges_n']==4
    assert a['actual_0BP_same_fraction']==0
    assert a['alternate_21ka_same_fraction']==0

def test_no_fourway_source_edges_is_none_not_zero_climate_effect():
    first=cloud([1,1]);second=cloud([1,1],[False,True])
    z=full_quad_comparison([0,1],[(0,1)],first,first,second,second,100)
    assert z['fourway_exact_common_original_edges_n']==0
    assert z['actual_minus_alternate_paired_delta'] is None

def test_prespecified_50_100_150_km_contact_thresholds_change_sample_not_response():
    first=cloud([1,1]);further=(np.array([30.,140.]),np.array([1,1]),np.array([GOOD,GOOD],object))
    z=full_quad_comparison([0,1],[(0,1)],first,first,further,further,100)
    q=full_quad_comparison([0,1],[(0,1)],first,first,further,further,150)
    assert z['fourway_exact_common_original_edges_n']==0
    assert q['fourway_exact_common_original_edges_n']==1
    assert q['actual_minus_alternate_paired_delta']==0

def test_duplicate_modern_locality_is_invalid():
    s=cloud([1,1])
    with pytest.raises(ValueError):
        full_quad_comparison([1,1],[(1,1)],s,s,s,s,100)
