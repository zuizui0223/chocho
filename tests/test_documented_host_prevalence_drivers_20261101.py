"""Source-free tests: never cherry-pick species within an empirical host stratum."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_documented_host_prevalence_drivers.py"


def load():
    sys.path.insert(0,str(ROOT/"scripts"))
    spec=importlib.util.spec_from_file_location("documented_host_bins",SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_complete_prevalence_bin_assignment_without_selection():
    m=load()
    f=np.array([1,2,3,5,6,10,25,50,99,101,150,201,284],int)
    mask=sum(1<<i for i in range(len(f)))
    bins,counts=m.split_by_prevalence(mask,f)
    assert sum(counts.values())==len(f)
    assert all(bins[b].bit_count()==counts[b] for b in bins)
    assert counts[0]==1
    assert counts[1]==1
    assert counts[2]==2
    assert counts[6]==1
    assert counts[7]==2
    assert counts[8]==2
    assert set(bins)==set(range(9))
    assert sum(bins.values())==mask


def test_pair_selector_filters_both_endpoints_consistently():
    m=load()
    original={
        "ia":np.array([0,0,1,2]),
        "ib":np.array([1,2,2,3]),
        "exposure":np.array([1.,2.,3.,4.]),
        "resource_delta":np.array([.1,.2,.3,.4]),
        "control_sets":{"colonial_vs_resource_beta_after_flora":[
            np.array([10.,20.,30.,40.]),
            np.array([3.,2.,4.,6.])
        ]},
    }
    include=np.array([True,False,True,False])
    subset=m.reindex_response_pairs(original,include)
    assert list(subset["ia"])==[0,1]
    assert list(subset["ib"])==[1,2]
    assert list(subset["resource_delta"])==[.1,.3]
    assert list(subset["control_sets"]["colonial_vs_resource_beta_after_flora"][0])==[10.,30.]
    assert len(original["ia"])==4
    with pytest.raises(ValueError,match="Pair mask"):
        m.reindex_response_pairs(original,np.array([True,False]))
