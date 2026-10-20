import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from aggregate_insect_host_climatic_graph_v02 import (
    host_model_status, union_known_host_potential,compare_component_metrics,SHAPE,BITS
)

def masks():
    a=np.zeros(BITS,bool);b=np.zeros(BITS,bool)
    a[0:3]=1;b[2:5]=1
    q=np.stack([np.packbits(a),np.packbits(b)]).astype("u1")
    return q,a,b

def test_or_of_two_complete_known_hosts_is_not_addition_of_probabilities():
    packed,a,b=masks()
    c,d=union_known_host_potential({"1","2"},{"1":0,"2":1},packed,packed)
    assert c.shape==d.shape==SHAPE
    assert np.sum(c)==5
    assert np.array_equal(c.reshape(-1),a|b)

def test_missing_host_model_must_not_be_assumed_absent():
    src={"1":0,"2":1}
    assert host_model_status({"1","2"},src)=="ALL_RECORDED_HOST_MODELS_AVAILABLE"
    assert host_model_status({"1","3"},src)=="PARTIAL_RECORDED_HOST_MODELS_NOT_EVALUATED"
    assert host_model_status({"3"},src)=="NO_RECORDED_HOST_MODEL_AVAILABLE"
    with pytest.raises(ValueError):
        union_known_host_potential({"1","3"},src,*([masks()[0]]*2))

def test_4_neighbor_potential_area_not_genetic_distance():
    a=np.zeros(SHAPE,bool);b=np.zeros(SHAPE,bool)
    a[1,1:4]=1;b[1,1:2]=1
    summary=compare_component_metrics(a,b,89.75-np.arange(360)*.5)
    assert summary["area_km2_delta_21ka_minus_0BP"]<0
    assert summary["component_count_0BP"]==1
    assert summary["component_count_21ka"]==1
    assert summary["sites_no_longer_suitable_21ka"]==2
    assert summary["sites_suitable_both"]==1
