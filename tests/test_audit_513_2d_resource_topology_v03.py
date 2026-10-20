"""Pure unit tests: source-independent hostset topology and scientific denominator guards."""
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_513_2d_resource_topology_v03 import (
    compute_masks,group_metrics,summarize,accepted,GRID,N
)

def test_4_and_8_neighbor_land_components_differ_without_confusing_area():
    a=np.zeros(GRID,dtype=bool);b=np.zeros(GRID,dtype=bool)
    a[50,50]=a[51,51]=True
    b[50,50]=b[51,50]=True
    lat=89.75-0.5*np.arange(360)
    out=group_metrics(a,b,lat)
    assert out["component_count_0BP_4"]==2
    assert out["component_count_0BP_8"]==1
    assert out["component_count_21ka_4"]==1
    assert out["component_count_21ka_8"]==1
    assert out["component_delta_4"]==-1
    assert out["component_delta_8"]==0

def test_same_plant_model_used_by_multiple_insects_yields_identical_climate_graph():
    mask_a=np.zeros(N,bool);mask_b=np.zeros(N,bool)
    mask_a[100:105]=True;mask_b[101:108]=True
    x=np.stack([np.packbits(mask_a),np.packbits(mask_b)])
    p=compute_masks([0,1],x,x)
    q=compute_masks([1,0],x,x)
    assert np.array_equal(p[0],q[0])
    assert np.array_equal(p[1],q[1])
    assert p[0].sum()==8
    assert p[0].shape==GRID

def test_summarize_cannot_confuse_loss_of_components_with_less_fragmentation():
    records=[dict(area_delta_km2=-2,component_delta_4=-2,component_delta_8=-1,
                  largest_fraction_delta_4=-.3,largest_fraction_delta_8=-.2,
                  largest_area_delta_4_km2=-2,largest_area_delta_8_km2=-2),
             dict(area_delta_km2=2,component_delta_4=1,component_delta_8=1,
                  largest_fraction_delta_4=.1,largest_fraction_delta_8=.2,
                  largest_area_delta_4_km2=2,largest_area_delta_8_km2=2)]
    d=summarize(records)
    assert d["n"]==2
    assert d["four_neighbor_components_decline"]==1
    assert d["four_neighbor_largest_fraction_declines"]==1
    assert d["four_neighbor_components_and_largest_fraction_both_decline"]==1

def test_wcvp_id_cannot_lose_precision():
    assert accepted("400000.0")=="400000"
    for val in ("0","NaN","-1","3.2"):
        with pytest.raises(ValueError):accepted(val)

def test_duplicate_host_models_rejected():
    a=np.zeros((2,32400),dtype=np.uint8)
    with pytest.raises(ValueError):compute_masks([0,0],a,a)
