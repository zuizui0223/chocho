"""Circuit effective resistances and LOIO training guardrails, no web calls."""
import pytest
for package in ("numpy","pandas"):pytest.importorskip(package)
import sys,importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(root))
sp=importlib.util.spec_from_file_location("island_circuit",root/"exploratory_island51_circuit_network_holdout.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def laplacian(n,edges):
    L=np.zeros((n,n))
    for a,b,r in edges:
        g=1./r
        L[a,a]+=g;L[b,b]+=g
        L[a,b]-=g;L[b,a]-=g
    return L

def test_simple_series_resistance_two_steps_equals_two():
    L=laplacian(3,[(0,1,1),(1,2,1)])
    r=m.resistance_from_laplacian(L)
    assert np.isclose(r[0,2],2.0,rtol=1e-9)
    assert np.isclose(r[0,1],1.0,rtol=1e-9)
    assert np.allclose(r,r.T) and np.allclose(np.diag(r),0)

def test_parallel_paths_four_node_ring_effective_resistance():
    L=laplacian(4,[(0,1,1),(1,2,1),(2,3,1),(3,0,1)])
    R=m.resistance_from_laplacian(L)
    assert np.isclose(R[0,2],1.0,atol=1e-10)
    assert np.isclose(R[0,1],0.75,atol=1e-10)

def test_graph_k_nearest_and_MST_bridge_and_distance_only():
    coords=np.array([0,1,2,20,21,22],dtype=float)
    d=np.abs(coords[:,None]-coords[None,:])*30+np.eye(len(coords))*0
    np.fill_diagonal(d,0)
    L,info=m.adjacency(d,1,1)
    assert info["mst_bridge_edges_added"]>=1
    assert info["n_total_edges"]>=5
    assert info["max_edge_geodesic_km"]>=500
    assert m.resistance_from_laplacian(L)[0,5]>0

def test_leave_one_island_out_never_uses_heldout_pairs_for_training():
    n=51
    ids=[str(i) for i in range(n)]
    rows=[]
    for a in range(n):
        for b in range(a+1,n):
            d=30*(b-a)+1
            rows.append({
                "island_a":str(a),"island_b":str(b),
                "i":a,"j":b,
                "port_distance_km":d,
                "host_jaccard":float((a+b)%13)/40,
                "earlier_butterfly_jaccard":float(a%5)/12,
                "new_butterfly_jaccard":float((a+b)%7)/25,
            })
    p=pd.DataFrame(rows)
    base=m.make_features(p,None,response="host_jaccard")
    d=np.array([[abs(i-j)*30+1 if i!=j else 0 for i in range(n)] for j in range(n)])
    L,_=m.adjacency(d,3,1)
    R=m.resistance_from_laplacian(L)
    aug=m.make_features(p,R,response="host_jaccard")
    z=m.holdout(p,base,aug,"host_jaccard",ids)
    assert z["n_islands_left_out"]==51
    assert z["n_predictions_heldout_island_pairs"]==2550
    assert z["n_unique_observed_pairs"]==1275
    assert all(r["train_pairs"]==1225 and r["test_pairs"]==50 for r in z["per_island_LOIO"])
    assert z["MSE_baseline"]>=0 and z["MSE_augmented"]>=0

def test_disconnected_graph_is_rejected():
    with pytest.raises(ValueError,match="disconnected"):
        m.resistance_from_laplacian(laplacian(4,[(0,1,1),(2,3,1)]))
