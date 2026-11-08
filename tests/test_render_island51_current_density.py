"""Circuit identity and superposition tests independent of full global sources."""
import importlib.util
from pathlib import Path

import numpy as np
from scipy import sparse

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"render_island51_current_density.py"
spec=importlib.util.spec_from_file_location("current_map",SCRIPT)
current=importlib.util.module_from_spec(spec)
spec.loader.exec_module(current)


def test_sorted_identity_matches_brute_force_and_gauge_shift():
    p=np.array([[0,0,0],[0,1,2],[0,0,0]],float)
    u,v=np.array([0,1]),np.array([1,2])
    expected=[]
    for a,b in zip(u,v):
        drops=p[a]-p[b]
        expected.append(np.mean([abs(drops[i]-drops[j])
                                 for i in range(3) for j in range(i+1,3)]))
    np.testing.assert_allclose(current.pair_mean_absolute_edge_current(p,u,v,np.ones(2)),
                               expected)
    np.testing.assert_allclose(current.pair_mean_absolute_edge_current(p+17,u,v,np.ones(2)),
                               expected)


def test_three_node_series_gives_two_thirds_mean_edge_current():
    u,v=np.array([0,1]),np.array([1,2])
    p,res=current.solve_island_potentials(3,u,v,np.ones(2),
                                           sparse.eye(3,format="csr"))
    edge=current.pair_mean_absolute_edge_current(p,u,v,np.ones(2))
    np.testing.assert_allclose(edge,[2/3,2/3],rtol=1e-12)
    np.testing.assert_allclose(current.node_through_current(edge,u,v,3),
                               [1/3,2/3,1/3],rtol=1e-12)
    assert res < 1e-10


def test_parallel_cycle_and_distributed_port_sources():
    u,v=np.array([0,1,2,0]),np.array([1,2,3,3])
    g=np.array([2.,1.,3.,1.])
    weights=sparse.csr_matrix([[.5,0,0,.5],[0,1,0,0],[0,0,1,0]])
    p,res=current.solve_island_potentials(4,u,v,g,weights)
    observed=current.pair_mean_absolute_edge_current(p,u,v,g)
    brute=np.zeros(len(g))
    for i in range(3):
        for j in range(i+1,3):
            d=p[:,i]-p[:,j]
            brute+=g*np.abs(d[u]-d[v])
    np.testing.assert_allclose(observed,brute/3,rtol=1e-12,atol=1e-12)
    assert res < 1e-10
