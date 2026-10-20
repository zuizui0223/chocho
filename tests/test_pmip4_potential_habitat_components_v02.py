import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from pmip4_potential_habitat_components_v02 import (
    components,graph_metrics,regular_cell_areas)

def test_four_connected_discards_diagonal_without_throwing_away_island():
    m=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0]],bool)
    a,n=components(m,4);b,n8=components(m,8)
    assert n==3 and n8==1 and np.count_nonzero(a)==3
    assert b[0,0]==b[1,1]==b[2,2]>0

def test_dateline_east_west_neighbors_join_only_when_land_is_suitable():
    m=np.zeros((5,6),bool);m[1,0]=m[1,-1]=1
    l,n=components(m,4)
    assert n==1 and l[1,0]==l[1,-1]
    m[1,-1]=False;m[1,-2]=True
    assert components(m,4)[1]==2

def test_periodic_eight_neighbor_can_connect_across_row_at_dateline():
    m=np.zeros((4,8),bool)
    m[1,0]=1;m[2,-1]=1
    assert components(m,4)[1]==2
    assert components(m,8)[1]==1

def test_zero_habitat_component_returns_explicit_nonapplicability():
    m=np.zeros((3,6),bool)
    d=graph_metrics(m,[30,0,-30],connectivity=4)
    assert d["components"]==0 and d["suitable_km2"]==0
    assert d["fraction_area_in_largest_component"] is None

def test_area_weighted_largest_component_and_fragmentation():
    m=np.zeros((3,8),bool)
    m[0,0:2]=True;m[1,4:7]=True
    d=graph_metrics(m,[80,0,-80],4)
    assert d["components"]==2
    assert d["suitable_grid_cells"]==5
    assert d["largest_component_cells"]==3
    expected=3/(3+2*np.cos(np.deg2rad(80)))
    assert d["fraction_area_in_largest_component"]==pytest.approx(expected,abs=1e-4)
    assert d["fraction_cells_in_largest_component"]==pytest.approx(3/5)

def test_equator_weight_bigger_than_pole_and_invalid_configs_rejected():
    areas=regular_cell_areas([0,50,80],0.5)
    assert areas[0]>areas[1]>areas[2]>0
    with pytest.raises(ValueError):components(np.ones((3,5)),connectivity=6)
    with pytest.raises(ValueError):graph_metrics(np.ones((3,5)),[1,2])
