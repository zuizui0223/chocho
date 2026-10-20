"""Regression tests for real PMIP4 0 and 21ka land/ice source interpretation."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_pmip4_original_landice_0_21ka_v01 import (
    icefree,component_ids,cell_area,nearest_grid_indices,area)

def test_only_ice_free_land_is_terrestrial_passability():
    land=np.array([[1,1,0,0],[0,1,1,1]],bool)
    ice=np.array([[0,1,1,0],[0,0,0,1]],bool)
    assert icefree(land,ice).astype(int).tolist()==[[1,0,0,0],[0,1,1,0]]

def test_east_west_dateline_seam_merges_4_neighbor_components():
    mask=np.zeros((4,8),bool)
    mask[1,0]=mask[1,-1]=True
    comp,n=component_ids(mask)
    assert n==1 and comp[1,0]==comp[1,-1]>0

def test_diagonal_is_not_traversable_land_connection():
    mask=np.array([[1,0,0],[0,1,0],[0,0,1]],bool)
    _,n=component_ids(mask)
    assert n==3

def test_ice_bridge_blocks_4_neighbor_land_connectivity():
    # Avoid the longitude seam: first and last columns are neighbors.
    mask=np.array([[0,0,0,0,0],[0,1,1,1,0],[0,0,0,0,0]],bool)
    assert component_ids(mask)[1]==1
    mask[1,2]=False
    assert component_ids(mask)[1]==2

def test_area_weighting_decreases_toward_pole_not_pixel_counts():
    lat=np.array([0,60,80])
    a=cell_area(lat)
    assert np.all(a>0) and a[0]>a[1]>a[2]
    mask=np.ones((3,4),bool)
    assert area(mask,a)==pytest.approx(4*sum(a))

def test_both_near_dateline_and_negative_longitude_index_correctly():
    lat=np.array([0,0,0],float)
    lon=np.array([0,-0.16666667,359.83333333])
    r,c=nearest_grid_indices(lat,lon)
    assert np.array_equal(c,[0,2159,2159])
    assert np.all(r==540)

def test_out_of_range_coords_not_silently_wrapped():
    with pytest.raises(ValueError):
        nearest_grid_indices([100.],[-10.])
    with pytest.raises(ValueError):
        nearest_grid_indices([0.],[np.nan])
