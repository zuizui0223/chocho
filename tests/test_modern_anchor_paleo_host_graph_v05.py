"""Unit tests for modern geographic anchors onto source-verified glacial climate graphs."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_modern_anchor_paleo_host_graph_v05 import (
    latlon_to_cell,grid_cell_centers,contact_modern_anchors,
    contact_summary,build_terrain_contact_index,SHAPE,R
)
from pmip4_potential_habitat_components_v02 import components

def test_0p5_grid_coordinate_seam_and_hemi():
    xy=[[89.75,-179.75],[0.25,0.25],[-89.75,179.75]]
    cell=latlon_to_cell(xy)
    assert cell.tolist()==[0,179*720+360,359*720+719]
    centers=grid_cell_centers(cell)
    assert np.allclose(centers,xy,atol=1e-9)

def test_physical_ice_or_sea_anchor_cannot_snap_across_close_water():
    physical=np.zeros(SHAPE,bool)
    potential=np.zeros(SHAPE,bool)
    # Single physical land cell on first island and distinct nearby island.
    physical[170,350]=True;physical[170,352]=True
    potential[170,352]=True
    plabel,_=components(physical,4)
    source=build_terrain_contact_index(potential,plabel)
    # A locality in first physically disconnected island must not snap to host
    # on second even when both centers <150 km on the spherical surface.
    xyz=[grid_cell_centers([170*720+350])[0]]
    d,comp,status=contact_modern_anchors(xyz,plabel,source)
    assert status[0]=="NO_HOST_POTENTIAL_ON_SAME_LANDMASS"
    assert np.isinf(d[0]) and comp[0]==0

def test_ocean_ice_covered_anchor_does_not_become_genetic_barrier():
    physical=np.zeros(SHAPE,bool);physical[200,410]=True
    plabel,_=components(physical,4)
    model=np.zeros(SHAPE,bool);model[200,410]=True
    source=build_terrain_contact_index(model,plabel)
    coast=grid_cell_centers([201*720+410])[0]
    d,comp,status=contact_modern_anchors([coast],plabel,source)
    assert status.tolist()==["PHYSICAL_ANCHOR_BLOCKED"]
    assert not np.isfinite(d[0])

def test_nearest_valid_land_potential_and_edge_decomposition():
    physical=np.zeros(SHAPE,bool);physical[190,360:365]=True
    potential=np.zeros(SHAPE,bool)
    potential[190,360]=True;potential[190,362]=True;potential[190,364]=True
    plabel,_=components(physical,4)
    source=build_terrain_contact_index(potential,plabel)
    xyz=grid_cell_centers([190*720+360,190*720+361,190*720+364])
    dist,patch,status=contact_modern_anchors(xyz,plabel,source)
    assert np.all(np.isfinite(dist))
    assert patch[0]!=patch[2]
    r=contact_summary([0,1,2],xyz,[(0,1),(0,2),(1,2)],dist,patch,status,100)
    assert r["edges_n"]==3 and r["edge_both_anchor_contacts_n"]==3
    assert r["edge_different_host_climate_components_n"]>=1
    assert r["edge_neither_anchor_contact_n"]==0
    assert r["edge_same_conditional_fraction"]+r["edge_different_host_climate_components_n"]/3==pytest.approx(1)

def test_no_eligible_anchor_does_not_fabricate_disconnected_edges():
    p=np.full(3,np.inf)
    components=np.zeros(3,int)
    status=np.array(["PHYSICAL_ANCHOR_BLOCKED"]*3)
    result=contact_summary([10,20,30],[[0,0],[0,1],[0,2]],
                           [(10,20),(20,30)],p,components,status,100)
    assert result["edge_both_anchor_contacts_n"]==0
    assert result["edge_neither_anchor_contact_n"]==2
    assert result["edge_same_conditional_fraction"] is None
    assert result["edge_same_all_original_fraction"]==0.

def test_nearest_land_assignment_radius_is_enforced_after_nearest_query():
    physical=np.zeros(SHAPE,bool);physical[180,360:365]=True
    potential=np.zeros(SHAPE,bool);potential[180,364]=True
    plabel,_=components(physical,4)
    source=build_terrain_contact_index(potential,plabel)
    pos=grid_cell_centers([180*720+360,180*720+364])
    d,patch,status=contact_modern_anchors(pos,plabel,source)
    assert d[0]>100
    r=contact_summary([0,1],pos,[(0,1)],d,patch,status,100)
    assert r["edge_one_anchor_contact_n"]==1
    assert r["edge_both_anchor_contacts_n"]==0
