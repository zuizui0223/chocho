"""Tests for exact LGM resource-climate vs land/ice domain spatial accounting."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from land_vs_climate_common_domain_accounting_v07 import (
    partition_area,area_by_flat_cell,describe,source_plan,CELL_N
)

def test_per_cell_signed_identity_changes_with_land_substrate_and_climate():
    p0=np.zeros(CELL_N,bool);p21=np.zeros(CELL_N,bool)
    h0=np.zeros(CELL_N,bool);h21=np.zeros(CELL_N,bool)
    p0[:2]=True;p21[1:3]=True
    h0[0:2]=True;h21[2]=True
    w=np.ones(CELL_N,float)
    a=partition_area(h0,h21,p0,p21,w)
    assert a["potential_area_0BP_km2"]==2
    assert a["potential_area_21ka_km2"]==1
    assert a["unchanged_physically_passable_domain_climate_change_km2"]==-1
    assert a["21ka_only_passable_land_suitable_area_km2"]==1
    assert a["0BP_only_passable_land_suitable_area_km2"]==1
    assert a["paired_total_change_km2"]==-1
    assert a["area_identity_remainder_km2"]==0

def test_climate_only_cannot_be_falsely_attributed_to_land_changes():
    p=np.zeros(CELL_N,bool);p[5:8]=True
    h0=np.zeros(CELL_N,bool);h0[5:7]=True
    h21=np.zeros(CELL_N,bool);h21[5]=True
    q=partition_area(h0,h21,p,p,np.ones(CELL_N))
    assert q["paired_total_change_km2"]==-1
    assert q["unchanged_physically_passable_domain_climate_change_km2"]==-1
    assert q["net_passable_terrain_replacement_km2"]==0

def test_geography_only_cannot_be_falsely_attributed_to_climate_change():
    p0=np.zeros(CELL_N,bool);p21=np.zeros(CELL_N,bool)
    p0[10:12]=True;p21[11:13]=True
    h0=p0.copy();h21=p21.copy()
    q=partition_area(h0,h21,p0,p21,np.ones(CELL_N))
    assert q["paired_total_change_km2"]==0
    assert q["unchanged_physically_passable_domain_climate_change_km2"]==0
    assert q["21ka_only_passable_land_suitable_area_km2"]==1
    assert q["0BP_only_passable_land_suitable_area_km2"]==1

def test_passable_land_source_constraint_violations_fail_closed():
    z=np.zeros(CELL_N,bool);w=np.ones(CELL_N)
    h0=z.copy();h0[10]=True
    with pytest.raises(ValueError):
        partition_area(h0,z,z,z,w)
    with pytest.raises(ValueError):
        partition_area(z,np.zeros(CELL_N-1),z,z,w)

def test_spherical_area_weights_downweight_high_latitude():
    w=area_by_flat_cell()
    assert w.shape==(CELL_N,) and np.all(np.isfinite(w))
    assert w[180*720]>w[0]>0
    assert w[180*720]>w[359*720]>0

def test_median_of_components_is_not_forced_to_equal_median_of_totals():
    q=partition_area(np.zeros(CELL_N),np.zeros(CELL_N),np.ones(CELL_N),np.ones(CELL_N),np.ones(CELL_N))
    a={**q,"paired_total_change_km2":-1,
          "unchanged_physically_passable_domain_climate_change_km2":-2,
          "net_passable_terrain_replacement_km2":1}
    b={**q,"paired_total_change_km2":2,
          "unchanged_physically_passable_domain_climate_change_km2":3,
          "net_passable_terrain_replacement_km2":-1}
    result=describe([a,b])
    assert result["n"]==2
    assert result["full_suitable_area_decrease"]==1
    assert result["both_climate_and_terrain_negative"]==0
