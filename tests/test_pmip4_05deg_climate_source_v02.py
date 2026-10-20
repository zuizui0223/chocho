"""Source-only 0.5-degree ICE-6G_C / CHELSA compatibility checks."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_pmip4_icefree_05deg_climate_coverage_v02 import (
    source_metrics, coarse_pmip_icefree, climate_05deg)

def test_mask_and_climate_count_intersections_are_not_interchangeable():
    a=np.zeros((360,720),dtype=bool)
    b=np.zeros_like(a);present=np.zeros_like(a);ancient=np.zeros_like(a)
    a[11,11]=1
    b[11,11]=1
    b[11,12]=1
    present[11,11]=1
    ancient[11,12]=1
    out=source_metrics(a,b,present,ancient)
    assert out["source"]["PMIP_0ka_ice_free_land_cells"]==1
    assert out["source"]["PMIP_21ka_ice_free_land_cells"]==2
    assert out["source"]["newly_ice_free_land_at21ka_05deg_cells"]==1
    assert out["coverage"]["21ka_climate_available_at_newly_icefree_land"]==1
    assert out["coverage"]["21ka_climate_missing_at_retained_icefree_land"]==1
    assert out["coverage"]["0ka_climate_missing_0ka_icefree_land"]==0

def test_pmip_longitude_roll_and_north_up_transform():
    # Synthetic 10arcmin longitude 0..360, latitude south->north.
    free=np.zeros((1080,2160),dtype=bool)
    # Source rows 1077-1079 northmost band. Source longitudes 180-180.5
    # are placed in westernmost output 0.5deg longitude bin.
    free[1077:1080,1080:1083]=True
    x=coarse_pmip_icefree(free)
    assert x.shape==(360,720)
    assert x.dtype.kind in "iu"
    assert x[0,0]==9
    assert x.sum()==9
    # 0 longitude should arrive at output col 360.
    free=np.zeros((1080,2160),bool)
    free[1077:1080,0:3]=True
    y=coarse_pmip_icefree(free)
    assert y[0,360]==9 and y.sum()==9

def test_majority_any_and_all_masks_are_ordered():
    g=np.zeros((1080,2160),bool)
    g[0,0:1]=1; g[3,0:5]=1; g[6:9,0:3]=1
    x=coarse_pmip_icefree(g)
    assert np.count_nonzero(x>=1)==3
    assert np.count_nonzero(x>=5)==2
    assert np.count_nonzero(x==9)==1

def test_bad_physical_grid_rejected():
    with pytest.raises(ValueError):
        coarse_pmip_icefree(np.ones((359,720),bool))
    with pytest.raises(ValueError):
        source_metrics(np.ones((10,10),bool),
                       *[np.ones((10,10),bool) for _ in range(3)])

def test_missing_CHELSA_is_a_missing_data_status_not_unsuitable():
    arr=np.ones((360,720),bool)
    observed=np.ones_like(arr)
    observed[0,0]=0
    values=source_metrics(arr,arr,observed,observed)
    assert values["coverage"]["21ka_climate_missing_at_retained_icefree_land"]==1
    assert values["coverage"]["21ka_climate_available_at_retained_icefree_land"]==360*720-1
