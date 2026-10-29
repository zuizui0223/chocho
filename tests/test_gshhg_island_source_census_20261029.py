"""GSHHG source unit dissolution and 5/20 km² threshold gates."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import pytest
# Keep original GEB paper CI independent of heavyweight exploratory GIS deps;
# dedicated geographic source workflow installs geopandas and runs all tests.
gpd=pytest.importorskip("geopandas")
from shapely.geometry import box
from audit_gshhg_island_source_census_20261029 import dissolve_original,PRIMARY_MIN_KM2,SENSITIVITY_MIN_KM2,ISLAND_MAINLAND_THRESHOLD_KM2

def test_area_threshold_order_independent():
    assert PRIMARY_MIN_KM2==5
    assert SENSITIVITY_MIN_KM2==20
    assert ISLAND_MAINLAND_THRESHOLD_KM2==7_000_000

def test_sibling_split_same_landmass_dissolves():
    g=gpd.GeoDataFrame({"sibling_id":[4,4,9]},geometry=[
        box(100,10,100.02,10.02),box(100.02,10,100.04,10.02),
        box(150,-20,150.01,-19.98)],crs=4326)
    out=dissolve_original(g)
    assert len(out)==2
    assert sorted(out.n_original_parts.tolist())==[1,2]

def test_missing_original_id_fails_closed():
    g=gpd.GeoDataFrame({"x":[1]},geometry=[box(0,0,.1,.1)],crs=4326)
    with pytest.raises(ValueError,match="ID field missing"):
        dissolve_original(g)
