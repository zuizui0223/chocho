"""Independent NOAA WPI/GSHHG exact source reconstruction guards."""
import pytest
for pkg in ("geopandas","shapely","numpy","pandas","pyproj"):pytest.importorskip(pkg)
import importlib.util, json
from pathlib import Path
import geopandas as gpd
from shapely.geometry import Point, box
base=Path(__file__).resolve().parents[1]/"scripts"
import sys
sys.path.insert(0,str(base))
src=base/"reconstruct_island51_from_original_NOAA_WPI_GSHHG.py"
sp=importlib.util.spec_from_file_location("original_wpi_gshhg",src)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_WPI_original_coordinate_integrity_rejects_wrong_source_size(tmp_path):
    p=tmp_path/"not_original_3630.geojson"
    p.write_text(json.dumps({"type":"FeatureCollection","features":[
      {"type":"Feature","properties":{"OBJECTID":2},"geometry":{"type":"Point","coordinates":[0,0]}}
    ]}))
    with pytest.raises(ValueError,match="3630"):
        m.load_WPI(p)

def test_local_geodesic_buffer_intersects_true_land_only():
    g=gpd.GeoDataFrame({"sibling_id":[10,15],
        "landmass_id":["10","15"],
        "geometry":[box(-.03,-.03,-.004,.03),box(.3,.3,.5,.5)]},crs=4326)
    assert m.region_overlaps(g,0,0,5000)==["10"]
    assert m.region_overlaps(g,1,1,5000)==[]
    assert m.region_overlaps(g,.35,.35,5000)==["15"]

def test_near_another_island_results_in_ambiguous_not_fabricated_single():
    g=gpd.GeoDataFrame({"sibling_id":[10,15],"landmass_id":["10","15"],
       "geometry":[box(-.03,-.03,-.004,.03),box(.004,-.03,.03,.03)]},crs=4326)
    assert m.region_overlaps(g,0,0,5000)==["10","15"]
