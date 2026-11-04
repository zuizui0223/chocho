"""Source-free GIS vector/raster math, missingness, noncausal N grid QA."""
import pytest
for dep in ("numpy","pandas","statsmodels","scipy","rasterio","pyproj","pycountry"):
    pytest.importorskip(dep)
import importlib.util, sys, numpy as np
from pathlib import Path
import rasterio
from rasterio.transform import from_origin

S=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(S))
sp=importlib.util.spec_from_file_location("port_cirad_grid",S/"exploratory_cirad_n_excretion_ports.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)


def test_pixel_circle_keeps_zero_separate_from_nodata(tmp_path):
    p=tmp_path/"synthetic.tif"
    vals=np.full((60,60),10,dtype=np.float32)
    vals[27:33,27:33]=0
    vals[:4,:]=-9999
    with rasterio.open(p,"w",driver="GTiff",crs="ESRI:54035",width=60,height=60,
        count=1,dtype="float32",nodata=-9999,
        transform=from_origin(-90000,90000,3000,3000)) as dst:
        dst.write(vals,1)
    with rasterio.open(p) as src:
        result=m.circle_sample(src,0,0,30000)
    assert result["circle_cells"]>0
    assert result["valid_cells"]>10
    assert 0<result["mean"]<10
    assert 0<result["zero_share_valid"]<1


def test_non_extrapolation_when_outside_source(tmp_path):
    p=tmp_path/"mini.tif"
    with rasterio.open(p,"w",driver="GTiff",crs="ESRI:54035",width=10,height=10,
        count=1,dtype="float32",nodata=-9999,
        transform=from_origin(0,30000,3000,3000)) as dst:
        dst.write(np.ones((10,10),dtype=np.float32),1)
    with rasterio.open(p) as src:
        v=m.circle_sample(src,-500000,500000,5000)
        assert v["valid_cells"]==0
        assert v["mean"] is None


def test_source_layer_metadata_fixed_to_N_and_not_P():
    assert m.FILES["human"][1]=="global_grid_human_N_excretion_2022"
    assert m.FILES["crop"][1]=="global_grid_N_supply_to_cropland_2022"
    assert m.RADIUS_MAIN==25000
