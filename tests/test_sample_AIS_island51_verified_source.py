"""Pixel-level AIS georeference tests independent of raster download."""
from __future__ import annotations
import pytest
for module in ("rasterio","numpy"):pytest.importorskip(module)
import sys,importlib.util
from pathlib import Path
import numpy as np,rasterio
from rasterio.transform import Affine
root=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(root))
sp=importlib.util.spec_from_file_location("AIS_real",root/"sample_AIS_island51_verified_source.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_worldfile_source_0p005_degree_overview_doubles_cell_size():
    meta={"original_source_world_file":"0.005\n0.0\n0.0\n-0.005\n-179.9975\n84.9975\n"}
    t,x=m.worldfile_transform(meta,36003,16999)
    assert abs(x["overview_dlon_deg"]-.01)<1e-12
    assert abs(x["overview_dlat_deg"]+.01)<1e-12
    assert abs(x["left"]+180)<1e-12
    assert abs(x["top"]-85)<1e-12
    assert abs(x["bottom"]+84.99)<.02

def test_raster_masked_nodata_is_not_zero_shipping(tmp_path):
    p=tmp_path/"small.tif"
    arr=np.zeros((50,50),dtype=np.int32)
    arr[23:28,23:28]=125
    arr[25,25]=2147483647
    tr=Affine(.01,0,0,0,-.01,5)
    with rasterio.open(p,"w",driver="GTiff",height=50,width=50,count=1,crs="EPSG:4326",
         transform=tr,nodata=2147483647,dtype="int32") as f:f.write(arr,1)
    with rasterio.open(p) as d:
        a=m.sample_port_area(d,tr,.255,4.745,5)
        assert a["status"]=="AIS_SOURCE_VALID"
        assert a["n_valid"]>0
        assert a["n_positive"]>0
        assert a["n_positive"]<a["n_valid"]
        assert a["max_AIS_positions"]==125
        assert a["mean_log1p_AIS_positions"]>0
