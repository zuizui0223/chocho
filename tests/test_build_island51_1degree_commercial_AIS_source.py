"""Original World Bank AIS ZIP must be source verified and missing not zero."""
import pytest
rasterio=pytest.importorskip("rasterio")
import numpy as np
import hashlib,importlib.util
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"build_island51_1degree_commercial_AIS_source.py"
sp=importlib.util.spec_from_file_location("aisgrid",f)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_source_hash_is_binary_and_not_filename(tmp_path):
    p=tmp_path/"sourcedata.zip";p.write_bytes(b"geography-not-traffic")
    assert m.hash_file(p)==hashlib.sha256(b"geography-not-traffic").hexdigest()
    assert m.hash_file(p)!=m.ORIGINAL_SHA

def test_original_raster_averaging_never_relabels_nodata_as_0(tmp_path):
    p=tmp_path/"ais_test.tif"
    from rasterio.transform import from_origin
    arr=np.array([[0,2,-9999,2],[10,30,-9999,4],
                  [0,0,8,0],[6,6,10,10]],dtype=np.float32)
    with rasterio.open(p,"w",driver="GTiff",width=4,height=4,count=1,
        crs="EPSG:4326",transform=from_origin(-2,2,1,1),
        dtype="float32",nodata=-9999) as w:
        w.write(arr,1)
    from rasterio.enums import Resampling
    with rasterio.open(p) as r:
        v=r.read(1,out_shape=(2,2),masked=True,resampling=Resampling.average,out_dtype="float32")
    assert np.isclose(v[0,0],10.5)
    assert np.isclose(v[0,1],3.0)
    assert not np.ma.getmaskarray(v).all()
