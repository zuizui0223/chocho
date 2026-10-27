"""Offline climate-consistency and missingness tests, no GBIF/WorldClim download."""
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest
from rasterio.transform import from_origin

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_butterfly_73_worldclim_point_niche_20261027 import (
    weighted_quantile, native_reference, extract_climate, outside_rectangle,
    load_climate_from_zip,
)


def fake_climate():
    shape=(8,8)
    return {
        "bio1":np.arange(64,dtype=np.float32).reshape(shape)/10.,
        "bio12":(np.arange(64,dtype=np.float32).reshape(shape)*10.+100),
        "valid":np.ones(shape,dtype=bool),
        "transform":from_origin(-180,90,1/6,1/6),
    }


def test_same_product_native_host_and_point_climate_percentiles():
    x=fake_climate()
    # Ten original native resource polygons, each with several pixels.
    grid=np.tile(np.arange(1,9,dtype=np.uint16),(8,1))
    # 8 distinct native source codes is below predeclared 10-region gate.
    result=native_reference("Sp a",{"Sp a":{"native":sum(1<<i for i in range(8))}},
                            grid,x)
    assert result["status"]=="NATIVE_HOST_REGIONAL_SUPPORT_LT10"
    # Extend to 10 documented source IDs (two unused by the small fake raster).
    native=sum(1<<i for i in range(10))
    result=native_reference("Sp a",{"Sp a":{"native":native}},grid,x)
    assert result["status"]=="VALID_NATIVE_HOST_RESOURCE_WORLDCLIM_REFERENCE"
    assert result["source_native_worldclim_valid_pixels"]==64
    assert 0 <= result["bounds"]["area_weighted_5_95"]["bio1"][0] < result["bounds"]["area_weighted_5_95"]["bio1"][1] <= 6.3
    p=extract_climate(x,89.91,-179.91)
    assert p==(0.0,100.0)
    assert outside_rectangle(p, result["bounds"]["area_weighted_5_95"])
    inside=(result["bounds"]["area_weighted_5_95"]["bio1"][0]+.1,
            result["bounds"]["area_weighted_5_95"]["bio12"][0]+1)
    assert not outside_rectangle(inside,result["bounds"]["area_weighted_5_95"])


def test_weighted_quantile_must_have_valid_and_positive_area_weights():
    assert weighted_quantile([0,1,2],[1,1,1],[.5])==[1.]
    assert weighted_quantile([0,1,2],[1,3,1],[.5])==[1.]
    with pytest.raises(ValueError):
        weighted_quantile([0,1],[0,1],[.5])
    with pytest.raises(ValueError):
        weighted_quantile([0],[1,2],[.5])


def test_missing_worldclim_pixels_cannot_be_interpreted_as_niche_novel():
    x=fake_climate()
    x["valid"][0,0]=False
    assert extract_climate(x,89.91,-179.91) is None
    assert extract_climate(x,98,-179) is None


def test_invalid_bundle_never_falls_back_to_other_climate_product(tmp_path):
    path=tmp_path/"bad.zip"
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("some-weather.txt","WRONG DATA")
    with pytest.raises(ValueError,match="Wrong original WorldClim"):
        load_climate_from_zip(path,tmp_path)
