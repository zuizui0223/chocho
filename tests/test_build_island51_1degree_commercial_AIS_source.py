"""Original World Bank AIS ZIP must be source verified and missing not zero."""
import pytest
rasterio=pytest.importorskip("rasterio")
import numpy as np
import hashlib,importlib.util,zipfile
from pathlib import Path
from affine import Affine
f=Path(__file__).resolve().parents[1]/"scripts"/"build_island51_1degree_commercial_AIS_source.py"
sp=importlib.util.spec_from_file_location("aisgrid",f)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_source_hash_is_binary_and_not_filename(tmp_path):
    p=tmp_path/"sourcedata.zip";p.write_bytes(b"geography-not-traffic")
    assert m.hash_file(p)==hashlib.sha256(b"geography-not-traffic").hexdigest()
    assert m.hash_file(p)!=m.ORIGINAL_SHA

def write_overview(path, values, *, nodata=-9999):
    """An external TIFF has pixel values but no geographic CRS of its own."""
    with rasterio.open(path,"w",driver="GTiff",width=values.shape[1],
        height=values.shape[0],count=1,dtype="float32",nodata=nodata,
        transform=Affine(1,0,37,0,-1,20)) as dst:
        dst.write(values.astype("float32"),1)


def test_offset_source_uses_fixed_geographic_cells_not_extent_resize(tmp_path):
    # A plain 2-column resize gives [3.2,12.8]; geographic overlap gives [4,12].
    p=tmp_path/"offset.tif.ovr"
    write_overview(p,np.array([[0,4,8,12,16],[0,4,8,12,16]]))
    values,valid=m.average_overview_on_grid(p,Affine(.5,0,-.25,0,-.5,1),
        target_shape=(1,2),target_transform=Affine(1,0,0,0,-1,1))
    np.testing.assert_allclose(values,[[4.,12.]],rtol=0,atol=1e-6)
    assert valid.tolist()==[[True,True]]


def test_averaging_retains_observed_zero_and_excludes_source_nodata(tmp_path):
    # Missing-only cell stays missing; the partly sampled cell has mean 4, not 2.
    p=tmp_path/"missing_and_zero.tif.ovr"
    write_overview(p,np.array([
        [0,0,-9999,-9999,4,-9999],
        [0,0,-9999,-9999,4,-9999],
        [2,6,8,12,0,0],
        [2,6,8,12,0,0],
    ]))
    values,valid=m.average_overview_on_grid(p,Affine(.5,0,0,0,-.5,2),
        target_shape=(2,3),target_transform=Affine(1,0,0,0,-1,2))
    np.testing.assert_allclose(values,[[0,np.nan,4],[4,10,0]],equal_nan=True)
    assert valid.tolist()==[[True,False,True],[True,True,True]]


def test_overview_without_nodata_is_rejected_before_sampling(tmp_path):
    p=tmp_path/"no_missing_value_contract.tif.ovr"
    write_overview(p,np.array([[0,4],[8,12]]),nodata=None)
    with pytest.raises(ValueError,match="NoData"):
        m.average_overview_on_grid(p,Affine(.5,0,0,0,-.5,1),
            target_shape=(1,1),target_transform=Affine(1,0,0,0,-1,1))


def test_unverified_archive_never_creates_derived_outputs(tmp_path):
    p=tmp_path/"not_official.zip";p.write_bytes(b"unverified source")
    grid=tmp_path/"output.npz";receipt=tmp_path/"receipt.json"
    with pytest.raises(ValueError,match="source SHA"):
        m.derive_raster(p,grid,receipt)
    assert not grid.exists()
    assert not receipt.exists()


def test_extracted_overview_uses_worldfile_and_never_opens_native_tiff(tmp_path):
    # The native TIFF is deliberately unreadable: opening it instead of .ovr fails.
    original=tmp_path/"original_overview.tif"
    with rasterio.open(original,"w",driver="GTiff",width=2,height=2,count=1,
        dtype="int32",nodata=2147483647,transform=Affine(1,0,37,0,-1,20)) as dst:
        dst.write(np.array([[0,4],[8,12]],dtype="int32"),1)
    archive=tmp_path/"overview_fixture.zip"
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("ShipDensity_Commercial1.tif",b"native TIFF must not be opened")
        z.writestr("ShipDensity_Commercial1.tfw",
            "0.005\n0\n0\n-0.005\n-180.012811275\n85.000147937\n")
        z.write(original,"ShipDensity_Commercial1.tif.ovr")
    extracted=tmp_path/"extracted.tif.ovr"
    transform,metadata=m.extract_external_overview(archive,extracted)
    assert extracted.read_bytes()==original.read_bytes()
    np.testing.assert_allclose(tuple(transform)[:6],
        [.01,0,-180.015311275,0,-.01,85.002647937],rtol=0,atol=1e-12)
    assert metadata["overview_width"]==2
    assert metadata["overview_height"]==2
    assert metadata["source_base_raster_read"] is False


def test_unsupported_worldfile_fails_before_overview_extraction(tmp_path):
    archive=tmp_path/"wrong_resolution.zip"
    with zipfile.ZipFile(archive,"w") as z:
        z.writestr("ShipDensity_Commercial1.tif",b"not read")
        z.writestr("ShipDensity_Commercial1.tif.ovr",b"not read either")
        z.writestr("ShipDensity_Commercial1.tfw","0.1\n0\n0\n-0.1\n-180\n85\n")
    extracted=tmp_path/"must_not_exist.tif.ovr"
    with pytest.raises(ValueError,match="worldfile"):
        m.extract_external_overview(archive,extracted)
    assert not extracted.exists()


def test_halfdegree_grid_retains_source_contrast_lost_at_one_degree(tmp_path):
    # Two adjacent half-degree ocean areas have distinct source values 0 and 40.
    # Upsampling their existing 1-degree mean (20) cannot recover that contrast.
    p=tmp_path/"halfdegree_source_contrast.tif.ovr"
    write_overview(p,np.array([[0,0,40,40],[0,0,40,40]]))
    transform=Affine(.25,0,0,0,-.25,1)
    coarse,coarse_valid=m.average_overview_on_grid(p,transform)
    fine,fine_valid=m.average_overview_on_grid(p,transform,resolution_degrees=.5)
    assert coarse.shape==(170,360)
    assert coarse[84,180]==20
    assert coarse_valid.sum()==1
    assert fine.shape==(340,720)
    assert fine[168,360]==0
    assert fine[168,361]==40
    assert fine_valid.sum()==2


def test_coverage_gate_uses_fraction_when_grid_resolution_changes():
    # 100,000 cells exceed the old absolute40,000 gate but cover only41% here.
    valid=np.zeros((340,720),dtype=bool)
    valid.ravel()[:100000]=True
    with pytest.raises(ValueError,match="fraction"):
        m.validate_source_coverage(valid)
    valid.ravel()[:160000]=True
    count,fraction=m.validate_source_coverage(valid)
    assert count==160000
    assert .65<fraction<.66


def test_nested_nearest_pyramids_do_not_replace_first_overview_average(tmp_path):
    # Nested nearest-neighbor pyramids contain only0; original first-level mean=1.
    from rasterio.enums import Resampling
    p=tmp_path/"misleading_nested_pyramids.tif.ovr"
    write_overview(p,np.tile(np.array([[0,0],[0,4]]),(4,4)))
    with rasterio.open(p,"r+") as dst:
        dst.build_overviews([2,4],Resampling.nearest)
    value,valid=m.average_overview_on_grid(p,Affine(.125,0,0,0,-.125,1),
        target_shape=(1,1),target_transform=Affine(1,0,0,0,-1,1))
    assert value[0,0]==1
    assert valid[0,0]
