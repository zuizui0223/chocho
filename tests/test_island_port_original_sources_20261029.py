"""No-network source identity, port point and geoshape fail closed checks."""
import io,sys,zipfile
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_island_port_original_sources_20261029 import validate_ne_land,check_wpi_metadata

def test_ne_zip_four_geoshape_parts():
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:
        for suffix in (".shp",".shx",".dbf",".prj"):
            z.writestr("ne_10m_land"+suffix,"stub")
    out=validate_ne_land(b.getvalue())
    assert out["island_connected_components_not_yet_validated"]
    assert not any("polygon" in k for k in out if k!="source")

def test_ne_missing_shapefile_fails():
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:z.writestr("ne_10m_land.shp","stub")
    with pytest.raises(ValueError):
        validate_ne_land(b.getvalue())

def test_wpi_requires_actual_wgs84_coordinates():
    metadata={"geometryType":"esriGeometryPoint","name":"WPI","fields":[{"name":x} for x in "ABCDE"]}
    ids={"objectIdFieldName":"OBJECTID","objectIds":list(range(1,3201))}
    point={"features":[{"geometry":{"x":120.2,"y":12.3}} for _ in range(10)]}
    src=check_wpi_metadata(metadata,ids,point)
    assert src["source_feature_count"]==3200
    assert src["real_commodity_flow_or_ship_calls_not_recorded"]
    point["features"][0]["geometry"]={"x":400,"y":12}
    with pytest.raises(ValueError,match="geometry"):
        check_wpi_metadata(metadata,ids,point)

def test_no_future_ports_means_no_establishment_inference():
    with pytest.raises(ValueError):
        check_wpi_metadata({"geometryType":"esriGeometryPoint","fields":[]}, {}, {})
