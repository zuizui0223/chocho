import pytest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_portwatch_maritime_connections_metadata import (
    ORIGINAL,ORIGINAL_GIT_BLOB,SERVICE_BASE,PORTS_LAYER,SAMPLE_MAX,
    git_blob_sha,fixed_wpi,service_candidates,compact_fields,safe_sample
)

def test_original_51_islands_194_ports_182_marine_nodes_source_locked():
    if not ORIGINAL.is_file():
        pytest.skip("Original full WPI geography omitted from anonymous review snapshot")
    ports=fixed_wpi()
    assert git_blob_sha(ORIGINAL)==ORIGINAL_GIT_BLOB
    assert len(ports)==194
    assert len({x["island_id"] for x in ports})==51
    assert sum(x["included"].lower()=="true" for x in ports)==182

def test_official_portwatch_registry_and_nonbulk_source_bound():
    assert SERVICE_BASE.startswith("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/")
    assert PORTS_LAYER.endswith("PortWatch_ports_database/FeatureServer/0")
    assert SAMPLE_MAX<=5

def test_arcgis_public_service_directory_detects_connection_models():
    directory={"services":[
      {"name":"Maritime_Connections","type":"FeatureServer"},
      {"name":"PortWatch_ports_database","type":"FeatureServer"},
      {"name":"Unrelated","type":"FeatureServer"}]}
    matches=service_candidates(directory,{})
    assert len(matches)==1
    assert matches[0]["url"].endswith("Maritime_Connections/FeatureServer")

def test_field_schema_and_sample_exclude_vessel_identifier_fields():
    source={"fields":[
      {"name":"portid_from","type":"string"},
      {"name":"portid_to","type":"string"},
    ]}
    assert [x["name"] for x in compact_fields(source)]==["portid_from","portid_to"]
    features={"features":[{"attributes":{
        "mmsi":"123456789","origin_portid":"PORT0","dest_portid":"PORT1",
        "capacity":111}},{"attributes":{"vessel_name":"name","origin_portid":"PORT2"}}]}
    result=safe_sample(features)
    assert len(result)==2
    assert result[0]["origin_portid"]=="PORT0"
    assert "mmsi" not in result[0]
    assert "vessel_name" not in result[1]
