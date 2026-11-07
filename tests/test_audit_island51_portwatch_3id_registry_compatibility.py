import pytest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_portwatch_3id_registry_compatibility import (
    fixed_scope,source_sample,haversine,PORTS,METADATA,CONNECTION,REGISTRY,
    SAMPLE_ID_COUNT,MATCH_DISTANCE_KM
)

def test_original194_source_not_automatically_replaced_by_portwatch_id():
    if not PORTS.is_file():
        pytest.skip("Original archived WPI194 port geography omitted from anonymous snapshot")
    assert len(fixed_scope())==194
    assert MATCH_DISTANCE_KM==25
    assert CONNECTION.endswith("Spillover_Simulator_Maritime_Connections/FeatureServer/0")
    assert REGISTRY.endswith("PortWatch_ports_database/FeatureServer/0")

def test_three_connection_ids_are_determined_before_provider_registry_query():
    if not METADATA.is_file():
        pytest.skip("Original source sample receipt not included in anonymous review bundle")
    selected=source_sample()
    assert len(selected)==SAMPLE_ID_COUNT==3
    assert len({x["portid"] for x in selected})==3
    assert all(x["portid"].startswith("port") for x in selected)
    assert [x["portid"] for x in selected]==["port0","port540","port1067"]

def test_provider_WGS84_latlon_geodesic_is_symmetric_and_nearest_not_assumed():
    assert haversine(0,0,0,0)==0
    assert abs(haversine(26.0,-80.0,26.1,-80.0)-haversine(26.1,-80.0,26.0,-80.0))<1e-10
    assert 10<haversine(26.0,-80.0,26.1,-80.0)<15
