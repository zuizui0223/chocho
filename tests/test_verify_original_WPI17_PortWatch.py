import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from verify_original_WPI17_PortWatch_geography import (
  WPI,IDS,frozen_seventeen,normalize,haversine,verify_provider_box,
  MAX_PROVIDER_ROWS_PER_BOX,REGISTRY,PRIOR_VALIDATED
)

def test_exact_original_seventeen_fixed_from_GL_UNCTAD_and_verified_2020_AIS():
    if not WPI.is_file():
        pytest.skip("Original WPI source omitted by anonymous review snapshot")
    selected=frozen_seventeen()
    assert len(selected)==17 and set(selected)==set(IDS)
    assert len({v[0]["island_id"] for v in selected.values()})==14
    assert len(PRIOR_VALIDATED)==4 and all(key in IDS for key in PRIOR_VALIDATED)
    assert MAX_PROVIDER_ROWS_PER_BOX==15
    assert REGISTRY.endswith("PortWatch_ports_database/FeatureServer/0")

def test_source_only_name_country_and_geodesic_are_independent():
    assert normalize(" Santa Cruz de Tenerife ")=="SANTA CRUZ DE TENERIFE"
    assert normalize("Port-Louis")=="PORT LOUIS"
    assert normalize("PUERTO DEL ROSARIO")=="PUERTO DEL ROSARIO"
    assert haversine(0,0,0,0)==0
    assert 20<haversine(0,0,.2,0)<25

def test_ambiguous_or_truncated_source_registry_never_promoted():
    w={"original_latitude":"28.47","original_longitude":"-16.23",
       "port_name":"SANTA CRUZ DE TENERIFE","country_code":"ES"}
    match={"attributes":{"portid":"portTEST","portname":"Santa Cruz de Tenerife",
        "country":"Spain","ISO3":"ESP","lat":28.471,"lon":-16.231,"LOCODE":"ES SCT"}}
    rows,port,status=verify_provider_box("428",w,"ESP",[match],1)
    assert status=="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH"
    assert port["portid"]=="portTEST"
    assert port["source_geodesic_km"]<1
    rows,port,status=verify_provider_box("428",w,"ESP",[match,match],2)
    assert status=="AMBIGUOUS_MULTIPLE_PORTWATCH_NAMES_IN_BOX" and port is None
    rows,port,status=verify_provider_box("428",w,"ESP",[match],16)
    assert status=="INCOMPLETE_OR_FAILED_REGISTRY_BOX" and port is None
    wrong={"attributes":dict(match["attributes"],ISO3="GBR")}
    rows,port,status=verify_provider_box("428",w,"ESP",[wrong],1)
    assert status=="NO_EXACT_SOURCE_PORT_MATCH_IN_COMPLETE_BOX" and port is None
    far={"attributes":dict(match["attributes"],lat=29.5)}
    rows,port,status=verify_provider_box("428",w,"ESP",[far],1)
    assert status=="NO_EXACT_SOURCE_PORT_MATCH_IN_COMPLETE_BOX" and port is None
