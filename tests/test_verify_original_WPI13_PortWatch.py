import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from verify_original_WPI13_PortWatch_geography import (
    WPI, IDS, frozen_thirteen, normalize, haversine, verify_provider_box,
    MAX_PROVIDER_ROWS_PER_BOX, REGISTRY, PRIOR_VALIDATED
)

def test_frozen_all_13_source_selected_GLSN_UNCTAD_but_not_AIS2020():
    if not WPI.is_file():
        pytest.skip("Original complete WPI source omitted from anonymous review snapshot")
    selected=frozen_thirteen()
    assert len(selected)==13 and set(selected)==set(IDS)
    assert len({record[0]["island_id"] for record in selected.values()})==5
    assert not PRIOR_VALIDATED
    assert all(record[1]["Zenodo2020_actual_WPI_port_identity_verified"]!="True"
               for record in selected.values())
    assert MAX_PROVIDER_ROWS_PER_BOX==15
    assert REGISTRY.endswith("PortWatch_ports_database/FeatureServer/0")

def test_exact_name_ISO3_25km_geodesic_and_fail_closed_registry():
    assert normalize(" KHOLMSK ")=="KHOLMSK"
    assert normalize("Benoa")=="BENOA"
    assert normalize("Port Louis")=="PORT LOUIS"
    assert haversine(0,0,0,0)==0
    assert 20<haversine(0,0,.2,0)<25
    w={"original_latitude":"1.23","original_longitude":"109.4","port_name":"KUCHING","country_code":"MY"}
    port={"attributes":{"portid":"portEXAMPLE","portname":"Kuching",
        "country":"Malaysia","ISO3":"MYS","lat":1.231,"lon":109.401,"LOCODE":"MY KCH"}}
    data,match,status=verify_provider_box("415",w,"MYS",[port],1)
    assert status=="ONE_UNIQUE_SOURCE_PORT_NAME_ISO3_GEOGRAPHY_MATCH"
    assert match["portid"]=="portEXAMPLE"
    assert 0<match["source_geodesic_km"]<1
    data,match,status=verify_provider_box("415",w,"MYS",[port,port],2)
    assert status=="AMBIGUOUS_MULTIPLE_PORTWATCH_NAMES_IN_BOX" and match is None
    data,match,status=verify_provider_box("415",w,"MYS",[port],17)
    assert status=="INCOMPLETE_OR_FAILED_REGISTRY_BOX" and match is None
    wrong={"attributes":dict(port["attributes"],ISO3="RUS")}
    data,match,status=verify_provider_box("415",w,"MYS",[wrong],1)
    assert status=="NO_EXACT_SOURCE_PORT_MATCH_IN_COMPLETE_BOX" and match is None
