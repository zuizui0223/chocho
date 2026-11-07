import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from reconcile_island51_zenodo2020_WPI_port_voyages import (
    original,name_key,km,map_node,source_endpoint,PRIMARY_MAX_DISTANCE_KM,
    GEO_DIAGNOSTIC_CAP_KM
)

def test_original_WPI_194_port_51_island_source_and_marine_scope():
    if not (Path(__file__).resolve().parents[1]/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv").is_file():
        pytest.skip("Original WPI source omitted intentionally from anonymous review snapshot")
    ports,by_name=original()
    assert len(ports)==194
    assert len({r["island_id"] for r in ports.values()})==51
    assert sum(r["marine_eligible"] for r in ports.values())==182
    assert PRIMARY_MAX_DISTANCE_KM==25
    assert GEO_DIAGNOSTIC_CAP_KM==10
    assert all(p["name_key"] in by_name for p in ports.values())

def test_port_names_stable_and_do_not_guess_from_unknown():
    assert name_key("Port Everglades")=="EVERGLADES"
    assert name_key("PORT OF EVERGLADES")=="EVERGLADES"
    assert name_key("Everglades Harbor")=="EVERGLADES"
    assert name_key("Bintulu Port")=="BINTULU"
    assert name_key("nan")==""
    assert name_key("None")==""
    assert name_key("Port")==""
    assert not name_key("???")

def test_actual_port_coordinate_rule_requires_name_and_unique_close_WPI_match():
    dummy={
        "10":{"port_id":"10","island_id":"i","port_name":"Port Everglades","country":"US",
              "name_key":"EVERGLADES","latitude":26.07,"longitude":-80.11,"marine_eligible":True},
        "20":{"port_id":"20","island_id":"j","port_name":"Miri","country":"MY",
              "name_key":"MIRI","latitude":4.38,"longitude":114.0,"marine_eligible":True},
    }
    names={"EVERGLADES":[dummy["10"]],"MIRI":[dummy["20"]]}
    good=source_endpoint("16041","Port Everglades","-80.115","26.070")
    mapped=map_node(good,dummy,names)
    assert mapped["status"]=="PRIMARY_NAME_AND_COORD_VERIFIED"
    assert mapped["verified_wpi_id"]=="10"
    wrong=source_endpoint("16041","nan","-80.115","26.070")
    assert map_node(wrong,dummy,names)["verified_wpi_id"] is None
    other=source_endpoint("16041","Port Everglades","114.0","4.38")
    assert map_node(other,dummy,names)["verified_wpi_id"] is None
    assert 0<km(26.070,-80.115,26.07,-80.11)<25

def test_same_name_two_WPI_ports_within_cap_is_ambiguous_not_invented():
    d1={"port_id":"1","name_key":"MIRI","latitude":4.4,"longitude":114.0}
    d2={"port_id":"2","name_key":"MIRI","latitude":4.41,"longitude":114.01}
    item=source_endpoint("5","Miri","114.005","4.405")
    out=map_node(item,{"1":d1,"2":d2},{"MIRI":[d1,d2]})
    assert out["status"]=="AMBIGUOUS_MULTIPLE_MATCHED_WPI_NAMES_NEARBY"
    assert out["verified_wpi_id"] is None
