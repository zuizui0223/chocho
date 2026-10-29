"""Source-independent WPI all-original-port coordinate integrity and no false shipping inference."""
import pytest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_world_port_index_complete_20261029 import validate

def test_complete_source_coordinates_not_traffic():
    rows=[{"attributes":{"objectid":i,"latitude":float(i/100),"longitude":float(i/100),
              "country":"SYNTHETIC"}} for i in range(1,31)]
    x=validate(rows,list(range(1,31)))
    assert x["valid_original_WGS84_attribute_coordinates"]==30
    assert x["raw_port_names_and_coordinates_or_IDs_exported"] is False
    assert x["original_source_2025_ports_are_NOT_trade_or_ship_arrival_counts"]

def test_source_identity_missing_or_duplicate_fails():
    rows=[{"attributes":{"objectid":1,"latitude":12,"longitude":13}}]
    with pytest.raises(ValueError,match="incomplete"):
        validate(rows,[1,2])
    with pytest.raises(ValueError,match="duplicated"):
        validate(rows*2,[1])

def test_obvious_invalid_coordinates_fail():
    rows=[{"attributes":{"objectid":i,"latitude":900,"longitude":13}} for i in range(1,31)]
    with pytest.raises(ValueError,match="invalid"):
        validate(rows,list(range(1,31)))
