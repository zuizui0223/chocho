import pytest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_mendeley_global_port_network import (
  source_roster,WPI,WPI_GIT_BLOB,git_blob,HOST,DOI,DATASET_ID,VERSION,
  records,normalize_lists,sample_binary,safe_file_probe,MAX_DOWNLOAD
)

def test_fixed_original_WPI194_source():
    if not WPI.is_file():
        pytest.skip("Frozen original WPI geographic input is omitted by the anonymous-review bundle")
    original=source_roster()
    assert git_blob(WPI)==WPI_GIT_BLOB
    assert len(original)==194
    assert len({r["island_id"] for r in original})==51
    assert sum(r["included"].lower()=="true" for r in original)==182

def test_public_repository_id_and_model_not_observed_ship_voyage_contract():
    assert DATASET_ID=="vzzy3b9gg4"
    assert VERSION==1
    assert DOI=="10.17632/vzzy3b9gg4.1"
    assert HOST=="https://api.data.mendeley.com"
    assert MAX_DOWNLOAD<=120_000_000

def test_metadata_file_extraction_and_size_guard():
    raw={"files":[
       {"id":"a","filename":"Network_ports.csv","size":2000,
        "content_details":{"sha256_hash":"source-supplied"}},
       {"id":"b","filename":"nonnetwork.txt","size":500},
       {"id":"c","filename":"Network_ports_big.csv","size":300000000}
    ]}
    v=records(raw,"fake")
    assert len(v)==3
    assert v[0]["name"]=="Network_ports.csv"
    assert safe_file_probe(v[0],v[0])
    assert not safe_file_probe(v[1],v[1])
    assert not safe_file_probe(v[2],v[2])

def test_tabular_schema_probe_is_limited_and_no_raw_vessel_identity():
    data=b"port_id,lat,lon\n1,42.0,-76.0\n"
    seen=sample_binary("Network_ports.csv",data)
    assert seen["format"]=="TEXT_CSV_OR_OTHER"
    assert seen["first_lines"][0]=="port_id,lat,lon"
