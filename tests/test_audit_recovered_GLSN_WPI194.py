import pytest
pytest.importorskip("pycountry", reason="Dedicated original-source CI installs pycountry")
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_original_GLSN2018_UNCTAD_WPI194_join import (
    SOURCE,WPI,WPI_BLOB,GLSNPORT,GLSNEDGE,GLSNRECEIPT,UNCTAD,AIS,
    verify_source,gitblob,canonical,iso3_to_iso2
)

def test_original_GLSN_2018_source_and_all_original_194_WPI_port_records():
    if not WPI.is_file() or not GLSNEDGE.is_file():
        pytest.skip("Full original research source inputs omitted in anonymous review bundle")
    ports,edges,original,unctad,ais,meta=verify_source()
    assert gitblob(WPI)==WPI_BLOB
    assert len(ports)==887 and len(edges)==28104 and len(original)==194
    assert len(unctad)==194 and len(ais)==194
    assert len({row["island_id"] for row in original})==51
    assert sum(row["included"].lower()=="true" for row in original)==182
    assert meta["sources"][0]["original_publisher_git_blob"]=="5486d8c204a7c49b375a14762fe5cb36c7877d96"
    assert meta["sources"][1]["original_publisher_git_blob"]=="18098cab432ddbcd7e0a7091fc1c4cedbda40fb7"

def test_source_country_and_port_canonical_names_not_raw_country_guessing():
    assert iso3_to_iso2("USA")=="US"
    assert iso3_to_iso2("GBR")=="GB"
    assert iso3_to_iso2("JPN")=="JP"
    assert canonical("PORT OF SOUTHAMPTON")=="SOUTHAMPTON"
    assert canonical("Southampton Harbour")=="SOUTHAMPTON"
    assert canonical("BINTULU PORT")=="BINTULU"
    assert canonical("Port Louis")=="LOUIS"
    assert canonical("nan")==""
