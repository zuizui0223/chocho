import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_zenodo2020_voyage_source import (
  original_port_roster,blob,PORTS,PORTS_GIT_BLOB,ZENODO,ZENODO_MD5,
  MAX_ARCHIVE_BYTES
)

def test_original_194_source_port_roster_reproduces_previous_study():
    if not PORTS.is_file():
        pytest.skip("Original WPI source omitted intentionally from anonymous review snapshot")
    r=original_port_roster()
    assert len(r)==194
    assert len({x["island_id"] for x in r})==51
    assert len({int(x["port_id"]) for x in r})==194
    assert sum(x["included"].lower()=="true" for x in r)==182
    assert blob(PORTS)==PORTS_GIT_BLOB

def test_external_2020_source_filename_checksum_and_expected_scope_fixed():
    assert "zenodo.org/records/10783698/files/triplet_2020.zip" in ZENODO
    assert ZENODO_MD5=="af27e1bef1f706e2e29ed53aee47ff10"
    assert MAX_ARCHIVE_BYTES>=32_453_859
