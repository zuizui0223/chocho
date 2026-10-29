"""No-network science/data integrity tests for BSBI Plant Atlas source."""
import io,sys,stat,zipfile
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_bsbi_historical_hectad_source_20261029 as m

def test_nonexistent_source_fails_closed(tmp_path):
    with pytest.raises(ValueError,match="missing"):
        m.audit_zip(tmp_path/"missing.zip")

def test_wrong_source_digest_fails_closed(tmp_path):
    p=tmp_path/"example.zip"
    with zipfile.ZipFile(p,"w") as z:
        z.writestr("PlantAtlas.csv","taxon,hectad,period\nA,A1,1950-1969\n")
    with pytest.raises(ValueError,match="MD5 differs"):
        m.audit_zip(p)

def test_original_source_header_and_no_original_site_export(tmp_path,monkeypatch):
    p=tmp_path/"source.zip"
    with zipfile.ZipFile(p,"w") as z:
        z.writestr("historical/hectads.csv","taxon,hectad,period\nSecretPlant,Site01,1950-1969\n")
    import hashlib
    monkeypatch.setattr(m,"MD5",hashlib.md5(p.read_bytes()).hexdigest())
    r=m.audit_zip(p)
    assert r["source_members"][0]["first_row_field_names"]==["taxon","hectad","period"]
    assert "SecretPlant" not in str(r) and "Site01" not in str(r)
    assert r["raw_plant_hectad_records_exported"] is False
