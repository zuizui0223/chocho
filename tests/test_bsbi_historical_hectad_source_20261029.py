"""Offline source integrity, privacy and date-class header checks."""
import io,sys,zipfile
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_historical_hectad_source_20261029 import inspect_archive, classify_header

def test_public_synthetic_hectad_source_only_schema():
    f=io.BytesIO()
    with zipfile.ZipFile(f,"w") as z:
        z.writestr("atlas/atlas.csv","scientific_name,hectad,date_class\nA a,SP12,DC2\nB b,SP13,DC3\n")
    result=inspect_archive(f.getvalue())
    assert result["source_rows_scanned_total"]==2
    m=result["original_source_members_schema"][0]
    assert m["plausibly_geographic"] and m["plausibly_taxonomic"] and m["plausibly_dated"]
    assert "SP12" not in str(result)
    assert "A a" not in str(result)

def test_reject_nonzip_and_unsafe_paths():
    with pytest.raises(ValueError,match="genuine ZIP"):
        inspect_archive(b"<html>Login</html>")
    f=io.BytesIO()
    with zipfile.ZipFile(f,"w") as z:z.writestr("../bad.csv","species,hectad\nx,y\n")
    with pytest.raises(ValueError,match="Unsafe"):
        inspect_archive(f.getvalue())

def test_no_hectad_records_from_columns():
    x=classify_header(["taxon","hectad","dateclass"])
    assert x["plausibly_geographic"] and x["plausibly_taxonomic"] and x["plausibly_dated"]
    assert "2000–19" not in str(x)
