"""Source-only archive verification, no original plant observations inspected."""
import io,sys,zipfile
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_atlas_and_clarke_source_20261029 import catalogue_members, SOURCES

def test_untrusted_zenodo_zip_rejects_html():
    with pytest.raises(ValueError):catalogue_members(b"<html>Zenodo requires login</html>")

def test_zip_csv_headers_only_no_species_or_site_values():
    buffer=io.BytesIO()
    with zipfile.ZipFile(buffer,"w") as z:
        z.writestr("PA2020/dateclass.csv","hectad,plant,date_period\nSP10,Urtica dioica,1970-1986\nSP11,Salix caprea,2010-2019\n")
    out=catalogue_members(buffer.getvalue())
    assert out[0]["header_fields"]==["hectad","plant","date_period"]
    assert "Urtica dioica" not in str(out)
    assert "SP10" not in str(out)
    assert "source_rows" not in str(out)

def test_zip_traversal_fails_closed():
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:z.writestr("../original.csv","a,b\n1,2\n")
    with pytest.raises(ValueError,match="unsafe path"):
        catalogue_members(b.getvalue())

def test_two_source_published_md5_fixed():
    assert SOURCES[0][2]=="PA2020_10by10kmSquareData.zip"
    assert SOURCES[0][3]=="66a8bb4305df484839587b5fe8acf1d2"
    assert SOURCES[1][3]=="a5f60929e6118b77094110d6ca77bb48"
