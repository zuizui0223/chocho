"""No network, never confuse OLE file magic or ZIP data with survey nondetections."""
import io,sys,zipfile
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_uk_public_transect_source_20261028 import describe

def test_ole_excel_only_format_not_zeros():
    d=describe(bytes.fromhex("D0CF11E0A1B11AE1")+b"hello","xls")
    assert d["format_detected"]=="LEGACY_OLE_XLS"
    assert "absent" not in d

def test_source_zip_csv_only_headers_rows():
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:
        z.writestr("visits.csv","site,date,species,count\na,2000-04-01,A,0\nb,2000-04-08,B,3\n")
    out=describe(b.getvalue(),"zipcsv")
    assert out["format_detected"]=="ZIP"
    assert out["total_csv_source_rows"]==2
    assert out["csv_schema_metadata"][0]["columns"]==["site","date","species","count"]
    assert "2000-04-01" not in str(out)
    assert "a,2000" not in str(out)

def test_html_download_rejected_by_stage():
    assert describe(b"<html>login required</html>","zipcsv")["format_detected"]=="HTML_OR_JSON_NOT_TABLE"

def test_bad_zip_fails_closed():
    with pytest.raises(zipfile.BadZipFile):
        describe(b"PK\x03\x04broken","zipcsv")
