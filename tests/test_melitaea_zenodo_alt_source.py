"""No silent transfer of data or conclusions between distinct Melitaea archives."""
from __future__ import annotations
import importlib.util
import io
import zipfile
from pathlib import Path

MODULE=Path(__file__).resolve().parents[1]/"scripts/audit_melitaea_zenodo_alt_source.py"
s=importlib.util.spec_from_file_location("alt_source",MODULE)
assert s and s.loader
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

def zipped(filename="study/panel.csv", content=b"patch,year,Pl,Vs,occupancy\np,2010,1,3,1\n"):
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,"w") as z:
        z.writestr(filename,content)
    return buf.getvalue()

def test_inventory_csv_and_provenance_scope():
    x=m.inventory(zipped())
    assert len(x)==1
    assert x[0]["columns"]==["patch","year","Pl","Vs","occupancy"]
    assert x[0]["approx_line_count"]==2
    assert m.RECORD=="4987060"
    assert m.EXPECTED_MD5=="69122a1d82b1fb970fb6638b02da3db4"

def test_reject_path_traversal_member():
    try:
        m.inventory(zipped("../fake.csv"))
    except ValueError as err:
        assert "unsafe" in str(err)
    else:
        raise AssertionError("ZIP traversal undetected")

def test_reject_html_blockpage():
    try:
        m.inventory(b"<html>Forbidden</html>")
    except ValueError as err:
        assert "ZIP" in str(err)
    else:
        raise AssertionError("HTML treated as ZIP")

def test_source_failure_never_claims_result():
    original=m.fetch
    try:
        m.fetch=lambda: (None,None,[{"problem":"HTTP denied"}])
        r=m.audit()
    finally:
        m.fetch=original
    assert r["status"]=="SOURCE_ACCESS_BLOCKED"
    assert r["host_fallback_effect_estimated"] is False
    assert r["separate_from_DiLeo_2024"] is True
