"""Only archive/source validation; no live field experimental data embedded in tests."""
from __future__ import annotations
import importlib.util
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/audit_kyoto_butterfly_regrowth_source.py"
spec=importlib.util.spec_from_file_location("kyoto",p)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_header_dates_and_cage_columns():
    raw=b"plot.no,treat.no,s.density,a.density,date,total.s,total.a\n1,1,2,2,2020-05-01,2,2\n"
    r=m.csv_metadata(raw)
    assert r["records"]==1
    assert "plot.no" in r["columns"]
    assert "date" in r["columns"]

def test_inconsistent_source_fails():
    raw=b"plot.no,treat.no,date\n1,1,x,extra\n"
    try:m.csv_metadata(raw)
    except ValueError as e:
        assert "inconsistent" in str(e)
    else:raise AssertionError("malformed source accepted")

def test_archive_identity_and_hypothesis_separation():
    assert m.ARTICLE_ID==23170898
    assert "larvaldata.csv" in m.EXPECTED[0]
    assert "plantdata.csv" in m.EXPECTED[1]
