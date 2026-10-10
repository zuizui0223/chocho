"""Synthetic provenance/source tests; no scientific inference from test fixtures."""
from __future__ import annotations
import importlib.util
import hashlib
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/audit_nettle_parasitoid_figshare_source.py"
spec=importlib.util.spec_from_file_location("prior_audit",SCRIPT)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def source():
    raw=b"butterfly,sampled,parasitised,site,year\nAglais urticae,10,3,A,2017\n"
    meta={"id":14260211,"title":"Rewiring of interactions in a changing environment",
          "doi":"10.17045/STHLMUNI.14260211",
          "files":[{"name":"batch_monitoring.csv","id":10,"size":len(raw),
                    "download_url":"https://example.com/test.csv",
                    "computed_md5":hashlib.md5(raw).hexdigest()}]}
    return meta,raw

def test_source_audit_has_no_biological_result():
    meta,raw=source()
    result=m.gate_inventory(meta,fetch=lambda url,limit:raw)
    assert result["status"]=="RAW_SOURCE_SCHEMA_VERIFIED_ONLY"
    assert result["raw_butterfly_records_verified"] is True
    assert result["source_files"][0]["schema"]["records"]==1
    assert result["new_ecological_effect_estimated"] is False
    assert result["botanical_introduction_tested"] is False
    assert result["GEB_PR38_unchanged"] is True

def test_integrity_mismatch_stops_schema_promotion():
    meta,raw=source()
    meta["files"][0]["computed_md5"]="00000000000000000000000000000000"
    result=m.gate_inventory(meta,fetch=lambda url,limit:raw)
    assert result["status"]=="SOURCE_INVENTORY_ONLY"
    assert result["source_files"][0]["status"]=="SOURCE_CSV_ACCESS_OR_SCHEMA_FAILURE"

def test_bad_data_row_width_rejected():
    try:
        m.parse_csv_preview(b"site,butterfly,parasitised\nA,Aglais\n","batch_monitoring.csv")
    except ValueError as err:
        assert "width mismatch" in str(err)
    else:
        raise AssertionError("malformed original CSV silently accepted")

def test_wrong_article_id_rejected():
    meta,raw=source()
    meta["id"]=42
    try:
        m.gate_inventory(meta,fetch=lambda url,limit:raw)
    except ValueError as err:
        assert "wrong original" in str(err)
    else:
        raise AssertionError("wrong archive accepted")
