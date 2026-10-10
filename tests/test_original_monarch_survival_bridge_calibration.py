from __future__ import annotations
import importlib.util
import csv
import io
from pathlib import Path

path=Path(__file__).resolve().parents[1]/"scripts/audit_original_monarch_survival_bridge_calibration.py"
spec=importlib.util.spec_from_file_location("auditmonarch",path)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def synthetic():
    buf=io.StringIO()
    wr=csv.DictWriter(buf,fieldnames=["ID","Plant_ID","Temp","Milkweed","OE_treatment","Surv_pupa","Surv_adult","Notes"])
    wr.writeheader()
    for temp in mod.GROUPS:
      for plant in mod.PLANTS:
       for infection in mod.TREATMENTS:
        for i in range(15):
         pid=f"{temp}-{plant}-{infection}-{i}"
         for j in (0,1):
          wr.writerow({"ID":f"{pid}-{j}","Plant_ID":pid,"Temp":temp,"Milkweed":plant,
                       "OE_treatment":infection,"Surv_pupa":"1","Surv_adult":"1","Notes":""})
    return buf.getvalue().encode("utf-8")

def test_audit_exactly_eight_cells_two_larvae_per_plant():
    raw=synthetic()
    original=mod.git_sha1
    try:
        mod.git_sha1=lambda _:mod.ORIGINAL_GIT_BLOB_SHA1
        r=mod.audit(raw)
    finally:
        mod.git_sha1=original
    assert r["n_assigned"]==240
    assert r["unique_plants"]==120
    assert r["adult_eclosions"]==240
    assert len(r["groups"])==8
    assert r["cross_species_pathogen_bridge_tested"] is False
    assert r["GEB_submission_PR38_untouched"]

def test_source_blob_is_verified():
    try:
        mod.audit(synthetic())
    except ValueError as e:
        assert "blob mismatch" in str(e)
    else:
        raise AssertionError("unaltered-scientific-source constraint bypassed")

def test_accidental_unknown_fate_not_counted_as_death():
    raw=synthetic().decode("utf-8")
    raw=raw.replace(",1,1,\r\n",",NA,NA,accidental injury\r\n",1)
    original=mod.git_sha1
    try:
        mod.git_sha1=lambda _:mod.ORIGINAL_GIT_BLOB_SHA1
        r=mod.audit(raw.encode("utf-8"))
    finally:
        mod.git_sha1=original
    assert r["n_known_survival"]==239
    assert r["adult_eclosions"]==239
    assert len(r["uncertain_fates"])==1
