"""Fixtures only validate original-data chronology; no simulated ecological outcome."""
import hashlib
import importlib.util
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"scripts/audit_kyoto_stage_demand_structure.py"
spec=importlib.util.spec_from_file_location("stages",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def synthetic_source():
    records=[",".join(m.COLUMNS)]
    for cage in range(1,31):
        a_n=4
        s_n=4
        for day in range(23):
            early_a=4 if day<4 else 0
            late_a=4 if day>=4 and day<12 else 0
            pup_a=4 if day>=12 else 0
            d={k:"0" for k in m.COLUMNS}
            d.update({"plot.no":str(cage),"treat.no":"1","s.density":str(s_n),
                      "a.density":str(a_n),"date":"2012-09-01",
                      "day":str(day),"total.s":"4",
                      "total.a":str(early_a+late_a),
                      "a1":str(early_a),"a3":str(late_a),"cumul.ap":str(pup_a)})
            records.append(",".join(d[k] for k in m.COLUMNS))
    return ("\n".join(records)+"\n").encode("utf-8")

def test_exact_original_schema_and_stage_onset():
    raw=synthetic_source()
    old=m.MD5
    try:
        m.MD5=hashlib.md5(raw).hexdigest()
        r=m.analyze(raw)
    finally:
        m.MD5=old
    assert r["n_source_cages"]==30
    assert r["n_original_cage_date_rows"]==690
    assert r["instar_support"]["cages_with_late_A_observed"]==30
    assert r["per_cage_stage_onsets"][0]["first_late_A_instar_observed_day"]==4
    assert r["new_butterfly_survival_effect_estimated"] is False
    assert r["resource_access_or_quality_observed_at_stage"] is False

def test_original_source_checksum_rejects_replaced_rows():
    try:m.analyze(synthetic_source())
    except ValueError as ex:
        assert "checksum mismatch" in str(ex)
    else:
        raise AssertionError("arbitrary fixture passed original source checksum")
