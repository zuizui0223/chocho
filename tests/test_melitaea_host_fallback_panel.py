from __future__ import annotations
import importlib.util
from pathlib import Path

PATH=Path(__file__).resolve().parents[1]/"scripts/audit_melitaea_host_fallback_panel.py"
spec=importlib.util.spec_from_file_location("fallback",PATH)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_source_rejects_html():
    class Fake:
        status=200
        def __enter__(self):return self
        def __exit__(self,*x):return False
        def read(self,n):return b"<html>forbidden</html>"
    data,source,errs=m.acquire(opener=lambda *a,**kw:Fake())
    assert data is None and source is None and len(errs)==2

def test_prospective_exclusion_of_missing_years_and_unoccupied_patches():
    rows={}
    for p,pl,vs_list,occ_list in [
        ("a",3,[3,2,1,1],[1,1,1,0]),
        ("b",0,[3,2,1,1],[1,1,1,1]),
        ("c",3,[3,2,1,1],[0,0,0,0]),
    ]:
        for i in range(4):
            rows[(p,2005+i)]={"patch":p,"year":2005+i,"network":"n1",
            "occupancy":occ_list[i],"pl":pl,"vs":vs_list[i]}
    res=m.count_events(rows)
    assert res["risk_set_occupied_t_triples"]==4
    assert res["event_support"]["decline_high"]["patch_years"]==2
    assert res["event_support"]["decline_low"]["patch_years"]==2
    assert res["independent_support_gate_passed"] is False

def test_duplicate_patch_year_fails_loudly():
    raw=("Patch,Year,Network,Area,Occupancy,Nest_count,Pl,Vs\n"
         "p1,2006,n,1,1,2,3,1\n"
         "p1,2006,n,1,0,0,2,1\n").encode()
    try:
        m.decode(raw)
    except ValueError as err:
        assert "duplicate" in str(err)
    else:
        raise AssertionError("duplicate key not rejected")
