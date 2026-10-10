from __future__ import annotations
import importlib.util,io,csv,hashlib
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/audit_kyoto_aai_paired_leaf_schema.py"
spec=importlib.util.spec_from_file_location("aai",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def mock():
    fa=io.StringIO();fb=io.StringIO()
    x=csv.DictWriter(fa,fieldnames=["h.id","strain","species","initial.mass","day1.mass","l.id","l.or.r","treatment","loss"])
    y=csv.DictWriter(fb,fieldnames=["l.id","l.or.r","treatment","initial.leaf.area","h.id","day1.leaf.area","consumed.leaf.area","note"])
    x.writeheader();y.writeheader()
    for i in range(120):
        leaf=f"L{i//2}"
        row={"h.id":f"h{i}","strain":"x","species":"A","initial.mass":"10","day1.mass":"11",
             "l.id":leaf,"l.or.r":"l" if i%2 else "r","treatment":"AAI" if i%2 else "control","loss":"0"}
        x.writerow(row)
        y.writerow({"l.id":leaf,"l.or.r":row["l.or.r"],"treatment":row["treatment"],
                    "initial.leaf.area":"100","h.id":row["h.id"],
                    "day1.leaf.area":"80","consumed.leaf.area":"20","note":""})
    return fa.getvalue().encode(),fb.getvalue().encode()

def test_frozen_source_lengths_and_matching_rows():
    a,b=mock()
    out=m.summarize(a,b)
    assert out["larvae"]["n"]==120
    assert out["matching_row_keys"]==120
    assert out["source_effect_estimated"] is False
    assert out["GEB_PR38_unchanged"] is True

def test_missing_pair_duplicate_rejected():
    a,b=mock()
    b=b.replace(b"h1,",b"h0,",1)
    try:
        m.summarize(a,b)
    except ValueError as ex:
        assert "compound IDs" in str(ex)
    else:
        raise AssertionError("incorrectly accepted duplicate original row identities")
