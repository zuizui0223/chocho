"""Synthetic split leaves test algorithm ONLY; no synthetic data support scientific effect."""
import csv,io,importlib.util
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/analyze_kyoto_original_aai_paired_effects.py"
spec=importlib.util.spec_from_file_location("aai_effect",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def sample():
    left=io.StringIO();right=io.StringIO()
    a=csv.DictWriter(left,fieldnames=["h.id","strain","species","initial.mass","day1.mass","l.id","l.or.r","treatment","loss"])
    b=csv.DictWriter(right,fieldnames=["l.id","l.or.r","treatment","initial.leaf.area","h.id","day1.leaf.area","consumed.leaf.area","note"])
    a.writeheader();b.writeheader()
    for sp in ("a","s"):
      for ix in range(30):
        leaf=f"{sp}-{ix}"
        for trt in ("a","c"):
          h=f"{sp}-{ix}-{trt}"
          lost=int(sp=="a" and ix==0 and trt=="a")
          a.writerow({"h.id":h,"strain":"stock","species":sp,"initial.mass":"10","day1.mass":"12" if trt=="a" else "11","l.id":leaf,"l.or.r":"r" if trt=="a" else "l","treatment":trt,"loss":lost})
          b.writerow({"l.id":leaf,"l.or.r":"r" if trt=="a" else "l","treatment":trt,
                      "initial.leaf.area":"40","h.id":h,"day1.leaf.area":"30","consumed.leaf.area":"10","note":""})
    return left.getvalue().encode(),right.getvalue().encode()

def test_paired_source_identity_and_losses():
    leaves,rows,mismatch,quality=m.parse(*sample())
    assert len(rows)==120 and len(leaves)==60 and mismatch==[]
    pairs=[val for (sp,lid),val in leaves.items() if sp=="a"]
    assert sum(p["a"]["loss"]==0 and p["c"]["loss"]==0 for p in pairs)==29
    assert abs(m.estimate(pairs,"growth_ratio",True)-.1)<1e-9
    assert m.estimate(pairs,"consumption_ratio",True)==0
    lohi=m.bootstrap(pairs,"growth_ratio",True,20261010,draws=150)
    assert lohi[0]<=.1<=lohi[1]

def test_leaf_area_join_identity_failure():
    a,b=sample()
    b=b.replace(b"a-0-a",b"bogus",1)
    try:m.parse(a,b)
    except ValueError as e:
        assert "ID sets" in str(e)
    else:raise AssertionError("unmatched original insect accepted")

def test_original_side_typo_flagged_but_identity_kept():
    a,b=sample()
    b=b.replace(b"a-0,r,a",b"a-0,l,a",1)
    leaves,rows,mismatch,quality=m.parse(a,b)
    assert len(mismatch)==1
    assert mismatch[0]["larval_half_label"]=="r"
    assert len(leaves)==60
