"""Procedural synthetic fixtures only; never present as butterfly research evidence."""
from __future__ import annotations
import csv,io,importlib.util
from pathlib import Path

PATH=Path(__file__).resolve().parents[1]/"scripts/audit_kyoto_aai_pretreatment_joint_loss.py"
spec=importlib.util.spec_from_file_location("leafbaseline",PATH)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def source():
    lf=io.StringIO()
    pl=io.StringIO()
    la=csv.DictWriter(lf,fieldnames=["h.id","species","l.id","treatment","strain",
                                     "loss","initial.mass","l.or.r"])
    aa=csv.DictWriter(pl,fieldnames=["h.id","l.id","treatment","l.or.r","initial.leaf.area"])
    la.writeheader()
    aa.writeheader()
    for sp in ("a","s"):
        for i in range(30):
            lid=f"{sp}.{i}"
            if i==3:lid+=",composite"
            for t in ("a","c"):
                hid=f"{sp}{i}" if t=="a" else f"{sp}-{i}-control"
                if sp=="a" and i==3 and t=="c":hid="a4"
                if hid=="a4" and not(sp=="a" and i==3 and t=="c"):
                    hid=f"fixture-{hid}"
                half="r" if t=="a" else "l"
                if hid=="a4":half="r,l"
                loss=int(sp=="a" and i<7)
                la.writerow({"h.id":hid,"species":sp,"l.id":lid,"treatment":t,
                             "strain":f"batch_{i//10}","loss":loss,
                             "initial.mass":".02","l.or.r":half})
                aa.writerow({"h.id":hid,"l.id":lid,"treatment":t,
                             "l.or.r":"r,r" if hid=="a4" else half,
                             "initial.leaf.area":str(6+i/10)})
    return lf.getvalue().encode(),pl.getvalue().encode()

def test_original_shape_leaf_unit_and_quality_mismatch():
    rows,ids,half=m.parse(*source())
    assert len(rows)==30
    assert sum(r["both_lost"] for r in rows)==7
    assert ids["composite"]==1
    assert len(half)==1
    assert half[0]["h.id"]=="a4"

def test_constrained_randomization_reproducible_and_holm_valid():
    rows,_,_=m.parse(*source())
    one=m.randomization_test(rows,20261010,draws=999)
    two=m.randomization_test(rows,20261010,draws=999)
    assert one==two
    assert one["source_strain_counts"]["batch_0"]=={"n_pairs":10,"both_lost":7}
    assert {r["pretreatment_feature"] for r in one["statistics"]}==set(m.STATS)
    for test in one["statistics"]:
        assert 0<test["conditional_strain_fixed_permutation_two_sided_p"]<=1
        assert test["holms_multiplicity_adjusted_p"]>=test["conditional_strain_fixed_permutation_two_sided_p"]
    assert one["plant_chemical_or_handling_cause_identified"] is False

def test_source_half_labels_are_not_silently_reconciled():
    a,b=source()
    b=b.replace(b"s0,s.0,a,r,6.0",b"s0,s.0,a,WRONG,6.0",1) # a new, second mismatch
    try:m.parse(a,b)
    except ValueError as e:
        assert "discrepancy" in str(e)
    else:raise AssertionError("original half-label deviation accepted")

def test_original_pretreatment_positive_only():
    a,b=source()
    b=b.replace(b",6.0\r\n",b",0.0\r\n",1)
    try:m.parse(a,b)
    except ValueError as e:
        assert "PRETREATMENT" in str(e)
    else:raise AssertionError("nonsensical baseline leaf area accepted")
