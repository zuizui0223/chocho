"""Synthetic leaf-identity and exact-pair statistics tests; not ecological results."""
from __future__ import annotations
import csv,io
import importlib.util
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/analyze_kyoto_aai_paired_loss.py"
spec=importlib.util.spec_from_file_location("aai_loss",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def fake():
    buf=io.StringIO()
    keys=["h.id","l.id","species","treatment","loss","strain"]
    w=csv.DictWriter(buf,fieldnames=keys)
    w.writeheader()
    for s in ["a","s"]:
        for i in range(30):
            for t in ["a","c"]:
                w.writerow({"h.id":f"{s}_{i}_{t}","l.id":f"{s}_{i}","species":s,
                             "treatment":t, "strain":f"family_{s}_{i}",
                             "loss":1 if s=="a" and i<6 and t=="a" else 0})
    return buf.getvalue().encode()

def test_exact_discordance_p_not_unpaired():
    assert m.exact_mcnemar(6,1)==0.125
    assert m.exact_mcnemar(0,0)==1
    assert m.exact_mcnemar(4,0)==0.125
    assert m.exact_mcnemar(2,2)==1

def test_paired_loss_source_structure():
    pairs=m.load_pairs(fake())
    assert len(pairs["a"])==30
    assert len(pairs["s"])==30
    out=m.summarize_pair_data(fake())
    native=out["species"]["Atrophaneura alcinous"]
    assert native["pair_status"]=={"neither_lost":24,"only_AAI_lost":6,
                                   "only_control_lost":0,"both_lost":0}
    assert native["paired_loss_risk_difference_AAI_minus_control"]==0.2
    assert out["source_loss_is_haphazard_death_not_proven_AAI_toxicity"] is True
    assert out["GEB_PR38_untouched"] is True

def test_one_missing_half_fails_closed():
    raw=fake().decode().splitlines()
    raw.pop()
    try:m.load_pairs(("\n".join(raw)+"\n").encode())
    except ValueError as exc:assert "cohort length" in str(exc)
    else:raise AssertionError("missing randomized larva accepted")

def test_repeated_leaf_treatment_fails():
    raw=fake().decode().splitlines()
    raw[2]=raw[2].replace(",c,",",a,")
    try:m.load_pairs(("\n".join(raw)+"\n").encode())
    except ValueError as exc:assert "duplicate treatment" in str(exc)
    else:raise AssertionError("duplicate leaf-treatment accepted")

def test_paired_bootstrap_is_reproducible():
    pairs=[(1,0)]*6+[(0,1)]+[(0,0)]*23
    a=m.bootstrap_paired_risk_difference(pairs,20261010,draws=599)
    b=m.bootstrap_paired_risk_difference(pairs,20261010,draws=599)
    assert a==b
    assert a[0] <= (6-1)/30 <= a[1]

def test_fisher_exact_paired_loss_dependence_not_treatment_effect():
    assert abs(m.fisher_two_sided_pair_dependence(30,13,8,7)-0.009357543450496972)<1e-12
    assert abs(m.fisher_two_sided_pair_dependence(30,0,0,0)-1)<1e-12

def test_strain_concordance_is_a_recorded_confounder_not_leaf_mechanism():
    out=m.summarize_pair_data(fake())
    assert out["strain_pair_concordance"]["Atrophaneura alcinous"]["leaf_pairs_sharing_same_recorded_strain"]==30
    assert out["exploratory_after_margin_inspection"] is True
    assert out["species"]["Atrophaneura alcinous"]["biological_cause_of_loss_identified"] is False
