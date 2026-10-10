"""Synthetic tests of source chronology, no outcome models."""
from __future__ import annotations
import importlib.util
import csv
import io
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/audit_nettle_lagged_enemy_support.py"
spec=importlib.util.spec_from_file_location("nettle_support",P)
assert spec and spec.loader
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)

def synthetic():
    rows=[]
    for w,new in [(20,True),(21,False),(22,True)]:
        for sp,n,k in [("au",10,1),("aio",8,2),
                        ("alev",12 if new else 0,0)]:
            if sp=="alev" and not new:
                continue
            rows.append(dict(BMS_id="siteA",year="2017",week_ISO=str(w),butterfly_species=sp,
                             nb_larvae_in_lab=str(n),larvae_Sturmia_bella=str(k),
                             pupae_Sturmia_bella="0",county="southern",presence_alev=str(int(new)),
                             lat4326="59.7",instar_rond="2",total_larvae="30"))
    return rows

def test_previous_not_current_butterfly_is_used():
    q=s.source_report(synthetic())
    assert q["previous_by_current_support"]["previous_1_current_0"]["n_batches"]==2
    assert q["previous_by_current_support"]["previous_0_current_1"]["n_batches"]==2
    assert q["effect_estimated"] is False
    assert not q["data_sufficiency_gate_pass"]

def test_missing_outcome_is_not_assumed_zero():
    rows=synthetic()
    rows[0]["larvae_Sturmia_bella"]=""
    q=s.source_report(rows)
    assert q["count_integrity"]["resident_invalid_outcome"]==1
    assert q["checks"]["no_invalid_resident_outcomes"] is False

def test_source_year_separation_and_nonfuture_exposure():
    rows=synthetic()
    rows[2]["year"]="2018" # newcomer previous week removed in 2017
    q=s.source_report(rows)
    assert q["previous_by_current_support"].get("previous_1_current_0",{}).get("n_batches",0)==0
