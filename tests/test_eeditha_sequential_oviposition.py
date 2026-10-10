from __future__ import annotations

import importlib.util
import io
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("eedi", ROOT/"scripts/analyze_eeditha_sequential_oviposition.py")
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_source_schema_clutches_never_used_as_resources():
    csvtext=("matriline,date,pot,choice.number,CAHI.cl,CALE.cl,PLLA.cl,other.cl,CAHI.eggs,CALE.eggs,PLLA.eggs,other.eggs,notes\n"
             "A,4/26/2017,1,1,0,0,3,0,0,0,104,0,\n"
             "A,4/27/2017,2,2,1,1,0,0,12,8,0,0,\n")
    rows=mod.parse_rows(csvtext.encode())
    trials,_=mod.sequence_dataset(rows)
    assert len(trials)==1
    assert trials[0]["base"] == [1,2]
    assert trials[0]["lag"]==1
    assert trials[0]["y"]==0
    assert rows[0]["PLLA_cl"]==3


def test_all_host_plants_assumed_offered_but_clutches_not_offer_counts():
    lines=["matriline,date,pot,choice.number,CAHI.cl,CALE.cl,PLLA.cl,other.cl,CAHI.eggs,CALE.eggs,PLLA.eggs,other.eggs,notes"]
    for k in range(11):
        group=f"F{k:02d}"
        for j in range(1,5):
            p=35 if (k+j)%2 else 5
            n=40-p
            lines.append(f"{group},4/{24+j}/2017,1,{j},0,0,0,0,{n},0,{p},0,")
    source="\n".join(lines)+"\n"
    rows=mod.parse_rows(source.encode())
    trial,_=mod.sequence_dataset(rows)
    assert len(trial)==33
    assert len(set(r["group"] for r in trial))==11
    r=mod.run(rows)
    assert r["same_trial_clutch_count_used_as_predictor"] is False
    assert r["n_original_trials"]==44
    assert r["decision"] in ("NO_ROBUST_PREDICTIVE_CARRYOVER_GATE_PASS",
                             "PREDICTIVE_CARRYOVER_PILOT_PASSES_NONCAUSAL_GATE")
    assert r["causal_evolutionary_hysteresis_estimated"] is False
    assert len(r["matriline_bootstrap_95ci"])==2


def test_outcome_gate_stops_when_only_one_matriline():
    csvtext=("matriline,date,pot,choice.number,CAHI.cl,CALE.cl,PLLA.cl,other.cl,CAHI.eggs,CALE.eggs,PLLA.eggs,other.eggs,notes\n"
             "A,4/26/2017,1,1,1,0,0,0,20,0,0,0,\n"
             "A,4/27/2017,2,2,0,0,1,0,0,0,30,0,\n")
    out=mod.run(mod.parse_rows(csvtext.encode()))
    assert out["decision"]=="NOT_ESTIMABLE_FROZEN_FEASIBILITY_GATE"
