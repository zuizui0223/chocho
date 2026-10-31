"""Verify complete near-island botanical acquisition gate; never make host-absence inference."""
from pathlib import Path
import importlib.util, json
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_host_controls_20261031"
spec=importlib.util.spec_from_file_location("historic_host_gate",D/"build_gate.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
def test_source_pinned_14_pair_gate():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory host/port source excluded from anonymous GEB manuscript")
    rows,s=m.audit(ROOT)
    assert m.render_csv(rows)==(D/"pair_priority.csv").read_text(encoding="utf8")
    assert s==json.loads((D/"pair_priority_summary.json").read_text(encoding="utf8"))
    assert (s["near_island_pairs"],s["unique_islands"],s["shared_new_butterfly_pair_records"])==(14,17,41)
    assert s["gate_counts"]["PRIORITY_UNCLASSIFIED"]=={"pairs":5,"new":9,"unknown":9}
    assert s["gate_counts"]["FOCAL_HOST_OBSERVED"]["new"]==3
    assert s["gate_counts"]["LOW_BOTANICAL_RECORDS"]["new"]==29
    assert s["no_confirmed_host_absence_controls"] is True
