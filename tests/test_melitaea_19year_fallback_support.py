"""Feasibility/chronology tests; never infer ecological efficacy from synthetic rows."""
from __future__ import annotations
import importlib.util
from pathlib import Path

path=Path(__file__).resolve().parents[1]/"scripts/audit_melitaea_19year_fallback_support.py"
s=importlib.util.spec_from_file_location("mel",path)
assert s and s.loader
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

def sample():
    lines=["year\tpatch\tpopulation\tplantago\tveronica\tprevious_population"]
    for i in range(3):
        patch=f"p{i}"
        for yr in range(2005,2009):
            v=3 if yr==2005 else 2
            lines.append(f"{yr}\t{patch}\t{2 if yr<2008 else 0}\t{3 if i else 0}\t{v}\t2")
    return ("\n".join(lines)+"\n").encode()

def test_read_source_columns_and_three_year_outcomes():
    rows,meta=m.ingest(sample())
    assert meta["unique_patches"]==3
    assert meta["unique_patch_years"]==12
    stat=m.support(rows)
    assert stat["host_scores_ordinal_0to3"] is True
    assert stat["three_year_consecutive_occupied_risk_set"]==6
    assert stat["decline_event_counts"]["decline_high"]["patch_years"]==2
    assert stat["decline_event_counts"]["decline_low"]["patch_years"]==1
    assert stat["data_support_pass"] is False

def test_invalid_host_categories_prevent_effect_analysis():
    rows,_=m.ingest(sample().replace(b"\t3\t2\t2\n",b"\t5\t2\t2\n",1))
    assert m.support(rows)["data_support_pass"] is False

def test_duplicated_patch_year_does_not_inflate_counts():
    raw=sample()+b"2005\tp0\t1\t3\t3\t1\n"
    rows,meta=m.ingest(raw)
    assert len(rows)==12
    assert meta["duplicates"]==1

def test_future_outcome_is_not_used_for_host_decline_definition():
    rows,_=m.ingest(sample())
    a=m.support(rows)["decline_event_counts"]
    for row in rows.values():
        if row["year"]==2008:
            row["population"]=7
    b=m.support(rows)["decline_event_counts"]
    assert a["decline_low"]["patch_years"]==b["decline_low"]["patch_years"]
    assert a["decline_high"]["patch_years"]==b["decline_high"]["patch_years"]
