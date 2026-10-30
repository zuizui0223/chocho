"""Synthetic-only validation for original-year resolver; source raw parquet stays local."""
from __future__ import annotations

import importlib.util
import pathlib
import pytest

ROOT=pathlib.Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/resolve_original_island_chain_years.py"
spec=importlib.util.spec_from_file_location("original_island_chain_years",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def make_screen():
    links=[]
    for i in range(5):
        links.append({
            "butterfly":f"Species {i}",
            "catalogue_compatible_prior_host":f"Plant {i}",
            "port_representative_great_circle_km":100.0*(i+1),
            "endpoints":[
                {"island":str(i*2+1),"first_year_in_2011_2017_followup_table":2014},
                {"island":str(i*2+2),"first_year_in_2011_2017_followup_table":None}
            ]
        })
    return {"schema":"chocho_prior_host_island_record_chain_screen_20261030_v1","links":links}


def test_missing_raw_year_never_imputed():
    screen=make_screen()
    keys={f"Species {i}":str(100+i) for i in range(5)}
    raw=[{"speciesKey":str(100+i),"landmass_id":str(i*2+1),"year":2014} for i in range(5)]
    r=mod.full_year_resolution(screen,keys,raw)
    assert len(r["rows"])==5
    assert r["fully_resolved_pair_count"]==0
    assert r["distinct_order_and_prior_re_record_count"]==0
    for row in r["rows"]:
        assert row["endpoints"][1]["year_state"]=="NO_ROW_IN_PROVIDED_YEAR_TABLE"
        assert row["upstream_re_recorded_before_later_island_first_record"] is None


def test_complete_years_resolve_order_and_re_record_before_downstream():
    screen=make_screen()
    keys={f"Species {i}":str(100+i) for i in range(5)}
    raw=[]
    for i in range(5):
        raw.extend([
            {"speciesKey":str(100+i),"landmass_id":str(i*2+1),"year":2014},
            {"speciesKey":str(100+i),"landmass_id":str(i*2+1),"year":2015},
            {"speciesKey":str(100+i),"landmass_id":str(i*2+2),"year":2019},
        ])
    r=mod.full_year_resolution(screen,keys,raw)
    assert r["fully_resolved_pair_count"]==5
    assert r["distinct_order_and_prior_re_record_count"]==5
    assert r["rows"][0]["earlier_indexed_island"]=="1"
    assert r["rows"][0]["other_island_recorded_after_years"]==5
    assert not r["not_biological_reproduction_or_actual_island_relay"] is False


def test_conflict_with_audited_2011_2017_year_fails_closed():
    screen=make_screen()
    keys={f"Species {i}":str(100+i) for i in range(5)}
    bad=[{"speciesKey":"100","landmass_id":"1","year":2015}]
    with pytest.raises(ValueError,match="Source-year disagreement"):
        mod.full_year_resolution(screen,keys,bad)
