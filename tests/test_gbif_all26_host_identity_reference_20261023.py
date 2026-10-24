"""Exact combinatorial bipartite degree/reference checks, no web."""
from __future__ import annotations

import importlib.util
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_gbif_all26_host_identity_reference_20261023.py"
sp=importlib.util.spec_from_file_location("host_identity_exact",SCRIPT)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

BASE=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"

def test_two_by_two_fixed_margin_exact_reference():
    b=[{"butterfly":"A a","candidate_k":1,"gadm":{"X.1":1},"N":1,"observed_count":1},
       {"butterfly":"B b","candidate_k":1,"gadm":{"Y.1":1},"N":1,"observed_count":1}]
    p=[{"accepted_id":"h1","links":frozenset({"A a"}),"gadm":frozenset({"X.1"})},
       {"accepted_id":"h2","links":frozenset({"B b"}),"gadm":frozenset({"Y.1"})}]
    v=m.enumerate_graphs(b,p)
    assert v["states"]==2
    assert v["joint"]==Counter({(2,0):1,(0,0):1})
    assert v["observed_all"]==2 and v["observed_n5"]==0
    d,E=m.exact_pmf([v],0)
    assert E==1 and sum(d.values())==2
    assert Fraction(sum(x for k,x in d.items() if k>=2),sum(d.values()))==Fraction(1,2)

def test_preserve_shared_plant_column_degrees():
    b=[{"butterfly":"A a","candidate_k":1,"gadm":{"X.1":2},"N":2,"observed_count":2},
       {"butterfly":"B b","candidate_k":1,"gadm":{"X.1":1},"N":1,"observed_count":1}]
    p=[{"accepted_id":"h1","links":frozenset({"A a","B b"}),"gadm":frozenset({"X.1"})}]
    v=m.enumerate_graphs(b,p)
    assert v["states"]==1
    assert v["joint"]==Counter({(3,0):1})
    assert v["observed_all"]==3

def test_reject_invalid_bipartite_margins():
    with pytest.raises(ValueError,match="degrees"):
        m.enumerate_graphs([{"butterfly":"A a","candidate_k":2,"gadm":{"X.1":1},"N":1,"observed_count":1}],
         [{"accepted_id":"h1","links":frozenset({"A a"}),"gadm":frozenset({"X.1"})}])

def test_original_26_source_reconstructs_exactly():
    bot=json.loads((BASE/"gbif_all26_host_gadm_source_20261023.json").read_text(encoding="utf8"))
    bf=json.loads((BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json").read_text(encoding="utf8"))
    cases,plants=m.build_countries(bot,bf)
    assert len(cases)==len(plants)==11
    assert sum(map(len,cases.values()))==26
    assert sum(len(v) for v in plants.values())==44
    assert sum(x["N"] for v in cases.values() for x in v)==453
    assert sum(x["observed_count"] for v in cases.values() for x in v)==231
