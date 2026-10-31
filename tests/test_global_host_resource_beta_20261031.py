"""Guard against structural claims in global native/current interaction geography."""
import importlib.util
from collections import Counter
from pathlib import Path
import random

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_global_host_resource_beta.py"


def load():
    spec=importlib.util.spec_from_file_location("global_resource_beta",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_adding_plant_can_decrease_or_increase_pairwise_similarity():
    m=load()
    assert m.jaccard_bits(1,1)==1
    assert m.jaccard_bits(1|2,1)==.5
    assert m.jaccard_bits(1|2,1|2)==1
    assert m.jaccard_bits(1|2,1) < m.jaccard_bits(1|2,1|2)
    assert m.jaccard_bits(0,0) is None


def test_host_edge_rewiring_preserves_both_degree_sequences_and_family():
    m=load()
    original=[(0,"x1"),(0,"y1"),(1,"x2"),(1,"y2"),(2,"x3"),
              (2,"x4"),(3,"x1"),(3,"x4"),(4,"y1"),(4,"y2")]
    fam={f"x{i}":"X" for i in range(1,5)}
    fam.update({f"y{i}":"Y" for i in range(1,3)})
    rewired,accepted=m.swap_family_edges(original,fam,random.Random(17),1200)
    assert accepted>0
    assert len(rewired)==len(original)
    assert len(set(rewired))==len(rewired)
    assert Counter(b for b,p in rewired)==Counter(b for b,p in original)
    assert Counter(p for b,p in rewired)==Counter(p for b,p in original)
    assert Counter((b,fam[p]) for b,p in rewired)==Counter(
        (b,fam[p]) for b,p in original)


def test_resource_availability_is_not_confirmed_butterfly_occurrence():
    m=load()
    regions={"native":{"a":1,"b":2},
             "current":{"a":1|2,"b":2}}
    z=m.build_butterfly_by_region([(0,"a"),(1,"b")],regions,2,2)
    assert z["native"]==[1,2]
    assert z["current"]==[1,3]
    assert m.mean_resource_delta(z)==(1/3,1)
    vals=m.describe_delta([-.1,0,.2,None])
    assert vals["n_valid"]==3
    assert vals["negative_fraction"]==pytest.approx(1/3)
    assert vals["positive_fraction"]==pytest.approx(1/3)


def test_exact_geb_gate_only_accepts_239_species():
    m=load()
    assert m.EXPECTED_CONSUMERS==239
    assert m.EXPECTED_NATIVE_UNITS==26530
    assert m.EXPECTED_CURRENT_UNITS==41083
