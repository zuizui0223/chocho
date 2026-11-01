"""Noncausal historical-port study source and numerical comparability tests."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/test_island_port_colonial_association.py"


def mod():
    pytest.importorskip("numpy")
    sys.path.insert(0,str(SOURCE.parent))
    spec=importlib.util.spec_from_file_location("historic_islands",SOURCE)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_source_42_island_history_support_and_country_shipping_are_distinct():
    m=mod()
    if not (m.BASE/"island_port_colonial_10km.islands.csv").exists():
        pytest.skip("Source colonial WPI geocodes not in anonymous GEB review package")
    ids,pairs,audit=m.joined_pairs(10)
    assert len(ids)==42
    assert len(pairs)==861
    assert audit["eligible_history_islands"]==42
    assert audit["source_2010_country_shipping_pairs"]>500
    assert audit["shared_focal_empire_pairs"]>0
    assert all(r["colonial_shared_duration_years"]>=0 for r in pairs)
    assert all(r["colonial_binary"] in (0,1) for r in pairs)
    assert all(r["host_x_history"]>=0 for r in pairs)
    assert any(r["shipping_country_2010"] is None for r in pairs)
    assert any(r["shipping_country_2010"] is not None for r in pairs)
    assert any(r["colonial_shared_duration_years"]>0 for r in pairs)
    with pytest.raises(ValueError,match="fixed geocoding"):
        m.joined_pairs(7)


def test_rank_partial_is_not_same_as_invented_independent_pair_pvalue():
    m=mod()
    rs=[{"x":float(i%5),"y":float(i%5+(.1 if i%2 else 0)),
         "z":float(i//7)} for i in range(100)]
    v=m.partial_rank(rs,"x","y",("z",))
    assert v is not None and .9<v<=1.
    rs=[{"x":3.0,"y":float(i),"z":float(i%4)} for i in range(100)]
    assert m.partial_rank(rs,"x","y",("z",)) is None


def test_colonial_history_does_not_overwrite_frozen_51_island_evidence():
    m=mod()
    assert "B1_plus_shared_hosts" in m.MODELS
    assert "H1_history_plus_hosts" in m.MODELS
    assert "S0_2010_country_shipping" in m.MARITIME_MODELS
    assert "host" in m.MODELS["B1_plus_shared_hosts"]
    assert "shipping_country_2010" not in m.MODELS["B1_plus_shared_hosts"]
    assert "shipping_country_2010" in m.MARITIME_MODELS["S0_2010_country_shipping"]
