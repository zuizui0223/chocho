"""Keep real historical VOC Europe–island journeys distinct from Asian coast hops."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/derive_das_colonial_island_voyage_endpoints.py"


def mod():
    spec=importlib.util.spec_from_file_location("huygens_das_island",SRC)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_precise_city_broad_island_label_and_dutch_metropole_not_conflated():
    m=mod()
    precise,broad,home=m.validate_lock()
    assert m.classify("Batavia",precise,broad,home)==("ISLAND_PRECISE","19")
    assert m.classify("Atjeh",precise,broad,home)==("ISLAND_PRECISE","12")
    assert m.classify("Colombo",precise,broad,home)==("ISLAND_PRECISE","31")
    assert m.classify("Macassar",precise,broad,home)==("ISLAND_PRECISE","17")
    assert m.classify("Texel",precise,broad,home)==("DUTCH_HISTORIC_METROPOLE",None)
    assert m.classify("Ceylon",precise,broad,home)==("UNMAPPED_HISTORIC_PORT_OR_PLACE",None)
    assert m.classify("Ceylon",precise,broad,home,True)==("ISLAND_BROAD","31")
    assert m.classify("Mauritius",precise,broad,home,True)==("ISLAND_BROAD","218")
    assert m.classify("Bali",precise,broad,home,True)[0]=="UNMAPPED_HISTORIC_PORT_OR_PLACE"
    assert m.classify("China",precise,broad,home,True)[0]=="UNMAPPED_HISTORIC_PORT_OR_PLACE"


def test_historical_voyage_category_marks_europe_to_island_not_direct_cargo():
    m=mod()
    home=("DUTCH_HISTORIC_METROPOLE",None)
    java=("ISLAND_PRECISE","19")
    ceylon=("ISLAND_BROAD","31")
    other=("UNMAPPED_HISTORIC_PORT_OR_PLACE",None)
    assert m.voyage_category(home,java)=="DUTCH_METROPOLE_TO_OR_FROM_WPI51_ISLAND"
    assert m.voyage_category(java,home)=="DUTCH_METROPOLE_TO_OR_FROM_WPI51_ISLAND"
    assert m.voyage_category(java,ceylon)=="BETWEEN_WPI51_ISLAND_ENDPOINTS"
    assert m.voyage_category(java,java)=="WITHIN_SAME_WPI51_ISLAND"
    assert m.voyage_category(java,other)=="ISLAND_ENDPOINT_TO_UNMAPPED_OTHER_SOURCE_PORT"
    assert m.voyage_category(home,other)=="NO_WPI51_ISLAND_ENDPOINT"
    assert m.year("1654-02-03")==1654
    assert m.year("1766")==1766
    assert m.year("2020-01-01") is None
    assert m.year("unknown") is None
    assert m.date_bin(1699)=="1650_1699"
    assert m.date_bin(1700)=="1700_1749"


def test_source_VOC_8178_is_not_IISG_57669_or_2020_AIS():
    m=mod()
    assert m.SOURCE_SHA256=="c875a0a59568713347b1d50e4a2f62708006be524133d30ee5ff4fe23f05b77f"
    assert "DAS" in m.SOURCE.name or "voyages" in m.SOURCE.name
