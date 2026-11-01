"""Synthetic security and source-grain tests of original 17th century VOC legs."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/derive_dataversenl_17c_VOC_island_routes.py"


def load():
    spec=importlib.util.spec_from_file_location("realVOC17c",SRC)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_real_source_DataverseNL_field_columns_and_precommitted_harbors():
    m=load()
    assert m.EXPECTED_ROWS==57669
    assert m.EXPECTED_SHA1=="eed2479019a8570e392f37f3d781d5b0c67448b1"
    assert m.PORT_FROM=="VERPLAATS,C,52"
    assert m.PORT_TO=="APLAATS,C,51"
    assert m.DATE_FROM=="VERTREK,D"
    assert m.DATE_TO=="AANKOMST,D"
    strict,primary,broad=m.load_aliases(include_broad=False)
    assert strict["batavia"]=="19"
    assert strict["japara"]=="19"
    assert strict["atjeh"]=="12"
    assert strict["gale"]=="31"
    assert strict["macassar"]=="17"
    assert strict["tayouan"]=="44"
    assert strict["deshima"]=="43"
    assert "ceylon" not in strict
    assert "taiwan" not in strict
    relaxed,_,_=m.load_aliases(include_broad=True)
    assert relaxed["ceylon"]=="31"
    assert relaxed["taiwan"]=="44"
    assert relaxed["mauritius"]=="218"
    assert len(primary)==15
    assert len(broad)==4
    assert "amboina" not in relaxed
    assert "piscadores" not in relaxed


def test_historical_date_parser_and_no_false_20th_century_voyages():
    m=load()
    assert m.source_year("1652-02-03")==1652
    assert m.source_year("16520301")==1652
    assert m.source_year("1707/01/31")==1707
    assert m.source_year("1795")==1795
    assert m.source_year("1594") is None
    assert m.source_year("2020-01-01") is None
    assert m.source_year("") is None
    assert m.year_bin(1595)=="1595_1649"
    assert m.year_bin(1687)=="1650_1699"
    assert m.year_bin(1700)=="1700_1749"
    assert m.year_bin(1795)=="1750_1795"


def test_historic_17c_places_are_not_2019_wpi_port_ids_or_cargo():
    m=load()
    strict,_,_=m.load_aliases(include_broad=False)
    assert m.normalize_port("  BATAVIA   ")=="batavia"
    assert m.normalize_port("  BANTAM ")=="bantam"
    assert strict["batavia"]==strict["bantam"]=="19"
    assert m.normalize_port("Tayouan")=="tayouan"
    assert m.LOCK.exists()
    # The biological outcome is original 51-island indexed host/GBIF record
    # Jaccard, not the 1600s voyage count.
    assert "island_host_network_20261030" in str(m.ORIGINAL)
