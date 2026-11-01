"""Source-free regression checks for dated 17th-century island voyages and cargo."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031"


def load(module_file,module_name):
    sys.path.insert(0,str(SRC))
    spec=importlib.util.spec_from_file_location(module_name,SRC/module_file)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_real_dataverse_date_format_is_month_day_year_not_year_first():
    x=load("derive_dataversenl_17c_VOC_island_routes.py","historic_voc_year")
    assert x.source_year("01/31/1675")==1675
    assert x.source_year("12/17/1645")==1645
    assert x.source_year("02/29/1600")==1600
    assert x.source_year("02/29/1700") is None
    assert x.source_year("13/01/1675") is None
    assert x.source_year("1652-02-03")==1652
    assert x.source_year("16520301")==1652
    assert x.source_year("1795")==1795
    assert x.source_year("2020-01-01") is None


def test_source_locked_commodity_categories_are_not_species_identities():
    mod=load("join_VOC_17c_cargo_to_historic_island_routes.py",
             "historic_voc_plant_commodity")
    lock=json.loads(mod.LOCK.read_text(encoding="utf8"))
    assert lock["schema"]=="chocho_17c_VOC_archived_cargo_commodity_category_lock_v1"
    c=lock["categories"]
    assert set(mod.classify("PEPER en KANEEL",c))=={"pepper","cinnamon"}
    assert set(mod.classify("rijst en sandelhout",c))=={"rice","sandalwood"}
    assert mod.classify("silver and silk only",c)==[]
    assert "generic_seeds_plants" in c
    assert "propagule" in " ".join(lock["classification_downgrade"]).casefold()
