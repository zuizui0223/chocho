"""Strict source-free tests of original VOC commodity × modern naturalized host.

No raw historic voyage ID, individual cargo description or GBIF coordinate
is used in these source-fixture tests.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/audit_VOC_cargo_exact_HOSTS_GloNAF_local_plant_species.py"


def load():
    spec=importlib.util.spec_from_file_location("voc_three_source_exact_plants",SCRIPT)
    x=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(x)
    return x


def test_exact_host_species_and_known_naturalized_region_needed_not_just_genus():
    m=load()
    recorded={"eurema hecabe","vanessa cardui"}
    hosts={
        "eurema hecabe":{"indigofera tinctoria","indigofera hirsuta"},
        "vanessa cardui":{"nicotiana tabacum"},
    }
    cargo={"indigo":{"indigofera"},"tobacco":{"nicotiana"}}
    aliases={
        "indigofera tinctoria":(2,),
        "indigofera hirsuta":(3,),
        "nicotiana tabacum":(4,),
    }
    aligned=m.local_triads(recorded,hosts,cargo,aliases,(1<<2)|(1<<4))
    assert aligned=={
        ("eurema hecabe","indigofera tinctoria","indigo"),
        ("vanessa cardui","nicotiana tabacum","tobacco"),
    }
    assert not m.local_triads(recorded,hosts,{"rice":{"oryza"}},aliases,(1<<2)|(1<<4))
    # A published genus on the island without exact species is not enough.
    assert not m.local_triads({"eurema hecabe"},hosts,cargo,aliases,1<<9)
    assert m.normalized_host_binomial("Indigofera","tinctoria")=="indigofera tinctoria"
    assert m.normalized_host_binomial("Indigofera","Indigofera tinctoria")=="indigofera tinctoria"
    assert m.normalized_host_binomial("Indigofera","sp.") is None


def test_pair_requires_same_butterfly_and_same_exact_naturalized_host_species_at_both():
    m=load()
    first={"12":{"eurema hecabe"},"19":{"eurema hecabe"},
           "31":{"vanessa cardui"}}
    host={
        "eurema hecabe":{"indigofera tinctoria"},
        "vanessa cardui":{"nicotiana tabacum"},
    }
    mapping={"indigo":{"indigofera"},"tobacco":{"nicotiana"}}
    aliases={"indigofera tinctoria":(1,),"nicotiana tabacum":(2,)}
    masks={"SUM":1<<1,"JAW":1<<1,"SRL":1<<2}
    source=m.aggregate_pair("12","19",{"indigo"},first,host,mapping,aliases,masks)
    assert source["n_both_exact_triplets"]==1
    assert source["identical_butterfly_host_species_category_both"]==[
        ("eurema hecabe","indigofera tinctoria","indigo")]
    missing=m.aggregate_pair("12","31",{"indigo","tobacco"},
             first,host,mapping,aliases,masks)
    assert missing["n_both_exact_triplets"]==0
    assert missing["n_left"]==1 and missing["n_right"]==1


def test_six_original_islands_and_nine_historical_source_links_remain_frozen():
    m=load()
    historic,plant_categories,first=m.load_sources()
    assert len(historic)==9
    assert len(first)==6
    assert sum(len(x) for x in first.values())==49
    assert sum(int(x["source_cargo_rows_joined_by_unrepeated_REISN"]>0)
               for x in historic)==7
    assert m.SIX=={
        "12":"SUM","17":"SUL","19":"JAW","31":"SRL","43":"JAP","44":"TAI"}
    assert m.CUTOFFS==(.01,.05,.1)
    assert plant_categories["indigo"]=={"indigofera"}
