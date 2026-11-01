"""Exact pinned HOSTS odd larval-plant records remain a SOURCE question."""
import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/audit_VOC_two_exact_HOSTS_relationship_sources.py"


def load():
    spec=importlib.util.spec_from_file_location("voc_host_original",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_only_exact_frozen_Hosts_rows_are_included():
    x=load()
    rows=[
        {"Insect Genus":"Zizina","Insect Species":"otis",
         "Hostplant Genus":"Indigofera","Hostplant Species":"zollingeriana",
         "Source":"Historic checklist","Life Stage":"UNKNOWN"},
        {"Insect Genus":"Pieris","Insect Species":"rapae",
         "Hostplant Genus":"Nicotiana","Hostplant Species":"tabacum",
         "Reference":"Original HOSTS row"},
        {"Insect Genus":"Pieris","Insect Species":"rapae",
         "Hostplant Genus":"Brassica","Hostplant Species":"rapa"},
    ]
    out=x.classify_rows(rows)
    assert out["source_rows"]==3
    assert len(out["matched"])==2
    assert out["matched"]["pieris rapae nicotiana tabacum"]["original_exact_match_count"]==1
    assert out["matched"]["zizina otis indigofera zollingeriana"]["source_annotations"][0]["Life Stage"]=="UNKNOWN"
    assert "Source" in out["source_provenance_field_names"]
