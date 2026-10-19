"""Exact source host-quality filter never equates unknown evidence to a bad host."""
from __future__ import annotations
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
# The official double-anonymous submission snapshot intentionally excludes
# exploratory independently sourced host-quality audit scripts.
if not (ROOT/"scripts/analyze_monarch_host_quality_geographic_novelty.py").is_file():
    pytest.skip("Exploratory host-quality audit excluded from frozen original GEB snapshot", allow_module_level=True)
import sys
sys.path.insert(0,str(ROOT/"scripts"))
P=ROOT/"scripts/analyze_monarch_host_quality_geographic_novelty.py"
spec=importlib.util.spec_from_file_location("monarch_quality",P)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def miniature_40_HOSTS_case():
    n={str(i):set() for i in range(40)}
    c={str(i):set() for i in range(40)}
    # synthetic 181 native, 263 contemporary: 82 introduced additions
    n["0"]={f"R{i:03}" for i in range(181)}
    c["0"]=n["0"]|{f"R{i:03}" for i in range(181,263)}
    for i in range(1,40):
        n[str(i)]={"R000"}
        c[str(i)]={"R000"}
    links=[{
        "insect_species":"Danaus plexippus",
        "accepted_plant_name_id":str(i),
        "input_host_name":"Asclepias curassavica" if i==0 else f"Auscepa taxon{i}",
        "accepted_name":"Asclepias curassavica" if i==0 else f"Auscepa taxon{i}",
    } for i in range(40)]
    return n,c,links


def test_exact_source_classes_preserve_unmatched_study_hosts():
    n,c,rows=miniature_40_HOSTS_case()
    out=m.reconstruct_quality_pilot(
        {"asclepias curassavica":"H3"}, rows,n,c
    )
    assert out["original_frozen_geography"] == {
        "accepted_HOSTS_species":40,"native_regions":181,
        "contemporary_regions":263,"distinct_added_regions":82
    }
    assert out["classes"]["literature_high_H1_H2_H3"]["accepted_host_species"]==1
    assert out["classes"]["literature_confirmed_high_H3_only"]["source_quality_supported_fraction_of_original_82"]==1
    assert out["classes"]["literature_low_L1_L2_L3"]["source_quality_supported_fraction_of_original_82"]==0
    assert out["classes"]["not_exactly_matched_in_published_127_host_review"]["accepted_host_species"]==39
    assert len(out["taxa_requiring_manual_synonym_review"])==39
    assert out["status"]=="EXPLORATORY_SOURCE_QUALITY_POTENTIAL_GEOGRAPHY_ONLY"
    assert "larval survival" in out["classes"]["all_original_HOSTS"]["note"]


def test_quality_source_coverage_cannot_be_fabricated_from_missing_data():
    n,c,rows=miniature_40_HOSTS_case()
    out=m.reconstruct_quality_pilot({},rows,n,c)
    assert out["classes"]["literature_high_H1_H2_H3"]["accepted_host_species"]==0
    assert out["classes"]["literature_low_L1_L2_L3"]["accepted_host_species"]==0
    assert out["classes"]["not_exactly_matched_in_published_127_host_review"]["accepted_host_species"]==40
    assert out["classes"]["literature_high_H1_H2_H3"]["original_82_cells_supported_by_any_source_quality_class_introduced_host"]==0


def test_true_count_drift_must_fail():
    n,c,rows=miniature_40_HOSTS_case()
    c["0"].remove("R262")
    with pytest.raises(ValueError,match="40/181/263/82"):
        m.reconstruct_quality_pilot({"asclepias curassavica":"H3"},rows,n,c)


def test_conflicting_evidence_classes_fail_closed():
    n,c,rows=miniature_40_HOSTS_case()
    rows[0]["input_host_name"]="Asclepias curassavica"
    rows[0]["accepted_name"]="Araujia sericifera"
    with pytest.raises(ValueError,match="conflicting quality classes"):
        m.reconstruct_quality_pilot({
            "asclepias curassavica":"H3",
            "araujia sericifera":"L3"
        },rows,n,c)


def test_botanical_name_extractor_discards_author_string_but_not_taxon():
    assert m.identity("Asclepias curassavica L.")=="asclepias curassavica"
    assert m.identity("Araujia sericifera (Brot.)")=="araujia sericifera"
    assert m.identity("gibberish")==""


def test_original_82_cells_need_quality_covered_native_counterfactual():
    n,c,rows=miniature_40_HOSTS_case()
    # Make documented plant quality unknown while retaining its native and
    # introduced mapped geography: native and introduced remain guaranteed
    # 181/263, but high-source regional coverage is now zero.
    out=m.reconstruct_quality_pilot({},rows,n,c)
    assert out["original_frozen_geography"]["distinct_added_regions"]==82
    assert out["classes"]["literature_high_H1_H2_H3"]["new_regions_under_quality_filter"]==0
    assert out["classes"]["all_original_HOSTS"]["new_regions_under_quality_filter"]==82
