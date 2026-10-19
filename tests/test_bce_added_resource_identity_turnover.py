"""Source-identity turnover is different from net new-area and temporal range change."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/summarize_bce_added_resource_identity_turnover.py"
spec=importlib.util.spec_from_file_location("bce_turnover",SCRIPT)
assert spec and spec.loader
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def stub(name:str,old:int,new:int,reclassified:int,found:int) -> dict:
    return {
        "species":name,"base_added":old,"augmented_added":new,
        "old_added_regions_reclassified_as_native":reclassified,
        "previously_unadded_regions_newly_added":found,
        "BCE_candidate_host_links_added":1,
    }


def test_equal_area_but_different_places_not_mislabelled_stable() -> None:
    row=m.species_turnover(stub("Thymelicus lineola",65,65,27,27))
    assert row["net_change"]==0
    assert row["retained_new_resource_cells"]==38
    assert row["union_new_resource_cells"]==92
    assert row["jaccard"]==pytest.approx(38/92)
    assert row["is_turnover_despite_unchanged_count"]


def test_adding_native_host_reclassifies_previous_novel_regions() -> None:
    row=m.species_turnover(stub("Minois dryas",157,47,125,15))
    assert row["net_change"]==-110
    assert row["retained_new_resource_cells"]==32
    assert row["union_new_resource_cells"]==172
    assert row["jaccard"]==pytest.approx(32/172)
    assert row["both_reclassification_and_discovery"]


def test_source_without_new_regions_has_undefined_overlap() -> None:
    row=m.species_turnover(stub("No new regions",0,0,0,0))
    assert row["union_new_resource_cells"]==0
    assert row["jaccard"] is None
    with pytest.raises(ValueError,match="cannot reconcile"):
        m.species_turnover(stub("Invalid",1,1,2,2))


def test_full_panel_and_affected_panel_jaccard_are_not_identical() -> None:
    retained_81=4628-1403
    union_81=4628+1817
    retained_all=14553-1403
    union_all=14553+1817
    assert retained_81==3225 and union_81==6445
    assert retained_all==13150 and union_all==16370
    assert retained_81/union_81==pytest.approx(0.5003878975950349)
    assert retained_all/union_all==pytest.approx(0.8032987171655467)
    assert (1403+1817)/union_81==pytest.approx(0.4996121024049651)
