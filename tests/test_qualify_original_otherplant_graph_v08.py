"""Source-only WCVP model eligibility for preselected real OTHER plants."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from qualify_original_otherplant_graph_v08 import eligible_actual_and_other_plants

def entry(actual,groups):
    return {"known_host_taxa":actual,"alternative_groups":[
        {"donor_taxa":[donor],"group_index":i,"mean_plant_signature_gap":.1+i*.05}
        for i,donor in enumerate(groups)]}

def test_four_preselected_fully_modelled_plant_groups_required_without_outcomes():
    model={"10","20","30","40","50","60"}
    alt={"A a":entry(["10"],["20","30","40","50","60"]),
         "B b":entry(["10"],["20","30","40"]),
         "C c":entry(["10"],["20","30","40","50"])}
    kept,rows=eligible_actual_and_other_plants(alt,model)
    assert sorted(kept)==["A a","C c"]
    assert len(kept["A a"]["alternative_groups"])==5
    assert next(x for x in rows if x["insect_species"]=="B b")["status"]=="SOURCE_TOO_FEW_MODELED_ALTERNATIVE_PLANT_GROUPS"

def test_missing_actual_host_invalidates_even_when_many_controls_exist():
    model={"10","20","30","40","50"}
    kept,rows=eligible_actual_and_other_plants({"A a":entry(["99"],["10","20","30","40","50"])},model)
    assert not kept
    assert rows[0]["status"]=="SOURCE_ACTUAL_KNOWN_HOST_MODEL_INCOMPLETE"

def test_recorded_host_cannot_become_negative_control():
    model={"10","20","30","40"}
    with pytest.raises(ValueError):
        eligible_actual_and_other_plants({"A a":entry(["10"],["10","20","30","40"])},model)

def test_no_postoutcome_group_threshold_change():
    with pytest.raises(ValueError):
        eligible_actual_and_other_plants({},set(),min_groups=1)
