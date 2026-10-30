"""Crosswalk checks: genus is NOT species-level shipping evidence."""
import importlib.util
import json
import pathlib
import pytest

ROOT=pathlib.Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/audit_container_seed_focal_crosswalk.py"
spec=importlib.util.spec_from_file_location("container_seed_crosswalk",SCRIPT)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_port_seed_crosswalk_does_not_equate_congeners():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island source is intentionally excluded from blinded GEB export")
    x=audit.crosswalk(ROOT)
    stored=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/container_seed_host_crosswalk.json"
    assert json.loads(stored.read_text(encoding="utf-8")) == x
    assert x["verified_species_level_identity_for_container_plant_vs_port_host"]==0
    assert x["genus_level_overlapping_focal_groups"]==2
    focal={r["container_study_taxon"]:r for r in x["four_focal_taxon_crosswalk"]}
    assert focal["Phragmites spp."]["port_species_with_same_genus"]==[
        {"candidate_host_species":"Phragmites australis","ports_with_current_alien_candidate":3}]
    assert focal["Saccharum spontaneum"]["port_species_with_same_genus"]==[
        {"candidate_host_species":"Saccharum officinarum","ports_with_current_alien_candidate":9}]
    assert focal["Saccharum spontaneum"]["exact_species_level_overlap"] is None
    assert focal["Typha spp."]["port_species_with_same_genus"] == []
    assert focal["Andropogon spp."]["port_species_with_same_genus"] == []
    hist={r["source_historically_early_host"]:
          r["ports_in_2001_2024_retrospective_alien_candidate_inventory"]
          for r in x["early_host_species_port_inventory_counts"]}
    assert hist == {"Alnus incana":57,"Lepidium didymum":151,
                    "Medicago sativa":202,"Zea mays":53}
    assert all(not r["paper_named_focal_transport_study_confirmation"]
               for r in x["early_host_species_port_inventory_counts"])
    assert x["no_actual_vessel_transport_or_host_plant_origin_assigned"]
