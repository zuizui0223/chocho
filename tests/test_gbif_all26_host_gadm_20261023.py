"""28/24/44 source-only GBIF botanical GADM tests, offline."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/audit_gbif_all26_host_gadm_20261023.py"
spec=importlib.util.spec_from_file_location("all26_plant_facets",SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
BASE=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"

def test_all_26_are_source_selected_before_botanical_gbif_outcomes():
    data,all48,old,names,pairs,prior,meta=m.load_inputs(
        BASE/"gbif_all26_source_host_country_inventory_20261023.json",
        BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json",
        BASE/"gbif_five_host_plant_gadm_source_20261021.json")
    assert len(data["all_48_positive_country_cells"])==48
    assert len([r for r in data["all_48_positive_country_cells"] if r["added_host_resource_country"]])==26
    assert len(names)==24 and len(pairs)==44 and len(prior)==20
    assert sum(len(v["butterflies_linked"]) for v in pairs.values())==46
    assert meta["source_inventory_sha256"]

def test_fail_closed_on_incomplete_taxon_data():
    with pytest.raises(ValueError):
        m.load_inputs(BASE/"gbif_five_host_plant_gadm_source_20261021.json",
           BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json",
           BASE/"gbif_five_host_plant_gadm_source_20261021.json")

def test_no_butterfly_causal_claim_when_plant_facets_missing():
    data,all48,old,names,pairs,prior,meta=m.load_inputs(
        BASE/"gbif_all26_source_host_country_inventory_20261023.json",
        BASE/"gbif_positive48_gadm_lifestage_facets_20261020.json",
        BASE/"gbif_five_host_plant_gadm_source_20261021.json")
    taxa={pid:{"status":"UNRESOLVED_OR_NOT_EXACT_ACCEPTED_SPECIES"} for pid in names}
    facets={key:{"status":"SOURCE_PLANT_TAXON_UNRESOLVED_NOT_ABSENCE"} for key in pairs}
    z=m.audit(data,all48,old,names,pairs,prior,taxa,facets,meta)
    assert z["status"]=="PARTIAL_GBIF_OR_TAXONOMY_NOT_PLANT_ABSENCE"
    assert z["complete_butterfly_country_cases"]==0
    assert z["plant_country_GADM_facets_complete"]==0
    assert z["source_only_host_geography_not_local_host_use"] is True
