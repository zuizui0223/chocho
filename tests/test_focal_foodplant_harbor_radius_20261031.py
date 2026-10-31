"""Guard historical butterfly foodplant 1km versus 5km harbor-radius inference.

The 3 focal butterfly species share only 4 independent foodplant-island
contexts; no 11 biologically independent vessel introductions are claimed.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/focal_host_harbor_radius_20261031"
SOURCES=ROOT/"provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031"

def test_earlier_host_years_vs_1km_coastal_indices_are_source_reproducible():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Independent GEB anonymous bundle redacts exploratory island botanicals and changes pinned source SHA")
    if not (D/"focal_harbor_host_1km_vs_5km_results.json").exists():
        pytest.skip("Exploratory focal island sources excluded from independent blinded GEB bundle")
    spec=importlib.util.spec_from_file_location("focal_host_radius",D/"reconcile_focal_host_radius.py")
    assert spec is not None and spec.loader is not None
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    computed=mod.run()
    expected=json.loads((D/"focal_harbor_host_1km_vs_5km_results.json").read_text(encoding="utf8"))
    assert computed==expected
    assert computed["independent_host_island_units"]==4
    assert computed["butterfly_species_island_mentions_not_independent"]==6
    assert computed["butterfly_species"]==3 and computed["plant_species"]==2
    assert computed["butterfly_species_matched_to_published_author_selected_migrants"]==3
    assert all(x["migrant_proxy"]=="REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"
               for x in computed["butterfly_endpoints"])
    assert computed["total_one_km_port_times_host_island_combinations"]==11
    assert computed["total_one_km_offport_times_host_island_combinations"]==22
    assert computed["one_km_port_plant_indexes_per_host_island"]==0
    assert computed["one_km_offport_plant_indexes_per_host_island"]==2
    assert set(x["plant"] for x in computed["host_island"])=={"Medicago sativa","Alnus incana"}
    assert all(x["one_km_all_sites_have_positive_botanical_index"] for x in computed["host_island"])

def test_fresh_5km_source_requery_discriminates_outer_ring_with_no_port_center_index():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Independent GEB anonymous bundle redacts exploratory island botanicals and changes pinned source SHA")
    if not (D/"focal_nested_harbor_queries.csv").exists():
        pytest.skip("Separate island study sources withheld from independent GEB anonymous bundle")
    with (D/"focal_nested_harbor_queries.csv").open(newline="",encoding="utf8") as f:
        data=list(csv.DictReader(f))
    x=json.loads((D/"focal_nested_harbor_query_results.json").read_text(encoding="utf8"))
    assert len(data)==11
    assert x["focal_host_x_port_samples"]==11
    assert x["original_5km_source_focal_combinations"]==4
    assert all(r["gbif_1km_two_source_query_concordant"]=="True" for r in data)
    assert all(int(r["gbif_1km_exact_accepted_speciesKey_records"])==0 for r in data)
    assert all(int(r["gbif_5km_exact_accepted_speciesKey_records"])>0 for r in data)
    assert all(int(r["gbif_indexed_ring_1to5km_record_count_difference"])>0 for r in data)
    assert all(r["historic_plant_absence_inferred"]=="False" and r["ship_transport_inferred"]=="False" for r in data)
    assert {r["key"] for r in data}=={"9151957","9148577"}
    assert all(v["one_to_five_ring_sites_with_index"]==v["port_sites"] for v in x["host_island_summary"])
    assert all(v["one_km_port_sites_with_index"]==0 for v in x["host_island_summary"])

def test_no_false_harbor_introduction_or_butterfly_rescue_claim():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Independent GEB anonymous bundle redacts exploratory island botanicals and changes pinned source SHA")
    path=D/"README.md"
    if not path.exists():
        pytest.skip("Focal butterfly harbor comparison excluded from independent GEB review")
    t=path.read_text(encoding="utf8")
    assert "11" in t and "5km" in t and "1km" in t
    assert "島" in t and "非港" in t
    assert "不在" in t
    assert "船舶" in t
    context=json.loads((D/"powo_denmark_introduced_context_20261031.json").read_text(encoding="utf8"))
    assert context["country"]=="Denmark"
    assert {x["plant"] for x in context["plant_records"]}=={"Medicago sativa","Alnus incana"}
    assert all(x["powo_country_status"]=="introduced"
               and x["island_specific_status_not_confirmed"] for x in context["plant_records"])
    assert "POWO" in context["status"]
