"""Source-consistency and toy checks for indexed local host portfolio overlap."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_gbif_all26_local_host_redundancy_20261024.py"
SPEC=importlib.util.spec_from_file_location("gbif_local_redundancy",SCRIPT)
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)
DIR=ROOT/"provenance/reviewer_defenses/results/gbif_country_20261020"


def source():
    return [json.loads((DIR/p).read_text(encoding="utf-8")) for p in (
        "gbif_all26_host_gadm_source_20261023.json",
        "gbif_positive48_gadm_lifestage_facets_20261020.json",
        "gbif_juvenile_gadm1_conditional_20261021.json",
        "gbif_source_introduced_hosts_five_cases_20261021.json")]


def test_original_country_and_plant_network_margins_reproduce():
    d=m.summary_for(*source())
    assert d["schema"]==m.SCHEMA
    assert d["original_26_cases"]==26
    assert d["original_44_plant_country_queries"]==44
    assert d["original_46_butterfly_country_plant_links"]==46
    assert d["source_original_mapped_butterfly_record_count"]==453
    assert d["source_original_mapped_butterfly_gadm_units"]==106
    r=d["minimum_plant_records_per_species_per_GADM1_sensitivity"]["1"]
    assert r["butterfly_record_degree"]=={"0":222,"1":226,"2plus":5}
    assert r["butterfly_gadm_unit_degree"]=={"0":66,"1":38,"2plus":2}
    assert sum(r["butterfly_record_degree"].values())==453
    assert d["case_counts"]=={"single_source_candidate_host":19,"multiple_source_candidate_hosts":7}
    assert d["five_source_juvenile_cases"]["stage_label_candidate_plant_species_degree_records"]=={"0":8,"1":14,"2plus":2}


def test_threshold_increasing_cannot_create_new_host_record_support():
    d=m.summary_for(*source())
    sensitivities=d["minimum_plant_records_per_species_per_GADM1_sensitivity"]
    previous=10**9
    for k in m.THRESHOLDS:
        r=sensitivities[str(k)]["butterfly_record_degree"]
        covered=r["1"]+r["2plus"]
        assert covered <= previous
        previous=covered
        assert sum(r.values())==453


def test_single_and_multi_host_indexed_support_preserves_source_sample():
    d=m.summary_for(*source())
    r=d["minimum_plant_records_per_species_per_GADM1_sensitivity"]["1"]
    multi=r["multi_candidate_record_degree"]
    assert sum(multi.values()) < 453
    assert sum(multi.values()) > 0
    assert multi["2plus"] == 5
    assert multi["1"] == 124
    assert all(row["largest_index_coverage_loss_after_removing_one_source_plant"][
        "indexed_butterfly_coverage_lost_if_only_this_plant_index_removed"] <=
        row["indexed_butterfly_records_in_at_least_one_candidate_plant_GADM1"]
        for row in d["case_plant_recording_portfolios"]
        if row["catalogue_candidate_host_species"] >= 2)


def test_fail_closed_on_host_network_source_drift():
    p,b,j,f=source()
    p["butterfly_country_results"][0]["source_candidate_host_species"]+=1
    with pytest.raises(ValueError,match="Original butterfly host-network"):
        m.summary_for(p,b,j,f)


def test_degree_class_clamps_2_or_more_to_one_category():
    assert [m.degree_bucket(i) for i in range(5)]==["0","1","2plus","2plus","2plus"]
