"""Replay original VOC 8178 ↔ 2020 AIS ↔ pre2011 harbor botany coverage.

The tests distinguish evidence *coverage* from species absence or routes.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"provenance/reviewer_defenses/results/island_observational_robustness_20261031/audit_cross_era_island_historical_routes_plant_coverage.py"


def mod():
    spec=importlib.util.spec_from_file_location("island_cross_era",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_source_verified_historical_current_ports_coverage_not_interpreted_as_route_absence():
    m=mod()
    if not m.VOYAGE.exists() or not m.BOTANY.exists():
        pytest.skip("Historical 51-island source artifacts are not in anonymized GEB review")
    table,result=m.audit()
    assert len(table)==result["original_island_nodes"]==51
    assert result["source_historic_DAS_city_endpoint_islands"]==4
    assert set(result["source_historic_DAS_city_endpoint_island_ids"])=={"12","17","19","31"}
    assert set(result["historic_city_endpoint_and_early_plant_observation_island_ids"])=={"17","31"}
    assert result["historic_VOC_ports_with_2001_2010_GBIF_index_zero"]==2
    # A wider one-to-one geocoding radius need not monotonically increase
    # matches: an extra second WPI port makes an otherwise valid site ambiguous.
    assert result["AIS_match_radius_sensitivity"]["2"]["historic_four_source_island_ids_with_ais_endpoint"]==["19"]
    assert set(result["AIS_match_radius_sensitivity"]["5"]["historic_four_source_island_ids_with_ais_endpoint"])=={"19","31"}
    assert set(result["AIS_match_radius_sensitivity"]["10"]["historic_four_source_island_ids_with_ais_endpoint"])=={"12","31"}
    assert result["AIS_match_radius_sensitivity"]["5"]["triple_source_island_ids"]==["31"]
    assert result["AIS_match_radius_sensitivity"]["10"]["triple_source_island_ids"]==["31"]
    assert not any("proof" in str(x).lower() for x in table)


def test_no_missing_port_botany_as_certain_absence():
    m=mod()
    assert "observation" in m.__doc__.lower()
    assert m.VOYAGE.name.startswith("das_1595_")
