"""CI: lock original island-port site support before claiming gateway mechanism.

Existing coastal WorldCover is 2021 / 1km; indexed GBIF botanics are historical
5km WPI circles. These are NOT 79 sampled port-vs-off-port larval-host contrasts.
"""
from pathlib import Path
import importlib.util
import json
import math
import pytest

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"provenance/reviewer_defenses/results/island_port_first_paper_20261031"
CODE40=DIR/"audit_primary40_control_support.py"
CODE79=DIR/"audit_joint79_port_plant_landcover_support.py"
RESULT40=DIR/"primary_40island_landcover_and_host_support.json"
RESULT79=DIR/"joint_79port_WPI_plant_landcover_support.json"

def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 module=importlib.util.module_from_spec(spec)
 spec.loader.exec_module(module)
 return module

def test_40_island_cover_and_richness_source_replay():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Exploratory island-port evidence omitted from independent GEB review pack")
 d=json.loads(RESULT40.read_text(encoding="utf8"))
 module=load(CODE40,"audit40")
 actual=module.summarize(ROOT)
 assert actual==d
 assert d["preexisting_eligible_islands"]==51
 assert d["islands_with_1km_matched_port_nonport_2021_landcover"]==40
 assert d["sampled_1km_port_sites_in_eligible_islands"]==273
 assert d["historically_dated_plant_islands_with_1km_controls"]==12
 assert len(d["matched_1km_port_aggregate_island_rows"])==40
 assert sum(x["matched_port_sites"] for x in d["matched_1km_port_aggregate_island_rows"])==273
 r=d["island_weighted_exploratory_spearman_correlations"]
 assert math.isclose(r["built_fraction_delta_vs_aliens"],-.05647169673485953,abs_tol=1e-12)
 assert math.isclose(r["built_fraction_delta_vs_later_butterflies"],.08899158829303697,abs_tol=1e-12)
 assert math.isclose(r["built_fraction_delta_vs_prior_butterflies"],.30126638175413084,abs_tol=1e-12)
 assert math.isclose(r["total_island_aliens_vs_later_butterflies"],-.26571022894549956,abs_tol=1e-12)
 assert d["limits"]["not_botanical_host_matched_port_vs_off_port"]
 assert d["limits"]["no_climate_habitat_reporting_effort_adjustment"]

def test_79_identical_ports_plant_5km_and_urban_1km_source_gates():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Exploratory port observational source gate omitted from blinded GEB reviewer pack")
 d=json.loads(RESULT79.read_text(encoding="utf8"))
 module=load(CODE79,"audit79")
 actual=module.summarize(ROOT)
 assert actual==d
 assert d["original_eligible_islands"]==51
 assert d["botanical_2001_2024_indexed_port_catalogue_ids"]==194
 assert d["same_island_2021_landcover_matched_1km_port_ids"]==273
 assert d["identical_WPI_port_and_landmass_ID_intersection"]==79
 assert d["identical_WPI_ports_in_distinct_islands"]==36
 assert d["historical_species_key_facets_from_2001_2010_available_port_IDs"]==15
 assert d["historical_species_key_facets_from_2001_2010_available_islands"]==12
 bins=d["pre2011_indexed_botanical_coverage_thresholds"]
 assert [x["minimum_2001_2010_indexed_botanical_query_records"] for x in bins]==[0,1,10,30,100,300,1000]
 assert [(x["same_port_ID_plant_and_landcover_sites"],x["distinct_islands"]) for x in bins][:5]==[
  (79,36),(67,32),(52,24),(44,19),(37,14)]
 assert len(d["exact_linked_source_2021_landcover_and_2001_2024_botanical_port_IDs"])==79
 limits=d["limits"]
 assert limits["original_botanical_port_radius_km"]==5
 assert limits["same_port_2021_worldcover_coastal_probe_radius_km"]==1
 assert limits["source_port_ID_match_does_not_mean_same_spatial_footprint"]
 assert limits["nonport_site_plant_specieskey_records_not_acquired"]
 assert limits["unknown_island_local_plant_introduction_year_and_shipping_route"]

def test_island_first_submission_fencing():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Independently submitted GEB manuscript unchanged by island first-paper plan")
 paper=(ROOT/"manuscript/island_port_paper_track_20261031/ISLAND_PORT_JBI_SUBMISSION_DECISION.md").read_text(encoding="utf8")
 evidence=(DIR/"ISLAND_PORT_PAPER_GATE_RESULTS_ja.md").read_text(encoding="utf8")
 assert "Journal of Biogeography" in paper
 assert "Biological Invasions" in paper
 assert "Diversity and Distributions" in paper
 for x in ("51島","40島","79","36島","12島","15港","5km","1km"):
  assert x in paper+evidence
 assert "非港沿岸" in paper and "植物" in paper
 assert "## 今回再計算した" in evidence
 assert "ρ=+0.300" in evidence
 assert "−0.056" in evidence
 assert "未取得" in evidence
 assert "投稿操作はしていない" in paper
