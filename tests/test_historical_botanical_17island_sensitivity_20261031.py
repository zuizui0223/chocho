"""Validate historical plant source chronology vs observed near-island butterfly similarity."""
from pathlib import Path
import importlib.util,json,math
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/historical_botanical_17island_20261031"
spec=importlib.util.spec_from_file_location("historical17",D/"validate_historical_butterfly_similarity_sensitivity.py")
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
def test_historical_source_and_quality_gates_are_not_caused_by_future_butterflies():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Exploratory source island data absent from independent GEB reviewer pack")
 expected=json.loads((D/"historical_butterfly_similarity_sensitivity.json").read_text())
 actual=mod.summarize(ROOT)
 assert math.isclose(actual["all"]["rho_historical_broad_plants_vs_new_butterflies"],expected["all"]["rho_historical_broad_plants_vs_new_butterflies"],abs_tol=1e-12)
 assert actual["all"]["n"]==14
 assert actual["acquisition_primary_quality"]["n"]==13
 assert actual["exploratory_stricter_quality"]["n"]==9
 assert abs(actual["exploratory_stricter_quality"]["rho"])<0.1
 assert actual["nonzero_new_pairs"]["n"]==12
 assert sum(x["n_pairs"] for x in actual["regional_components"])==14
 assert len(actual["regional_components"])==6
 assert len(actual["leave_one_geography_group_out"])==6
 assert expected["no_independent_pairwise_significance_tests"]
 assert expected["not_species_specific_larval_host"]
