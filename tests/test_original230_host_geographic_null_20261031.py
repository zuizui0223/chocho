"""Validate positive overlap attenuation under geo-matched host-source nulls.

Public snapshot only; no redistributing original author-only HOSTS source
or GBIF raw butterfly/plant occurrence coordinates.
"""
from __future__ import annotations
import ast,hashlib,json,math
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/host_geographic_prevalence_null_20261031"
F=D/"host_geo_frequency_null_results.json"
SCRIPT=D/"audit_geographic_host_null.py"
BASE=ROOT/"provenance/reviewer_defenses/results/original230_dated_host_join_20261031"
def blob(p):
 raw=p.read_bytes()
 return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def test_source_geo_null_population_and_structural_integrity():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Original source host and new island observational analyses are excluded from GEB anonymous bundle")
 x=json.loads(F.read_text(encoding="utf8"))
 ast.parse(SCRIPT.read_text(encoding="utf8"))
 assert blob(BASE/"run_dated_host_recurrence.py")==x["source_original_algo_git_blob_sha"]
 assert blob(BASE/"original230_prior_host_recurrence_audit.json")==x["source_original_frozen_audit_git_blob_sha"]
 assert (x["source_species_island_events"],x["source_geographic_islands"],x["source_geographic_near_island_pairs"])==(113,17,14)
 assert x["resolved_butterfly_host_edges"]==2396 and x["host_source_original_edges"]==2524
 assert x["observed_compatible_preindexed_host_events"]==34
 assert (x["observed_any_five_year_reindex_host_yes"],x["observed_any_five_year_reindex_host_no_index"])==(26,70)
 assert x["taxon_mapping_gzip_sha256"]=="d8d21e6686e3f069d4374d4be1c4e47744788f9244c46a1bbbbb9015e5cd7fb9"
 n=x["prior_family_and_degree_only_original_null_replayed"]
 assert (n["n_draws"],n["draws_ge_observed"],n["draws_ge_observed_recorder_delta"])==(299,10,194)
 assert len(x["source_species_and_geographic_regions"])==6
 assert sum(len(r["islands"]) for r in x["source_species_and_geographic_regions"])==17
 assert x["original_host_edges_and_original_GBIF_points_unpublished"]
 assert x["original_host_non_detection_never_treated_as_host_absence"]

def test_plant_geographic_sampling_explains_unconditioned_host_overlap_excess():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Exploratory historical host network evidence excluded from GEB blinded bundle")
 x=json.loads(F.read_text(encoding="utf8"))
 out=x["more_constrained_structural_nulls"]
 assert set(out)=={"family_plus_island_prevalence_band","family_plus_region_signature_and_prevalence_band"}
 a=out["family_plus_island_prevalence_band"]
 b=out["family_plus_region_signature_and_prevalence_band"]
 assert (a["draws"],b["draws"])==(499,499)
 assert (a["null_draws_ge_observed_fraction"],b["null_draws_ge_observed_fraction"])==(293,356)
 assert (a["null_draws_ge_observed_reindex_delta"],b["null_draws_ge_observed_reindex_delta"])==(485,381)
 assert (a["swaps_accepted"]>1000000 and b["swaps_accepted"]>1000000)
 assert (a["edges_in_swappable_strata"]>2100 and b["edges_in_swappable_strata"]>2100)
 for k in (a,b):
  assert math.isclose(k["observed_prior_host_index_fraction"],34/113,abs_tol=1e-12)
  assert k["null_fraction_q025_q975"][0]<34/113<k["null_fraction_q025_q975"][1]
  assert abs(k["null_mean_fraction"]-34/113)<0.004
  assert k["each_plant_degree_and_each_butterfly_plant_family_count_preserved"]
  assert k["plant_presence_pre2011_17islands_preserved_by_band"]
  assert k["no_climate_or_observer_effort_control"]
 assert b["plant_geographic_region_signature_exact_preserved"]
 assert not a["plant_geographic_region_signature_exact_preserved"]
 assert x["not_tested_direct_mainland_source_or_butterfly_breeding"]
