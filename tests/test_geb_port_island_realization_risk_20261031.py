"""Public aggregate checks for private original-year source GEB-to-island gateway risk.

This test deliberately does NOT claim a private original Parquet source replay
and does NOT treat plants not indexed at selected WPI port circles as absent.
"""
from __future__ import annotations
from pathlib import Path
import ast,hashlib,json,math
import pytest

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/geb_port_island_realization_20261031"
RAW=D/"GEB_port_island_realization_risk_20261031.json"
SCRIPT=D/"recompute_geb_port_activation_risk.py"
PUBLIC_GEB=ROOT/"data/frozen/figure_sources/host_contribution_metrics.csv"
PUBLIC_225=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv"

def gitblob(path):
 x=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(x)).encode()+b"\0"+x).hexdigest()

def test_private_original_year_hash_and_source_subpopulation():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Original island source and exploratory gateway excluded from anonymous GEB review bundle")
 z=json.loads(RAW.read_text(encoding="utf8"))
 ast.parse(SCRIPT.read_text(encoding="utf8"))
 assert z["original_private_species_island_year_parquet_SHA256"]=="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
 assert gitblob(PUBLIC_GEB)==z["source_GEB_contribution_metrics_git_blob"]
 assert gitblob(PUBLIC_225)==z["source_225_butterfly_name_GBIF_key_git_blob"]
 assert (z["all_island_first_2011_2017_before_original230_name_screen"]["events"],
         z["all_island_first_2011_2017_before_original230_name_screen"]["species"])==(121,64)
 assert z["selected_name_join"]["source_events"]==121
 assert (z["selected_name_join"]["resolved_event_count"],z["selected_name_join"]["resolved_larval_host_species"],z["selected_name_join"]["unresolved_source_events"])==(113,60,8)
 assert len(z["geographic_components"])==6
 assert len({x for comp in z["geographic_components"] for x in comp})==17
 assert not (ROOT/"work/barriers/global_original_species_island_year.parquet").exists()

def test_at_risk_vs_nonindex_not_absent_and_holdout_null():
 if "_anonymous_stage" in ROOT.parts:
  pytest.skip("Exploratory private-year risk panel not part of independent GEB manuscript")
 d=json.loads(RAW.read_text(encoding="utf8"))
 early=d["periods"]["2011_2017_first_index_among_no_pre2011_butterfly_index"]
 later=d["periods"]["2018_2022_first_index_among_still_unindexed_after_2017"]
 assert (early["recording_risk"]["cells"],early["recording_risk"]["indexed_first_events"])==(905,113)
 assert (later["recording_risk"]["cells"],later["recording_risk"]["indexed_first_events"])==(792,33)
 assert (early["recording_risk"]["host_early_recorded"]["risk_cells"],
         early["recording_risk"]["host_early_recorded"]["new_indexed_events"])==(207,34)
 assert (early["recording_risk"]["host_early_NOT_indexed"]["risk_cells"],
         early["recording_risk"]["host_early_NOT_indexed"]["new_indexed_events"])==(698,79)
 assert (later["recording_risk"]["host_early_recorded"]["risk_cells"],
         later["recording_risk"]["host_early_recorded"]["new_indexed_events"])==(173,20)
 assert (later["recording_risk"]["host_early_NOT_indexed"]["risk_cells"],
         later["recording_risk"]["host_early_NOT_indexed"]["new_indexed_events"])==(619,13)
 assert early["recording_risk"]["cells"]-early["recording_risk"]["indexed_first_events"]==later["recording_risk"]["cells"]
 for z,observed,exceeds in [(early,34,1047),(later,20,770)]:
  r=z["recording_risk"];n=z["geographic_margin_null"]
  assert r["no_index_is_not_botanical_absence"]
  assert sum(x["risk_cells"] for x in r["by_geographic_component"])==r["cells"]
  assert sum(x["first_index_events"] for x in r["by_geographic_component"])==r["indexed_first_events"]
  assert sum(x["first_index_events_at_preindexed_compatible_host"] for x in r["by_geographic_component"])==observed
  assert len(n["chains"])==3 and n["total_draws"]==1497
  assert n["observed_matched_events"]==observed and n["reference_tail_draws"]==exceeds
  assert sum(v["draws_at_least_observed"] for v in n["chains"])==exceeds
  assert math.isclose(n["reference_tail_fraction"],exceeds/1497,rel_tol=0,abs_tol=1e-12)
  assert n["no_independent_mc_pvalue_or_causal_interpretation"]
 assert math.isclose(early["geographic_margin_null"]["null_mean_matched_events"],33.91983967935872,abs_tol=1e-9)
 assert math.isclose(later["geographic_margin_null"]["null_mean_matched_events"],19.50367401469606,abs_tol=1e-9)
 assert d["study_limits"]["no_port_vs_nonport_flora_control"]
 assert d["study_limits"]["no_proof_of_positive_or_negative_butterfly_colonization_effect"]
