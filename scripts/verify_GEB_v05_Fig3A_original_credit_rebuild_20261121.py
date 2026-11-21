#!/usr/bin/env python3
"""Verify a completed Fig3A full-input HOSTS/WCVP source CI audit RECEIPT.

This offline test does NOT redo the full R 7-minute rebuild. It verifies the
permanently archived source receipt against original source hashes and source
numeric report, and the current GEB source/rights ledger. An original end-to-end
rebuild was executed in GitHub Actions #44847332630, not by this offline test.
"""
from __future__ import annotations
import hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"results/geb_v05_private_source_recovery_20261120"
REBUILT=DATA/"ORIGINAL_FIG3A_WCVP_HOSTS_CREDIT_REBUILT_SOURCE_QA_20261121.json"
LEDGER=DATA/"THREE_LIMITED_GEB_REPLAY_COHORT_AND_RELEASE_GATE_20261121.json"
ORIGINAL=ROOT/"provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json"
SHA="af5d7327bcb8883c85be94250a218348b40906e8e35b5f8636a87cea21960eff"

def blob(path):
    raw=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def check():
    raw=REBUILT.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==SHA
    x=json.loads(raw)
    original=json.loads(ORIGINAL.read_text(encoding="utf8"))
    ledger=json.loads(LEDGER.read_text(encoding="utf8"))
    assert x["status"]=="ORIGINAL_FIXED_SOURCE_FIG3A_FRACTIONAL_PLANT_CREDIT_RECOMPUTED_AGREES_WITH_2026_ORIGINAL_RECEIPT"
    assert x["schema"]=="GEB_Fig3A_original_WCVP_HOSTS_host_credit_REBUILT_20261121_v1"
    assert x["upstream_source_commits"]=={
      "HOSTS":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
      "WCVP_v13":"65bed76bae9d644ccb6ad200c05f9f5071d89e05"}
    assert blob(ROOT/"scripts/build_wcvp_hosts_null_sidecars.R")==x["original_source_code_pins"]["R_sidecar_git_blob"]
    assert blob(ROOT/"scripts/analyze_butterfly_host_contribution_concentration.py")==x["original_source_code_pins"]["Python_credits_git_blob"]
    assert hashlib.sha256((ROOT/"data/frozen/figure_sources/s1_resource_descriptors.csv").read_bytes()).hexdigest()==x["original_source_code_pins"]["S1_descriptors_sha256"]
    assert x["original_butterflies_with_new_opportunities"]==original["reconstruction"]["expanded_butterflies"]==206
    assert x["introduced_added_butterfly_region_units"]==original["reconstruction"]["added_butterfly_x_wgsrpd3_units"]==14553
    assert x["contributors_plant_species"]==original["plant_species"]["contributors"]==670
    assert x["contributors_plant_species_to_half"]==original["plant_species"]["contributors_for_50pct"]==38
    assert x["contributing_genera"]==original["genera"]["contributors"]==431
    assert x["contributing_genera_to_half"]==original["genera"]["contributors_for_50pct"]==25
    assert math.isclose(x["plant_credit_gini"],original["plant_species"]["gini"],abs_tol=1e-12)
    assert math.isclose(x["genus_credit_gini"],original["genera"]["gini"],abs_tol=1e-12)
    assert math.isclose(x["plant_credit_top10_share"],original["plant_species"]["top_k_share"]["10"],abs_tol=1e-12)
    assert x["all_original_report_stat_values_match_1e9"]
    assert len(x["source_intermediate_and_output_sha256"])==8
    assert x["source_hosts_WCVP_tables_NOT_REUPLOADED"]
    assert x["source_species_level_credit_tables_NOT_RELEASED"]
    assert x["review_deposit_rights_NOT_CLEARED"]
    z=ledger["fig3A_original_upstream_credit_reconstruction"]
    assert z["source_actions_run"]==44847332630
    assert z["receipt_archive_actions_run"]==44848263272
    assert z["receipt_sha256"]==SHA
    assert z["fractional_credited_plant_species"]==670
    assert z["plant_species_covering_50_percent_credit"]==38
    assert z["added_butterfly_x_region_units"]==14553
    assert z["not_in_earlier_seven_file_fig1_fig3BC_ZIP"]
    assert not z["source_rights_verified"] and not z["final_anonymous_reviewer_data_deposit"]
    assert len(ledger["components"])==3
    print("ORIGINAL_FIG3A_PINNED_WCVP_HOSTS_670_PLANT_CREDIT_FULL_REBUILD_RECEIPT_AUDIT_PASS")
    print("FIG3A_SEPARATE_7FILE_OLD_ZIP_NOT_REPLAYING_PLANT_EDGES_AND_NO_RIGHTS_APPROVAL_PASS")
if __name__=="__main__":
    check()
