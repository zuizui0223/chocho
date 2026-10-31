"""Protect the observation-first island/port interpretation from causal overclaiming.

The established descriptive correlations and source-specific time/space forecast
improvements may be presented as empirical associations, irrespective of whether
the independent ship-transport or larval-establishment mechanisms are identified.
"""
from pathlib import Path
import json
import pytest

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/"provenance/reviewer_defenses/results"
REPORT=RESULTS/"geb_port_island_realization_20261031/GEB_port_island_realization_results_ja.md"
LEDGER=ROOT/"provenance/reviewer_defenses/EVIDENCE_STATUS_20261031.md"
GEB_PORT=RESULTS/"geb_port_island_realization_20261031/GEB_to_port_species_bridge.json"
GATEWAY=RESULTS/"geb_port_island_realization_20261031/GEB_port_island_realization_risk_20261031.json"

def test_positive_observational_island_port_claim_preserved_without_shipping_causality():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory harbor/island synthesis excluded from standalone GEB blinded paper")
    report=REPORT.read_text(encoding="utf8")
    ledger=LEDGER.read_text(encoding="utf8")
    required_in_report=[
        "## 現在の生態学的結論：港・島は寄主再分配",
        "54.9%", "155種", "4.12%", "3.00%",
        "51島", "1,275", "ρ=+0.300", "+0.215", "ρ=+0.652",
        "港・島は、人為的に再配置された寄主植物",
        "因果経路ではない", "ρ=+0.418",
    ]
    for phrase in required_in_report:
        assert phrase in report, f"Source-grounded observable result disappeared: {phrase}"
    required_in_ledger=[
        "## 採用する観測的結論", "+54.9%", "4.12%",
        "ρ=+0.300", "+0.215", "ρ=+0.652",
        "地理的接点", "因果機構の実証ではない",
    ]
    for phrase in required_in_ledger:
        assert phrase in ledger, f"Central evidence ladder regressed: {phrase}"
    bridge=json.loads(GEB_PORT.read_text(encoding="utf8"))
    assert bridge["GEB_resource_opportunity"]["expanded_species"]==206
    assert bridge["GEB_resource_opportunity"]["butterfly_species"]==239
    assert abs(bridge["aligned_port_candidate_comparison"][
        "spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"]-0.6524282493600642)<1e-12
    assert bridge["primary_limits"]["port_index_does_not_establish_transport_or_larval_use"]
    risk=json.loads(GATEWAY.read_text(encoding="utf8"))
    assert risk["periods"]["2011_2017_first_index_among_no_pre2011_butterfly_index"][
        "geographic_margin_null"]["observed_matched_events"]==34
    assert risk["study_limits"]["no_proof_of_positive_or_negative_butterfly_colonization_effect"]
