import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs/exploratory/geb_occurrence_realm_protocol_v0.1.json"


def load():
    return json.loads(PROTOCOL.read_text(encoding="utf-8"))


def test_realm_generality_panel_and_transport_are_frozen():
    p = load()
    assert p["schema"] == "chocho_geb_occurrence_realm_protocol_v0.1"
    assert p["panel"]["species_count"] == 239
    assert p["panel"]["source_sha256"] == (
        "b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa"
    )
    assert p["gbif"]["period"] == "2010-2026"
    assert p["gbif"]["ordinal_windows_per_species"] == 6
    assert p["gbif"]["window_size"] == 300
    assert p["gbif"]["scientific_filter_change_from_existing_climate_route"] is False


def test_realm_evaluability_gate_is_fail_closed():
    p = load()
    g = p["evaluability_gate"]
    assert g["minimum_realm_informative_species"] == 150
    assert g["minimum_interpretable_realms"] == 4
    assert g["interpretable_realm_minimum_n"] == 10
    assert g["no_backfill_after_result"] is True
    assert g["no_threshold_change_after_result"] is True


def test_realm_assignment_and_sampling_sensitivity_are_predeclared():
    p = load()
    assert p["species_realm_summary"]["minimum_mapped_occurrences"] == 30
    assert "primary_realm_share >= 0.60" in p["species_realm_summary"]["core_realm_sensitivity"]
    assert p["species_realm_summary"]["tie_handling"]["realm_replication"].startswith(
        "TIE is not a realm"
    )
    pooled = p["generality_tests"]["pooled_within_realm"]
    assert pooled["inference"]["permutations"] == 9999
    assert pooled["inference"]["p_value"].startswith("two-sided Monte Carlo")
    assert "1-degree-cell-weighted primary realm assignment" in pooled["sensitivities"]


def test_holt_realm_source_is_explicit_and_external():
    p = load()
    r = p["realm_source"]
    assert r["doi"] == "10.1126/science.1228282"
    assert r["realm_field"] == "Realm"
    assert len(r["expected_realms"]) == 11
    assert "butterfly-derived" in p["claim_boundary"][2]
