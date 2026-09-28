from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_repository_layout_uses_independent_ecology_namespace() -> None:
    assert not (ROOT / "src" / "ttf").exists()
    package = ROOT / "src" / "butterfly_specialization_ecology"
    assert package.is_dir()
    for name in (
        "__init__.py",
        "butterfly_climate_release.py",
        "butterfly_resource_envelope.py",
        "checkpointed_gbif_occurrence.py",
        "lepidoptera_host_resource.py",
        "resource_envelope_climate.py",
    ):
        assert (package / name).is_file()


def test_migration_history_is_not_at_repository_root() -> None:
    assert not (ROOT / "MIGRATION_PROVENANCE.md").exists()
    assert not (ROOT / "SOURCE_SNAPSHOT.txt").exists()
    assert (ROOT / "provenance" / "migration" / "MIGRATION_PROVENANCE.md").is_file()
    assert (ROOT / "provenance" / "migration" / "SOURCE_SNAPSHOT.txt").is_file()


def test_frozen_s1_manifest_is_a_data_input() -> None:
    assert (ROOT / "data" / "frozen" / "s1_species_manifest.json").is_file()
    assert not (
        ROOT
        / "benchmarks"
        / "frozen"
        / "relational_prior_S1_species_exclusion_v0.1.json"
    ).exists()


def test_navigation_indexes_exist() -> None:
    for rel in (
        "manuscript/README.md",
        "benchmarks/README.md",
        "benchmarks/exploratory/README.md",
        "docs/README.md",
        "docs/exploratory/README.md",
        "scripts/README.md",
        "tests/README.md",
        "data/README.md",
        "provenance/README.md",
        "provenance/climate/README.md",
        "provenance/archive/climate/README.md",
        "provenance/archive/climate/operations/README.md",
        "provenance/archive/climate/pilot/README.md",
    ):
        assert (ROOT / rel).is_file()


def test_superseded_climate_contracts_are_archived_not_active() -> None:
    active_protocols = ROOT / "docs" / "exploratory"
    active_receipts = ROOT / "benchmarks" / "exploratory"
    archive = ROOT / "provenance" / "archive" / "climate"

    for name in (
        "butterfly_climate_release_independent_test_v0.1.json",
        "butterfly_climate_release_independent_test_v0.2.json",
    ):
        assert not (active_protocols / name).exists()
        assert (archive / "protocols" / name).is_file()

    for name in (
        "butterfly_climate_release_execution_binding_v0.1.json",
        "butterfly_climate_release_execution_binding_v0.2.json",
    ):
        assert not (active_receipts / name).exists()
        assert (archive / "bindings" / name).is_file()

    assert (
        active_protocols
        / "butterfly_climate_release_independent_test_v0.2.1.json"
    ).is_file()
    assert not (
        active_receipts
        / "butterfly_climate_release_execution_binding_v0.2.1.json"
    ).exists()
    assert (
        ROOT
        / "provenance"
        / "climate"
        / "butterfly_climate_release_execution_binding_v0.2.1.json"
    ).is_file()


def test_audit_and_resource_history_are_archived_outside_active_benchmarks() -> None:
    active = ROOT / "benchmarks" / "exploratory"
    climate_archive = ROOT / "provenance" / "archive" / "climate" / "audits"
    resource_archive = ROOT / "provenance" / "archive" / "resource"

    climate_audits = (
        "butterfly_climate_release_panel_balance_audit_v0.1.json",
        "butterfly_climate_release_contemporary_resource_balance_v0.1.json",
        "butterfly_climate_release_effect_score_audit_v0.1.json",
        "butterfly_climate_release_transport_blocker_audit_v0.1.json",
    )
    for name in climate_audits:
        assert not (active / name).exists()
        assert (climate_archive / name).is_file()

    resource_history = (
        "butterfly_resource_envelope_pilot_reconstruction_v0.1.json",
        "butterfly_host_breadth_geography_result_v0.1.json",
    )
    for name in resource_history:
        assert not (active / name).exists()
        assert (resource_archive / name).is_file()

    assert (resource_archive / "README.md").is_file()


def test_active_gbif_acquisition_uses_current_repository_identity() -> None:
    text = (
        ROOT / "scripts" / "acquire_butterfly_resource_envelope_occurrences.py"
    ).read_text(encoding="utf-8")
    assert "chocho-butterfly-resource-envelope/1.0" in text
    assert "github.com/zuizui0223/chocho" in text
    assert "github.com/zuizui0223/TTF" not in text


def test_climate_execution_history_is_outside_active_protocols_and_results() -> None:
    active_protocols = ROOT / "docs" / "exploratory"
    active_results = ROOT / "benchmarks" / "exploratory"
    operations = ROOT / "provenance" / "archive" / "climate" / "operations"
    pilot = ROOT / "provenance" / "archive" / "climate" / "pilot"
    current = ROOT / "provenance" / "climate"

    for name in (
        "butterfly_climate_release_postgate_transport_completion_v0.1.json",
        "butterfly_climate_release_transport_completion_rule_v0.1.json",
        "butterfly_climate_release_transport_recovery_v0.1.json",
    ):
        assert not (active_protocols / name).exists()
        assert (operations / name).is_file()

    pilot_gate = "butterfly_resource_envelope_climate_pilot_gate_v0.1.json"
    assert not (active_protocols / pilot_gate).exists()
    assert (pilot / pilot_gate).is_file()

    for name in (
        "butterfly_climate_release_independent_panel_v0.1.json",
        "butterfly_climate_release_execution_binding_v0.2.1.json",
        "butterfly_climate_release_preclimate_result_v0.2.1.json",
    ):
        assert not (active_results / name).exists()
        assert (current / name).is_file()

    pilot_result = "butterfly_resource_envelope_climate_pilot_result_v0.1.json"
    assert not (active_results / pilot_result).exists()
    assert (pilot / pilot_result).is_file()


def test_active_protocol_and_result_surfaces_are_whitelisted() -> None:
    protocol_names = {
        p.name
        for p in (ROOT / "docs" / "exploratory").iterdir()
        if p.is_file()
    }
    assert protocol_names == {
        "README.md",
        "BUTTERFLY_SPECIALIZATION_ECOLOGY_SYNTHESIS_V0_1.md",
        "butterfly_anthropogenic_resource_expansion_protocol_v0.1.json",
        "butterfly_climate_release_independent_test_v0.2.1.json",
        "butterfly_host_specialization_hierarchy_protocol_v0.1.json",
        "butterfly_resource_envelope_protocol_v0.1.json",
        "butterfly_resource_expansion_mechanism_protocol_v0.1.json",
    }

    result_names = {
        p.name
        for p in (ROOT / "benchmarks" / "exploratory").iterdir()
        if p.is_file()
    }
    assert result_names == {
        "README.md",
        "butterfly_anthropogenic_resource_expansion_result_v0.1.json",
        "butterfly_climate_release_postgate_independent_result_v0.1.json",
        "butterfly_host_specialization_hierarchy_result_v0.1.json",
        "butterfly_resource_expansion_mechanism_result_v0.1.json",
        "butterfly_specialization_dimensionality_result_v0.1.json",
        "butterfly_specialization_ecology_synthesis_v0.1.json",
    }
