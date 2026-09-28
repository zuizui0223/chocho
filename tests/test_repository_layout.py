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
        "provenance/archive/climate/README.md",
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
    assert (
        active_receipts
        / "butterfly_climate_release_execution_binding_v0.2.1.json"
    ).is_file()
