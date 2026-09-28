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
    ):
        assert (ROOT / rel).is_file()
