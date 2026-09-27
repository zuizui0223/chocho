from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_LEPTRAITS_SHA256 = (
    "6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b"
)


def test_frozen_leptraits_snapshot_is_present_and_hash_pinned() -> None:
    path = ROOT / "data/external/leptraits_consensus_v1.0.csv"
    data = path.read_bytes()
    assert len(data) > 2_000_000
    assert hashlib.sha256(data).hexdigest() == EXPECTED_LEPTRAITS_SHA256


def test_historical_figure_workflow_is_preserved_in_provenance() -> None:
    path = (
        ROOT
        / "provenance/workflows/butterfly-specialization-manuscript-figures-v01.yml"
    )
    text = path.read_text(encoding="utf-8")
    for run_id in ("36222306369", "36243549495", "36270981110", "36270581743"):
        assert f"run-id: {run_id}" in text


def test_data_code_availability_uses_ecology_repository_structure() -> None:
    text = (
        ROOT / "manuscript/butterfly_specialization_ecology_v0.1.md"
    ).read_text(encoding="utf-8")
    assert "versioned in the `chocho` ecology repository" in text
    assert (
        "provenance/workflows/"
        "butterfly-specialization-manuscript-figures-v01.yml"
    ) in text
    assert "versioned in the TTF repository" not in text
    assert ".github/workflows/butterfly-specialization-manuscript-figures-v01.yml" not in text


def test_r_sidecar_rebuild_uses_base_r_and_is_documented() -> None:
    scripts = [
        ROOT / "scripts/build_wcvp_hosts_sidecar.R",
        ROOT / "scripts/build_wcvp_hosts_contemporary_sidecar.R",
    ]
    for path in scripts:
        text = path.read_text(encoding="utf-8")
        lowered = text.lower()
        assert "library(" not in lowered
        assert "require(" not in lowered
        assert "install.packages" not in lowered

    doc = (ROOT / "docs/REPRODUCIBILITY.md").read_text(encoding="utf-8")
    assert 'r-version: "release"' in doc
    assert "base R only" in doc
    assert "65bed76bae9d644ccb6ad200c05f9f5071d89e05" in doc
    assert "808e0b869f9ec1adf8efff87cf6a395adda103e0" in doc


def test_shapely_is_explicit_in_clean_test_environment() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'test = ["pytest>=8", "shapely>=2,<3"]' in pyproject
    workflow = (ROOT / ".github/workflows/paper-ci.yml").read_text(encoding="utf-8")
    assert 'python -m pip install -e ".[test]"' in workflow
