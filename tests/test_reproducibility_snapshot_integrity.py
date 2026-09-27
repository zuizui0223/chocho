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
