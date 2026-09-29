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


def test_data_code_availability_uses_current_v02_submission_surface() -> None:
    text = (
        ROOT / "manuscript/butterfly_specialization_ecology_v0.2.md"
    ).read_text(encoding="utf-8")
    data_start = text.index("## Data and Code Availability")
    figure_start = text.index("## Figure legends")
    data_section = text[data_start:figure_start]
    assert "Analysis code and the inputs required to reproduce" in data_section
    assert "versioned in the study repository" in data_section
    assert "TTF repository" not in data_section
    assert ".github/workflows/" not in data_section


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


FROZEN_FIGURE_SOURCE_SHA256 = {
    "s1_resource_descriptors.csv": "894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc",
    "anthropogenic_species_metrics.csv": "7b2a20387d656dbfcd7b3c38ec6a7fe2505474ad51f783dd1f201eb6b7eabd30",
    "host_contribution_metrics.csv": "b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa",
    "independent_climate_primary_result.json": "a73dca6e8b669f721d8f2745f27198d9b47a5cbd05dde6146cc8d4f1ddfbf79b",
}


def test_frozen_figure_source_inputs_are_vendored_byte_exact() -> None:
    root = ROOT / "data/frozen/figure_sources"
    for name, expected in FROZEN_FIGURE_SOURCE_SHA256.items():
        data = (root / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == expected
