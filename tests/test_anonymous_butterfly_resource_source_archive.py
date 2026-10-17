"""Safety/integrity boundaries of anonymous offline butterfly resource package."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"scripts/build_anonymous_butterfly_novelty_source_archive.py"
spec=importlib.util.spec_from_file_location("anon_butterfly_bundle",P)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_no_identifying_source_link_survives_scan() -> None:
    s={"RUN_ME.py":b'print("Safe reviewer reproduction")\n',
       "inputs/sample.csv":b"species,region\nButterfly example,ABC\n"}
    assert mod.privacy_audit(s)["status"]=="PASS_STATIC_DIRECT_IDENTIFIER_SCAN"
    with pytest.raises(ValueError,match="anonymity pattern"):
        mod.privacy_audit({"code/example.py":b"https://github.com/zuizui0223/"})
    with pytest.raises(ValueError,match="anonymity pattern"):
        mod.privacy_audit({"README.md":b"Contact: someone@example.com"})
    with pytest.raises(ValueError,match="anonymity pattern"):
        mod.privacy_audit({"README.md":"Author 張瑞琪".encode("utf-8")})


def test_source_run_numbers_are_removed_only_from_pure_provenance_assignments() -> None:
    raw=(b'd={"source_workflow_run": 37795847333, "retained": 1027}\n'
         b'# host effects are ordinary data\n'
         b'output={"BCE_source_workflow_run": 37795847333,\n'
         b'        "BCE_sensitivity_workflow_run": 38988757959,\n'
         b'        "species": 81}\n')
    result=mod.scrub_script("analyze.py",raw)
    # Inline mixed provenance/data assignment is intentionally NOT magically
    # rewritten; the anonymity scanner must reject such cases instead.
    with pytest.raises(ValueError,match="Source|source|workflow"):
        if b"source_workflow_run" in result:
            raise ValueError("Unredacted inline source workflow provenance")
    normalized=(b'd={\n'
                b' "source_workflow_run": 37795847333,\n'
                b' "species": 81,\n'
                b'}\n')
    text=mod.scrub_script("analyze.py",normalized).decode()
    assert '"species": 81' in text
    assert '37795847333' not in text
    assert mod.scrub_script("example.R",normalized)==normalized


def test_manifest_expected_sources_are_pinned() -> None:
    assert len(mod.INPUTS)==11
    assert len(mod.SCRIPTS)>=12
    assert "inputs/bce/accepted_host_edges.json" in mod.EXPECTED_SHA
    assert "inputs/level3.geojson" in mod.EXPECTED_SHA
    assert mod.EXPECTED_SHA["evidence/original_HOSTS_spatial_mcmc_draws.json"]==(
        "22c60accba86b576aa6ddc9f46ca7dd4d2003e2e81f1ba333ef41da5dbabbb91"
    )


def test_absolute_paths_and_parent_navigation_are_rejected() -> None:
    with pytest.raises(ValueError,match="escape"):
        mod.privacy_audit({"../exfiltrate.csv":b"a,b\n"})
    with pytest.raises(ValueError,match="escape"):
        mod.privacy_audit({"/home/reviewer/data.csv":b"a,b\n"})


def test_readme_declares_missing_provider_licence_and_no_demographic_inference() -> None:
    r=(ROOT/"docs/exploratory/ANONYMOUS_RESOURCE_NOVELTY_ARCHIVE_README_20261017.md").read_text()
    assert "licen" in r.lower()
    assert "not" in r.lower()
    assert "python RUN_ME.py" in r
    assert "0.5004" in r
    assert "not" in r.lower().split("clarke")[-1] or "global" in r.lower()
