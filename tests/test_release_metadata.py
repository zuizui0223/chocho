from __future__ import annotations

import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VERSION = "1.0.0"
EXPECTED_TAG = "v1.0.0-butterfly"


def test_release_version_identifiers_are_consistent() -> None:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert pyproject["project"]["version"] == EXPECTED_VERSION

    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r'^version:\s*["\']?([^"\'\n]+)["\']?\s*$', cff, re.MULTILINE)
    assert match is not None
    assert match.group(1).strip() == EXPECTED_VERSION

    manifest = (ROOT / "provenance/RELEASE_MANIFEST.md").read_text(encoding="utf-8")
    assert f"version `{EXPECTED_VERSION}`" in manifest
    assert f"tag `{EXPECTED_TAG}`" in manifest


def test_release_preflight_uses_current_v02_title() -> None:
    manuscript = (
        ROOT / "manuscript" / "butterfly_specialization_ecology_v0.2.md"
    ).read_text(encoding="utf-8")
    current_title = manuscript.splitlines()[0].removeprefix("# ").strip()

    preflight = (ROOT / "scripts" / "release_preflight.py").read_text(
        encoding="utf-8"
    )
    expected = (
        "Anthropogenic host redistribution expands butterfly resource geography "
        "across the specialization spectrum"
    )
    assert current_title == expected
    assert expected in preflight
    assert "largely independently of diet breadth" not in preflight
