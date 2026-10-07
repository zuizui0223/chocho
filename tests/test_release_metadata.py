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
        "Plant globalization expands and homogenizes "
        "butterfly larval-resource geography"
    )
    assert current_title == expected
    match = re.search(
        r'EXPECTED_TITLE\s*=\s*\(\s*"([^"]*)"\s*"([^"]*)"\s*\)',
        preflight,
        re.DOTALL,
    )
    assert match is not None
    assert "".join(match.groups()) == expected
    assert "largely independently of diet breadth" not in preflight


def test_release_preflight_does_not_require_external_review_url() -> None:
    preflight = (ROOT / "scripts" / "release_preflight.py").read_text(
        encoding="utf-8"
    ).lower()
    assert "anonymized reviewer-access url is not yet inserted" not in preflight
    assert "obsolete anonymized-reviewer-link placeholder remains" in preflight

    release_doc = (ROOT / "docs" / "RELEASE_PROCEDURE.md").read_text(
        encoding="utf-8"
    ).lower()
    assert "no external anonymous reviewer url is required" in release_doc
    assert "uploaded directly as supplementary review material" in release_doc
