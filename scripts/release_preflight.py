#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VERSION = "1.0.0"
EXPECTED_TAG = "v1.0.0-butterfly"
CURRENT_MANUSCRIPT = "manuscript/butterfly_specialization_ecology_v0.2.md"
CURRENT_TITLE_PAGE = "manuscript/butterfly_specialization_geb_title_page_template_v0.2.md"
CURRENT_BLINDED = "manuscript/generated/butterfly_specialization_ecology_blinded_v0.2.md"
CURRENT_CLAIM_MAP = "manuscript/butterfly_specialization_claim_map_v0.2.json"
EXPECTED_TITLE = (
    "Anthropogenic host redistribution expands butterfly resource geography "
    "across the specialization spectrum"
)


def cff_version(text: str) -> str | None:
    m = re.search(r'^version:\s*["\']?([^"\'\n]+)["\']?\s*$', text, re.MULTILINE)
    return m.group(1).strip() if m else None


def inspect_release(root: Path = ROOT) -> dict[str, object]:
    blockers: list[str] = []
    complete: list[str] = []

    pyproject = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    pkg_version = pyproject["project"]["version"]
    if pkg_version == EXPECTED_VERSION:
        complete.append(f"package version = {EXPECTED_VERSION}")
    else:
        blockers.append(f"package version is {pkg_version!r}; expected {EXPECTED_VERSION!r}")

    cff = (root / "CITATION.cff").read_text(encoding="utf-8")
    cv = cff_version(cff)
    if cv == EXPECTED_VERSION:
        complete.append(f"CITATION.cff version = {EXPECTED_VERSION}")
    else:
        blockers.append(f"CITATION.cff version is {cv!r}; expected {EXPECTED_VERSION!r}")

    if (root / "LICENSE").exists():
        complete.append("LICENSE exists")
    else:
        blockers.append("LICENSE is not chosen/added")

    manuscript = (root / CURRENT_MANUSCRIPT).read_text(encoding="utf-8")
    title_line = manuscript.splitlines()[0].removeprefix("# ").strip()
    if title_line == EXPECTED_TITLE:
        complete.append("current v0.2 manuscript title is synchronized")
    else:
        blockers.append(
            f"current manuscript title is {title_line!r}; expected {EXPECTED_TITLE!r}"
        )

    abstract_start = manuscript.index("## Abstract")
    abstract_end = manuscript.index("---", abstract_start)
    abstract_words = len(
        re.findall(r"\b[\w–'-]+\b", manuscript[abstract_start:abstract_end])
    )
    if abstract_words <= 300:
        complete.append(f"structured abstract is within 300-word limit ({abstract_words})")
    else:
        blockers.append(f"structured abstract exceeds 300 words ({abstract_words})")

    if (root / CURRENT_CLAIM_MAP).exists():
        complete.append("current v0.2 claim map exists")
    else:
        blockers.append("current v0.2 claim map is missing")

    title = (root / CURRENT_TITLE_PAGE).read_text(encoding="utf-8")
    author_tokens = (
        "[Author 1 full name]",
        "[Author 2 full name]",
        "[Author 3 full name]",
        "[Department / Graduate School]",
        "[Institution]",
        "[Full name]",
        "[Email]",
        "[ORCID]",
    )
    unresolved_authors = [token for token in author_tokens if token in title]
    if unresolved_authors:
        blockers.append("title-page author/affiliation/contact placeholders remain")
    else:
        complete.append("title-page authors/affiliations/contact metadata filled")

    if "- Conceptualization: [ ]" in title or "- Writing — original draft: [ ]" in title:
        blockers.append("CRediT contribution placeholders remain")
    else:
        complete.append("CRediT contributions filled")

    if "[Insert acknowledgements" in title:
        blockers.append("acknowledgements placeholder remains")
    else:
        complete.append("acknowledgements finalized")

    if "[Insert funding sources" in title:
        blockers.append("funding placeholder remains")
    else:
        complete.append("funding statement finalized")

    if "[Insert the final disclosure" in title:
        blockers.append("conflict-of-interest placeholder remains")
    else:
        complete.append("conflict-of-interest statement finalized")

    has_cff_doi = bool(re.search(r"^doi:\s*\S+", cff, re.MULTILINE))
    if has_cff_doi and "[DOI PLACEHOLDER]" not in title:
        complete.append("public archive DOI recorded in CFF/title-page metadata")
    else:
        blockers.append("public archive DOI is not yet minted/recorded")

    blinded = (root / CURRENT_BLINDED).read_text(encoding="utf-8")
    if "[ANONYMIZED REVIEW LINK]" in blinded:
        blockers.append("obsolete anonymized-reviewer-link placeholder remains in blinded manuscript")
    else:
        complete.append("blinded manuscript is free of obsolete reviewer-link placeholders")

    manifest = (root / "provenance/RELEASE_MANIFEST.md").read_text(encoding="utf-8")
    if EXPECTED_TAG in manifest:
        complete.append(f"intended GitHub tag recorded: {EXPECTED_TAG}")
    else:
        blockers.append(f"release manifest does not record tag {EXPECTED_TAG}")

    return {
        "schema": "chocho_release_preflight_v1_1",
        "manuscript": CURRENT_MANUSCRIPT,
        "expected_version": EXPECTED_VERSION,
        "expected_tag": EXPECTED_TAG,
        "ready": not blockers,
        "complete": complete,
        "blockers": blockers,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    result = inspect_release()
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        print(f"Release preflight: {'READY' if result['ready'] else 'PENDING'}")
        for item in result["complete"]:
            print(f"  [x] {item}")
        for item in result["blockers"]:
            print(f"  [ ] {item}")
    return 0 if result["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
