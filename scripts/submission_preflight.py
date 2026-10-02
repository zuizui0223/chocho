#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = "manuscript/butterfly_specialization_ecology_v0.2.md"
BLINDED = "manuscript/generated/butterfly_specialization_ecology_blinded_v0.2.md"
TITLE_PAGE = "manuscript/butterfly_specialization_geb_title_page_template_v0.2.md"
COVER_LETTER = "manuscript/butterfly_specialization_geb_cover_letter_v0.2.md"
SUPPLEMENT = "manuscript/butterfly_specialization_supplement_v0.2.md"
ANON_BUNDLE_WORKFLOW = ".github/workflows/build-anonymous-review-bundle.yml"
ANON_BUNDLE_BUILDER = "scripts/build_anonymous_review_bundle.py"

EXPECTED_TITLE = (
    "Anthropogenic host redistribution expands butterfly resource geography "
    "across the specialization spectrum"
)
EXPECTED_RUNNING_TITLE = "Host redistribution and resource gain"


def words(text: str) -> int:
    return len(re.findall(r"\b[\w–'-]+\b", text))


def inspect_initial_submission(root: Path = ROOT) -> dict[str, object]:
    blockers: list[str] = []
    complete: list[str] = []

    manuscript = (root / MANUSCRIPT).read_text(encoding="utf-8")
    title = manuscript.splitlines()[0].removeprefix("# ").strip()
    if title == EXPECTED_TITLE:
        complete.append("manuscript title synchronized")
    else:
        blockers.append(f"manuscript title mismatch: {title!r}")

    abstract_start = manuscript.index("## Abstract")
    abstract_end = manuscript.index("---", abstract_start)
    abstract_words = words(manuscript[abstract_start:abstract_end])
    if abstract_words <= 300:
        complete.append(f"structured abstract within 300 words ({abstract_words})")
    else:
        blockers.append(f"structured abstract exceeds 300 words ({abstract_words})")

    main_start = manuscript.index("## 1. Introduction")
    refs_start = manuscript.index("## References")
    main_words = words(manuscript[main_start:refs_start])
    if main_words <= 5000:
        complete.append(f"main text within approximately 5,000 words ({main_words})")
    else:
        blockers.append(f"main text exceeds approximately 5,000 words ({main_words})")

    blinded = (root / BLINDED).read_text(encoding="utf-8")
    running = re.search(r"^\*\*Running title:\*\*\s*(.+)$", blinded, re.MULTILINE)
    if running and running.group(1).strip() == EXPECTED_RUNNING_TITLE:
        complete.append("running title synchronized")
    else:
        blockers.append("running title is missing or unsynchronized")

    if len(EXPECTED_RUNNING_TITLE) < 40:
        complete.append(f"running title under 40 characters ({len(EXPECTED_RUNNING_TITLE)})")
    else:
        blockers.append("running title is 40 characters or longer")

    keywords = re.search(r"^\*\*Keywords:\*\*\s*(.+)$", blinded, re.MULTILINE)
    if keywords:
        items = [x.strip() for x in keywords.group(1).split(",") if x.strip()]
        if 6 <= len(items) <= 10:
            complete.append(f"keyword count within GEB range ({len(items)})")
        else:
            blockers.append(f"keyword count outside 6-10 ({len(items)})")
    else:
        blockers.append("keywords missing from blinded manuscript")

    if "[ANONYMIZED REVIEW LINK]" in blinded:
        blockers.append("obsolete anonymous-reviewer-link placeholder remains")
    elif "an anonymized supplementary review archive" in blinded.lower():
        complete.append("blinded data/code statement uses attached anonymous review archive")
    else:
        blockers.append("blinded data/code statement does not describe attached review archive")

    for forbidden in ("zuizui0223", "ruiqi", "zhang.ruiqi", "repository provenance"):
        if forbidden in blinded.lower():
            blockers.append(f"blinded manuscript contains identifying/internal token: {forbidden}")

    gbif_doi_re = re.compile(r"https://doi\.org/10\.15468/dl\.[A-Za-z0-9]+")
    gbif_placeholder = "GBIF_OCCURRENCE_DOWNLOAD_DOI_PLACEHOLDER"
    if gbif_placeholder in manuscript or gbif_placeholder in blinded:
        blockers.append("GBIF occurrence-download DOI placeholder remains")
    elif gbif_doi_re.search(manuscript) and gbif_doi_re.search(blinded):
        complete.append("GBIF occurrence-download DOI recorded in manuscript and blinded manuscript")
    else:
        blockers.append("GBIF occurrence-download DOI is missing from manuscript or blinded manuscript")

    title_page = (root / TITLE_PAGE).read_text(encoding="utf-8")
    title_page_tokens = (
        "[Author 1 full name]",
        "[Email]",
        "[ORCID]",
        "[Department / Graduate School]",
        "[Institution]",
        "[Full name]",
        "[Postal address]",
    )
    unresolved = sorted({x for x in title_page_tokens if x in title_page})
    if unresolved:
        blockers.append("title-page author/affiliation/contact placeholders remain")
    else:
        complete.append("title-page author/affiliation/contact metadata filled")

    if "- Conceptualization: [ ]" in title_page or "- Writing — original draft: [ ]" in title_page:
        blockers.append("CRediT contribution placeholders remain")
    else:
        complete.append("CRediT contributions filled")

    if "[Insert acknowledgements" in title_page:
        blockers.append("acknowledgements placeholder remains")
    else:
        complete.append("acknowledgements finalized")

    if "[Insert funding sources" in title_page:
        blockers.append("funding placeholder remains")
    else:
        complete.append("funding statement finalized")

    if "[Insert the final disclosure" in title_page:
        blockers.append("conflict-of-interest placeholder remains")
    else:
        complete.append("conflict-of-interest statement finalized")

    if "supplementary review material" in title_page.lower():
        complete.append("title-page data/code statement supports attached review archive")
    else:
        blockers.append("title-page data/code statement does not describe review archive")

    cover = (root / COVER_LETTER).read_text(encoding="utf-8")
    paragraphs = [p.strip() for p in re.split(r"\n\s*\n", cover) if p.strip()]
    interest = next((p for p in paragraphs if p.startswith("This manuscript addresses")), "")
    interest_words = words(interest)
    if interest and interest_words < 250:
        complete.append(f"cover-letter GEB-interest paragraph under 250 words ({interest_words})")
    else:
        blockers.append(f"cover-letter interest paragraph missing or too long ({interest_words})")

    if any(token in cover for token in ("[Corresponding author]", "[Affiliation]", "[Email]", "[ORCID]")):
        blockers.append("cover-letter corresponding-author placeholders remain")
    else:
        complete.append("cover-letter corresponding-author signature filled")

    for path in (SUPPLEMENT, ANON_BUNDLE_WORKFLOW, ANON_BUNDLE_BUILDER):
        if (root / path).exists():
            complete.append(f"submission asset exists: {path}")
        else:
            blockers.append(f"submission asset missing: {path}")

    # A public DOI, public repository release/tag and external anonymous reviewer URL
    # are intentionally NOT initial-submission blockers. GEB permits data/code access
    # during peer review via supplementary materials; stable public archiving is
    # required for publication.
    deferred_until_publication = [
        "choose/apply public archive licensing",
        "deposit stable public data/code archive",
        "record persistent public archive DOI",
        "tag public v1.0.0-butterfly release",
    ]

    return {
        "schema": "chocho_geb_initial_submission_preflight_v0.1",
        "journal": "Global Ecology and Biogeography",
        "ready_for_initial_submission": not blockers,
        "complete": complete,
        "blockers": blockers,
        "deferred_until_publication": deferred_until_publication,
        "counts": {
            "abstract_words": abstract_words,
            "main_text_words": main_words,
            "cover_letter_interest_words": interest_words,
        },
        "peer_review_data_code_mode": "attach anonymized review archive as supplementary review material",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    result = inspect_initial_submission()
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        state = "READY" if result["ready_for_initial_submission"] else "PENDING"
        print(f"GEB initial submission preflight: {state}")
        for item in result["complete"]:
            print(f"  [x] {item}")
        for item in result["blockers"]:
            print(f"  [ ] {item}")
        if result["deferred_until_publication"]:
            print("  Deferred until publication/archive:")
            for item in result["deferred_until_publication"]:
                print(f"    - {item}")
    return 0 if result["ready_for_initial_submission"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
