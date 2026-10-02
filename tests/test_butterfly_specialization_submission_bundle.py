from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "render_butterfly_specialization_blinded_manuscript.py"
SPEC = importlib.util.spec_from_file_location(
    "render_blinded_butterfly_manuscript", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_generated_blinded_v02_manuscript_matches_renderer():
    source = (
        ROOT / "manuscript" / "butterfly_specialization_ecology_v0.2.md"
    ).read_text(encoding="utf-8")
    generated = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.2.md"
    ).read_text(encoding="utf-8")
    assert generated == mod.render_blinded(source)


def test_blinded_v02_manuscript_has_geb_front_matter_and_limits():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.2.md"
    ).read_text(encoding="utf-8")

    running = re.search(r"^\*\*Running title:\*\*\s*(.+)$", text, re.MULTILINE)
    assert running is not None
    assert running.group(1).strip() == "Host redistribution and resource gain"
    assert len(running.group(1).strip()) < 40

    abstract_start = text.index("## Abstract")
    intro_start = text.index("## 1. Introduction")
    abstract = text[abstract_start:intro_start]
    for heading in (
        "Aim",
        "Location",
        "Time period",
        "Major taxa studied",
        "Methods",
        "Results",
        "Main conclusions",
    ):
        assert f"**{heading}:**" in abstract

    abstract_words = len(re.findall(r"\b[\w–'-]+\b", abstract.split("---")[0]))
    assert abstract_words <= 300

    keywords = re.search(r"^\*\*Keywords:\*\*\s*(.+)$", abstract, re.MULTILINE)
    assert keywords is not None
    items = [item.strip() for item in keywords.group(1).split(",") if item.strip()]
    assert 6 <= len(items) <= 10
    assert items == sorted(items, key=str.casefold)


def test_blinded_v02_manuscript_removes_identity_and_internal_history():
    text = (
        ROOT
        / "manuscript"
        / "generated"
        / "butterfly_specialization_ecology_blinded_v0.2.md"
    ).read_text(encoding="utf-8").lower()
    for token in (
        "ttf repository",
        ".github/workflows/",
        "repository provenance",
        "zuizui0223",
        "original ttf program",
        "added after manuscript review",
        "review_response_map",
        "geb_v02_decision_memo",
    ):
        assert token not in text

    assert "[anonymized review link]" not in text
    assert "an anonymized supplementary review archive" in text
    assert "supplied with this submission" in text
    assert "stable public archival snapshot" in text


def test_v02_title_page_cover_letter_and_checklist_are_synchronized():
    expected_title = (
        "Anthropogenic host redistribution expands butterfly resource geography "
        "across the specialization spectrum"
    )
    title_page = (
        ROOT / "manuscript" / "butterfly_specialization_geb_title_page_template_v0.2.md"
    ).read_text(encoding="utf-8")
    cover = (
        ROOT / "manuscript" / "butterfly_specialization_geb_cover_letter_v0.2.md"
    ).read_text(encoding="utf-8")
    checklist = json.loads(
        (
            ROOT
            / "manuscript"
            / "butterfly_specialization_geb_submission_checklist_v0.2.json"
        ).read_text(encoding="utf-8")
    )

    assert expected_title in title_page
    assert expected_title in cover
    for required in (
        "Authors",
        "Affiliations",
        "Corresponding author",
        "Author contributions",
        "Acknowledgements",
        "Funding",
        "Conflict of interest",
        "Data and code availability",
        "Double-anonymous review note",
    ):
        assert required in title_page

    assert checklist["journal"] == "Global Ecology and Biogeography"
    assert checklist["author_guideline_snapshot"]["double_anonymous_review"] is True
    assert checklist["author_guideline_snapshot"]["separate_title_page_required"] is True
    assert checklist["manuscript"]["keyword_count"] == 7
    assert checklist["figures"]["main_count"] == 3
    assert checklist["figures"]["supplementary_count"] == 1
    assert checklist["figures"]["total_count"] == 4


def test_anonymous_bundle_excludes_internal_response_memos():
    source = (ROOT / "scripts" / "build_anonymous_review_bundle.py").read_text(
        encoding="utf-8"
    )
    assert 'Path("provenance/reviewer_defenses/GEB_V02_DECISION_MEMO.md")' in source
    assert 'Path("provenance/reviewer_defenses/REVIEW_RESPONSE_MAP.md")' in source


def test_cover_letter_interest_paragraph_and_title_page_contacts_are_submission_ready():
    cover = (
        ROOT / "manuscript" / "butterfly_specialization_geb_cover_letter_v0.2.md"
    ).read_text(encoding="utf-8")
    title_page = (
        ROOT / "manuscript" / "butterfly_specialization_geb_title_page_template_v0.2.md"
    ).read_text(encoding="utf-8")
    checklist = json.loads(
        (
            ROOT
            / "manuscript"
            / "butterfly_specialization_geb_submission_checklist_v0.2.json"
        ).read_text(encoding="utf-8")
    )

    paragraphs = [p.strip() for p in re.split(r"\n\s*\n", cover) if p.strip()]
    interest = next(p for p in paragraphs if p.startswith("This manuscript addresses"))
    interest_words = len(re.findall(r"\b[\w–'-]+\b", interest))
    assert interest_words < 250
    assert checklist["manuscript"]["cover_letter_interest_paragraph_words"] == interest_words
    assert "Submission draft" not in cover
    assert "Points to customize before submission" not in cover

    author_lines = [
        line for line in title_page.splitlines()
        if re.match(r"^[1-3]\. \[Author ", line)
    ]
    assert len(author_lines) == 3
    assert all("[Email]" in line and "[ORCID]" in line for line in author_lines)
    assert title_page.count("## Corresponding author") == 1


def test_initial_submission_uses_attached_anonymous_review_archive():
    title_page = (
        ROOT / "manuscript" / "butterfly_specialization_geb_title_page_template_v0.2.md"
    ).read_text(encoding="utf-8").lower()
    cover = (
        ROOT / "manuscript" / "butterfly_specialization_geb_cover_letter_v0.2.md"
    ).read_text(encoding="utf-8").lower()
    renderer = (
        ROOT / "scripts" / "render_butterfly_specialization_blinded_manuscript.py"
    ).read_text(encoding="utf-8").lower()
    bundle = (ROOT / "scripts" / "build_anonymous_review_bundle.py").read_text(
        encoding="utf-8"
    ).lower()

    assert "supplementary review material" in title_page
    assert "before publication" in title_page
    assert "supplementary review material" in cover
    assert "[anonymized review link]" not in renderer
    assert "anonymized supplementary review archive" in renderer
    assert "uploaded directly with the manuscript" in bundle
