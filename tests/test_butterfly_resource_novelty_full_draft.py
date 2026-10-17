"""Protect the separate full butterfly-resource-novelty manuscript against claim drift.

A textual and numeric audit is not a peer review or a test of causality;
it ensures the edited manuscript is consistent with frozen source receipts.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs/exploratory/BUTTERFLY_RESOURCE_NOVELTY_FULL_DRAFT_20261017.md"
RESULTS = ROOT / "provenance/reviewer_defenses/results"


def read_result(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def manuscript() -> str:
    return MANUSCRIPT.read_text(encoding="utf-8")


def test_journal_article_scope_and_sections() -> None:
    text = manuscript()
    assert text.startswith("# When is an introduced host a new resource?")
    assert "PR #38" in text
    for heading in (
        "## Structured abstract",
        "## 1. Introduction",
        "## 2. Methods",
        "## 3. Results",
        "## 4. Discussion",
        "## 5. Limitations",
        "## 6. Conclusions",
        "## References",
        "## Data and Code Availability",
        "## Candidate Figure Legends",
    ):
        assert text.count(heading) == 1, heading
    short = re.search(r"\*\*Running title:\*\*\s*(.+)", text)
    assert short and len(short.group(1)) < 40


def test_geb_structured_abstract_word_count() -> None:
    text = manuscript()
    abstract = text.split("## Structured abstract", 1)[1].split("**Keywords:**", 1)[0]
    headings = [
        "Aim", "Location", "Time period", "Major taxa studied",
        "Methods", "Results", "Main conclusions",
    ]
    seen = [re.search(r"\*\*" + re.escape(x) + r":\*\*", abstract).start() for x in headings]
    assert seen == sorted(seen)
    words = len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b", abstract))
    assert words <= 300, words
    keys = text.split("**Keywords:**", 1)[1].split("\n", 1)[0]
    items = [x.strip() for x in keys.split(",")]
    assert 6 <= len(items) <= 10
    assert items == sorted(items, key=str.casefold)


def test_full_manuscript_main_body_and_reference_ceiling() -> None:
    text = manuscript()
    body = text.split("## 1. Introduction", 1)[1].split("## References", 1)[0]
    # GEB estimates approximately 5000 main-text words, excluding
    # abstract, references and captions.
    count = len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b", body))
    assert count <= 5000, count
    references = text.split("## References", 1)[1].split("## Data and Code Availability", 1)[0]
    ref_lines = [x for x in references.splitlines() if x.startswith("- ")]
    assert 12 <= len(ref_lines) <= 50
    assert "https://doi.org/10.1002/ece3.10834" in references
    assert "https://doi.org/10.1038/s41586-020-2640-y" in references


def test_stated_source_budget_is_exact_and_denominators_explicit() -> None:
    text = manuscript()
    b = read_result("butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json")
    orig = b["original_HOSTS"]
    aug = b["BCE_Clarke_partial_augmented"]
    for src in (orig, aug):
        H = src["candidate_incidences"]
        A = src["native_absorbed_incidences"]
        R = src["redundant_new_cell_incidences"]
        U = src["distinct_added_cells"]
        assert H == A + R + U
        assert f"{H:,}" in text
        assert f"{A:,}" in text
        assert f"{U:,}" in text
    assert "61.6%" in text and "65.5%" in text
    assert "19.3%" in text and "14.5%" in text
    assert "incidence" in text
    assert "not probabilities" in text or "not probabilities of" in text


def test_turnover_grouping_does_not_hide_158_unchanged_species() -> None:
    text = manuscript()
    d = read_result("butterfly_bce_added_region_identity_turnover_v0.1.json")
    info = d["affected_81"]
    assert info["original_new_resource_cells"] == 4628
    assert info["augmented_new_resource_cells"] == 5042
    assert info["retained_in_both"] == 3225
    assert info["reclassified_as_native"] == 1403
    assert info["newly_added_by_augmentation"] == 1817
    assert info["union_of_either_source"] == 6445
    assert 3225 / 6445 == info["intersection_over_union_jaccard"]
    assert 158 == d["full_239"]["unchanged_host_inventory_species"]
    assert info["species_increase_net"] == 36
    assert info["species_decrease_net"] == 35
    assert info["species_equal_net"] == 10
    assert all(v in text for v in (
        "4,628", "5,042", "3,225", "1,403", "1,817", "6,445",
        "0.500", "81", "158",
        "Lycaena dispar", "Minois dryas", "Thymelicus lineola",
    ))


def test_spatial_null_qualified_but_never_promoted_to_source_robustness() -> None:
    text = manuscript()
    n = read_result("butterfly_spatial_null_longchain_qualified_v0.2.json")
    assert n["status"] == "QUALIFIED_EXPLORATORY_MCMC_SENSITIVITY"
    assert n["qualification"]["status"] == n["status"]
    assert len(n["chains"]) == 2
    assert all(x["pass"] for x in n["qualification"]["checks"])
    assert [x["observed_absorption_minus_null"] for x in n["chains"]] == [271, 276]
    assert all(k in text for k in ("1,500 km", "271", "276", "0.36", "post-hoc"))
    assert "not" in text.lower()
    assert "not" in text.lower().split("BCE-augmented", 1)[-1]
    assert "does not" in text.lower()


def test_original_submitted_manuscript_remains_separately_named() -> None:
    original = ROOT / "manuscript/butterfly_specialization_ecology_v0.2.md"
    assert original.is_file()
    assert original.resolve() != MANUSCRIPT.resolve()
    assert "Plant globalization expands and homogenizes" in original.read_text(encoding="utf-8")
