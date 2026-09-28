from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _json(path: str):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _words(text: str) -> int:
    return len(re.findall(r"\b[\w–'-]+\b", text))


def test_butterfly_specialization_manuscript_bundle_is_internally_consistent():
    manuscript_path = ROOT / "manuscript/butterfly_specialization_ecology_v0.1.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.1.json")
    readiness = _json(
        "manuscript/butterfly_specialization_submission_readiness_v0.1.json"
    )
    figures = _json("manuscript/butterfly_specialization_figures_v0.1.json")

    dimensionality = _json(
        "benchmarks/exploratory/butterfly_specialization_dimensionality_result_v0.1.json"
    )
    anthropogenic = _json(
        "benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json"
    )
    mechanism = _json(
        "benchmarks/exploratory/butterfly_resource_expansion_mechanism_result_v0.1.json"
    )
    hierarchy = _json(
        "benchmarks/exploratory/butterfly_host_specialization_hierarchy_result_v0.1.json"
    )
    climate = _json(
        "benchmarks/exploratory/butterfly_climate_release_postgate_independent_result_v0.1.json"
    )

    assert claim_map["status"] == "WORKING_SUBMISSION_CLAIM_MAP"
    assert readiness["target_journal"]["primary"] == "Global Ecology and Biogeography"
    assert readiness["target_journal"]["article_type"] == "Research Article"
    assert figures["status"] == "RENDERED_AND_VISUALLY_AUDITED"
    assert figures["workflow"]["run_id"] == 36292128188
    assert len(figures["figures"]) == 5

    # Structured GEB abstract should remain below 300 words.
    abstract_start = manuscript.index("## Abstract")
    abstract_end = manuscript.index("---", abstract_start)
    abstract = manuscript[abstract_start:abstract_end]
    abstract_words = _words(abstract)
    assert abstract_words <= 300

    # Approximate main-text budget remains under the GEB typical 5,000 words.
    intro_start = manuscript.index("## 1. Introduction")
    legend_start = manuscript.index("## Figure legends")
    main_text_words = _words(manuscript[intro_start:legend_start])
    assert main_text_words <= 5000

    # Readiness metadata must stay synchronized with the actual manuscript.
    assert readiness["manuscript"]["approx_total_words"] == _words(manuscript)
    assert readiness["manuscript"]["approx_abstract_words"] == abstract_words
    assert readiness["manuscript"]["approx_main_text_words"] == main_text_words

    # GEB reference list is alphabetical by first author / corporate author.
    refs_start = manuscript.index("## References (working)")
    refs_end = manuscript.index("## Repository provenance", refs_start)
    reference_lines = [
        line[2:]
        for line in manuscript[refs_start:refs_end].splitlines()
        if line.startswith("- ")
    ]
    reference_keys = [
        re.split(r",|\. \d{4}\.", line, maxsplit=1)[0].strip()
        for line in reference_lines
    ]
    assert reference_keys == sorted(reference_keys, key=str.casefold)

    # Submission sections expected in the working bundle.
    for heading in (
        "## 1. Introduction",
        "## 2. Methods",
        "## 3. Results",
        "## 4. Discussion",
        "## 5. Limitations",
        "## 6. Conclusions",
        "## Figure legends",
        "## Data and Code Availability",
        "## References (working)",
    ):
        assert heading in manuscript

    # Core source values used in the paper remain identical to frozen receipts.
    assert (
        dimensionality["host_taxonomy_lower_bound_adequate_subset"][
            "spearman_host_family_vs_geographic_resource_breadth"
        ]
        == 0.398238154462201
    )

    full = anthropogenic["full_resource_eligible_panel"]
    assert full["species"] == 239
    assert full["species_expanded_by_introduced_host_ranges"] == 206
    assert full["total_native_species_units"] == 26530
    assert full["total_contemporary_species_units"] == 41083
    assert full["total_introduced_added_species_units"] == 14553
    assert (
        full["spearman"]["host_family_count_vs_log_resource_expansion"]
        == 0.00798541214928876
    )

    adequate_mech = mechanism["host_taxonomy_lower_bound_adequate_subset"]
    assert (
        adequate_mech["spearman"][
            "host_family_count_vs_effective_contributor_number"
        ]
        == 0.4919670615030713
    )
    assert (
        adequate_mech["spearman"][
            "host_family_count_vs_maximum_single_host_fractional_share"
        ]
        == -0.4861161332504768
    )
    assert (
        adequate_mech["host_breadth_strata"]["1_family"][
            "median_maximum_single_host_fractional_share"
        ]
        == 0.7465747904577693
    )
    assert (
        adequate_mech["host_breadth_strata"]["6plus_families"][
            "median_maximum_single_host_fractional_share"
        ]
        == 0.313289241622575
    )

    one_family = hierarchy["primary_focus"]
    assert one_family["species"] == 82
    assert (
        one_family["spearman"][
            "resolved_host_species_vs_effective_contributor_number"
        ]
        == 0.7338967715221224
    )
    assert (
        one_family["spearman"][
            "resolved_host_species_vs_maximum_single_host_fractional_share"
        ]
        == -0.7065837667974557
    )

    assert climate["climate_crossfit"]["climate_informative_species"] == 24
    assert climate["climate_crossfit"]["species_with_score_above_neutral_0_5"] == 23
    assert climate["primary_test"]["observed_partial_spearman"] == -0.16643611123731772
    assert climate["primary_test"]["one_sided_p_value"] == 0.2237
    assert (
        climate["primary_test"]["decision"]
        == "PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED"
    )

    # Abstract should lead with biological magnitude and architecture, not only correlations.
    for literal in (
        "206/239",
        "26,530",
        "41,083",
        "54.9% increase",
        "74.7%",
        "31.3%",
        "rho = 0.008",
        "median score = 0.801",
        "partial rho = -0.166",
        "p = 0.2237",
    ):
        assert literal in abstract
    assert abstract.index("206/239") < abstract.index("rho = 0.276")

    # Novelty framing should foreground the anthropogenic resource-portfolio result.
    title_line = manuscript.splitlines()[0]
    assert title_line == (
        "# Anthropogenic host redistribution expands butterfly resource geography "
        "through contrasting host portfolios in specialists and generalists"
    )
    assert "climate" not in title_line.lower()
    assert (
        "**Aim:** To determine how anthropogenic host redistribution reshapes "
        "butterfly resource geography"
    ) in abstract
    assert (
        "**Main conclusions:** Anthropogenic host redistribution expands butterfly "
        "resource geography without a proportional advantage"
    ) in abstract
    assert "A complementary resource-side problem remains unresolved" in manuscript

    discussion_start = manuscript.index("## 4. Discussion")
    discussion_end = manuscript.index("### 4.1", discussion_start)
    discussion_lead = manuscript[discussion_start:discussion_end]
    for literal in ("206/239", "54.9%", "74.7%", "31.3%"):
        assert literal in discussion_lead
    assert (
        "The key result is therefore not simply that introduced hosts can add opportunity"
        in discussion_lead
    )

    cover_letter = (
        ROOT / "manuscript/butterfly_specialization_geb_cover_letter_v0.1.md"
    ).read_text(encoding="utf-8")
    assert title_line[2:] in cover_letter
    assert "Our contribution is the complementary resource-side reconstruction" in cover_letter

    # Figure legends should expose sample size and ecological magnitude without requiring the text.
    for literal in (
        "aggregate species × region units increased from 26,530 to 41,083",
        "191 host-taxonomy-adequate butterflies",
        "74.7% to 31.3%",
        "82 host-taxonomy-adequate butterflies",
        "resolved host-species richness spans 1–37 species",
        "median filtering score was 0.801 and 23/24 species exceeded 0.5",
    ):
        assert literal in manuscript

    renderer = (
        ROOT / "scripts" / "render_butterfly_specialization_manuscript_figures.py"
    ).read_text(encoding="utf-8")
    for literal in (
        "Aggregate species × region units",
        "median effective contributors:",
        "median share:",
        "resolved host species =",
        "species > 0.5",
        "Contemporary resource breadth",
    ):
        assert literal in renderer

    # Detailed portfolio correlations remain available in the full manuscript.
    for literal in ("rho = 0.734", "rho = -0.707"):
        assert literal in manuscript

    # A non-supported association must not be rewritten as statistical independence.
    lowered = manuscript.lower()
    for overclaim in (
        "independently of family-level diet breadth",
        "largely separate filter",
        "climate strongly filters realized distributions",
        "climate is a strong filter",
        "rejects the simple interpretation",
        "climate helps determine which portions",
    ):
        assert overclaim not in lowered

    # Portfolio-concentration results must expose their structural upper-bound caveat.
    assert "host-contribution concentration metrics are structurally bounded by portfolio size" in lowered
    assert "associations between host richness and these concentration metrics" in lowered


def test_manuscript_claim_map_preserves_inference_boundaries():
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.1.json")
    claims = {row["id"]: row for row in claim_map["claims"]}

    assert claims["C1"]["status"] == "exploratory descriptive"
    assert claims["C2"]["status"] == "exploratory descriptive"
    assert claims["C3"]["status"] == "exploratory post-result mechanism analysis"
    assert claims["C4"]["status"] == "exploratory post-result follow-up"
    assert claims["C5"]["status"] == "independent cross-fit descriptive result"
    assert claims["C6"]["status"] == "independent hypothesis test"
    assert (
        claims["C6"]["decision"]
        == "PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED"
    )

    forbidden = "\n".join(claim_map["forbidden_overclaims"])
    assert "caused butterfly range expansion" in forbidden
    assert "true absence" in forbidden
    assert "confirmatory tests" in forbidden
