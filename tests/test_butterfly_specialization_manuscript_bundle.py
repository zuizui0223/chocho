from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _json(path: str):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _words(text: str) -> int:
    return len(re.findall(r"\b[\w–'-]+\b", text))


def test_current_v02_manuscript_bundle_is_internally_consistent():
    manuscript_path = ROOT / "manuscript/butterfly_specialization_ecology_v0.2.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.2.json")
    readiness = _json(
        "manuscript/butterfly_specialization_submission_readiness_v0.2.json"
    )
    concentration = _json(
        "provenance/reviewer_defenses/results/"
        "butterfly_host_contribution_concentration_v0.1.json"
    )
    equivalence = _json(
        "provenance/reviewer_defenses/results/"
        "butterfly_expansion_equivalence_v0.1.json"
    )
    occurrence = _json(
        "provenance/reviewer_defenses/results/"
        "butterfly_occurrence_species_robustness_v0.1.json"
    )
    pgls = _json(
        "provenance/reviewer_defenses/results/butterfly_kawahara_pgls_v0.2.json"
    )
    poaceae = _json(
        "provenance/reviewer_defenses/results/butterfly_poaceae_sensitivity_v0.1.json"
    )
    family_generality = _json(
        "provenance/reviewer_defenses/results/butterfly_taxonomic_family_expansion_v0.1.json"
    )

    expected_title = (
        "Plant globalization expands and homogenizes "
        "butterfly larval-resource geography"
    )
    assert manuscript.splitlines()[0] == f"# {expected_title}"
    assert claim_map["title"] == expected_title
    assert readiness["manuscript"]["working_title"] == expected_title

    # Submission-scale format checks must apply to the current v0.2 manuscript.
    abstract_start = manuscript.index("## Abstract")
    abstract_end = manuscript.index("---", abstract_start)
    abstract = manuscript[abstract_start:abstract_end]
    abstract_words = _words(abstract)
    assert abstract_words <= 300

    intro_start = manuscript.index("## 1. Introduction")
    refs_start = manuscript.index("## References")
    main_text_words = _words(manuscript[intro_start:refs_start])
    assert main_text_words <= 5000

    assert readiness["manuscript"]["approximate_abstract_words"] == abstract_words
    assert readiness["manuscript"]["approximate_main_text_words"] == main_text_words
    assert readiness["manuscript"]["keyword_count"] == 7

    # Current manuscript references remain alphabetical by first/corporate author.
    refs_start = manuscript.index("## References")
    refs_end = manuscript.index("## Data and Code Availability", refs_start)
    reference_lines = [
        line[2:]
        for line in manuscript[refs_start:refs_end].splitlines()
        if line.startswith("- ")
    ]
    reference_keys = [
        ("van Kleunen" if line.startswith("van Kleunen") else
         "GBIF" if line.startswith("GBIF.org") else
         re.split(r",|\. \d{4}\.", line, maxsplit=1)[0].strip())
        for line in reference_lines
    ]
    assert reference_keys == sorted(reference_keys, key=str.casefold)

    # The ecological contribution must remain explicit in the current manuscript.
    for literal in (
        "resource release and resource homogenization",
        "Plant globalization homogenizes regional resource assemblages",
        "potential interspecific resource-sharing exposure",
        "potential resource co-use, not realized competition",
        "Exact shared-host butterfly-pair × region units increased",
        "58.9%",
        "resource filtering",
        "interspecific resource competition",
        "coarse global stress test",
        "Nakadai et al. 2018",
        "Braga 2023",
        "Yoon & Read 2016",
    ):
        assert literal in manuscript

    # References remain alphabetized by first-author/corporate-author key.
    refs_start = manuscript.index("## References")
    refs_end = manuscript.index("## Data and Code Availability", refs_start)
    ref_lines = [
        line[2:]
        for line in manuscript[refs_start:refs_end].splitlines()
        if line.startswith("- ")
    ]
    ref_keys = [
        re.split(r",|\. \d{4}\.", line, maxsplit=1)[0].strip()
        for line in ref_lines
    ]
    assert ref_keys == sorted(ref_keys, key=str.casefold)

    # Core ecological magnitudes and the new positive plant-side result.
    assert concentration["reconstruction"]["expanded_butterflies"] == 206
    assert concentration["reconstruction"]["added_butterfly_x_wgsrpd3_units"] == 14553
    species = concentration["plant_species"]
    genera = concentration["genera"]
    assert species["contributors"] == 670
    assert species["contributors_for_50pct"] == 38
    assert species["top_k_share"]["10"] == 0.2510105136113853
    assert species["top_k_share"]["50"] == 0.5700407340691861
    assert species["top10_poaceae"] == 7
    assert genera["contributors"] == 431
    assert genera["contributors_for_50pct"] == 25
    assert genera["top_k_share"]["10"] == 0.28686383837607116
    assert genera["top_k_share"]["50"] == 0.6703956605963941

    # Near-zero point estimate is not allowed to become an exact-null claim.
    assert equivalence["n"] == 239
    assert equivalence["rho"] == 0.00798541214928876
    assert equivalence["bootstrap"]["replicates"] == 49999
    assert equivalence["bootstrap"]["ci95"] == [
        -0.1113648840403796,
        0.12777832073032636,
    ]
    assert equivalence["tost_fisher_z_approx"]["margin"] == 0.1
    assert equivalence["tost_fisher_z_approx"]["equivalent_at_alpha_0_05"] is False
    assert equivalence["tost_fisher_z_approx"]["p_tost"] == 0.0779926620144849

    # Phylogenetic coverage and model boundary are explicit.
    assert pgls["exact_tree_match_species"] == 124
    assert pgls["pagel_rank_pgls"]["lambda"] == 0.0327902289475728
    assert pgls["pagel_rank_pgls"]["coefficient"] == -0.05108744710537585
    assert pgls["brownian_rank_pgls"]["coefficient"] == 0.1401214490702778
    assert pgls["brownian_rank_pgls"]["ci95"][1] > 0.3

    # Grass-feeding guilds do not generate the near-zero diet-breadth slope.
    assert poaceae["poaceae_users"]["species"] == 58
    assert poaceae["excluding_any_poaceae_user"]["species"] == 181
    assert poaceae["excluding_any_poaceae_user"]["rho"] == 0.025860902287550575
    assert poaceae["family_level_poaceae_specialists"]["species"] == 34
    assert (
        poaceae["excluding_family_level_poaceae_specialists"]["rho"]
        == 0.01324328000717657
    )

    # Expansion is widespread across the five major butterfly families.
    major = {row["family"]: row for row in family_generality["major_families"]}
    assert major["Nymphalidae"]["expanded_species"] == 86
    assert major["Nymphalidae"]["species"] == 103
    assert major["Hesperiidae"]["expanded_species"] == 46
    assert major["Pieridae"]["expanded_species"] == 37
    assert major["Lycaenidae"]["expanded_species"] == 22
    assert major["Papilionidae"]["expanded_species"] == 14
    assert min(row["expanded_fraction"] for row in major.values()) > 0.83
    assert major["Hesperiidae"]["median_host_family_count"] == 1.0
    assert major["Pieridae"]["median_host_family_count"] == 1.0

    # Secondary occurrence validation remains species-robust.
    assert occurrence["species"] == 23
    assert occurrence["observed"]["recovered_units"] == 66
    assert occurrence["observed"]["outside_native_units"] == 115
    assert occurrence["leave_one_out"]["pyrgus_communis_excluded"]["recovered_units"] == 44
    assert occurrence["leave_one_out"]["pyrgus_communis_excluded"]["outside_units"] == 93

    # Current main Results order: expansion -> homogenization -> shared exposure -> occurrence.
    result_headings = [
        "### 3.1 Plant globalization expands butterfly resource geography through uneven host contributions",
        "### 3.2 Plant globalization homogenizes regional resource assemblages and butterfly resource geography",
        "### 3.3 Introduced hosts expand shared-resource exposure and create a low-redundancy resource network",
        "### 3.4 Added resource geography aligns with contemporary butterfly occurrence",
    ]
    positions = [manuscript.index(h) for h in result_headings]
    assert positions == sorted(positions)
    assert "### 3.5 Secondary climate analysis" not in manuscript
    assert "### 4.3 Resource opportunity is filtered before realization" not in manuscript

    supplement = (
        ROOT / "manuscript/butterfly_specialization_supplement_v0.2.md"
    ).read_text(encoding="utf-8")
    assert "## Supplementary Methods S1. Climate filtering within contemporary resource opportunity" in supplement
    assert "## Supplementary Table S5. Climate-distance sensitivity and effect-size precision" in supplement
    assert "## Supplementary Table S7. Crop-host exclusion sensitivity" in supplement
    assert "## Supplementary Table S8. Resource homogenization, shared-resource exposure and host-removal stress" in supplement
    assert "47.0%" in supplement
    assert "46.0%" in supplement
    assert "| 0.052 | 38 |" in supplement
    assert "| 0.057 | 37 |" in supplement
    assert "broader diets retained a more distributed host-contribution architecture" not in supplement.lower()

    # Three main figures plus one supplementary figure, with Figure 1 carrying C1 + C1b.
    for literal in (
        "**Figure 1. Anthropogenic host redistribution expands butterfly resource geography across specialization classes and through concentrated host contributions.**",
        "**Figure 2. Introduced host geography recovers butterfly occurrences beyond structural overlap expectations.**",
        "**Figure 3. Plant globalization homogenizes butterfly resource geography and increases shared-resource exposure.**",
        "**Supplementary Figure S1. Climate-associated filtering persists after geographic controls, whereas the predicted host-breadth release is unsupported.**",
    ):
        assert literal in manuscript
    assert claim_map["figure_claim_mapping"]["Figure_1"] == ["C1", "C1b"]
    assert claim_map["figure_claim_mapping"]["Figure_3"] == ["C0", "C0b"]

    # Submission-facing text must not retain development/review-history language.
    lowered = manuscript.lower()
    for forbidden in (
        "added after manuscript review",
        "the original descriptive analysis",
        "systematically misses",
        "most importantly",
        "no broad-generalist proportional advantage",
        "the absence of a broad-generalist advantage",
    ):
        assert forbidden not in lowered

    # Strong null and causal plant-trait overclaims are explicitly forbidden.
    forbidden_claims = "\n".join(claim_map["forbidden_overclaims"]).lower()
    assert "exactly zero effect" in forbidden_claims
    assert "demonstrated equivalent to a zero effect within ±0.10" in forbidden_claims
    assert "causal plant traits" in forbidden_claims
    assert "poaceae membership explains" in forbidden_claims


def test_current_v02_claim_map_preserves_layered_inference_boundaries():
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.2.json")
    claims = {row["id"]: row for row in claim_map["claims"]}

    assert "homogenizes regional butterfly resource assemblages" in claims["C0"]["claim"]
    assert "exact shared resources" in claims["C0b"]["claim"]
    assert "little relationship" in claims["C1"]["claim"]
    assert claims["C1b"]["status"].startswith("post-hoc descriptive decomposition")
    assert claims["C1b"]["allowed_use"].startswith(
        "Describe which resources generate reconstructed added opportunity"
    )
    assert claims["C2"]["status"].startswith("secondary occurrence validation")
    assert "Supplementary Information" in claims["C3"]["status"]
    assert claims["C5"]["decision"] == "PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED"

    language = claim_map["required_language"]
    assert "equivalence" in language["diet_breadth_equivalence"].lower()
    assert "unit-preserving" in language["host_contribution_concentration"]
    assert "fixed Brownian" in language["phylogeny"]
    assert "Discussion/SI" in language["host_prominence"]
    assert "fixed-margin" in language["resource_homogenization"]
    assert "realized competition" in language["competition_exposure"]
    assert "network concentration" in language["removal_stress"]


def test_invasion_and_competition_interpretation_is_bounded():
    manuscript = (
        ROOT / "manuscript" / "butterfly_specialization_ecology_v0.2.md"
    ).read_text(encoding="utf-8")
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.2.json")

    assert "do **not** demonstrate stronger competition" in manuscript
    assert "resource filtering" in manuscript
    assert "local abundance and performance" in manuscript.lower()
    assert "ecological traps" in manuscript
    assert "not an argument for retaining invasive plants" in manuscript
    assert "site-level management" in manuscript
    assert "realized competition" in claim_map["required_language"]["competition_exposure"]
    assert "network concentration and redundancy" in claim_map["required_language"]["removal_stress"]

def test_submission_citations_cover_figures_supplement_and_external_data():
    manuscript = (
        ROOT / "manuscript" / "butterfly_specialization_ecology_v0.2.md"
    ).read_text(encoding="utf-8")
    supplement = (
        ROOT / "manuscript" / "butterfly_specialization_supplement_v0.2.md"
    ).read_text(encoding="utf-8")

    for literal in (
        "Fig. 1a",
        "Fig. 1b",
        "Fig. 1c",
        "Fig. 2a",
        "Fig. 2b",
        "Fig. 3a",
        "Fig. 3b",
        "Fig. 3c",
        "Fig. S1",
        "Supplementary Table S1",
        "Supplementary Table S2",
        "Supplementary Table S3",
        "Supplementary Table S4",
        "Supplementary Table S5",
        "Supplementary Table S6",
        "Supplementary Table S7",
        "Supplementary Table S8",
        "(GBIF.org 2026)",
        "(Karger et al. 2017, 2021)",
    ):
        assert literal in manuscript

    assert "## References (working)" not in manuscript
    assert "## References" in manuscript
    assert "secondary in v0.2" not in supplement
    assert "shown as Fig. S1" in supplement
