#!/usr/bin/env python3
"""Verify one integrated butterfly GEB manuscript against source-locked raw receipts.

Independent file-level crosschecks; no model refit, no inference from private
GBIF annual Parquet, no network/API request, no journal submission. Fail if
numbers, sampling grain, data source identity, abstract or figures drift.
"""
from __future__ import annotations
import json
import hashlib
import math
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / "manuscript/integrated_butterfly_resource_access_20261109"
MANUSCRIPT = D / "INTEGRATED_RESEARCH_ARTICLE_V01.md"
LEDGER = D / "EVIDENCE_LEDGER_V01.json"
FIGURES = D / "FIGURE_PLAN_V01.md"
DECISION = D / "ONE_PAPER_DECISION_JA.md"

def read_json(path):
    return json.loads((ROOT/path).read_text(encoding="utf8"))

def approx(actual, expected, tol=1e-7):
    if not math.isclose(float(actual), float(expected), rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"Expected {expected}, got {actual}")

def check() -> dict:
    doc = MANUSCRIPT.read_text(encoding="utf8")
    figures=FIGURES.read_text(encoding="utf8")
    decision=DECISION.read_text(encoding="utf8")
    ledger=json.loads(LEDGER.read_text(encoding="utf8"))
    if ledger.get("canonical_manuscript") != str(MANUSCRIPT.relative_to(ROOT)):
        raise AssertionError("Integrated canonical manuscript pointer changed")
    if len(ledger["claims"])<18 or len({x["id"] for x in ledger["claims"]})!=len(ledger["claims"]):
        raise AssertionError("Unclear or duplicated integrated evidence ledger")
    if not all(x.get("population") and x.get("interpretation") and x.get("source")
               for x in ledger["claims"]):
        raise AssertionError("All claims require a population, interpretation and source")
    assert "唯一の文章正本" in decision
    assert "GEB" in decision and "同時投稿" in decision
    assert "JBI" in decision
    assert "one" in doc.splitlines()[2].lower()
    assert all(f"**Figure {i}." in doc for i in range(1,6))
    assert all(f"## Figure {i}" in figures for i in range(1,6))
    for needed in ("239 butterflies","155-butterfly","225-butterfly","51 landmass"):
        if needed not in doc:
            raise AssertionError("Missing butterfly-source population disclosure: "+needed)
    if "not independently" not in doc.lower() and "not a validated" not in doc.lower():
        raise AssertionError("Biological identification boundary not included")
    if "114 unlisted" not in doc and "114 species" not in doc and "114 unlisted" not in figures:
        raise AssertionError("Positive-only butterfly migratory roster lacks UNKNOWN clarification")
    start=doc.index("## Structured abstract")
    end=doc.index("\n---",start)
    abstract=doc[start:end]
    abstract_words=re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",abstract)
    if len(abstract_words)>300:
        raise AssertionError(f"GEB structured abstract is too long: {len(abstract_words)}")
    for label in ("Aim.", "Location.", "Time period.", "Major taxa studied.", "Methods.",
                  "Results.", "Main conclusions."):
        if label not in abstract: raise AssertionError("Missing GEB abstract section "+label)
    if doc.count("\n## 3. Results")!=1 or doc.count("\n## 4. Discussion")!=1:
        raise AssertionError("Not a single standard Research Article")
    claims={x["id"]:x for x in ledger["claims"]}

    # Global butterfly resource estimates.
    x=read_json("benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json")
    full=x["full_resource_eligible_panel"]
    assert (full["species"],full["species_expanded_by_introduced_host_ranges"])==(239,206)
    assert (full["total_native_species_units"],full["total_contemporary_species_units"],
            full["total_introduced_added_species_units"])==(26530,41083,14553)
    approx(100*full["total_introduced_added_species_units"]/
           full["total_native_species_units"],claims["R1"]["estimate"]["relative_pct"],tol=.1)
    approx(full["spearman"]["host_family_count_vs_log_resource_expansion"],
           claims["R2"]["estimate"]["rho"],tol=.001)
    assert full["host_breadth_strata"]["1_family"]["species_expanded"]==82
    assert full["host_breadth_strata"]["1_family"]["species"]==97

    # Nonidentical WGSRPD3 class denominators, no hidden 14pp alternative.
    x=read_json("results/island-mainland-resource-q1/island_mainland_summary.json")
    isl=x["groups"]["island"];main=x["groups"]["mainland"]
    assert isl["all_regions"]==126 and main["all_regions"]==237
    assert (isl["native_butterfly_region_units"],isl["contemporary_butterfly_region_units"])==(4829,10570)
    assert (main["native_butterfly_region_units"],main["contemporary_butterfly_region_units"])==(21215,29901)
    approx(isl["pooled_increase_pct_if_native_positive"],claims["I1"]["estimate"]["island_gain_pct"])
    approx(main["pooled_increase_pct_if_native_positive"],claims["I1"]["estimate"]["mainland_gain_pct"])
    m=read_json("results/island-mainland-resource-area-v02/area_adjusted_summary.json")["stratum_test"]
    assert len(m["strata_with_both_classes"])==15
    approx(100*m["weighted_island_minus_mainland_share"],
           claims["I2"]["estimate"]["matched_island_added_share_excess_pp"],tol=1e-5)

    # Port count models are DISTINCT and same effort-covered n=98.
    flora=read_json("results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json")["models"]
    n1=flora["all_flora_with_gbif_total_effort"];n2=flora["host_taxa_with_effort"]
    conditional=read_json("results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json")["model_results"]["M1eff_port_host_taxa_given_alien_flora"]
    assert n1["n"]==n2["n"]==conditional["n"]==98
    approx(n1["regression_beta_per_log2_port_unit"],claims["P1"]["estimate"]["whole_introduced_flora_beta"])
    approx(n2["regression_beta_per_log2_port_unit"],claims["P1"]["estimate"]["distinct_introduced_host_beta"])
    approx(conditional["coefficient"],claims["P1"]["estimate"]["host_given_whole_flora_beta"])

    # Independent 155-species butterfly recording comparison.
    forecasts={}
    for period in ("temporal","spatial"):
        r=read_json("provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_"+period+"_results.json")
        assert r["test_species"]==155 and r["test_rows"]==208482
        old=r["results"]["M3_both"]["species_equal_mse"]
        add=r["results"]["M5_host_composition"]["species_equal_mse"]
        forecasts[period]=100*(old-add)/old
    approx(forecasts["temporal"],claims["O2"]["estimate"]["time_held_species_equal_MSE_reduction_pct"])
    approx(forecasts["spatial"],claims["O2"]["estimate"]["spatial_held_species_equal_MSE_reduction_pct"])

    # Original 51 nonindependent pair correlation (NOT same as new 51 index snapshot).
    pair=read_json("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json")
    assert pair["eligible_islands"]==51 and pair["all_pairs"]==1275
    for k,field in [("spearman_host_new","later_host_butterfly_rho"),
                    ("spearman_host_earlier","earlier_host_butterfly_rho"),
                    ("partial_rank_host_new","partial_later")]:
        approx(pair[k],claims["A1"]["estimate"][field])

    # New independent GBIF index source 225 × 51 and event eligibility.
    receipt=read_json("data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json")
    source=ROOT/"data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
    assert receipt["merged_table_sha256"]==hashlib.sha256(source.read_bytes()).hexdigest()
    assert receipt["cells"]==11475 and receipt["source_islands"]==51 and receipt["source_butterflies"]==225
    assert receipt["totals"]["first_indexed_in_2011_2017"]==221
    assert receipt["totals"]["first_indexed_in_2018_2022"]==96
    assert len(receipt["shard_audits"])==6
    summar=read_json("results/public51_butterfly_trait_access_20261109/public51_trait_summary.json")
    assert summar["source_cells"]==11475
    p1=summar["periods"]["2011_2017"]["cohort"];p2=summar["periods"]["2018_2022"]["cohort"]
    assert (p1["prior_source_and_trait_rows"],p1["prior_source_and_trait_events"])==(4374,187)
    assert (p2["prior_source_and_trait_rows"],p2["prior_source_and_trait_events"])==(4933,83)
    approx(summar["periods"]["2018_2022"]["heldout"]["geographic_cluster"]["sea_over_geodesic"]["brier"],
           claims["A6"]["estimate"]["sea_beyond_distance_Brier_reduction_pct"],tol=.001)

    # Source-verified conditional butterfly case-count null (not causal p-values).
    null=read_json("results/public51_butterfly_trait_access_20261109/public51_butterfly_conditional_trait_marine.json")
    assert null["source_sha256"]==receipt["merged_table_sha256"]
    late=null["periods"]["2018_2022"]["within_geographic_cluster_two_margin"]["tests"]
    approx(late["geo_access_mean"]["p_two_sided_monte_carlo"],
           claims["A5"]["estimate"]["region_fixed_geo_proximity_p"])
    approx(late["sea_residual_access_mean"]["p_two_sided_monte_carlo"],
           claims["A5"]["estimate"]["region_fixed_ocean_access_residual_p"])
    approx(late["wing_proxy_x_geo"]["p_two_sided_monte_carlo"],
           claims["A7"]["estimate"]["within_region_wing_geo_p"],tol=.001)
    approx(late["wing_proxy_x_sea_residual"]["p_two_sided_monte_carlo"],
           claims["A7"]["estimate"]["within_region_wing_ocean_p"],tol=.001)

    # Positive-only indirect migrant roster and sensitivity to removing two
    # heavy-port landmass groups both as targets and candidate sources.
    roster=read_json("data/exploratory/public51_author_migrant_20261109/migration_gbif_predicate_audit.json")
    assert roster["candidate_225_exact_matched"]==111
    assert roster["candidate_225_not_matched_or_unresolved_UNKNOWN"]==114
    pos=read_json("results/public51_author_selected_migrant_access_20261109/author_selected_migrant_two_margin_index_null.json")
    remove=read_json("results/public51_author_selected_migrant_access_20261109/author_selected_migrant_no9_15_source_null.json")
    assert remove["sensitivity_excluded_original_landmass_ids"]==["9","15"]
    assert remove["source_island_groups_remaining"]==49
    assert pos["source_sha256"]==remove["source_sha256"]==receipt["merged_table_sha256"]
    assert (pos["periods"]["2011_2017"]["indexed_first_cases"],
            pos["periods"]["2018_2022"]["indexed_first_cases"])==(145,68)
    late=pos["periods"]["2018_2022"]["within_original_geographic_region"]["tests"]
    newlate=remove["periods"]["2018_2022"]["within_original_geographic_region"]["tests"]
    approx(late["nearest_prior_conspecific_geo_access"]["p_two_sided_monte_carlo"],
           claims["A8"]["estimate"]["late_geo_regionfixed_p"])
    approx(late["residual_uniform_ocean_access"]["p_two_sided_monte_carlo"],
           claims["A8"]["estimate"]["late_sea_regionfixed_p"])
    approx(newlate["residual_uniform_ocean_access"]["p_two_sided_monte_carlo"],
           claims["A8"]["estimate"]["late_sea_after_removing_landmass9_15_p"],tol=.001)

    # Integrated report remains an observational study.
    must_state=("not observed", "not verified", "not established", "not a")
    if sum(t in doc.lower() for t in must_state)<2:
        raise AssertionError("Butterfly inference caveats lost")
    return {
      "status":"PASS_SOURCE_GROUNDED_SINGLE_BUTTERFLY_RESEARCH_ARTICLE_V01",
      "structured_abstract_words":len(abstract_words),
      "claim_count":len(ledger["claims"]),
      "panels":5,
      "source225x51_sha256":receipt["merged_table_sha256"],
      "checked_source_scopes":["239 butterfly global", "126/237 insular-mainland",
        "15 matched strata", "98 effort-covered island botanical regions",
        "155 butterfly forecast", "225 × 51 island GBIF index",
        "199 butterfly traits", "111 indirect candidate migrant + 114 unknown",
        "51-landmass and 49-landmass constrained observational checks"],
      "scope_limit":"No model refitting, source licensing, actual rendered composite figure or bibliography verification performed by this script",
    }

if __name__ == "__main__":
    print(json.dumps(check(),ensure_ascii=False,indent=2))
