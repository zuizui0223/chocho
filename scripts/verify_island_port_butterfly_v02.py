#!/usr/bin/env python3
"""Verify island/port-first butterfly integrated v0.2 scientific claims.

Source-integrity/claim audit only: independently reads the archived source
JSON/CSV receipts, DOES NOT reproduce models, re-query restricted historical
butterfly files, certify historical shipping, or generate final artwork.
"""
from __future__ import annotations
from pathlib import Path
import csv,hashlib,importlib.util,json,math,re

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"manuscript/integrated_butterfly_resource_access_20261109"
MS=DIR/"ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md"
FIG=DIR/"FIGURE_PLAN_V02.md"
LEDGER=DIR/"EVIDENCE_LEDGER_V02.json"

def load(path):
    return json.loads((ROOT/path).read_text(encoding="utf8"))

def same(actual,expected,tol=1e-6):
    if not math.isclose(float(actual),float(expected),rel_tol=tol,abs_tol=tol):
        raise ValueError(f"Source number differs: actual={actual}, reported={expected}")

def check():
    md=MS.read_text(encoding="utf8")
    figures=FIG.read_text(encoding="utf8")
    claims=load(str(LEDGER.relative_to(ROOT)))
    if not md.startswith("# Islands and ports"):
        raise AssertionError("Central manuscript does not lead with islands and ports")
    assert claims["canonical_manuscript"]==str(MS.relative_to(ROOT))
    assert len(claims["claims"])>=27
    if len(set(c["id"] for c in claims["claims"]))!=len(claims["claims"]):
        raise AssertionError("Duplicate claim provenance")
    assert all(c["source"] and c["interpretation"] and c["population"]
               for c in claims["claims"])
    c={z["id"]:z for z in claims["claims"]}
    for name in ("PORT_GLOBAL","BUTTERFLY_PORT_SPECIES","BUTTERFLY_HOST_PAIRS",
                 "HISTORICAL_PORT_PLANTS","PORT_LOCAL_BUTTERFLY",
                 "PRIOR_PORT_HOST_EVENT","PRIOR_HOST_RERECORD","ISLAND_HOST_FORECAST"):
        assert name in c, name
    a=md.index("## Structured abstract");b=md.index("\n---",a)
    abstract=md[a:b]
    wc=len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",abstract))
    if wc>300:raise AssertionError(f"Structured abstract >300 words: {wc}")
    assert abstract.index("Islands may")<abstract.index("butterfly")
    for section in ("## 1. Introduction","## 2. Materials and methods","## 3. Results",
                    "## 4. Discussion","## 5. Limitations","## 6. Conclusions"):
        if section not in md:raise AssertionError("Missing section "+section)
    results=md.split("## 3. Results")[1].split("## 4. Discussion")[0]
    order=["### 3.1. Islands", "### 3.2. Port", "### 3.3. Butterfly",
           "### 3.4. Butterfly", "### 3.5. Historically", "### 3.6. Marine"]
    indices=[results.index(x) for x in order]
    if indices!=sorted(indices):raise AssertionError("Island/port positive result structure lost")
    for x in ("98", "1,596", "142", "687", "25,742", "208", "1,163",
              "2001–2010","2011–2017","2018–2022","34/113","155",
              "225","239","−0.0415","+0.652","0.300",
              "0.215","0.418","4.12%","3.00%","−6.48%"):
        if x not in md:raise AssertionError("Missing essential butterfly-port evidence "+x)
    for i in range(1,6):
        if f"**Figure {i}" not in figures:
            raise AssertionError("No figure specification "+str(i))
    assert "not yet" in figures.lower() or "not yet" in md.lower()
    assert "not proof" in md.lower() or "not demonstrate" in md.lower()

    original=load("benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json")["full_resource_eligible_panel"]
    assert original["species"]==239
    assert original["species_expanded_by_introduced_host_ranges"]==206
    assert [original["total_native_species_units"],original["total_contemporary_species_units"]]==[26530,41083]
    same(100*(41083-26530)/26530,54.9,tol=.002)

    botanical=load("results/island-mainland-resource-q1/island_mainland_summary.json")["groups"]
    same(botanical["island"]["pooled_increase_pct_if_native_positive"],118.88589770138745)
    same(botanical["mainland"]["pooled_increase_pct_if_native_positive"],40.942729201037004)
    matched=load("results/island-mainland-resource-area-v02/area_adjusted_summary.json")["stratum_test"]
    assert len(matched["strata_with_both_classes"])==15
    same(matched["weighted_island_minus_mainland_share"]*100,10.702159757115921)
    same(matched["within_stratum_descriptive_p_two_sided"],.0002)

    broad=load("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_validation.json")
    x=c["PORT_GLOBAL"]["estimate"]
    for k in ("registered_ports","complete_ports","country_codes",
               "ports_with_alien_candidates"):
        target={"registered_ports":"registered_ports","complete_ports":"complete_ports",
                 "country_codes":"country_codes","ports_with_alien_candidates":"ports_with_introduced_host_candidates"}[k]
        assert broad[k]==x[target],k
    assert broad["alien_candidate_species"]==x["introduced_host_candidate_plant_species"]==687
    assert broad["port_alien_host_pairs"]==x["plant_x_port_combinations"]==25742
    assert broad["complete_ports"]<broad["registered_ports"],"Missing survey truncation"

    fm=load("results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json")["models"]
    flora=fm["all_flora_with_gbif_total_effort"];hosts=fm["host_taxa_with_effort"]
    residual=load("results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json")["model_results"]["M1eff_port_host_taxa_given_alien_flora"]
    assert flora["n"]==hosts["n"]==residual["n"]==98
    same(flora["regression_beta_per_log2_port_unit"],.24447852550694815)
    same(hosts["regression_beta_per_log2_port_unit"],.16390595789001258)
    same(residual["coefficient"],-.041526655753407464)

    bridge=load("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_to_port_species_bridge.json")
    assert bridge["aligned_port_candidate_comparison"]["species"]==225
    cp=bridge["aligned_port_candidate_comparison"]
    assert cp["with_global_geb_added_resource_units"]==197
    assert cp["species_with_any_2001_2024_island_port_alien_candidate_host_record"]==159
    same(cp["spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"],c["BUTTERFLY_PORT_SPECIES"]["estimate"]["rank_correlation_global_added_units_vs_port_host_count"])
    same(cp["partial_rank_assoc_after_adjusting_log1p_host_count_native_units_and_prior_recorded_islands"],c["BUTTERFLY_PORT_SPECIES"]["estimate"]["partial_rank_adjusted_hostcount_nativeprior"])

    plants=load("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json")
    assert plants["top10_species_with_plant_index_at_study_island_ports"]==8
    same(plants["GEB_top10_species_WITH_any_port_index_total_fraction_of_GLOBAL_added_units"],.20701372105789737)
    plant={r["plant_name"]:r for r in plants["source_contributor_plants"]}
    for p,n in (("Medicago sativa",7),("Poa pratensis",8),("Dactylis glomerata",9)):
        assert plant[p]["n_port_buffers_with_pre2011_record"]==n

    pair=load("provenance/reviewer_defenses/results/island_shared_host_match_20261030/island_shared_host_match_audit.json")
    family_null=load("provenance/reviewer_defenses/results/island_host_null_20261030/island_host_family_null_audit.json")
    assert (pair["matched_shared_new_species_pairs"],pair["shared_new_species_pairs"])==(208,1163)
    same(pair["fraction_matched"],208/1163)
    same(family_null["null_mean_fraction"],.11896109956489237)
    assert family_null["draws"]==199 and family_null["draws_at_least_observed"]==0
    assert family_null["family_degree_checks_all_draws"]
    local=load("provenance/reviewer_defenses/results/port_host_corecording_partial_20261030/port_alien_host_butterfly_corecording_audit.json")
    assert local["complete_plant_ports"]==17 and local["species_with_both_at_same_port"]==16
    assert local["species_port_pairs_with_both"]==59
    assert local["butterfly_unique_records_in_queries"]==8525

    pair51=load("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json")
    assert (pair51["eligible_islands"],pair51["all_pairs"])==(51,1275)
    same(pair51["spearman_host_new"],.29980576483377025)
    same(pair51["spearman_host_earlier"],.4178255833828514)
    same(pair51["partial_rank_host_new"],.2153535842700987)

    prior=load("provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_audit.json")
    assert prior["butterfly_event_cohort"]["target_with_resolved_source_catalogue"]==113
    earlier=prior["rates"]["all_source_eligible"]
    assert earlier["historically_recorded_compatible_host_events"]==34
    assert earlier["any_later_re_record"]["matched_yes"]==26
    assert earlier["any_later_re_record"]["no_index_yes"]==70
    geography=load("provenance/reviewer_defenses/results/host_geographic_prevalence_null_20261031/host_geo_frequency_null_results.json")
    assert geography["observed_compatible_preindexed_host_events"]==34
    geo=geography["more_constrained_structural_nulls"]["family_plus_region_signature_and_prevalence_band"]
    same(geo["null_mean_fraction"],.302587475836629)
    assert geo["null_draws_ge_observed_fraction"]==356 and geo["draws"]==499

    later=load("provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_temporal_results.json")
    other=load("provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_spatial_results.json")
    for p in (later,other):
        assert p["test_species"]==155 and p["test_rows"]==208482
    def improvement(q):
        a=q["results"]["M3_both"]["species_equal_mse"]
        b=q["results"]["M5_host_composition"]["species_equal_mse"]
        return 100*(a-b)/a
    same(improvement(later),4.12017560793182)
    same(improvement(other),2.9952059516988094)
    island_pred_csv=ROOT/"provenance/reviewer_defenses/results/island_forecast_overlap_20261030/island_forecast_overlap_summary.csv"
    with island_pred_csv.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    strict=next(r for r in rows if r["cohort"]=="target51_strict_coordinate_agreement")
    assert int(strict["rows"])==39753 and int(strict["species"])==95
    same(strict["relative_MSE_improvement_pct"],4.9825045631704645)

    pub51=load("data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json")
    assert (pub51["source_butterflies"],pub51["source_islands"],pub51["cells"])==(225,51,11475)
    assert pub51["totals"]["first_indexed_in_2011_2017"]==221
    assert pub51["totals"]["first_indexed_in_2018_2022"]==96
    null=load("results/public51_butterfly_trait_access_20261109/public51_butterfly_conditional_trait_marine.json")
    assert null["source_sha256"]==pub51["merged_table_sha256"]
    follow=null["periods"]["2018_2022"]["within_geographic_cluster_two_margin"]["tests"]
    same(follow["geo_access_mean"]["p_two_sided_monte_carlo"],.0003332222592469177)
    same(follow["sea_residual_access_mean"]["p_two_sided_monte_carlo"],.0023325558147284237)
    m=load("results/public51_butterfly_trait_access_20261109/public51_trait_summary.json")
    same(m["periods"]["2018_2022"]["heldout"]["geographic_cluster"]["sea_over_geodesic"]["brier"],-6.479351063178629)
    mw=load("results/public51_author_selected_migrant_access_20261109/author_selected_migrant_no9_15_source_null.json")
    assert mw["sensitivity_excluded_original_landmass_ids"]==["9","15"]
    same(mw["periods"]["2018_2022"]["within_original_geographic_region"]["tests"]["residual_uniform_ocean_access"]["p_two_sided_monte_carlo"],.22859046984338555)
    print(json.dumps({"status":"PASS_ISLAND_PORT_BUTTERFLY_ONE_PAPER_V02",
         "abstract_words":wc,"source_checked_claims":len(claims["claims"]),
         "primary_source_scales":["239 butterfly resource","126 insular/237 mainland botanical regions",
             "1596 completed ports/142 country codes","98 effort-covered island regions",
             "225 species x 51 landmass co-index","155 butterfly forecasting","113 dated butterfly events"],
         "positive_and_negative_conclusions_preserved":True,
         "source_limit":"No primary model refits, rendered figure QA, original Table S1, prior local introduced status, source licences or journal submission checked"},
         ensure_ascii=False,indent=2))

if __name__=="__main__":
    check()
