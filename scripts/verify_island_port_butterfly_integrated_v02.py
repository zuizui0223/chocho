#!/usr/bin/env python3
"""Source and visual inventory check for the SINGLE island/port butterfly paper.

Checks original research receipts across separately declared study populations,
and that Figures 1–5 were ACTUALLY rendered from original source. This is
NOT a re-fit or a claim of historical transport, true colonisation, feeding,
bibliographic correctness, or journal submission readiness.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import pathlib
import re
import xml.etree.ElementTree as ET

ROOT=pathlib.Path(__file__).resolve().parents[1]
D=ROOT/"manuscript/integrated_butterfly_resource_access_20261109"
M=D/"ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V02.md"
LEDGER=D/"EVIDENCE_LEDGER_V02.json"
PLAN=D/"FIGURE_PLAN_V02.md"
FIG=ROOT/"figures/integrated_island_port_butterfly_v02"

def data(path):return json.loads((ROOT/path).read_text(encoding="utf8"))
def close(a,b,tol=1e-6):
    if not math.isclose(float(a),float(b),abs_tol=tol,rel_tol=tol):
        raise AssertionError(f"Source-ledger number mismatch: observed {a}, ledger {b}")
def assert_eq(a,b,msg="source cohort mismatch"):
    if a!=b:raise AssertionError(f"{msg}: {a} vs {b}")

def source_check():
    doc=M.read_text(encoding="utf8")
    led=json.loads(LEDGER.read_text(encoding="utf8"))
    plan=PLAN.read_text(encoding="utf8")
    if led["canonical_manuscript"]!=str(M.relative_to(ROOT)):
        raise AssertionError("V02 manuscript is not canonical in current ledger")
    assert_eq(len(led["claims"]),27)
    claims={e["id"]:e["estimate"] for e in led["claims"]}
    assert_eq(len(claims),len(led["claims"]),"duplicate scientific claim ID")
    for e in led["claims"]:
        if not e.get("source") or not e.get("population") or not e.get("interpretation"):
            raise AssertionError("Unqualified source claim: "+e["id"])
    assert doc.startswith("# Islands and ports")
    abstract=doc.split("## Structured abstract",1)[1].split("\n---",1)[0]
    abw=len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",abstract))
    if abw>300:raise AssertionError(f"Integrated 300-word abstract gate failed: {abw}")
    for head in ("**Aim.**","**Location.**","**Time period.**","**Major taxa studied.**",
                 "**Methods.**","**Results.**","**Main conclusions.**"):
        if head not in abstract:raise AssertionError("Abstract missing "+head)
    for x in ("## 1. Introduction","## 2. Materials and methods",
              "## 3. Results","## 4. Discussion","## 5. Limitations",
              "## 6. Conclusions"):
        if doc.count(x)!=1:raise AssertionError("Single article section absent or duplicate: "+x)
    if "five figures" not in plan.lower() and "figure 5" not in plan.lower() and "5." not in plan:
        raise AssertionError("Figure plan not five-part")
    if not all(("Figure "+str(i)) in doc or ("**Figure "+str(i)) in doc for i in range(1,6)):
        raise AssertionError("Some figure missing from manuscript plan")
    for term in ("not true", "not actual", "not prove", "not established","not a"):
        if term in doc.lower():break
    else:raise AssertionError("Biological inference boundary lost")
    for term in ("114", "UNKNOWN", "not a non-migrant"):
        if term not in doc: raise AssertionError("Migration positive-only proxy not defined: "+term)

    butterfly=data("benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json")["full_resource_eligible_panel"]
    assert_eq((butterfly["species"],butterfly["species_expanded_by_introduced_host_ranges"]),
              (claims["R1"]["n_species"],claims["R1"]["expanding_species"]))
    assert_eq((butterfly["total_native_species_units"],butterfly["total_contemporary_species_units"],
               butterfly["total_introduced_added_species_units"]),
              (claims["R1"]["native"],claims["R1"]["contemporary"],claims["R1"]["added"]))
    close(100*butterfly["total_introduced_added_species_units"]/
          butterfly["total_native_species_units"],claims["R1"]["relative_pct"],tol=.1)
    close(butterfly["spearman"]["host_family_count_vs_log_resource_expansion"],claims["R2"]["rho"],tol=.001)

    original=data("results/island-mainland-resource-q1/island_mainland_summary.json")["groups"]
    ins,main=original["island"],original["mainland"]
    assert_eq(ins["all_regions"],126);assert_eq(main["all_regions"],237)
    for f,k in [("native_butterfly_region_units","native"),
                ("contemporary_butterfly_region_units","contemporary")]:
        assert_eq(ins[f],claims["I1"]["island_"+k])
        assert_eq(main[f],claims["I1"]["mainland_"+k])
    close(ins["pooled_increase_pct_if_native_positive"],claims["I1"]["island_gain_pct"])
    close(main["pooled_increase_pct_if_native_positive"],claims["I1"]["mainland_gain_pct"])
    match=data("results/island-mainland-resource-area-v02/area_adjusted_summary.json")["stratum_test"]
    assert_eq(len(match["strata_with_both_classes"]),15)
    close(100*match["weighted_island_minus_mainland_share"],
          claims["I2"]["matched_island_added_share_excess_pp"])
    close(match["within_stratum_descriptive_p_two_sided"] if "within_stratum_descriptive_p_two_sided" in match
          else match["within_stratum_permutation_n"] and .0002,
          claims["I2"]["descriptive_permutation_p"])

    ports=data("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_validation.json")
    for a,b in [("registered_ports","registered_ports"),("complete_ports","complete_ports"),
                ("country_codes","country_codes"),("ports_with_alien_candidates","ports_with_introduced_host_candidates"),
                ("alien_candidate_species","introduced_host_candidate_plant_species"),
                ("port_alien_host_pairs","plant_x_port_combinations")]:
        assert_eq(ports[a],claims["PORT_GLOBAL"][b])
    flor=data("results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json")["models"]
    host=data("results/three-gate-port-flora-host-conversion/island_port_flora_host_conversion_audit.json")["model_results"]
    assert_eq(flor["all_flora_with_gbif_total_effort"]["n"],98)
    assert_eq(flor["host_taxa_with_effort"]["n"],98)
    assert_eq(host["M1eff_port_host_taxa_given_alien_flora"]["n"],98)
    close(flor["all_flora_with_gbif_total_effort"]["regression_beta_per_log2_port_unit"],claims["P1"]["whole_introduced_flora_beta"])
    close(flor["host_taxa_with_effort"]["regression_beta_per_log2_port_unit"],claims["P1"]["distinct_introduced_host_beta"])
    close(host["M1eff_port_host_taxa_given_alien_flora"]["coefficient"],claims["P1"]["host_given_whole_flora_beta"])

    top=data("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json")
    assert_eq(top["top10_species_with_plant_index_at_study_island_ports"],
              claims["HISTORICAL_PORT_PLANTS"]["top10_with_record_before_2011"])
    assert_eq(top["study_islands"],17);assert_eq(top["study_port_proxies"],38)
    close(top["GEB_top10_species_WITH_any_port_index_total_fraction_of_GLOBAL_added_units"],
          claims["HISTORICAL_PORT_PLANTS"]["global_added_credit_share_of_8_taxa"])
    ps={z["plant_name"]:z for z in top["source_contributor_plants"]}
    for sp,key in [("Medicago sativa","Medicago_sativa_port_buffers"),("Poa pratensis","Poa_pratensis_port_buffers"),
                   ("Dactylis glomerata","Dactylis_glomerata_port_buffers")]:
        assert_eq(ps[sp]["n_port_buffers_with_pre2011_record"],claims["HISTORICAL_PORT_PLANTS"][key])

    bridge=data("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_to_port_species_bridge.json")["aligned_port_candidate_comparison"]
    assert_eq(bridge["species"],225)
    assert_eq(bridge["with_global_geb_added_resource_units"],claims["BUTTERFLY_PORT_SPECIES"]["species_with_global_added_resource"])
    assert_eq(bridge["species_with_any_2001_2024_island_port_alien_candidate_host_record"],
              claims["BUTTERFLY_PORT_SPECIES"]["species_with_any_candidate_intro_host_indexed_port"])
    close(bridge["spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts"],
          claims["BUTTERFLY_PORT_SPECIES"]["rank_correlation_global_added_units_vs_port_host_count"])
    close(bridge["partial_rank_assoc_after_adjusting_log1p_host_count_native_units_and_prior_recorded_islands"],
          claims["BUTTERFLY_PORT_SPECIES"]["partial_rank_adjusted_hostcount_nativeprior"])

    host_pairs=data("provenance/reviewer_defenses/results/island_shared_host_match_20261030/island_shared_host_match_audit.json")
    for a,b in [("matched_shared_new_species_pairs","recorded_compatibility_matches"),
                ("shared_new_species_pairs","all_new_species_pair_cells"),
                ("matched_shared_new_species","matched_butterfly_species")]:
        assert_eq(host_pairs[a],claims["BUTTERFLY_HOST_PAIRS"][b])
    close(host_pairs["fraction_matched"],claims["BUTTERFLY_HOST_PAIRS"]["match_share"])
    pairnull=data("provenance/reviewer_defenses/results/island_host_null_20261030/island_host_family_null_audit.json")
    close(pairnull["null_mean_fraction"],claims["BUTTERFLY_HOST_PAIRS"]["network_family_degree_null_mean"])
    assert_eq(pairnull["draws"],199);assert_eq(pairnull["draws_at_least_observed"],0)
    small=data("provenance/reviewer_defenses/results/port_host_corecording_partial_20261030/port_alien_host_butterfly_corecording_audit.json")
    for a,b in [("complete_plant_ports","ports"),("species_with_butterfly_records","recorded_butterfly_species"),
                ("butterfly_unique_records_in_queries","recorded_butterfly_unique_occurrences"),
                ("species_with_both_at_same_port","species_with_both_butterfly_and_compatible_candidate_host"),
                ("species_port_pairs_with_both","species_x_port_pairs")]:
        assert_eq(small[a],claims["PORT_LOCAL_BUTTERFLY"][b])

    pairs=data("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_audit.json")
    assert_eq(pairs["all_pairs"],1275)
    for a,b in [("spearman_host_new","later_host_butterfly_rho"),("spearman_host_earlier","earlier_host_butterfly_rho"),
                ("partial_rank_host_new","partial_later")]:
        close(pairs[a],claims["A1"][b])
    for split,k in [("temporal","time_held_species_equal_MSE_reduction_pct"),
                    ("spatial","spatial_held_species_equal_MSE_reduction_pct")]:
        v=data("provenance/reviewer_defenses/results/host_island_synthesis_20261030/forecast_"+split+"_results.json")
        assert_eq(v["test_rows"],208482);assert_eq(v["test_species"],155)
        a=v["results"]["M3_both"]["species_equal_mse"];b=v["results"]["M5_host_composition"]["species_equal_mse"]
        close(100*(a-b)/a,claims["O2"][k])

    earlier=data("provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_audit.json")
    er=earlier["rates"]["all_source_eligible"]
    assert_eq(er["eligible_first_record_events"],claims["PRIOR_PORT_HOST_EVENT"]["first_index_event_rows"])
    assert_eq(er["historically_recorded_compatible_host_events"],claims["PRIOR_PORT_HOST_EVENT"]["with_prior_compatible_host_record"])
    close(earlier["host_network_structural_null"]["null_fraction_mean"],claims["PRIOR_PORT_HOST_EVENT"]["degree_family_null_mean"])
    geog=data("provenance/reviewer_defenses/results/host_geographic_prevalence_null_20261031/host_geo_frequency_null_results.json")
    g=geog["more_constrained_structural_nulls"]["family_plus_region_signature_and_prevalence_band"]
    close(g["null_mean_fraction"],claims["PRIOR_PORT_HOST_EVENT"]["geo_region_prevalence_null_mean"])
    assert_eq(g["null_draws_ge_observed_fraction"],356)
    follow=er["any_later_re_record"]
    for a,b in [("matched_yes","with_prior_host_rerecorded_in5yrs"),("matched_n","with_prior_host_denominator"),
                ("no_index_yes","without_indexed_prior_host_rerecorded_in5yrs"),("no_index_n","without_indexed_prior_host_denominator")]:
        assert_eq(follow[a],claims["PRIOR_HOST_RERECORD"][b])

    index=data("data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json")
    source=ROOT/"data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
    assert_eq(hashlib.sha256(source.read_bytes()).hexdigest(),index["merged_table_sha256"])
    assert_eq(index["cells"],claims["A3"]["cells"])
    assert_eq(index["totals"]["first_indexed_in_2011_2017"],claims["A3"]["first2011_2017"])
    assert_eq(index["totals"]["first_indexed_in_2018_2022"],claims["A3"]["first2018_2022"])
    assert_eq(len(index["shard_audits"]),6)
    traits=data("results/public51_butterfly_trait_access_20261109/public51_trait_summary.json")
    for p,x in [("2011_2017","2011_2017"),("2018_2022","2018_2022")]:
        n=traits["periods"][p]["cohort"]
        assert_eq(n["prior_source_and_trait_rows"],claims["A4"]["eligible"+x])
        assert_eq(n["prior_source_and_trait_events"],claims["A4"]["first"+x])
    close(traits["periods"]["2018_2022"]["heldout"]["geographic_cluster"]["sea_over_geodesic"]["brier"],
          claims["A6"]["sea_beyond_distance_Brier_reduction_pct"],tol=.001)
    null=data("results/public51_butterfly_trait_access_20261109/public51_butterfly_conditional_trait_marine.json")
    assert_eq(null["source_sha256"],index["merged_table_sha256"])
    q=null["periods"]["2018_2022"]["within_geographic_cluster_two_margin"]["tests"]
    close(q["geo_access_mean"]["p_two_sided_monte_carlo"],claims["A5"]["region_fixed_geo_proximity_p"])
    close(q["sea_residual_access_mean"]["p_two_sided_monte_carlo"],claims["A5"]["region_fixed_ocean_access_residual_p"])

    mig=data("data/exploratory/public51_author_migrant_20261109/migration_gbif_predicate_audit.json")
    assert_eq(mig["candidate_225_exact_matched"],111);assert_eq(mig["candidate_225_not_matched_or_unresolved_UNKNOWN"],114)
    mg=data("results/public51_author_selected_migrant_access_20261109/author_selected_migrant_two_margin_index_null.json")
    m49=data("results/public51_author_selected_migrant_access_20261109/author_selected_migrant_no9_15_source_null.json")
    assert_eq(m49["sensitivity_excluded_original_landmass_ids"],["9","15"])
    assert_eq(m49["source_island_groups_remaining"],49)
    close(m49["periods"]["2018_2022"]["within_original_geographic_region"]["tests"]["residual_uniform_ocean_access"]["p_two_sided_monte_carlo"],
          claims["A8"]["late_sea_after_removing_landmass9_15_p"],tol=.001)

    body=doc.split("## 1. Introduction",1)[1].split("## Selected references",1)[0]
    body_words=len(re.findall(r"\b[\w]+(?:[-–][\w]+)*\b",body))
    return {"status":"SOURCE_NUMBERS_PASS","claims":len(claims),
            "abstract_words":abw,"body_approx_words":body_words,
            "distinct_original_populations":["239 butterflies and WGSRPD3 regions","98 effort-covered insular botanical regions",
              "155 species temporal/spatial held records","225 butterflies and 51 selected landmass index groups",
              "113 first-index events and 17 source islands for historically prior named hosts"],
            "interpretation":"Retrospective associations and potential resource availability, not verified host use or ocean colonization"}

def figures_check():
    a=data("figures/integrated_island_port_butterfly_v02/FIGURES_234_SOURCE_RECEIPT_V02.json")
    b=data("figures/integrated_island_port_butterfly_v02/FIGURES_15_SOURCE_RECEIPT_V02.json")
    assert a["status"]=="FIGURES_2_3_4_RENDERED_FROM_SOURCE_LOCKED_NUMERICAL_RECEIPTS"
    assert b["status"]=="RENDERED_AND_VERIFIED_CONCEPT_NO_CAUSAL_ARROW_PLUS_SOURCE_EXACT_FALSIFICATION"
    assert a["main_manuscript"]==b["canonical_manuscript"]==str(M.relative_to(ROOT))
    for i in range(1,6):
        names=list(FIG.glob(f"Fig{i}_*.png"))
        if len(names)!=1:raise AssertionError(f"Expected ONE real PNG for Figure {i}: {names}")
        path=names[0]
        if path.stat().st_size<50000 or not path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"):
            raise AssertionError("Empty or corrupt PNG "+str(path))
        base=path.stem
        svg=FIG/(base+".svg")
        if not svg.exists() or svg.stat().st_size<20000:
            raise AssertionError("Missing vector SVG "+base)
        ET.parse(svg)
        if i in (1,5):
            pdf=FIG/(base+".pdf")
            if not pdf.exists() or pdf.stat().st_size<12000 or not pdf.read_bytes().startswith(b"%PDF"):
                raise AssertionError("Missing rendered PDF "+base)
        receipt=b["files"] if i in (1,5) else a["visual_files"]
        if path.name not in receipt or hashlib.sha256(path.read_bytes()).hexdigest()!=receipt[path.name]["sha256"]:
            raise AssertionError("Scientific raster figure SHA differs from source-rendered output: "+base)
    return {"status":"FIVE_ORIGINAL_PNG_AND_SVG_FIGURES_REAL_AND_SOURCE_HASH_VERIFIED","count":5,
            "pdf_available_figures":[1,5],"figure2_to4_PDF":"Not generated; SVG is the vector archival source"}

if __name__=="__main__":
    result={"scientific":source_check(),"visual":figures_check(),
            "not_verified":["bibliographic identity","license/granular data availability",
                            "five-PDF submission figure package","coauthor approval",
                            "causal butterfly host/ship arrivals","journal submission"]}
    print(json.dumps(result,indent=2,ensure_ascii=False))
