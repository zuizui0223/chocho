#!/usr/bin/env python3
"""Check published source-receipt arithmetic and inference boundaries on main.

This is a CONSISTENCY verification of frozen results, not a rerun of upstream
HOSTS/WCVP, GBIF model fitting, AIS data, field records, or raw source acquisition.
Stdlib only; fail closed on denominator, model, or source-label drift.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "species": "data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "region": "results/island-mainland-resource-q1/island_mainland_summary.json",
    "pooled": "results/island_resource_gap_filling/pooled.csv",
    "decomposition": "results/island_resource_gap_filling/summary.json",
    "area": "results/island-mainland-resource-area-v02/area_adjusted_summary.json",
    "port126": "results/three-gate-port-all-flora/all_flora_port_association_audit.json",
    "port98": "results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json",
    "port_mod": "results/port_alien_flora_moderators_20261117/VERIFIED_SOURCE_CONTRASTS_20261117.json",
    "marine_flora": "results/WCVP126_crossed_native_alien_20261117/VERIFIED_SOURCE_NATIVE_ALIEN_CROSSED_20261117.json",
    "marine_butterfly": "results/island51_observed_ocean_circuit_20261106/model_comparison_summary.csv",
    "forecast": "provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_metric_validation.json",
    "forecast_pairs": "provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_paired_comparisons.csv",
    "literature": "results/hosts_island_local_use_20261116/VERIFIED_SUMMARY_ALL24_20261116.json",
    "family_control": "results/hosts_use_vs_flora_20261116/VERIFIED_FLORISTIC_AND_FAMILY_BASELINE_20261116.json",
    "endemic": "provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118/VERIFIED_ENDEMIC18_NONIDENTIFICATION_TRANSPORT_20261118.json",
    "endemic_8": "provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118/VERIFIED_INDEPENDENT8_HOST_GEOGRAPHY_AND_FIELD_GATES_20261118.json",
    "field": "provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118/VERIFIED_FOCAL_ALIEN_LARVAL_HOST_PLANT_SPECIES_20261118.json",
    "global_beta": "provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/gbif_global224_observed_beta_v1.json",
}
PIN_WCVP = "65bed76bae9d644ccb6ad200c05f9f5071d89e05"
PIN_HOSTS = "808e0b869f9ec1adf8efff87cf6a395adda103e0"


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def near(value: float, expected: float, message: str, atol: float = 1e-8) -> None:
    ensure(math.isfinite(float(value)) and abs(float(value) - float(expected)) <= atol,
           f"{message}: got {value!r}, wanted {expected!r}")


def load(name: str):
    path = ROOT / INPUTS[name]
    ensure(path.is_file(), f"Missing pinned evidence input: {path}")
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8-sig"))
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def main() -> dict:
    d = {k: load(k) for k in INPUTS}
    out = {"schema": "chocho_main_source_receipt_consistency_v1", "checks": {},
           "source_sha256": {k: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
                             for k, path in sorted(INPUTS.items())},
           "scope": "Frozen result arithmetic and source/estimand gates; not raw analyses or causal validation"}
    s = d["species"]
    ensure(len(s) == 239 and len({r["species"] for r in s}) == 239, "239 unique butterfly species")
    native = sum(int(r["native_resource_units"]) for r in s)
    current = sum(int(r["contemporary_resource_units"]) for r in s)
    added = sum(int(r["introduced_added_units"]) for r in s)
    ensure((native, current, added) == (26530, 41083, 14553), "Frozen 239-species totals")
    ensure(all(int(r["contemporary_resource_units"]) - int(r["native_resource_units"])
               == int(r["introduced_added_units"]) for r in s), "Source row addition identity")
    beneficiaries = sum(int(r["introduced_added_units"]) > 0 for r in s)
    ensure(beneficiaries == 206, "Beneficiaries out of original 239")
    near(100 * added / native, 54.855, "Global increase percent", atol=.03)
    out["checks"]["fixed_host_resources"] = {"n_butterflies": 239, "beneficiaries": beneficiaries,
         "native_pairs": native, "current_pairs": current, "new_pairs": added,
         "increase_pct": 100 * added / native, "unit": "butterfly x WGSRPD3 potential host-resource"}

    q = d["region"]
    ensure(q["source_commits"]["wcvp"] == PIN_WCVP and q["source_commits"]["hosts"] == PIN_HOSTS,
           "Original source commit identities")
    ensure(q["panel_total"]["panel_species"] == 239 and q["panel_total"]["n_panel_mismatches"] == 0,
           "Panel source identity and coverage")
    ensure(q["panel_total"]["native_total"] == native and q["panel_total"]["contemporary_total"] == current,
           "239-species total matches 126/237/mixed region arithmetic")
    groups = q["groups"]
    for name in ("island", "mainland", "mixed"):
        g = groups[name]
        ensure(g["native_butterfly_region_units"] + g["added_butterfly_region_units"]
               == g["contemporary_butterfly_region_units"], name + " cohort native/add identity")
    ensure(sum(g["native_butterfly_region_units"] for g in groups.values()) == native,
           "Geography native sum including mixed regions")
    ensure(sum(g["added_butterfly_region_units"] for g in groups.values()) == added,
           "Geography introduced sum including mixed regions")
    ensure(groups["island"]["all_regions"] == 126 and groups["mainland"]["all_regions"] == 237
           and groups["mixed"]["all_regions"] == 6, "Strict island/mainland/mixed source denominators")
    pooled = {r["geography"]: r for r in d["pooled"]}
    ensure(set(pooled) == {"island", "mainland"}, "Pooled baseline comparison cohort")
    for name in ("island", "mainland"):
        r, g = pooled[name], groups[name]
        ensure(int(float(r["native"])) == g["native_butterfly_region_units"], name + " pooled native")
        ensure(int(float(r["added"])) == g["added_butterfly_region_units"], name + " pooled added")
        near(float(r["coverage"]), float(r["native"]) / float(r["cells"]), name + " coverage fraction")
        near(float(r["fill"]), float(r["added"]) / float(r["empty"]), name + " gap filling")
        near(float(r["gain"]), float(r["added"]) / float(r["native"]), name + " gain")
    near(float(pooled["island"]["coverage"]) * 100, 16.035730889287, "Islands native starting coverage")
    near(float(pooled["mainland"]["coverage"]) * 100, 37.453877796021, "Mainland native starting coverage")
    near(float(pooled["island"]["fill"]) * 100, 22.7051611627447, "Island gap filling")
    near(float(pooled["mainland"]["fill"]) * 100, 24.5173309246923, "Mainland gap filling")
    gap = d["decomposition"]
    near(gap["baseline_component_pp"] + gap["fill_component_pp"],
         100 * (float(pooled["island"]["gain"]) - float(pooled["mainland"]["gain"])),
         "Baseline and fill decomposition of pooled gain")
    ensure(gap["equal_species_with_both_denominators"] == 239
           and sum(gap[k] for k in ("island_fill_higher", "mainland_fill_higher", "tied")) == 239,
           "Equal-species gap fill proportions total 239")
    area = d["area"]["stratum_test"]
    ensure(len(area["strata_with_both_classes"]) == 15, "15 matched Level1 x area strata")
    weights = sum(r["weight"] for r in area["strata_with_both_classes"])
    mean = sum(r["difference"] * r["weight"] for r in area["strata_with_both_classes"]) / weights
    near(mean, area["weighted_island_minus_mainland_share"], "Matched 15-stratum share")
    near(100 * mean, 10.70215975711592, "Contemporary introduced-share matched island excess, pp")
    ensure(abs(q["stratified_contrast"]["weighted_difference"] - mean) > .01,
           "Do not silently substitute unmatched 7-stratum pooled Q1 statistic")
    out["checks"]["island_amplification"] = {
        "island_gain_pct": 100 * float(pooled["island"]["gain"]),
        "mainland_gain_pct": 100 * float(pooled["mainland"]["gain"]),
        "island_coverage_pct": 100 * float(pooled["island"]["coverage"]),
        "mainland_coverage_pct": 100 * float(pooled["mainland"]["coverage"]),
        "island_gap_fill_pct": 100 * float(pooled["island"]["fill"]),
        "mainland_gap_fill_pct": 100 * float(pooled["mainland"]["fill"]),
        "matched_added_share_excess_pp": 100 * mean,
        "boundary": "relative gain is baseline-amplified; 15 matched share is not a colonisation effect"}

    port = d["port126"]
    ensure(port["input"]["regions"] == 126
           and port["input"]["all_flora_accepted_species"] == 360064, "All flora accepted source universe")
    primary = port["cohorts"]["all_126_islands"]["native_flora_adjusted"]
    near(primary["coefficient"], .22140536684115839, "126-island WPI vs introduced flora")
    # The locked 126 model's displayed source control is HOST RESOURCE units,
    # not WCVP native-plant richness. Do not silently rename this control.
    ensure("ln1p_native_host_resource_units" in primary["controls"],
           "126-island control is native host-resource units, not generic native-plant richness")
    effort = d["port98"]
    ensure(effort["scope"]["island_regions_with_successful_botanical_index_effort"] == 98,
           "Botanical effort source supported 98 regions")
    near(effort["models"]["all_flora_with_gbif_total_effort"]["regression_beta_per_log2_port_unit"],
         .24447852550694815, "98-island effort-controlled all-alien flora port slope")
    near(effort["models"]["host_taxa_with_effort"]["regression_beta_per_log2_port_unit"],
         .16390595789001258, "98-island effort-controlled introduced HOSTS candidate slope")
    moder = d["port_mod"]
    near(moder["source_provenance"]["original_port_flora_baseline_on_123"], .2334712198139267,
         "123-island corrected native-flora port baseline")
    ensure(moder["status"].endswith("NOT_HISTORICAL_SHIPPING_CAUSAL_EFFECT"),
           "Port moderation source causal boundary")
    out["checks"]["port_flora"] = {
        "whole_flora_126_beta": primary["coefficient"],
        "effort_98_whole_flora_beta": effort["models"]["all_flora_with_gbif_total_effort"]["regression_beta_per_log2_port_unit"],
        "effort_98_host_candidate_beta": effort["models"]["host_taxa_with_effort"]["regression_beta_per_log2_port_unit"],
        "native_flora_123_beta": moder["source_provenance"]["original_port_flora_baseline_on_123"],
        "control_warning": "126 coefficient conditions on native host RESOURCE counts, not native whole-flora richness; 123 model is separate"}

    water = d["marine_flora"]
    ensure(water["original_WCVP_accepted_species"] == 360064
           and water["original_GSHHG126_verified_model_coastal_islands"] == 104
           and water["original_exact_common_WCVP_native_alien_and_GSHHG_pair_sample"] == 5355,
           "Source-matched 104 WCVP coastal islands and 5355 dyads")
    for name in ("native", "alien"):
        delta = water["relative_MSE_positive_is_improvement"][name]
        ensure(delta["G_vs_B_pct"] > 0, name + " marine geometry should help in this model")
        ensure(delta["P_vs_B_pct"] < .1, name + " port-count proxy alone has small effect")
    near(water["relative_MSE_positive_is_improvement"]["native"]["G_vs_B_pct"], 2.234321123027717,
         "Native shortest water geometry gain")
    near(water["relative_MSE_positive_is_improvement"]["alien"]["G_vs_B_pct"], 1.7750986997451879,
         "Alien shortest water geometry gain")
    circ = [r for r in d["marine_butterfly"] if r["comparison"] == "base_vs_uniform"
            and r["target"] == "new_butterfly_jaccard"]
    ensure(len(circ) == 1, "Butterfly-specific independent 51-island model")
    z = circ[0]
    gain = 100 * (float(z["MSE_reference"]) - float(z["MSE_augmented"])) / float(z["MSE_reference"])
    near(gain, float(z["relative_MSE_reduction_percent"]), "51-island uniform sea circuit recomputed")
    near(gain, 5.0818868452225594, "Butterfly indexed-assemblage circuit increment")
    out["checks"]["marine_geography"] = {
        "native_water_shortest_geometry_pct": water["relative_MSE_positive_is_improvement"]["native"]["G_vs_B_pct"],
        "alien_water_shortest_geometry_pct": water["relative_MSE_positive_is_improvement"]["alien"]["G_vs_B_pct"],
        "butterfly_uniform_circuit_pct": gain,
        "warning": "Distinct flora/butterfly samples and models; no measured trajectories or plant transport"}

    forecast = d["forecast"]
    ensure(forecast["status"] == "PASS", "Original forecast metric audit status")
    vals = {(r["mode"], r["model"], r["metric"]): float(r["value"])
            for r in forecast["checks"]}
    paired = {(r["mode"], r["added"], r["reference"], r["metric"]): r
              for r in d["forecast_pairs"]}
    improvements = {}
    for mode, want in (("temporal", 4.12), ("spatial", 3.00),
                       ("species", -0.21)):
        key = (mode, "M5_host_composition", "M3_both", "mse")
        ref, new = vals[(mode, "M3_both", "mse")], vals[(mode, "M5_host_composition", "mse")]
        raw = paired[key]
        near(new - ref, float(raw["difference"]), "Forecast M5-M3 paired MSE " + mode)
        relative = 100 * (ref-new)/ref
        near(relative, want, "Forecast relative error gain " + mode, atol=.045)
        ensure(int(raw["species"]) == 155, "Forecast original species count " + mode)
        improvements[mode] = {"mse_relative_reduction_pct": relative,
                              "evaluated_cells": int(raw["rows"]),
                              "species_improved": int(raw["species_improved"])}
    # The >=100-km subset is in the paired-difference source, but not in the
    # full M3/M5 model-score receipt. Check only facts that are directly
    # reproducible from that source; do not invent a baseline to get 2.74%.
    remote = paired[("spatial_buffer100", "M5_host_composition", "M3_both", "mse")]
    ensure(int(remote["species"]) == 155 and int(remote["rows"]) == 102728
           and int(remote["species_improved"]) == 111
           and float(remote["difference"]) < 0,
           "Spatial >=100km M5-M3 paired source only, no invented percent")
    improvements["spatial_buffer100"] = {
        "paired_mse_difference": float(remote["difference"]),
        "evaluated_cells": int(remote["rows"]),
        "species_improved": int(remote["species_improved"]),
        "reported_relative_pct_requires_separate_model_score": 2.74}
    ensure(improvements["temporal"]["species_improved"] == 119
           and improvements["spatial"]["species_improved"] == 115,
           "Known-species forecast lift distribution")
    out["checks"]["independent_record_forecast"] = improvements

    literature = d["literature"]
    ensure(literature["status"] == "SOURCE_VERIFIED_DESCRIPTIVE_NOT_FIELD_CONSUMPTION",
           "Literature is not direct field larval consumption")
    l = literature["scopes"]["HOSTS_Papilionoidea__strict_exact_or_composite__all_hosts"]
    crop = literature["scopes"]["HOSTS_Papilionoidea__strict_exact_or_composite__exclude_FAO_exact"]
    ensure(l["island"]["localities"] == 13 and l["mainland"]["localities"] == 15,
           "Original mapped locality denominators")
    near(l["island"]["median_share"] * 100, 42.4242424242424, "Island reported alien larval-host fraction")
    near(l["mainland"]["median_share"] * 100, 9.09090909090909, "Mainland reported alien larval-host fraction")
    near(crop["island"]["median_share"]*100, 20.9005376344086, "Crop-removed island report fraction")
    near(crop["mainland"]["median_share"]*100, 8.11688311688311, "Crop-removed mainland report fraction")
    fam = d["family_control"]["scopes"]["HOSTS_Papilionoidea__exclude_FAO_exact__min10"]
    f = fam["contrasts"]["family_flora_adjusted_excess"]
    ensure(f["median_difference"] < 0 and f["descriptive_locality_bootstrap95"][1] > 0,
           "No robust family-adjusted extra island preference")
    out["checks"]["host_use_reporting"] = {
        "island_locality_median_pct": 100*l["island"]["median_share"],
        "mainland_locality_median_pct": 100*l["mainland"]["median_share"],
        "crop_excluded_island_pct": 100*crop["island"]["median_share"],
        "crop_excluded_mainland_pct": 100*crop["mainland"]["median_share"],
        "crop_and_family_adjusted_excess_pp": 100*f["median_difference"],
        "warning": "Same HOSTS source network; geographic reporting and family composition not wild feeding rate"}

    endemic = d["endemic"]
    scope = endemic["source_scope"]
    ensure(scope["living_butterflies"] == 41
           and scope["independently_origin_classified_endemic_extant_butterflies"] == 18
           and scope["living_nonendemic_butterflies"] == 23,
           "Independent 41 living butterflies population")
    ensure(endemic["original_selection"]["original_final239"]["endemic_n"] == 0
           and endemic["original_selection"]["original_final239"]["nonendemic_n"] == 13,
           "Original GEB239 excludes all 18 Macaronesian island endemics")
    effect = endemic["foodplant_source_selected_effect"]
    ensure(effect["exact_named_larval_host_plant_species_known_butterflies_n"] == 8
           and effect["source_list_observed_positive_added_WGSRPD3_n"] == 4
           and effect["source_list_observed_zero_added_WGSRPD3_n"] == 4,
           "8 selected endemics: 4 source-positive, 4 source-zero")
    e8 = d["endemic_8"]["source_endemic_resource"]["source_butterfly_rows"]
    ensure(len(e8) == 8 and sum(x["potential_global_added_region_units"] > 0 for x in e8) == 4,
           "Independent 8-endemic original-WCVP source rows")
    ensure(all(x["original_WCVP_native_plus_introduced_global_region_units"]
               - x["original_WCVP_native_global_region_units"]
               == x["potential_global_added_region_units"] for x in e8),
           "Per-endemic potential global plant-region arithmetic")
    ensure(d["endemic_8"]["origin_archipelago_not_individual_island"][
           "newly_host_present_origin_archipelago_region_units_created"] == 0,
           "Do not conflate extra global plant regions with new source-island host presence")
    field = d["field"]["primary_new_positive"]
    ensure(field["endemic_butterfly"] == "Pieris cheiranthi"
           and field["introduced_named_larval_host"] == "Tropaeolum majus",
           "Direct named introduced host field case")
    ensure(len(field["independently_observed_larvae"]) == 2
           and {x["island"] for x in field["independently_observed_larvae"]}
             == {"Tenerife","La Palma"},
           "Exact direct named-island records; garden is not wild")
    ensure(d["field"]["site_resident_versus_breeding"]["survival_and_recruitment_measured_in_current_audit"] is False,
           "Direct larval plant use must not imply successful recruitment")
    out["checks"]["endemic_larval_use"] = {
        "living_independent_butterflies": 41, "endemic": 18,
        "original_GEB239_endemic_overlap": 0,
        "named_host_species_endemics": 8, "source_positive_endemics": 4,
        "named_alien_host_larva_islands": ["Tenerife", "La Palma"],
        "warning": "Incomplete endemic host lists; living source populations unrepresentative; field larva is not colonization"}

    gbif = d["global_beta"]
    ensure(gbif["cohort"]["global_source_mapped_candidate_species"] == 224
           and gbif["cohort"]["region_count"] == 169
           and gbif["cohort"]["evaluated_region_pairs"] == 14196,
           "Independent larger GBIF assemblage endpoint denominator")
    models = gbif["heldout_region_model_scores"]
    m0, m1 = (models[k]["region_equal_MSE"]
              for k in ("M0_all_flora_effort_prior", "M1_known_host"))
    relative = 100 * (m0 - m1) / m0
    near(relative, gbif["heldout_region_model_comparisons"]["M1_known_host"]["MSE_reduction_vs_M0_percent"],
         "Global 224 indexed-butterfly assemblage gain arithmetic")
    near(relative, .1574407088114985, "Global butterfly community host composition gain")
    out["checks"]["independent_global_community"] = {
        "species": 224,"regions": 169,"pairs": 14196,
        "known_host_increment_pct": relative,
        "warning": "Contemporary plants, changing survey effort; NOT same outcome as 155-species forecast"}

    ensure(len(out["source_sha256"]) == len(INPUTS), "All evidence receipts available")
    out["status"] = "PASS"
    out["check_count"] = len(out["checks"])
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, help="Write deterministic source-receipt verification JSON")
    args = parser.parse_args()
    result = main()
    print(f"PASS {result['check_count']} independently scoped arithmetic/source-gate groups; "
          f"{len(result['source_sha256'])} original pinned file fingerprints.")
    for name, score in result["checks"].items():
        print(name + ": " + json.dumps(score, sort_keys=True, ensure_ascii=False))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
                               encoding="utf-8")
