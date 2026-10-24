#!/usr/bin/env python3
"""Source-locked posthoc host *recording* redundancy at butterfly GADM1 scale.

Counts accepted WCVP candidate plant *species* with indexed GBIF records,
not viable host plants, plant abundance, host use or recruitment.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

HOST_SCHEMA = "chocho_all26_gbif_plant_gadm1_feasibility_20261023_v1"
POS_SCHEMA = "chocho_subcountry_48_positive_GBIF_GADM1_lifestage_facet_v1"
STAGE_SCHEMA = "chocho_juvenile_conditioned_GADM1_source_facet_20261021_v1"
SOURCE_SCHEMA = "chocho_host_plant_candidates_5_juvenile_added_country_cases_20261021_v1"
SCHEMA = "gbif_all26_local_indexed_host_portfolio_redundancy_posthoc_v1"
THRESHOLDS = (1, 2, 3, 5)


def degree_bucket(n):
    return "2plus" if n >= 2 else str(n)


def _empty_hist():
    return {"0": 0, "1": 0, "2plus": 0}


def fixed_sources(plants, butterflies, stages, original_five):
    if plants.get("schema") != HOST_SCHEMA or plants.get("status") != "ALL26_BOTANICAL_GADM_QUERIES_COMPLETE":
        raise ValueError("Original completed 26-cell botanical source unavailable")
    if butterflies.get("schema") != POS_SCHEMA or butterflies.get("status") != "COMPLETE_48":
        raise ValueError("Original completed 48-cell butterfly GADM source unavailable")
    if stages.get("schema") != STAGE_SCHEMA or stages.get("status") != "COMPLETE_CONDITIONAL_JUVENILE_GADM1":
        raise ValueError("Original GADM-matched immature-stage source unavailable")
    if original_five.get("schema") != SOURCE_SCHEMA or original_five.get("status") != "SOURCE_VERIFIED_CANDIDATE_PLANTS_NOT_LOCAL_HOST_USE":
        raise ValueError("Five-case source-locked botanical selection missing")
    if (plants["added_butterfly_country_cases"],
        plants["plant_country_GADM_facets_complete"],
        plants["indexed_plant_source_records"],
        plants["taxa_exact_accepted_in_GBIF"]) != (26, 44, 3595, 24):
        raise ValueError("Original country×plant source margins differ")
    if (len(butterflies["observations"]), len(stages["cell_gadm1_juvenile_summary"]),
        len(original_five["focal_five_cases"])) != (48, 11, 5):
        raise ValueError("Original butterfly and juvenile cohorts differ")
    byplant = {}
    for p in plants["plant_country_results"]:
        key = (p["country"], p["plant_name"])
        if key in byplant or p["status"] != "COMPLETE_INDEXED_PLANT_FACET":
            raise ValueError("Repeated or incomplete plant-country GBIF facet")
        if sum(p["gadm1_counts"].values()) != p["gadm1_mapped_records"]:
            raise ValueError("Plant GADM1 mapped counts drifted")
        if p["gadm1_mapped_records"] + p["gadm1_unmapped_records"] != p["indexed_country_plant_records"]:
            raise ValueError("Plant country GADM and country totals drifted")
        byplant[key] = p
    if len(byplant) != 44:
        raise ValueError("Original 44 plant-country query source changed")
    bybutterfly = {(r["species"], r["country"]): r for r in butterflies["observations"]}
    if len(bybutterfly) != 48:
        raise ValueError("Repeated butterfly country source pair")
    byjuvenile = {(r["species"], r["country"]): r for r in stages["cell_gadm1_juvenile_summary"]}
    if len(byjuvenile) != 11:
        raise ValueError("Repeated juvenile country source pair")
    five = {(r["species"], r["country"]) for r in original_five["focal_five_cases"]}
    if len(five) != 5 or not five <= set(byjuvenile):
        raise ValueError("Original five preselected juvenile cases not reproduced")
    return byplant, bybutterfly, byjuvenile, five


def summary_for(plants, butterflies, stages, original_five):
    byplant, bybutterfly, byjuvenile, five = fixed_sources(
        plants, butterflies, stages, original_five)
    case_rows = []
    thresholds = {
        str(t): {"butterfly_record_degree": _empty_hist(),
                 "butterfly_gadm_unit_degree": _empty_hist(),
                 "multi_candidate_record_degree": _empty_hist(),
                 "multi_candidate_gadm_unit_degree": _empty_hist()}
        for t in THRESHOLDS
    }
    stage_hist = _empty_hist()
    stage_case_rows = []
    byspcountry = {}
    plant_links = 0
    for pcase in plants["butterfly_country_results"]:
        sp, country = pcase["butterfly"], pcase["country"]
        key = (sp, country)
        if key in byspcountry or key not in bybutterfly:
            raise ValueError("Missing or repeated source butterfly-country")
        original = bybutterfly[key]
        if original["added_resource_country"] != 1:
            raise ValueError("All26 plant source includes a non-added butterfly-country")
        bcounts = original["facets"]["gadm1"]["categories"]
        if sum(bcounts.values()) != pcase["all_butterfly_mapped_gadm1_records"]:
            raise ValueError("Mapped GBIF butterfly counts disagree between sources")
        candidates = [p for p in byplant.values()
                      if p["country"] == country and sp in p["linked_butterflies"]]
        if len(candidates) != pcase["source_candidate_host_species"] or len({
            p["plant_name"] for p in candidates}) != len(candidates):
            raise ValueError("Original butterfly host-network list drifted")
        plant_links += len(candidates)
        multi = len(candidates) >= 2
        local = {str(t): {"records": _empty_hist(), "units": _empty_hist()}
                 for t in THRESHOLDS}
        for t in THRESHOLDS:
            for region, bcount in bcounts.items():
                if bcount <= 0:
                    raise ValueError("Zero/negative indexed butterfly GADM record")
                degree = sum(p["gadm1_counts"].get(region, 0) >= t
                             for p in candidates)
                bucket = degree_bucket(degree)
                local[str(t)]["records"][bucket] += bcount
                local[str(t)]["units"][bucket] += 1
                thresholds[str(t)]["butterfly_record_degree"][bucket] += bcount
                thresholds[str(t)]["butterfly_gadm_unit_degree"][bucket] += 1
                if multi:
                    thresholds[str(t)]["multi_candidate_record_degree"][bucket] += bcount
                    thresholds[str(t)]["multi_candidate_gadm_unit_degree"][bucket] += 1
        base = local["1"]["records"]
        coverage = base["1"] + base["2plus"]
        if (coverage != pcase["butterfly_records_in_candidate_plant_positive_gadm1"]
            or base["0"] + coverage != pcase["all_butterfly_mapped_gadm1_records"]
            or local["1"]["units"]["1"] + local["1"]["units"]["2plus"]
               != pcase["butterfly_gadm1_units_with_candidate_plant_records"]):
            raise ValueError("Source 26-cell GBIF co-detection margins do not replay")
        removal = []
        for p in candidates:
            unique_loss = 0
            for region, count in bcounts.items():
                if p["gadm1_counts"].get(region, 0) <= 0:
                    continue
                if not any(other["gadm1_counts"].get(region, 0) > 0
                           for other in candidates if other["plant_name"] != p["plant_name"]):
                    unique_loss += count
            removal.append({"plant_name": p["plant_name"],
                            "indexed_butterfly_coverage_lost_if_only_this_plant_index_removed": unique_loss})
        removal.sort(key=lambda row: (-row["indexed_butterfly_coverage_lost_if_only_this_plant_index_removed"], row["plant_name"]))
        entry = {
            "species": sp, "country": country,
            "catalogue_candidate_host_species": len(candidates),
            "GBIF_plant_names_with_any_country_record":
                sum(p["indexed_country_plant_records"] > 0 for p in candidates),
            "butterfly_GADM1_record_count": sum(bcounts.values()),
            "source_GADM1_units": len(bcounts),
            "indexed_butterfly_records_in_at_least_one_candidate_plant_GADM1": coverage,
            "indexed_butterfly_records_in_two_or_more_candidate_species_GADM1": base["2plus"],
            "primary_gadm_plant_record_degrees": local["1"],
            "all_thresholds": local,
            "largest_index_coverage_loss_after_removing_one_source_plant":
                removal[0] if removal else None,
            "plant_name_removal_sensitivities": removal if multi else [],
            "single_plant_source_removal_is_not_botanical_extirpation": True,
        }
        case_rows.append(entry)
        byspcountry[key] = (candidates, bcounts)
    if len(case_rows) != 26 or plant_links != 46:
        raise ValueError("Expected 26 butterfly-countries and 46 source plant links")
    if sum(case["catalogue_candidate_host_species"] >= 2 for case in case_rows) != 7:
        raise ValueError("Source candidate host multiplicity changed")
    if sum(case["catalogue_candidate_host_species"] == 1 for case in case_rows) != 19:
        raise ValueError("Original nineteen single-host source cases changed")
    primary = thresholds["1"]
    if (primary["butterfly_record_degree"],
        primary["butterfly_gadm_unit_degree"]) != (
        {"0": 222, "1": 226, "2plus": 5},
        {"0": 66, "1": 38, "2plus": 2}):
        raise ValueError("Source 453/231/106/40 spatial marginals not reproduced")
    stage_total = 0
    for key in sorted(five):
        st = byjuvenile[key]
        candidates, all_butterfly = byspcountry[key]
        hist = _empty_hist()
        for gadm, n in st["juvenile_gadm1_counts"].items():
            if n > all_butterfly.get(gadm, 0):
                raise ValueError("Juvenile-label records exceed total focal butterfly record count")
            d = sum(p["gadm1_counts"].get(gadm, 0) >= 1 for p in candidates)
            hist[degree_bucket(d)] += n
            stage_hist[degree_bucket(d)] += n
            stage_total += n
        stage_case_rows.append({"species": key[0], "country": key[1],
                                "juvenile_source_label_count": sum(hist.values()),
                                "juvenile_label_candidate_plant_degrees": hist})
    if stage_total != 24 or stage_hist != {"0": 8, "1": 14, "2plus": 2}:
        raise ValueError("Source-locked five-case 24 juvenile-label GADM records changed")
    return {
        "schema": SCHEMA,
        "status": "EXACT_SOURCE_DESCRIPTIVE_POSTHOC_NOT_HOST_USE",
        "original_26_cases": 26,
        "original_44_plant_country_queries": 44,
        "original_46_butterfly_country_plant_links": plant_links,
        "source_original_total_plant_index_records": plants["indexed_plant_source_records"],
        "source_original_source_plant_taxa": plants["taxa_exact_accepted_in_GBIF"],
        "source_original_mapped_butterfly_record_count": 453,
        "source_original_mapped_butterfly_gadm_units": 106,
        "case_counts": {
            "single_source_candidate_host": 19, "multiple_source_candidate_hosts": 7
        },
        "minimum_plant_records_per_species_per_GADM1_sensitivity": thresholds,
        "five_source_juvenile_cases": {
            "cases": 5, "source_stage_label_count": stage_total,
            "stage_label_candidate_plant_species_degree_records": stage_hist,
            "per_case": stage_case_rows
        },
        "case_plant_recording_portfolios": case_rows,
        "inference": (
            "Plant labels count species with GBIF indexed GADM1 records, not real host "
            "abundance, feeding, larval survival, functional redundancy, or butterfly recruitment. "
            "All source-positive butterfly cases are selected on GBIF presence, and plant "
            "GBIF coverage is incomplete; no ecological p-values are computed."),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--all26-plants", type=Path, required=True)
    p.add_argument("--all48-butterflies", type=Path, required=True)
    p.add_argument("--juvenile-gadm", type=Path, required=True)
    p.add_argument("--original-five", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    result = summary_for(
        json.loads(a.all26_plants.read_text(encoding="utf-8")),
        json.loads(a.all48_butterflies.read_text(encoding="utf-8")),
        json.loads(a.juvenile_gadm.read_text(encoding="utf-8")),
        json.loads(a.original_five.read_text(encoding="utf-8")))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("GBIF_LOCAL_HOST_REDUNDANCY_POSTHOC=" + json.dumps({
        k:v for k,v in result.items() if k != "case_plant_recording_portfolios"},
        separators=(",", ":")), flush=True)


if __name__ == "__main__":
    main()
