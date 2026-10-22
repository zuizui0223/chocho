#!/usr/bin/env python3
"""Audit source-exact candidate-plant/juvenile-label co-detection within GADM1.

This *posthoc measurement diagnostic* never treats GBIF non-detection as
biological absence. Juvenile labels may be wrong, duplicated, or biased by
recording effort. Conditional hypergeometric assignment is a REFERENCE,
not a biologically valid p-value or population effect.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path

SCHEMA = "chocho_five_country_gadm_host_juvenile_reference_posthoc_v1"
INPUT_SCHEMA = "chocho_strict_plant_taxon_gadm_cocoverage_five_juvenile_country_cases_20261021_v1"
STAGE_SCHEMA = "chocho_juvenile_conditioned_GADM1_source_facet_20261021_v1"
EXPECTED_CASES = {
    ("Anteos clorinde", "CU"), ("Colobura dirce", "VE"),
    ("Papilio demoleus", "CU"), ("Papilio demoleus", "DO"),
    ("Papilio demoleus", "HT")
}


def exact_reference(rows):
    """Convolve conditional hypergeometrics with rational, exact masses."""
    pmf = [Fraction(1)]
    for r in rows:
        N, H, j = (r[k] for k in ("all_butterfly_gadm1_records",
                                 "all_butterfly_records_in_plant_gadm1",
                                 "immature_source_labels"))
        nxt = [Fraction(0)] * (len(pmf) + min(H, j))
        denom = comb(N, j)
        if denom <= 0 or j > N or H > N:
            raise ValueError("Invalid conditional source margins")
        for x in range(max(0, j - (N - H)), min(H, j) + 1):
            mass = Fraction(comb(H, x) * comb(N - H, j - x), denom)
            for z, prior in enumerate(pmf):
                nxt[z+x] += prior * mass
        pmf = nxt
    if sum(pmf) != 1:
        raise ValueError("Exact label-reassignment reference loses probability mass")
    return pmf


def summarize_reference(rows):
    j = sum(r["immature_source_labels"] for r in rows)
    o = sum(r["immature_labels_in_plant_detected_gadm1"] for r in rows)
    expectation = sum((
        Fraction(r["immature_source_labels"] * r["all_butterfly_records_in_plant_gadm1"],
                 r["all_butterfly_gadm1_records"])
        for r in rows
    ), Fraction(0))
    distribution = exact_reference(rows)
    return {
        "juvenile_label_count": j, "in_plant_detected_gadm": o,
        "null_expected_labels": round(float(expectation), 10),
        "observed_minus_reference_pp": round(100 * float((o-expectation)/j), 7) if j else None,
        "conditional_finite_label_reassignment_upper_tail_NOT_ecological_p": (
            round(float(sum(distribution[o:])), 10) if j else None),
        "reference_design": "Reassign juvenile source-stage labels without replacement within each species-country, keeping exact original all-butterfly GADM1 records and candidate-plant-positive GADM1 distribution; ignores observer heterogeneity, correlated records, and stage-annotation bias.",
    }


def analyse(source, stage, all_stage):
    if source.get("schema") != INPUT_SCHEMA or source.get("status") != "ALL_FIVE_CASES_SOURCE_ASSESSED":
        raise ValueError("Not complete fixed five source candidate host dataset")
    if (source.get("total_taxon_names"),source.get("taxon_resolved"),
        source.get("candidate_links"),source.get("plant_facet_responses")) != (8, 8, 20, 20):
        raise ValueError("Eight original accepted botanical species / 20 source host links not complete")
    if stage.get("schema") != STAGE_SCHEMA or stage.get("status") != "COMPLETE_CONDITIONAL_JUVENILE_GADM1":
        raise ValueError("Frozen juvenile GADM facets not complete")
    if all_stage.get("status") != "COMPLETE_48" or len(all_stage.get("observations", [])) != 48:
        raise ValueError("Original 48 source positive butterfly-country cells incomplete")
    source_keys = {}
    for k,t in source["exact_taxonomy"].items():
        if t["status"] != "EXACT_ACCEPTED_GBIF_PLANT" or t["plant_name"] != k or not isinstance(t["accepted_species_key"], int):
            raise ValueError("Ambiguous/changed original plant taxonomy")
        source_keys[k] = t["accepted_species_key"]
    if len(source_keys) != 8:
        raise ValueError("Eight source plant identities expected")
    source_pairs = {}
    positive_plant_links = 0
    total_plant_records = 0
    for item in source["per_source_host_country"]:
        key = item["species"], item["country"], item["plant_name"]
        if key in source_pairs or item["status"] != "COMPLETE_INDEXED_PLANT_FACET":
            raise ValueError("Duplicated/missing botanical GADM1 facet")
        if source_keys[item["plant_name"]] != item["accepted_species_key"]:
            raise ValueError("Plant taxonomy does not match frozen GBIF accepted key")
        counts = item["gadm1_counts"]
        n = item["indexed_country_plant_records"]
        if (n < 0 or sum(counts.values()) > n or
            len(counts) != item["gadm1_count"] or
            sum(counts.values()) != item["gadm1_mapped_records"]):
            raise ValueError("Invalid plant GBIF GADM1 counts")
        source_pairs[key] = item
        positive_plant_links += int(n > 0)
        total_plant_records += n
    if len(source_pairs) != 20:
        raise ValueError("Not all 20 fixed botanical country links were queried")
    stage_cases = {(r["species"], r["country"]): r for r in stage["cell_gadm1_juvenile_summary"]}
    all_cases = {(r["species"], r["country"]): r for r in all_stage["observations"]}
    if len(stage_cases) != 11 or len(all_cases) != 48:
        raise ValueError("Original positive-cell source identities drifted")
    cases, status_rows = [], []
    seen = set()
    for case in source["focal_case_results"]:
        key = case["species"], case["country"]
        if key not in EXPECTED_CASES or key in seen or key not in stage_cases or key not in all_cases:
            raise ValueError("Case selection drifted or outcome-selected")
        seen.add(key)
        jsource = stage_cases[key]["juvenile_gadm1_counts"]
        asource = all_cases[key]["facets"]["gadm1"]["categories"]
        plant_links = [v for k,v in source_pairs.items() if k[:2] == key]
        if len(plant_links) != case["source_accepted_host_candidates"]:
            raise ValueError("Candidate plants count mismatch against source")
        plant_gadm = Counter()
        plant_summary = []
        for link in plant_links:
            plant_gadm.update(link["gadm1_counts"])
            matched_juv_gadm = sorted(set(link["gadm1_counts"]) & set(jsource))
            plant_summary.append({
                "plant": link["plant_name"], "taxon_key": link["accepted_species_key"],
                "source_country_records": link["indexed_country_plant_records"],
                "country_gadm1_units_detected": link["gadm1_count"],
                "juvenile_source_gadm1_units_with_this_plant": matched_juv_gadm,
            })
        if stage_cases[key]["added_host_country"] != 1 or all_cases[key]["added_resource_country"] != 1:
            raise ValueError("Fixed cases no longer originally added-resource")
        if jsource != case["juvenile_labeled_focal_gadm1_counts"]:
            raise ValueError("GBIF juvenile data source mismatched")
        N = sum(asource.values())
        j = sum(jsource.values())
        H = sum(v for unit,v in asource.items() if plant_gadm.get(unit,0)>0)
        O = sum(v for unit,v in jsource.items() if plant_gadm.get(unit,0)>0)
        for unit, n in jsource.items():
            if n > asource.get(unit,0):
                raise ValueError("Juvenile GADM1 counts cannot exceed same all-stage butterfly observations")
        if (N != case["all_focal_butterfly_record_count"] or
            H != case["all_focal_records_in_candidate_plant_gadm"] or
            j != case["juvenile_records_all"] or
            O != case["juvenile_labeled_records_in_plant_detected_gadm_units"] or
            len(set(plant_gadm)&set(jsource)) != case["juvenile_labeled_gadm_units_with_plant_detection"]):
            raise ValueError("Exact source case/country co-detection was not reconstructed")
        cases.append({
            "butterfly":key[0], "country":key[1],
            "candidate_host_species":len(plant_links),
            "with_any_indexed_country_plant_observation":sum(v["indexed_country_plant_records"]>0 for v in plant_links),
            "all_butterfly_gadm1_records":N,
            "all_butterfly_records_in_plant_gadm1":H,
            "immature_source_labels":j,
            "immature_labels_in_plant_detected_gadm1":O,
            "source_immature_gadm1_units":len(jsource),
            "immature_gadm1_units_with_candidate_plant_observations":len(set(jsource)&set(plant_gadm)),
            "same_administrative_units_NOT_same_site_or_feeding":True,
            "candidate_plant_records":plant_summary,
        })
    if seen != EXPECTED_CASES:
        raise ValueError("Missing pre-existing selected host cases")
    if (positive_plant_links,total_plant_records)!=(9,38):
        raise ValueError("Original source plant-record totals unexpectedly changed")
    if (sum(r["immature_source_labels"] for r in cases),
        sum(r["immature_labels_in_plant_detected_gadm1"] for r in cases))!=(24,16):
        raise ValueError("Fixed juvenile labels and plant GADM overlap changed")
    all_ref = summarize_reference(cases)
    non_dom = [r for r in cases if (r["butterfly"],r["country"]) != ("Papilio demoleus","DO")]
    details = {"all_five": all_ref,
               "excluding_papilio_dominican_republic": summarize_reference(non_dom)}
    return {
        "schema":SCHEMA,
        "status":"COMPLETE_SOURCE_RECONSTRUCTION_INSUFFICIENT_BIOLOGICAL_EVIDENCE",
        "exact_gbif_accepted_botanical_species":len(source_keys),
        "source_candidate_host_country_links":len(source_pairs),
        "plant_country_links_with_at_least_one_indexed_record":positive_plant_links,
        "source_indexed_botanical_records_across_20_links":total_plant_records,
        "preselected_juvenile_positive_country_cases":len(cases),
        "positive_cases_with_administrative_plant_co_detection":sum(r["immature_labels_in_plant_detected_gadm1"]>0 for r in cases),
        "source_juvenile_label_count_in_five_cases":24,
        "source_juvenile_labels_in_candidate_plant_positive_gadm1":16,
        "juvenile_source_admin_units_in_five_cases":sum(r["source_immature_gadm1_units"] for r in cases),
        "juvenile_source_admin_units_with_plant_records":sum(r["immature_gadm1_units_with_candidate_plant_observations"] for r in cases),
        "source_status_warning": "Only 38 botanical GBIF records for 20 plant-country links; missing indexed plant records are not proof of plant absence.",
        "posthoc_label_allocation_reference": details,
        "cases": cases,
        "interpretation": "Administrative GADM1 co-detection only; source stage names are unverified and the null is descriptive, not a demographic host effect or valid independent-sample p-value.",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plant-source", required=True, type=Path)
    p.add_argument("--juvenile-source", required=True, type=Path)
    p.add_argument("--all-stage-source", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args=p.parse_args()
    plant_raw=args.plant_source.read_bytes()
    stage_raw=args.juvenile_source.read_bytes()
    all_raw=args.all_stage_source.read_bytes()
    result=analyse(json.loads(plant_raw),json.loads(stage_raw),json.loads(all_raw))
    result["source_sha256"]={
        "gbif_20_host_country_facets":hashlib.sha256(plant_raw).hexdigest(),
        "gbif_15_juvenile_stage_facets":hashlib.sha256(stage_raw).hexdigest(),
        "gbif_48_all_stage_positive_cells":hashlib.sha256(all_raw).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GBIF_FIVE_HOST_GADM_REFERENCE="+json.dumps(result,separators=(",",":")),flush=True)
    print("GBIF_FIVE_HOST_GADM_ARCHIVED_SOURCE="+json.dumps(json.loads(plant_raw),separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
