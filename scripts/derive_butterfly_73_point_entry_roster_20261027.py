#!/usr/bin/env python3
"""Rebuild exact 73 PR72 GBIF-positive novel/analog species-country cohort.

No new butterfly record bodies or point outcomes are opened here. Original
five Git blob identities and PR72 classifier are source locked.
"""
from __future__ import annotations

import argparse
import gzip
import json
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from analyze_butterfly_niche_entry_country_20261026 import (
    PRIMARY, check_source, load_areas, opportunity_class,
)
from analyze_gbif_country_host_opportunity_20261020 import load_observations

SCHEMA = "chocho_pr72_exact73_point_followup_roster_20261027_v01"


def select_roster(rows, counts, complete, source_keys, *, parameters=PRIMARY):
    """Return all prespecified positive cases, and all source strata margins."""
    basis, window = parameters["basis"], parameters["window"]
    allowed = set(complete[basis])
    indexed = counts[basis, window]
    effort = defaultdict(set)
    for (species, country), count in indexed.items():
        if species in allowed and count > 0:
            effort[country].add(species)

    strata = defaultdict(list)
    audit = Counter()
    for r in rows:
        species, country = r["species"], r["country"]
        if species not in allowed or r["native"]:
            continue
        kind = opportunity_class(r["added"], r["analog"],
                                 min_area=parameters["min_area"],
                                 purity=parameters["purity"])
        if kind not in {"analogue", "novel"}:
            continue
        if len(effort[country] - {species}) < parameters["effort_gate"]:
            continue
        original_count = indexed.get((species, country), 0)
        strata[species, r["level1"]].append({
            "species": species, "country": country, "level1": r["level1"],
            "resource_class": kind, "original_gbif_record_count": original_count,
            "added_host_country_area_fraction": r["added"],
            "analog_host_country_area_fraction": r["analog"],
            "novel_host_country_area_fraction": r["novel"],
        })
        audit["candidate_source_cells"] += 1
    paired = []
    source_strata = []
    for (species, level1), cells in sorted(strata.items()):
        n = len(cells)
        new = sum(c["resource_class"] == "novel" for c in cells)
        if not 0 < new < n:
            continue
        positives = [c for c in cells if c["original_gbif_record_count"] >=
                     parameters["min_records"]]
        k = len(positives)
        observed = sum(c["resource_class"] == "novel" for c in positives)
        source_strata.append({
            "species": species, "level1": level1,
            "N": n, "A": new, "K": k, "O": observed,
            "expected": new * k / n,
        })
        for c in positives:
            if species not in source_keys or not isinstance(source_keys[species], int):
                raise ValueError("Positive taxon name missing exact source accepted GBIF key")
            paired.append({
                **c,
                "gbif_accepted_species_key": source_keys[species],
            })
    paired.sort(key=lambda c: (c["species"], c["country"]))
    if len({(x["species"], x["country"]) for x in paired}) != len(paired):
        raise ValueError("Duplicate positive species-country source pair")
    return paired, source_strata, dict(audit)


def summarize(roster, strata):
    classes = Counter(r["resource_class"] for r in roster)
    by_species = Counter(r["species"] for r in roster)
    positive_count = sum(r["original_gbif_record_count"] for r in roster)
    maxima = sorted(roster, key=lambda r: -r["original_gbif_record_count"])
    obs = sum(r["O"] for r in strata)
    exp = sum(r["expected"] for r in strata)
    total_K = sum(r["K"] for r in strata)
    info = sum(0 < r["K"] < r["N"] for r in strata)
    result = {
        "positive_country_species_cases": len(roster),
        "novel_resource_country_cases": classes["novel"],
        "analog_resource_country_cases": classes["analogue"],
        "distinct_butterfly_species": len(by_species),
        "paired_strata": len(strata),
        "informative_strata": info,
        "paired_source_eligible_country_species_candidates": sum(r["N"] for r in strata),
        "paired_source_eligible_novel_candidates": sum(r["A"] for r in strata),
        "original_GBIF_indexed_record_total_in_positive_cases": positive_count,
        "source_positive_cases_above_cap_5000": sum(
            r["original_gbif_record_count"] > 5000 for r in roster),
        "source_positive_records_above_cap_5000": sum(
            r["original_gbif_record_count"] for r in roster
            if r["original_gbif_record_count"] > 5000),
        "positive_case_record_count_quantiles": {
            "min": min(r["original_gbif_record_count"] for r in roster),
            "median": sorted(r["original_gbif_record_count"] for r in roster)[len(roster)//2],
            "max": maxima[0]["original_gbif_record_count"],
        },
        "positive_cases_with_1_record": sum(r["original_gbif_record_count"] == 1
                                            for r in roster),
        "positive_cases_with_ge5_records": sum(r["original_gbif_record_count"] >= 5
                                              for r in roster),
        "original_country_novel_obs_minus_null_pp": 100 * (obs - exp) / total_K,
        "max_source_case": {
            "butterfly": maxima[0]["species"], "country": maxima[0]["country"],
            "index_records": maxima[0]["original_gbif_record_count"],
        },
        "by_resource_class": {
            key: {
                "positive_species_country_cases": classes[key],
                "original_index_records": sum(
                    r["original_gbif_record_count"] for r in roster
                    if r["resource_class"] == key),
                "distinct_species": len({
                    r["species"] for r in roster if r["resource_class"] == key
                })
            } for key in ("novel", "analogue")
        },
    }
    return result


def run(args):
    paths = {
        "areas": args.areas, "design": args.design, "counts": args.counts,
        "manifest": args.original_manifest, "receipt": args.country_receipt
    }
    for name, path in paths.items():
        check_source(path, name)
    areas, meta = load_areas(args.areas, args.design)
    with tempfile.TemporaryDirectory(prefix="chocho-pr72-cohort-") as temp:
        f = Path(temp) / "source_gbif_receipt.json"
        with gzip.open(args.country_receipt, "rt", encoding="utf-8") as src:
            f.write_text(src.read(), encoding="utf-8")
        counts, complete, original, full_receipt, excluded_alias = load_observations(
            args.counts, f, args.original_manifest
        )
    if sorted(excluded_alias) != ["Oeneis bore", "Oeneis polixenes"]:
        raise ValueError("Original accepted-key collision set changed")
    if len(complete["all_accepted_bases"]) != 229:
        raise ValueError("Country accepted-bases cohort changed")
    manifest = json.loads(args.original_manifest.read_text(encoding="utf-8"))
    taxa = {}
    for rec in manifest["species_receipts"]:
        if rec["status"] == "TAXON_UNRESOLVED":
            continue
        taxon = rec.get("taxon") or {}
        if not isinstance(taxon.get("key"), int):
            raise ValueError("Source accepted butterfly key missing")
        taxa[rec["species"]] = taxon["key"]
    roster, strata, extra = select_roster(areas, counts, complete, taxa)
    summary = summarize(roster, strata)
    assert (
        summary["positive_country_species_cases"],
        summary["novel_resource_country_cases"],
        summary["analog_resource_country_cases"],
        summary["distinct_butterfly_species"],
        summary["paired_source_eligible_country_species_candidates"],
        summary["paired_source_eligible_novel_candidates"],
        summary["informative_strata"]
    ) == (73, 32, 41, 20, 2646, 1167, 20), summary
    if abs(summary["original_country_novel_obs_minus_null_pp"]
           - (-3.1847519478462263)) > 1e-8:
        raise ValueError("PR72 point-roster original country-null result not recovered")
    out = {
        "schema": SCHEMA,
        "status": "SOURCE_IDENTITY_FROZEN_BEFORE_NEW_POINT_OBSERVATION",
        "source": meta,
        "GBIF_origin": "frozen original 2010-2017 and 2018-2025 country indexed facets; not a DOI download",
        "selection": PRIMARY,
        "summary": summary,
        "original_five_git_blobs": {
            k: __import__("analyze_butterfly_niche_entry_country_20261026").BLOBS[k]
            for k in paths
        },
        "roster": roster,
        "inferential_gates": [
            "73 means positive species-country source units, not 73 independent insects.",
            "Source country class describes added HOST-PLANT resource area, not each butterfly point climate.",
            "No source-only case can be classified as newly established wild population.",
            "Do not convert GBIF zero point queries to ecological absences or observed niche conservatism.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "source_status": out["status"], "source_summary": summary,
        "source_total_candidates": extra["candidate_source_cells"],
        "roster_cases": len(roster),
    }, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--areas", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--original-manifest", type=Path, required=True)
    parser.add_argument("--country-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
