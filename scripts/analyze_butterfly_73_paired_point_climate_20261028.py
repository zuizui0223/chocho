#!/usr/bin/env python3
"""Posthoc species-paired attrition audit of the source-frozen 73 GBIF point cases.

Only species x country aggregate point counts are used. The outcome is
WorldClim novelty relative to native larval-HOST-PLANT region climate, NOT
butterfly physiological niche evolution or establishment. This is not an
interchangeable-case randomisation, absence model, or ecological significance test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SOURCE_SCHEMA = "chocho_73_worldclim_case_counts_from_run_40435313884_v01"
OUTPUT_SCHEMA = "chocho_73_paired_point_climate_attrition_20261028_v01"
SOURCE_CSV_SHA256 = "e5d229bea32f93055d486d59a174bcca73ac8c7d6638e253d068238e5af764dc"
SOURCE_JSON_SHA256 = "5d3ca4137821e60aaca7e1bf965c509b1cddee76c643fb00ee15df44ae775937"
SOURCE_ARTIFACT_SHA256 = "a91cfd6ae4a1ef4c17c1a6e35b0d40eda078d312a69aa67cc6c1a17d164f9fa0"
ENDPOINTS = ("area95", "area90", "equal95")
DEFAULT_SOURCE = Path(__file__).resolve().parents[1] / "data/exploratory/butterfly_73_worldclim_case_counts_20261027.csv"


def read_cases(path: Path, verify_frozen: bool = True):
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if verify_frozen and digest != SOURCE_CSV_SHA256:
        raise ValueError(f"Frozen 73-case aggregate SHA256 mismatch: {digest}")
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    keys = set()
    cleaned = []
    integer_cols = ("n", "outside_area95", "outside_area90", "outside_equal95",
                    "strict_n", "strict_outside_area95", "strict_outside_area90",
                    "strict_outside_equal95", "source_added_region",
                    "source_not_added_region", "source_unmapped_region")
    for r in rows:
        key = (r["species"], r["country"])
        if not all(key) or key in keys or r["source_class"] not in ("novel", "analogue"):
            raise ValueError(f"Invalid or duplicate case: {key}")
        keys.add(key)
        d = {k: int(r[k]) for k in integer_cols}
        if d["n"] < 1 or d["strict_n"] < 0 or d["strict_n"] > d["n"]:
            raise ValueError(f"Invalid denominator: {key}")
        if sum(d[k] for k in ("source_added_region", "source_not_added_region", "source_unmapped_region")) != d["n"]:
            raise ValueError(f"Source geographic point categories do not partition case: {key}")
        for e in ENDPOINTS:
            if not 0 <= d[f"outside_{e}"] <= d["n"] or not 0 <= d[f"strict_outside_{e}"] <= d["strict_n"]:
                raise ValueError(f"Out-of-range climate count: {key}/{e}")
            if d[f"strict_outside_{e}"] > d[f"outside_{e}"]:
                raise ValueError(f"Strict outside subset exceeds full outside: {key}/{e}")
            if d["strict_n"]-d[f"strict_outside_{e}"] > d["n"]-d[f"outside_{e}"]:
                raise ValueError(f"Strict inside subset exceeds full inside: {key}/{e}")
        cleaned.append({"species": key[0], "country": key[1], "source_class": r["source_class"], **d})
    if verify_frozen:
        classes = Counter(r["source_class"] for r in cleaned)
        if len(cleaned) != 73 or classes != {"novel": 32, "analogue": 41} or \
                len({r["species"] for r in cleaned}) != 20 or sum(r["n"] for r in cleaned) != 6076:
            raise ValueError("Frozen source roster 73/32/41/20/6076 no longer matches")
    return cleaned, digest


def descriptive_bootstrap(differences, seed=20261028, iterations=20000):
    """An empirical-species resampling diagnostic, NOT a confidence interval for effects."""
    x = np.asarray(differences, dtype=float)
    if not len(x):
        return None
    gen = np.random.default_rng(seed)
    ix = gen.integers(0, len(x), size=(iterations, len(x)))
    return [float(v) for v in np.percentile(x[ix].mean(axis=1), [2.5, 97.5])]


def compare(cases, endpoint="area95", minimum_count=1, strict=False, boot_iterations=20000):
    if endpoint not in ENDPOINTS or minimum_count < 1:
        raise ValueError("Invalid comparison options")
    denom = "strict_n" if strict else "n"
    outside = f"strict_outside_{endpoint}" if strict else f"outside_{endpoint}"
    retained = [r for r in cases if r[denom] >= minimum_count]
    grouped = defaultdict(lambda: defaultdict(list))
    for r in retained:
        grouped[r["species"]][r["source_class"]].append(r)
    paired = []
    for species, classes in sorted(grouped.items()):
        if "novel" not in classes or "analogue" not in classes:
            continue
        rates = {}
        for group in ("novel", "analogue"):
            v = classes[group]
            rates[group] = {
                "n_countries": len(v),
                "n_point_records": sum(r[denom] for r in v),
                "equal_country_outside_share": sum(r[outside] / r[denom] for r in v) / len(v),
                "record_weighted_outside_share": sum(r[outside] for r in v) / sum(r[denom] for r in v),
            }
        paired.append({"species": species, "novel": rates["novel"], "analogue": rates["analogue"],
            "difference_equal_country": rates["novel"]["equal_country_outside_share"] - rates["analogue"]["equal_country_outside_share"],
            "difference_record_weighted": rates["novel"]["record_weighted_outside_share"] - rates["analogue"]["record_weighted_outside_share"]})
    d = [r["difference_equal_country"] for r in paired]
    positive = sum(z > 1e-12 for z in d)
    negative = sum(z < -1e-12 for z in d)
    ties = len(d) - positive - negative
    n = len(paired)
    by_class = {}
    for cl in ("novel", "analogue"):
        q = [r for r in retained if r["source_class"] == cl]
        by_class[cl] = {
            "cases": len(q), "points": sum(r[denom] for r in q),
            "equal_case_outside_share": sum(r[outside] / r[denom] for r in q) / len(q) if q else None,
            "record_weighted_outside_share": sum(r[outside] for r in q) / sum(r[denom] for r in q) if q else None,
        }
    return {
        "endpoint": endpoint, "minimum_points_per_species_country": minimum_count, "strict_coordinate_uncertainty_le1km": strict,
        "cases_retained": len(retained), "paired_species": n, "paired_species_countries_each_side_ge2": sum(r["novel"]["n_countries"] >= 2 and r["analogue"]["n_countries"] >= 2 for r in paired),
        "positive_species": positive, "negative_species": negative, "tied_species": ties,
        "mean_equal_species_difference": sum(d) / n if n else None,
        "median_equal_species_difference": float(np.median(d)) if n else None,
        "descriptive_species_resampling_interval": descriptive_bootstrap(d, iterations=boot_iterations),
        "by_class": by_class, "paired_species_details": paired,
        "inference_warning": "Posthoc and positive-records-only comparison. Class labels are constructed from native-HOST-PLANT climate geography; no exchangeability, independent butterfly-native climate baseline, true absence, or documented juvenile feeding. Intervals are descriptive empirical-species resampling ONLY; no causal/ecological p value.",
    }


def run(cases, digest):
    main = compare(cases, "area95", 1)
    sensitivity = [compare(cases, e, n, strict=s, boot_iterations=20000)
                   for e, n, s in (("area95", 3, False), ("area95", 5, False), ("area95", 10, False),
                                   ("area90", 1, False), ("equal95", 1, False),
                                   ("area95", 1, True), ("area95", 5, True))]
    largest = sorted(cases, key=lambda r: r["n"], reverse=True)[:2]
    assert [(z["species"], z["country"], z["n"]) for z in largest] == [
        ("Pseudozizeeria maha", "TW", 2220), ("Catopsilia pomona", "TW", 2198)]
    return {
        "schema": OUTPUT_SCHEMA,
        "source": {"input_csv_sha256": digest, "source_json_sha256": SOURCE_JSON_SHA256,
                   "source_artifact_zip_sha256": SOURCE_ARTIFACT_SHA256,
                   "github_run": "https://github.com/zuizui0223/chocho/actions/runs/40435313884",
                   "github_artifact_id": 12703686785, "original_source_schema": SOURCE_SCHEMA,
                   "records_are_gbif_indices_not_individuals": True},
        "frozen_cohort": {"species_country_cases": len(cases),
            "species": len({r["species"] for r in cases}),
            "novel_cases": sum(r["source_class"] == "novel" for r in cases),
            "analogue_cases": sum(r["source_class"] == "analogue" for r in cases),
            "indexed_point_records": sum(r["n"] for r in cases),
            "single_record_cases": sum(r["n"] == 1 for r in cases),
            "cases_ge5_records": sum(r["n"] >= 5 for r in cases),
            "top_two_record_counts": sum(r["n"] for r in largest),
            "top_two_record_share": sum(r["n"] for r in largest) / sum(r["n"] for r in cases)},
        "primary_posthoc": main, "sensitivity": sensitivity,
        "decision": "NOT_REALIZED_BUTTERFLY_NICHE_ENTRY: insufficient within-species point replication and no independent host-use, breeding, effort or geographic null. Do not update frozen GEB main manuscript.",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DEFAULT_SOURCE)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows, digest = read_cases(args.input)
    result = run(rows, digest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"cohort": result["frozen_cohort"], "main": {
        k: v for k, v in result["primary_posthoc"].items()
        if k in ("paired_species", "paired_species_countries_each_side_ge2", "positive_species", "negative_species", "tied_species", "mean_equal_species_difference")},
        "min5_paired_species": result["sensitivity"][1]["paired_species"],
        "decision": result["decision"]}, sort_keys=True))


if __name__ == "__main__":
    main()
