#!/usr/bin/env python3
"""Source-locked early botanical observation support for fixed 51-island pairs.

The number of all plant occurrence records is a coverage proxy only:
overlapping 5-km harbor polygons can repeat GBIF occurrences, and counts
are neither unique taxa nor evidence that any particular host is absent.
Outputs describe comparability, not biological host effects.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = Path("provenance/reviewer_defenses/results")
SOURCES = {
    "availability": ("island_dated_host_match_20261030/island_port_early_plant_availability.csv",
                     "0a01831056c2cd59a4bd41bd1dbd4c84d77cf922"),
    "pairs": ("island_host_network_20261030/island_host_network_pairs.csv",
              "5663d2ffe42419f5c4496b92c9f3bbf76175907b"),
    "complete": ("island_dated_host_match_20261030/island_port_early_inventory_islands.csv",
                 "143f9956ce4f871c7a97b4777ab69d416963c2bb"),
    "links": ("island_dated_host_match_20261030/island_dated_host_match_links.csv",
              "88db773c74dea9cd2d58bde0b67a394d93226869"),
    "network": ("island_host_network_20261030/island_host_network_islands.csv",
                "18cecde5bf3066f276380167f451f1fa5cc1caa1"),
}
THRESHOLDS = (0, 1, 10, 100)


def source_rows(root: Path, key: str):
    rel, expected = SOURCES[key]
    p = root / BASE / rel
    content = p.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()
    if blob != expected:
        raise ValueError(f"Source content changed from pinned snapshot: {key}")
    rows = list(csv.DictReader(io.StringIO(content.decode("utf-8"), newline="")))
    if not rows:
        raise ValueError(f"Empty source: {key}")
    return rows


def numeric_int(s, key: str) -> int:
    if not s.isdecimal():
        raise ValueError(f"{key} must be integer: {s}")
    return int(s)


def unordered_pair(a: str, b: str):
    if not a or not b or a == b:
        raise ValueError(f"Invalid distinct island pair: {a}, {b}")
    return tuple(sorted((a, b), key=lambda v: int(v)))


def coverage_by_island(availability):
    groups = {}
    seen = set()
    for row in availability:
        pid, lid = row["port_id"], row["landmass_id"]
        if pid in seen:
            raise ValueError(f"Duplicate port id: {pid}")
        seen.add(pid)
        early = numeric_int(row["early_plant_records"], "early")
        total = numeric_int(row["all_plant_records"], "all")
        if early > total:
            raise ValueError(f"Earlier count larger than all-period count: {pid}")
        val = groups.setdefault(lid, {
            "landmass_id": lid, "ports": 0, "early_positive_ports": 0,
            "early_max_per_port": 0, "early_port_record_sum": 0,
            "all_port_record_sum": 0,
        })
        val["ports"] += 1
        val["early_positive_ports"] += int(early > 0)
        val["early_max_per_port"] = max(val["early_max_per_port"], early)
        val["early_port_record_sum"] += early
        val["all_port_record_sum"] += total
    return groups


def staged_counts(analytical_pairs, threshold):
    selected = [row for row in analytical_pairs if row["min_of_island_early_max"] >= threshold]
    n = sum(row["shared_new_species"] for row in selected)
    m = sum(row["early_matched_butterfly_species"] for row in selected)
    if m > n:
        raise ValueError("Host matches exceed shared-new butterflies")
    return {
        "min_indexed_early_plant_records_per_island_port_max": threshold,
        "island_pairs": len(selected),
        "island_pairs_with_shared_new_butterfly_species": sum(row["shared_new_species"] > 0
                                                              for row in selected),
        "shared_new_butterfly_species_island_pairs": n,
        "historically_recorded_matching_host_species_island_pairs": m,
        "match_fraction_when_denominator_positive": m/n if n else None,
        "islands": len({a for row in selected for a in (row["island_a"], row["island_b"])}),
    }


def audit(root: Path = ROOT):
    availability = source_rows(root, "availability")
    pairs = source_rows(root, "pairs")
    complete = source_rows(root, "complete")
    links = source_rows(root, "links")
    network = source_rows(root, "network")
    if (len(availability), len(pairs), len(complete), len(links), len(network)) != (194, 1275, 51, 5, 103):
        raise ValueError("An input cohort changed")
    coverage = coverage_by_island(availability)
    complete_islands = {r["landmass_id"] for r in complete if r["all_ports_complete"] == "True"}
    eligible_islands = {r["landmass_id"] for r in network if r["eligible"] == "True"}
    if len(complete_islands) != 51 or complete_islands != eligible_islands or set(coverage) != eligible_islands:
        raise ValueError("Availability inventory does not span the original eligible 51 islands")
    for r in complete:
        g = coverage[r["landmass_id"]]
        if numeric_int(r["planned_ports"], "planned") != g["ports"] or (
                r["complete_ports"] != r["planned_ports"]):
            raise ValueError("Incomplete or changed island-port inventory")
    link_count = Counter(unordered_pair(r["island_a"], r["island_b"]) for r in links)
    seen_pairs = set()
    analytical = []
    for r in pairs:
        a, b = unordered_pair(r["island_a"], r["island_b"])
        if (a, b) in seen_pairs or a not in eligible_islands or b not in eligible_islands:
            raise ValueError("Duplicate or ineligible island pair")
        seen_pairs.add((a, b))
        shared = numeric_int(r["shared_new_species"], "shared new")
        known = link_count[(a, b)]
        if known > shared:
            raise ValueError("Known paired host links exceed shared-new butterflies")
        analytical.append({
            "island_a": a, "island_b": b,
            "min_of_island_early_max":
                min(coverage[a]["early_max_per_port"], coverage[b]["early_max_per_port"]),
            "shared_new_species": shared,
            "early_matched_butterfly_species": known,
        })
    if len(seen_pairs) != 1275 or sum(row["shared_new_species"] for row in analytical) != 1163:
        raise ValueError("Fixed 1,163 shared-new species-island-pair denominator not reproduced")
    if sum(row["early_matched_butterfly_species"] for row in analytical) != 5:
        raise ValueError("Five independently verified early-host matches not reconciled")

    island_rows = sorted(coverage.values(), key=lambda r: int(r["landmass_id"]))
    support = [staged_counts(analytical, threshold) for threshold in THRESHOLDS]
    if [r["island_pairs"] for r in support] != [1275, 946, 666, 231] or [
        r["shared_new_butterfly_species_island_pairs"] for r in support] != [1163, 800, 586, 189]:
        raise ValueError("Source-locked coverage support disagrees with independent check")

    total_early = sum(r["early_port_record_sum"] for r in island_rows)
    total_all = sum(r["all_port_record_sum"] for r in island_rows)
    most_sampled = max(island_rows, key=lambda r: r["early_port_record_sum"])
    positive_counts = {str(t): sum(r["early_max_per_port"] >= t for r in island_rows)
                       for t in THRESHOLDS if t > 0}
    if (len(island_rows), total_early, total_all, positive_counts) != (
            51, 155002, 638085, {"1": 44, "10": 37, "100": 22}):
        raise ValueError("Early botanical source coverage totals have drifted")
    links_output = []
    for r in links:
        a, b = unordered_pair(r["island_a"], r["island_b"])
        links_output.append({
            "butterfly": r["species"],
            "catalogue_host": r["host_name"],
            "island_a": a, "island_b": b,
            "early_max_per_port_island_a": coverage[a]["early_max_per_port"],
            "early_max_per_port_island_b": coverage[b]["early_max_per_port"],
            "not_historical_proof_of_plant_introduction": True,
        })
    return {
        "schema": "chocho_early_botanical_opportunity_51islands_20261030_v1",
        "role": "data_coverage_audit_not_host_facilitation_or_transport",
        "locked_git_blobs": {k: v[1] for k, v in SOURCES.items()},
        "original_islands": 51,
        "original_ports": 194,
        "original_island_pairs": 1275,
        "original_shared_new_species_pair_denominator": 1163,
        "ports_early_all_aggregate": {
            "early_indexed_port_query_sum_not_unique_records": total_early,
            "all_year_indexed_port_query_sum_not_unique_records": total_all,
            "early_as_fraction_of_all_year_port_query_sums": total_early / total_all,
            "single_largest_early_record_island": most_sampled["landmass_id"],
            "single_largest_early_record_island_port_query_sum": most_sampled["early_port_record_sum"],
            "largest_island_fraction_of_early_query_sum":
                most_sampled["early_port_record_sum"] / total_early,
            "islands_with_no_early_plant_indexed_records": 7,
            "islands_with_at_least_n_early_indexed_records_in_one_port": positive_counts,
            "port_radius_km": 5,
            "overlap_may_duplicate_original_plant_occurrences": True,
        },
        "fixed_thresholds": support,
        "five_early_host_link_coverage": links_output,
        "unobserved_plant_not_proven_absent": True,
        "high_general_plant_count_does_not_prove_complete_host_detection": True,
        "not_a_new_randomized_null_test": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = json.dumps(audit(args.root), indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    if args.output:
        dest = args.output if args.output.is_absolute() else args.root / args.output
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(output, encoding="utf-8")
        print(dest)
    else:
        print(output)


if __name__ == "__main__":
    main()
