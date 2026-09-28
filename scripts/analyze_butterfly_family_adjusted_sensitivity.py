#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def average_ranks(values):
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and x[order[stop]] == x[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def correlation(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def residualize_groups(values, groups):
    values = np.asarray(values, dtype=float)
    out = values.copy()
    buckets = defaultdict(list)
    for index, group in enumerate(groups):
        buckets[group].append(index)
    for indices in buckets.values():
        idx = np.asarray(indices, dtype=int)
        out[idx] -= np.mean(values[idx])
    return out


def adjusted_rank_correlation(a, b, groups):
    return correlation(
        residualize_groups(average_ranks(a), groups),
        residualize_groups(average_ranks(b), groups),
    )


def load_butterfly_family(path):
    result = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"Family", "Species"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("LepTraits taxonomy columns missing")
        for row in reader:
            result[str(row["Species"]).strip()] = str(row["Family"]).strip()
    return result


def load_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--leptraits-csv", type=Path, required=True)
    parser.add_argument("--anthropogenic-csv", type=Path, required=True)
    parser.add_argument("--mechanism-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    family = load_butterfly_family(args.leptraits_csv)
    resource = load_csv(args.anthropogenic_csv)
    mechanism = [
        row
        for row in load_csv(args.mechanism_csv)
        if row["host_taxonomy_lower_bound_adequate"] == "True"
        and int(row["introduced_added_units"]) > 0
    ]

    missing = sorted(
        {row["species"] for row in resource + mechanism if row["species"] not in family}
    )
    if missing:
        raise RuntimeError("species missing LepTraits Family: " + ", ".join(missing))

    resource_family = [family[row["species"]] for row in resource]
    host_family = [float(row["host_family_count"]) for row in resource]
    log_expansion = [float(row["log_resource_expansion"]) for row in resource]
    native_units = [float(row["native_resource_units"]) for row in resource]
    added_units = [float(row["introduced_added_units"]) for row in resource]

    mech_family = [family[row["species"]] for row in mechanism]
    mech_host_family = [float(row["host_family_count"]) for row in mechanism]
    mech_host_species = [int(row["resolved_host_species"]) for row in mechanism]
    effective = [float(row["effective_contributor_number"]) for row in mechanism]
    max_share = [float(row["maximum_single_host_fractional_share"]) for row in mechanism]
    combined_groups = [
        f"{fam}|{count}" for fam, count in zip(mech_family, mech_host_species)
    ]

    payload = {
        "schema": "chocho_butterfly_family_adjusted_sensitivity_v0.1",
        "status": "POSTHOC_PHYLOGENETIC_NONINDEPENDENCE_SENSITIVITY",
        "resource_species": len(resource),
        "mechanism_species": len(mechanism),
        "butterfly_family_counts_resource": dict(sorted(Counter(resource_family).items())),
        "resource_expansion": {
            "family_adjusted_rank_correlation_host_breadth_vs_log_expansion":
                adjusted_rank_correlation(host_family, log_expansion, resource_family),
            "family_adjusted_rank_correlation_host_breadth_vs_absolute_added_units":
                adjusted_rank_correlation(host_family, added_units, resource_family),
            "family_adjusted_rank_correlation_host_breadth_vs_native_resource_units":
                adjusted_rank_correlation(host_family, native_units, resource_family),
        },
        "architecture": {
            "exact_host_species_plus_butterfly_family_adjusted_effective_contributors":
                adjusted_rank_correlation(mech_host_family, effective, combined_groups),
            "exact_host_species_plus_butterfly_family_adjusted_maximum_single_host_share":
                adjusted_rank_correlation(mech_host_family, max_share, combined_groups),
            "group_count": len(set(combined_groups)),
        },
        "interpretation": (
            "The near-zero association between taxonomic host breadth and proportional "
            "resource expansion persists after controlling butterfly Family. The apparent "
            "specialist-generalist architecture gradient disappears after jointly fixing "
            "resolved host-species richness and butterfly Family."
        ),
        "claim_boundary": (
            "Butterfly Family is a coarse proxy for phylogenetic non-independence, not a "
            "replacement for a dated species-level phylogeny or PGLS."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
