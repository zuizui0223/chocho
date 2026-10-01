#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
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


def spearman(a, b):
    if len(a) < 3 or len(a) != len(b):
        return None
    x, y = average_ranks(a), average_ranks(b)
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def load_descriptors(path):
    out = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[row["species"]] = {
                "host_family_count": float(row["host_family_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
                "native_resource_units_descriptor": int(row["host_wgsrpd3_unit_count"]),
            }
    return out


def load_pairs(path):
    by_species = defaultdict(dict)
    pools = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id", "family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for row in reader:
            insect = row["insect_species"].strip()
            host = row["accepted_plant_name_id"].strip()
            family = row["family"].strip()
            if insect and host and family:
                by_species[insect][host] = family
                pools[family].add(host)
    return (
        {species: dict(hosts) for species, hosts in by_species.items()},
        {family: tuple(sorted(hosts)) for family, hosts in pools.items()},
    )


def load_units(path):
    out = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("distribution sidecar schema drift")
        for row in reader:
            host = row["accepted_plant_name_id"].strip()
            unit = row["area_code_l3"].strip()
            if host and unit:
                out[host].add(unit)
    return out


def portfolio_metrics(hosts, native_by_host, contemporary_by_host):
    native_union = set()
    contemporary_union = set()
    added_by_host = {}
    for host in hosts:
        native = set(native_by_host.get(host, set()))
        contemporary = set(contemporary_by_host.get(host, set()))
        if not native <= contemporary:
            raise RuntimeError(f"native not subset contemporary: {host}")
        native_union.update(native)
        contemporary_union.update(contemporary)
        added_by_host[host] = contemporary - native

    added_union = contemporary_union - native_union
    result = {
        "native_units": len(native_union),
        "contemporary_units": len(contemporary_union),
        "introduced_added_units": len(added_union),
        "log_expansion": float(
            math.log1p(len(contemporary_union)) - math.log1p(len(native_union))
        ),
    }
    if not added_union:
        result.update(
            maximum_single_host_share=None,
            effective_contributor_number=None,
        )
        return result

    credit = defaultdict(float)
    for unit in added_union:
        contributors = [host for host in hosts if unit in added_by_host[host]]
        share = 1.0 / len(contributors)
        for host in contributors:
            credit[host] += share
    proportions = [value / len(added_union) for value in credit.values()]
    result.update(
        maximum_single_host_share=float(max(proportions)),
        effective_contributor_number=float(
            1.0 / sum(value * value for value in proportions)
        ),
    )
    return result


def seed_int(tag, species, replicate, family):
    raw = f"{tag}|{species}|{replicate}|{family}".encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "big")


def sample_portfolio(composition, pools, tag, species, replicate):
    sampled = []
    for family, count in sorted(composition.items()):
        pool = pools[family]
        sampled.extend(
            random.Random(seed_int(tag, species, replicate, family)).sample(
                list(pool), count
            )
        )
    return tuple(sampled)


def upper_tail(null_values, observed):
    values = [float(x) for x in null_values if x is not None and math.isfinite(x)]
    return (1 + sum(value >= observed for value in values)) / (len(values) + 1)


def lower_tail(null_values, observed):
    values = [float(x) for x in null_values if x is not None and math.isfinite(x)]
    return (1 + sum(value <= observed for value in values)) / (len(values) + 1)


def two_sided_abs(null_values, observed):
    values = [float(x) for x in null_values if x is not None and math.isfinite(x)]
    return (1 + sum(abs(value) >= abs(observed) for value in values)) / (
        len(values) + 1
    )


def summarize_null(values):
    values = np.asarray([x for x in values if x is not None], dtype=float)
    return {
        "q025": float(np.quantile(values, 0.025)),
        "median": float(np.median(values)),
        "q975": float(np.quantile(values, 0.975)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol-json", type=Path, required=True)
    parser.add_argument("--descriptors-csv", type=Path, required=True)
    parser.add_argument("--insect-host-csv", type=Path, required=True)
    parser.add_argument("--native-distribution-csv", type=Path, required=True)
    parser.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "chocho_butterfly_family_composition_matched_null_protocol_v0.2":
        raise RuntimeError("unexpected matched-null protocol")
    null = protocol["null"]
    permutations = int(null["permutations"])
    seed_tag = str(null["seed_tag"])

    descriptors = load_descriptors(args.descriptors_csv)
    pairs, pools = load_pairs(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)
    pools = {
        family: tuple(
            host
            for host in hosts
            if native.get(host) and native[host] <= contemporary.get(host, set())
        )
        for family, hosts in pools.items()
    }

    focal = []
    excluded_insufficient_pool = []
    for species, descriptor in descriptors.items():
        if (
            descriptor["host_family_count"] <= 0
            or descriptor["native_resource_units_descriptor"] <= 0
            or descriptor["resolved_host_species"] < descriptor["host_family_count"]
        ):
            continue
        host_map = pairs.get(species, {})
        if not host_map:
            continue
        composition = Counter(host_map.values())
        if any(len(pools.get(family, ())) < count for family, count in composition.items()):
            excluded_insufficient_pool.append(species)
            continue
        observed = portfolio_metrics(tuple(sorted(host_map)), native, contemporary)
        focal.append(
            {
                "species": species,
                "descriptor": descriptor,
                "composition": composition,
                "observed": observed,
            }
        )

    if not focal:
        raise RuntimeError("no matched-null focal species")

    per_species = {
        row["species"]: {
            "log_expansion": [],
            "introduced_added_units": [],
            "effective_contributor_number": [],
            "maximum_single_host_share": [],
        }
        for row in focal
    }
    global_null = {
        "expanded_species": [],
        "mean_log_expansion": [],
        "median_log_expansion": [],
        "total_introduced_added_units": [],
        "rho_host_family_vs_log_expansion": [],
    }

    host_family_counts = [
        row["descriptor"]["host_family_count"] for row in focal
    ]

    for replicate in range(permutations):
        replicate_metrics = []
        for row in focal:
            species = row["species"]
            sampled = sample_portfolio(
                row["composition"], pools, seed_tag, species, replicate
            )
            metrics = portfolio_metrics(sampled, native, contemporary)
            replicate_metrics.append(metrics)
            per_species[species]["log_expansion"].append(metrics["log_expansion"])
            per_species[species]["introduced_added_units"].append(
                metrics["introduced_added_units"]
            )
            if metrics["introduced_added_units"] > 0:
                per_species[species]["effective_contributor_number"].append(
                    metrics["effective_contributor_number"]
                )
                per_species[species]["maximum_single_host_share"].append(
                    metrics["maximum_single_host_share"]
                )

        log_values = [m["log_expansion"] for m in replicate_metrics]
        added_values = [m["introduced_added_units"] for m in replicate_metrics]
        global_null["expanded_species"].append(
            sum(value > 0 for value in added_values)
        )
        global_null["mean_log_expansion"].append(float(np.mean(log_values)))
        global_null["median_log_expansion"].append(float(np.median(log_values)))
        global_null["total_introduced_added_units"].append(sum(added_values))
        global_null["rho_host_family_vs_log_expansion"].append(
            spearman(host_family_counts, log_values)
        )

    observed_log = [row["observed"]["log_expansion"] for row in focal]
    observed_added = [
        row["observed"]["introduced_added_units"] for row in focal
    ]
    observed_global = {
        "species": len(focal),
        "expanded_species": sum(value > 0 for value in observed_added),
        "mean_log_expansion": float(np.mean(observed_log)),
        "median_log_expansion": float(np.median(observed_log)),
        "total_introduced_added_units": int(sum(observed_added)),
        "rho_host_family_vs_log_expansion": spearman(
            host_family_counts, observed_log
        ),
    }

    species_rows = []
    for row in focal:
        species = row["species"]
        observed = row["observed"]
        null_values = per_species[species]
        effective_null = null_values["effective_contributor_number"]
        max_share_null = null_values["maximum_single_host_share"]
        species_rows.append(
            {
                "species": species,
                "host_family_count": row["descriptor"]["host_family_count"],
                "resolved_host_species": row["descriptor"]["resolved_host_species"],
                "observed_log_expansion": observed["log_expansion"],
                "null_log_expansion_median": float(
                    np.median(null_values["log_expansion"])
                ),
                "p_log_expansion_high": upper_tail(
                    null_values["log_expansion"], observed["log_expansion"]
                ),
                "observed_introduced_added_units": observed[
                    "introduced_added_units"
                ],
                "null_added_units_median": float(
                    np.median(null_values["introduced_added_units"])
                ),
                "p_added_units_high": upper_tail(
                    null_values["introduced_added_units"],
                    observed["introduced_added_units"],
                ),
                "observed_effective_contributor_number": observed[
                    "effective_contributor_number"
                ],
                "conditional_null_effective_median": (
                    None
                    if not effective_null
                    else float(np.median(effective_null))
                ),
                "conditional_p_effective_high": (
                    None
                    if observed["effective_contributor_number"] is None
                    or not effective_null
                    else upper_tail(
                        effective_null,
                        observed["effective_contributor_number"],
                    )
                ),
                "observed_maximum_single_host_share": observed[
                    "maximum_single_host_share"
                ],
                "conditional_null_max_share_median": (
                    None if not max_share_null else float(np.median(max_share_null))
                ),
                "conditional_p_max_share_high": (
                    None
                    if observed["maximum_single_host_share"] is None
                    or not max_share_null
                    else upper_tail(
                        max_share_null,
                        observed["maximum_single_host_share"],
                    )
                ),
                "conditional_p_max_share_low": (
                    None
                    if observed["maximum_single_host_share"] is None
                    or not max_share_null
                    else lower_tail(
                        max_share_null,
                        observed["maximum_single_host_share"],
                    )
                ),
                "family_composition": ";".join(
                    f"{family}:{count}"
                    for family, count in sorted(row["composition"].items())
                ),
            }
        )

    payload = {
        "schema": "chocho_butterfly_family_composition_matched_null_v0.2",
        "status": "POSTHOC_REVIEWER_DEFENSE_SENSITIVITY",
        "protocol": protocol["schema"],
        "permutations": permutations,
        "species": len(focal),
        "excluded_for_insufficient_family_pool": sorted(excluded_insufficient_pool),
        "null_definition": (
            "For each butterfly, preserve exact resolved host-species count and "
            "exact WCVP host-family composition; sample alternative HOSTS-WCVP "
            "host species without replacement within each family and reconstruct "
            "native and contemporary WGSRPD3 resource unions."
        ),
        "observed_global": observed_global,
        "global_null": {
            "expanded_species": {
                **summarize_null(global_null["expanded_species"]),
                "p_high": upper_tail(
                    global_null["expanded_species"],
                    observed_global["expanded_species"],
                ),
            },
            "mean_log_expansion": {
                **summarize_null(global_null["mean_log_expansion"]),
                "p_high": upper_tail(
                    global_null["mean_log_expansion"],
                    observed_global["mean_log_expansion"],
                ),
            },
            "median_log_expansion": {
                **summarize_null(global_null["median_log_expansion"]),
                "p_high": upper_tail(
                    global_null["median_log_expansion"],
                    observed_global["median_log_expansion"],
                ),
            },
            "total_introduced_added_units": {
                **summarize_null(global_null["total_introduced_added_units"]),
                "p_high": upper_tail(
                    global_null["total_introduced_added_units"],
                    observed_global["total_introduced_added_units"],
                ),
            },
            "rho_host_family_vs_log_expansion": {
                **summarize_null(
                    global_null["rho_host_family_vs_log_expansion"]
                ),
                "p_two_sided": two_sided_abs(
                    global_null["rho_host_family_vs_log_expansion"],
                    observed_global["rho_host_family_vs_log_expansion"],
                ),
            },
        },
        "species_with_log_expansion_above_matched_null_p05": sum(
            row["p_log_expansion_high"] <= 0.05 for row in species_rows
        ),
        "architecture_note": (
            "Architecture comparisons are conditional on a randomized portfolio "
            "producing at least one introduced-added unit. They are secondary "
            "because exact host-species-count diagnostics show that the marginal "
            "specialist-generalist architecture gradient is largely structural."
        ),
        "claim_boundary": protocol["claim_boundary"],
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(species_rows[0]))
        writer.writeheader()
        writer.writerows(sorted(species_rows, key=lambda row: row["species"]))
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
