#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
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
    if len(a) < 3:
        return None
    x, y = average_ranks(a), average_ranks(b)
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def load_descriptors(path):
    rows = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows[row["species"]] = {
                "host_family_count": float(row["host_family_count"]),
                "resolved_host_species": int(row["resolved_host_species"]),
                "native_resource_units": int(row["host_wgsrpd3_unit_count"]),
            }
    return rows


def load_pairs(path):
    by_species = defaultdict(dict)
    consumers = defaultdict(set)
    family_pool = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id", "family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("interaction sidecar schema drift")
        for row in reader:
            insect = str(row["insect_species"]).strip()
            host = str(row["accepted_plant_name_id"]).strip()
            family = str(row["family"]).strip()
            if not insect or not host or not family:
                continue
            by_species[insect][host] = family
            consumers[host].add(insect)
            family_pool[family].add(host)
    return (
        {species: dict(hosts) for species, hosts in by_species.items()},
        {host: frozenset(vals) for host, vals in consumers.items()},
        {family: tuple(sorted(hosts)) for family, hosts in family_pool.items()},
    )


def load_units(path):
    out = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            host = str(row["accepted_plant_name_id"]).strip()
            unit = str(row["area_code_l3"]).strip()
            if host and unit:
                out[host].add(unit)
    return {host: frozenset(units) for host, units in out.items()}


def portfolio_metrics(hosts, native, contemporary):
    native_union = set()
    contemporary_union = set()
    for host in hosts:
        native_union.update(native.get(host, ()))
        contemporary_union.update(contemporary.get(host, ()))
    if not native_union:
        return None
    added = contemporary_union - native_union
    return {
        "native_units": len(native_union),
        "contemporary_units": len(contemporary_union),
        "added_units": len(added),
        "log_expansion": float(
            math.log1p(len(contemporary_union)) - math.log1p(len(native_union))
        ),
    }


def stable_seed(tag, model, replicate, species):
    raw = f"{tag}|{model}|{replicate}|{species}".encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "big")


def weighted_choice(rng, candidates, weights):
    w = np.asarray(weights, dtype=float)
    if len(candidates) == 0:
        raise RuntimeError("empty candidate set")
    if not np.all(np.isfinite(w)) or np.sum(w) <= 0:
        raise RuntimeError("invalid candidate weights")
    w /= np.sum(w)
    return candidates[int(rng.choice(len(candidates), p=w))]


def sample_portfolio(
    species,
    observed_host_map,
    family_pool,
    native,
    consumers,
    rng,
    model,
    bandwidth,
    usage_exponent,
):
    observed_hosts = set(observed_host_map)
    selected = set()
    sampled = []
    native_differences = []
    sampled_usage = []

    by_family = defaultdict(list)
    for host, family in observed_host_map.items():
        by_family[family].append(host)

    for family, observed in sorted(by_family.items()):
        alternatives = [
            host
            for host in family_pool.get(family, ())
            if host not in observed_hosts
            and host in native
            and len(native[host]) > 0
        ]
        if len(alternatives) < len(observed):
            raise ValueError(f"insufficient alternatives in {family}")

        # Match the most range-extreme hosts first so they cannot be stranded
        # after more central hosts consume nearby alternatives.
        family_native = np.asarray(
            [math.log1p(len(native[h])) for h in alternatives], dtype=float
        )
        center = float(np.median(family_native))
        observed_order = sorted(
            observed,
            key=lambda h: abs(math.log1p(len(native[h])) - center),
            reverse=True,
        )

        for observed_host in observed_order:
            candidates = [h for h in alternatives if h not in selected]
            target = math.log1p(len(native[observed_host]))
            weights = []
            for candidate in candidates:
                native_distance = abs(math.log1p(len(native[candidate])) - target)
                weight = 1.0
                if model in {"native_range", "native_range_usage"}:
                    weight *= math.exp(-native_distance / bandwidth)
                if model in {"usage", "native_range_usage"}:
                    other_users = len(consumers.get(candidate, frozenset()) - {species})
                    weight *= (1.0 + other_users) ** usage_exponent
                weights.append(weight)

            chosen = weighted_choice(rng, candidates, weights)
            selected.add(chosen)
            sampled.append(chosen)
            native_differences.append(
                abs(math.log1p(len(native[chosen])) - target)
            )
            sampled_usage.append(
                len(consumers.get(chosen, frozenset()) - {species})
            )

    return tuple(sampled), native_differences, sampled_usage


def null_summary(values, observed):
    x = np.asarray(values, dtype=float)
    return {
        "q025": float(np.quantile(x, 0.025)),
        "median": float(np.median(x)),
        "q975": float(np.quantile(x, 0.975)),
        "p_high": float((1 + np.sum(x >= observed)) / (len(x) + 1)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--native-distribution-csv", type=Path, required=True)
    ap.add_argument("--contemporary-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol.get("schema") != "chocho_butterfly_hostbias_null_protocol_v0.1":
        raise RuntimeError("unexpected host-bias null protocol")
    permutations = int(protocol["permutations"])
    bandwidth = float(protocol["native_log_bandwidth"])
    usage_exponent = float(protocol["usage_exponent"])
    seed_tag = str(protocol["seed_tag"])

    descriptors = load_descriptors(args.descriptors_csv)
    pairs, consumers, family_pool = load_pairs(args.insect_host_csv)
    native = load_units(args.native_distribution_csv)
    contemporary = load_units(args.contemporary_distribution_csv)

    # All candidates must have internally consistent native/contemporary footprints.
    valid = {
        host
        for host, units in native.items()
        if units and units <= contemporary.get(host, frozenset())
    }
    family_pool = {
        family: tuple(host for host in hosts if host in valid)
        for family, hosts in family_pool.items()
    }

    focal = []
    excluded = {}
    for species, descriptor in descriptors.items():
        if (
            descriptor["host_family_count"] <= 0
            or descriptor["native_resource_units"] <= 0
            or descriptor["resolved_host_species"] < descriptor["host_family_count"]
        ):
            continue
        observed_host_map = pairs.get(species, {})
        if not observed_host_map:
            excluded[species] = "no_resolved_hosts"
            continue
        composition = Counter(observed_host_map.values())
        bad = [
            family
            for family, count in composition.items()
            if len(
                [
                    h for h in family_pool.get(family, ())
                    if h not in observed_host_map
                ]
            ) < count
        ]
        if bad:
            excluded[species] = "insufficient_alternative_pool:" + ",".join(sorted(bad))
            continue
        observed = portfolio_metrics(tuple(observed_host_map), native, contemporary)
        if observed is None:
            excluded[species] = "no_native_union"
            continue
        focal.append(
            {
                "species": species,
                "descriptor": descriptor,
                "host_map": observed_host_map,
                "observed": observed,
            }
        )

    models = ("native_range", "usage", "native_range_usage")
    observed_total = sum(row["observed"]["added_units"] for row in focal)
    observed_mean = float(np.mean([row["observed"]["log_expansion"] for row in focal]))
    observed_median = float(np.median([row["observed"]["log_expansion"] for row in focal]))
    observed_rho = spearman(
        [row["descriptor"]["host_family_count"] for row in focal],
        [row["observed"]["log_expansion"] for row in focal],
    )
    observed_usage = [
        len(consumers.get(host, frozenset()) - {row["species"]})
        for row in focal
        for host in row["host_map"]
    ]

    result_models = {}
    for model in models:
        total_added = []
        mean_log = []
        median_log = []
        rho = []
        median_native_match = []
        mean_sampled_usage = []

        for replicate in range(permutations):
            rep_metrics = []
            rep_native_diff = []
            rep_usage = []
            for row in focal:
                rng = np.random.default_rng(
                    stable_seed(seed_tag, model, replicate, row["species"])
                )
                sampled, native_diff, sampled_usage = sample_portfolio(
                    row["species"],
                    row["host_map"],
                    family_pool,
                    native,
                    consumers,
                    rng,
                    model,
                    bandwidth,
                    usage_exponent,
                )
                metrics = portfolio_metrics(sampled, native, contemporary)
                rep_metrics.append(metrics)
                rep_native_diff.extend(native_diff)
                rep_usage.extend(sampled_usage)

            logs = [m["log_expansion"] for m in rep_metrics]
            total_added.append(sum(m["added_units"] for m in rep_metrics))
            mean_log.append(float(np.mean(logs)))
            median_log.append(float(np.median(logs)))
            rho.append(
                spearman(
                    [row["descriptor"]["host_family_count"] for row in focal],
                    logs,
                )
            )
            median_native_match.append(float(np.median(rep_native_diff)))
            mean_sampled_usage.append(float(np.mean(rep_usage)))

        result_models[model] = {
            "total_introduced_added_units": null_summary(total_added, observed_total),
            "mean_log_expansion": null_summary(mean_log, observed_mean),
            "median_log_expansion": null_summary(median_log, observed_median),
            "rho_host_family_vs_log_expansion": {
                "q025": float(np.quantile(rho, 0.025)),
                "median": float(np.median(rho)),
                "q975": float(np.quantile(rho, 0.975)),
                "p_two_sided": float(
                    (1 + np.sum(np.abs(rho) >= abs(observed_rho)))
                    / (len(rho) + 1)
                ),
            },
            "matching_diagnostics": {
                "median_absolute_log1p_native_breadth_difference_across_replicates": float(
                    np.median(median_native_match)
                ),
                "q975_median_absolute_log1p_native_breadth_difference": float(
                    np.quantile(median_native_match, 0.975)
                ),
                "mean_other_lepidoptera_users_per_sampled_host_median": float(
                    np.median(mean_sampled_usage)
                ),
            },
        }

    payload = {
        "schema": "chocho_butterfly_hostbias_null_result_v0.1",
        "status": "POSTHOC_RECORDING_AND_HOST_COMMONNESS_SENSITIVITY",
        "protocol": protocol["schema"],
        "species": len(focal),
        "excluded_species": excluded,
        "observed": {
            "total_introduced_added_units": observed_total,
            "mean_log_expansion": observed_mean,
            "median_log_expansion": observed_median,
            "rho_host_family_vs_log_expansion": observed_rho,
            "mean_other_lepidoptera_users_per_observed_host": float(np.mean(observed_usage)),
            "median_other_lepidoptera_users_per_observed_host": float(np.median(observed_usage)),
        },
        "models": result_models,
        "interpretation_rule": (
            "A host-identity excess is robust to plant commonness/recording bias only if "
            "observed expansion remains above null expectation when candidate hosts are "
            "matched on native range breadth and/or weighted by use by other Lepidoptera."
        ),
        "claim_boundary": protocol["claim_boundary"],
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
