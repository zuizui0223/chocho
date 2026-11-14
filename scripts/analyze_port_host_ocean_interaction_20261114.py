#!/usr/bin/env python3
"""Exploratory time-ordered port-host x marine-access association for chocho.

All inputs are already committed original-source derivatives. No network access
or re-query is performed. First GBIF indexing is NOT establishment. An absent
host record is UNKNOWN, not a verified host absence. This analysis deliberately
does not alter the GEB manuscript or any previously frozen claims.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

FILES = {
    "butterflies": "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv",
    "hosts": "results/route_free_port_hosts/species_island_hosts.csv",
    "plant_effort": "results/route_free_port_hosts/port_coverage.csv",
    "circuits": "results/island51_observed_ocean_circuit_20261106/island51_observed_circuit_pairs.csv",
}
PERIODS = {
    "2011_2017": ("eligible_2011_2017", "first_indexed_in_2011_2017"),
    "2018_2022": ("eligible_2018_2022", "first_indexed_in_2018_2022"),
}


def source_table(root: Path, relative: str):
    path = root / relative
    raw = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return rows, hashlib.sha256(raw).hexdigest()


def source_join(root: Path):
    tables = {}
    hashes = {}
    for key, path in FILES.items():
        tables[key], hashes[key] = source_table(root, path)
    butterflies, hosts = tables["butterflies"], tables["hosts"]
    circuits, ports = tables["circuits"], tables["plant_effort"]

    if len(butterflies) != 225 * 51 or len(circuits) != 51 * 50 // 2:
        raise ValueError("Original 225 x 51 or 1,275 pair source changed")
    if len({(x["speciesKey"], x["landmass_id"]) for x in butterflies}) != len(butterflies):
        raise ValueError("Duplicate butterfly x island keys")
    if len({(x["key"], x["island"]) for x in hosts}) != len(hosts):
        raise ValueError("Duplicate host x island keys")

    host_by_key = {(h["key"], h["island"]): h for h in hosts}
    effort = defaultdict(int)
    for port in ports:
        effort[port["island"]] += int(port["records"])
    by_species = defaultdict(list)
    for b in butterflies:
        by_species[b["speciesKey"]].append(b)

    pair = {}
    q = []
    for z in circuits:
        a, b = z["island_a"], z["island_b"]
        distance = float(z["port_distance_km"])
        resistance = float(z["uniform_distributed_port_resistance"])
        if distance <= 0 or resistance < 0:
            raise ValueError("Invalid source marine distance/resistance")
        pair[a, b] = pair[b, a] = (distance, resistance)
        ld = np.log1p(distance)
        q.append((1., ld, ld * ld, np.log1p(resistance)))
    q = np.asarray(q)
    coefficients = np.linalg.lstsq(q[:, :3], q[:, 3], rcond=None)[0]

    rows = {period: [] for period in PERIODS}
    for period, (eligibility, event) in PERIODS.items():
        later = period == "2018_2022"
        for b in butterflies:
            if int(b[eligibility]) != 1:
                continue
            island, species = b["landmass_id"], b["speciesKey"]
            host = host_by_key.get((species, island))
            if host is None or effort[island] <= 0:
                # NOT "no hosts": no sufficiently indexed plant source.
                continue
            prior_sources = []
            for source in by_species[species]:
                other = source["landmass_id"]
                prior = int(source["pre2001_2010_indexed_records"]) > 0
                if later:
                    prior |= int(source["future2011_2017_indexed_records"]) > 0
                if other != island and prior and (island, other) in pair:
                    distance, resistance = pair[island, other]
                    prior_sources.append((distance, other, resistance))
            if not prior_sources:
                # The geographic access covariate is undefined, not a true zero.
                continue
            distance, earlier_island, resistance = min(prior_sources)
            log_distance = np.log1p(distance)
            expected_log_resistance = coefficients @ np.asarray(
                [1., log_distance, log_distance * log_distance]
            )
            access = -(np.log1p(resistance) - expected_log_resistance)
            rows[period].append({
                "species": species,
                "island": island,
                "earlier_island": earlier_island,
                "event": int(b[event]),
                "introduced": int(int(host["introduced"]) > 0),
                "native": int(int(host["native"]) > 0),
                "effort": effort[island],
                "distance_km": distance,
                "access_residual": float(access),
            })
    return rows, hashes, coefficients.tolist(), {
        "butterfly_cells": len(butterflies),
        "marine_pairs": len(circuits),
        "historical_port_host_cells": len(hosts),
        "host_crosswalk_islands": len({i for _, i in host_by_key if i in {
            b["landmass_id"] for b in butterflies
        }}),
    }


def tw_demean(values, islands, species):
    x = np.asarray(values, dtype=float).copy()
    groups = []
    for labels in (islands, species):
        _, group = np.unique(labels, return_inverse=True)
        groups.append((group, np.bincount(group).astype(float)))
    for _ in range(200):
        worst = 0.
        for group, count in groups:
            means = np.bincount(group, weights=x, minlength=len(count)) / count
            x -= means[group]
            worst = max(worst, float(np.max(np.abs(means))))
        if worst < 1e-11:
            break
    return x


def estimate(rows):
    if len(rows) < 50:
        return {"not_identifiable": "Too few joined records", "n": len(rows)}
    n = len(rows)
    island = np.asarray([r["island"] for r in rows])
    species = np.asarray([r["species"] for r in rows])
    host = np.asarray([r["introduced"] for r in rows], dtype=float)
    native = np.asarray([r["native"] for r in rows], dtype=float)
    access_raw = np.asarray([r["access_residual"] for r in rows], dtype=float)
    distance_raw = np.log1p([r["distance_km"] for r in rows])
    std_a, std_d = np.std(access_raw), np.std(distance_raw)
    if std_a <= 0 or std_d <= 0:
        return {"not_identifiable": "No sea/distance variation", "n": n}
    access = (access_raw - np.mean(access_raw)) / std_a
    distance = (distance_raw - np.mean(distance_raw)) / std_d
    x = np.column_stack((host, native, access, distance, host * access))
    y = np.asarray([r["event"] for r in rows], dtype=float)
    x_dm = np.column_stack([
        tw_demean(x[:, k], island, species) for k in range(x.shape[1])
    ])
    y_dm = tw_demean(y, island, species)
    rank = int(np.linalg.matrix_rank(x_dm))
    if rank < x.shape[1]:
        return {"not_identifiable": "Insufficient fixed-effect support",
                "rank": rank, "n": n, "events": int(sum(y))}
    beta = np.linalg.lstsq(x_dm, y_dm, rcond=None)[0]
    residual = y_dm - x_dm @ beta
    bread = np.linalg.pinv(x_dm.T @ x_dm)
    target = 4  # introduced x standardized distance-residual sea access
    influence = residual * (x_dm @ bread[target, :])

    def clustered_sum_of_squares(labels):
        _, group = np.unique(labels, return_inverse=True)
        score = np.bincount(group, weights=influence)
        return float(np.sum(score * score))

    two_way_var = (clustered_sum_of_squares(island) +
                   clustered_sum_of_squares(species) -
                   clustered_sum_of_squares(np.asarray(
                       [a + "|" + b for a, b in zip(island, species)])))
    g = len(np.unique(island))
    se = np.sqrt(max(0., two_way_var) * g / (g - 1))
    coefficient_pp = float(100 * beta[target])
    standard_error_pp = float(100 * se)
    return {
        "n": n, "events": int(y.sum()), "islands": g,
        "species": len(np.unique(species)),
        "introduced_record_cells": int(host.sum()),
        "introduced_event_cells": int(np.dot(host, y)),
        "access_residual_SD": float(std_a),
        "introduce_x_sea_access_effect_pp_per_SD": coefficient_pp,
        "two_way_cluster_se_pp": standard_error_pp,
        "approximate_95pct_interval_pp": [
            coefficient_pp - 2.025 * standard_error_pp,
            coefficient_pp + 2.025 * standard_error_pp,
        ],
        "native_main_effect_pp": float(100 * beta[1]),
        "introduced_main_effect_pp": float(100 * beta[0]),
        "sea_access_main_effect_pp": float(100 * beta[2]),
        "distance_main_effect_pp": float(100 * beta[3]),
        "method": "LPM with island and butterfly FE, native introduced access "
                  "and log-distance terms; approximate two-way clustered SE",
    }


def model_sets(all_rows):
    selected = {}
    for period, rows in all_rows.items():
        enough = [r for r in rows if r["effort"] >= 10]
        island_support = defaultdict(set)
        species_support = defaultdict(set)
        for r in enough:
            island_support[r["island"]].add(r["introduced"])
            species_support[r["species"]].add(r["introduced"])
        overlap = [r for r in enough
                   if len(island_support[r["island"]]) == 2
                   and len(species_support[r["species"]]) == 2]
        selected[period] = {
            "all_indexed_ports": estimate(rows),
            "plant_records_ge_10": estimate(enough),
            "two_sided_island_species_overlap": estimate(overlap),
        }
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path,
                        default=Path("results/port_host_ocean_interaction_20261114"))
    args = parser.parse_args()
    rows, hashes, coefficients, support = source_join(args.repo_root)
    report = {
        "analysis": "port_host_ocean_interaction_20261114",
        "source_sha256": hashes,
        "source_support": support,
        "marine_log_resistance_vs_log_distance_coefficients": coefficients,
        "periods": model_sets(rows),
        "interpretation": {
            "hypothesis": "early known introduced host near port interacts with "
                          "distance-residual uniform-sea access to predict later "
                          "species-island first GBIF indexing",
            "not_claimed": ["true colonization or breeding", "true host absence",
                            "historically dated plant introduction", "larval use",
                            "wind-driven travel", "ship transport", "causal mediation"],
            "unresolved": [
                "plant origin is regional contemporary WCVP status",
                "unindexed historical host is UNKNOWN, not absence",
                "source islands restricted to original 51 selected landmasses",
                "port sampling effort and butterfly GBIF effort remain incomplete",
                "species/island cluster interval is approximate, not confirmatory",
                "source plant status fixed in 2001-2010, but butterflies "
                "first indexed within broad time windows",
            ],
        },
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    path = args.out_dir / "source_join_pilot.json"
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    print(json.dumps({
        "output": str(path),
        "source_support": support,
        "periods": report["periods"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
