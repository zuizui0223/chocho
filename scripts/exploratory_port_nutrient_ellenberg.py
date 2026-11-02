#!/usr/bin/env python3
"""Exploratory, exact-name-only Ellenberg nutrient trait pilot for chocho port hosts.

This is a *descriptive* test among globally port-indexed candidate non-native
larval hosts; it does not estimate the effects of sewage, nitrogen deposition,
or ports relative to non-port coasts. Do not reuse it as confirmatory evidence.
External Ellenberg data are not vendored. See adjacent protocol README.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import statistics
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_HOSTS = REPO / (
    "provenance/reviewer_defenses/results/"
    "port_alien_broad_inventory_20261030/port_alien_host_broad_species.csv"
)


def parse_value(value: str, *, allow_zero: bool = False) -> int | None:
    s = value.strip()
    # Examples in the historic German Ellenberg table: 5, (5), 5~ and x.
    m = re.fullmatch(r"\(?([0-9])\)?~?", s)
    if not m:
        return None
    v = int(m.group(1))
    return v if (allow_zero or v >= 1) else None


def load_data(host_csv: Path, ellenberg_tsv: Path) -> tuple[list[dict], dict]:
    with ellenberg_tsv.open(encoding="utf-8-sig", newline="") as f:
        traits = list(csv.DictReader(f, delimiter="\t"))
    table = {}
    for row in traits:
        key = (row.get("Name") or "").strip().casefold()
        if key:
            table[key] = {
                "N": parse_value(row.get("N") or ""),
                "L": parse_value(row.get("L") or ""),
                "F": parse_value(row.get("F") or ""),
                "S": parse_value(row.get("S") or "", allow_zero=True),
            }
    with host_csv.open(encoding="utf-8-sig", newline="") as f:
        hosts = list(csv.DictReader(f))
    rows = []
    for row in hosts:
        name = row["host_name"].strip()
        trait = table.get(name.casefold())
        if trait is None or trait["N"] is None:
            continue
        rows.append(
            {"name": name, "genus": name.split(" ")[0], "ports": int(row["ports"]), **trait}
        )
    if len({x["name"] for x in rows}) != len(rows):
        raise ValueError("Duplicate exact-name host matches; must resolve before analysis.")
    return rows, {
        "n_port_candidate_host_species": len(hosts),
        "n_german_ellenberg_n_valued_names": sum(
            v["N"] is not None for v in table.values()
        ),
        "n_exact_species_matches_with_N": len(rows),
        "match_fraction": len(rows) / len(hosts) if hosts else None,
    }


def ranks(values: list[int | float]) -> list[float]:
    pairs = sorted(enumerate(values), key=lambda t: t[1])
    out = [math.nan] * len(values)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and pairs[j][1] == pairs[i][1]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            out[pairs[k][0]] = avg_rank
        i = j
    return out


def spearman(x: list[int | float], y: list[int | float]) -> float:
    if len(x) != len(y) or len(x) < 3:
        return math.nan
    a, b = ranks(x), ranks(y)
    ma, mb = statistics.mean(a), statistics.mean(b)
    ab = sum((z - ma) * (w - mb) for z, w in zip(a, b))
    aa = sum((z - ma) ** 2 for z in a)
    bb = sum((w - mb) ** 2 for w in b)
    return ab / math.sqrt(aa * bb) if aa > 0 and bb > 0 else math.nan


def corr(rows: list[dict], x: str = "N") -> float:
    usable = [r for r in rows if r.get(x) is not None]
    return spearman([r[x] for r in usable], [r["ports"] for r in usable])


def bootstrap_genus(rows: list[dict], *, seed: int, n_draws: int) -> list[float]:
    groups = defaultdict(list)
    for r in rows:
        groups[r["genus"]].append(r)
    genera = sorted(groups)
    rng = random.Random(seed)
    draws = []
    for _ in range(n_draws):
        re_sample = [
            record for _ in genera for record in groups[rng.choice(genera)]
        ]
        value = corr(re_sample)
        if math.isfinite(value):
            draws.append(value)
    return sorted(draws)


def analyze(rows: list[dict], meta: dict, *, n_draws: int, seed: int) -> dict:
    if len(rows) < 3:
        raise ValueError("Fewer than three matched species; no correlation is estimable.")
    upper = [r["ports"] for r in rows if r["N"] >= 7]
    lower = [r["ports"] for r in rows if r["N"] <= 3]
    genera = set(r["genus"] for r in rows)
    loo = [corr([r for r in rows if r["genus"] != g]) for g in genera]
    boot = bootstrap_genus(rows, seed=seed, n_draws=n_draws)
    q = lambda p: boot[max(0, min(len(boot) - 1, round((len(boot) - 1) * p)))]
    sizes = sorted(r["ports"] for r in rows)
    cutoff = sizes[int((len(sizes) - 1) * 0.9)]
    traits = {}
    for k in ("L", "F", "S"):
        a = [r for r in rows if r[k] is not None]
        traits[k] = {
            "n": len(a),
            "spearman_trait_vs_ports": corr(a, k),
            "same_subset_N_vs_ports": corr(a),
        }
    return {
        "status": "DESCRIPTIVE_TRAIT_ASSOCIATION_NOT_SEWAGE_EFFECT",
        "source": {
            "hosts": "zuizui0223/chocho:port_alien_host_broad_species.csv",
            "traits": "bathyscapher/EVA:Ellenberg_VascularPlants.csv",
            "trait_scope": "historic Germany-based indicator values; no global validity",
        },
        "counts": {**meta, "n_matched_genera": len(genera)},
        "spearman_n_vs_number_of_indexed_ports": corr(rows),
        "group_descriptives": {
            "N_1_to_3": {"species": len(lower), "median_ports": statistics.median(lower) if lower else None},
            "N_7_to_9": {"species": len(upper), "median_ports": statistics.median(upper) if upper else None},
        },
        "leave_one_genus_out_rho_range": [min(loo), max(loo)],
        "genus_resampled_descriptive_interval_95": [q(0.025), q(0.975)],
        "bootstrap_draws": n_draws,
        "bootstrap_seed": seed,
        "sensitivity_excluding_top_10_percent_ports": {
            "cutoff_inclusive": cutoff,
            "n": sum(r["ports"] <= cutoff for r in rows),
            "rho": corr([r for r in rows if r["ports"] <= cutoff]),
        },
        "other_ellenberg_traits": traits,
        "boundaries": [
            "No port-versus-nonport botanical comparison.",
            "Port index count mixes geographic prevalence and GBIF observation effort.",
            "Only exact binomials with Germany-derived N index are evaluated.",
            "N is a plant site-fertility affinity indicator, not direct tissue nitrogen.",
            "No measured wastewater nitrogen loading or site-specific soil nitrogen.",
            "Genus bootstrap measures subset resampling variability, not causal uncertainty.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ports", type=Path, default=DEFAULT_HOSTS)
    parser.add_argument("--ellenberg", type=Path, required=True,
                        help="Local external Ellenberg_VascularPlants.csv (tab-separated)")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=72819)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    rows, meta = load_data(args.ports, args.ellenberg)
    result = analyze(rows, meta, n_draws=args.bootstrap, seed=args.seed)
    out = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(out, encoding="utf-8")
    else:
        print(out, end="")


if __name__ == "__main__":
    main()
