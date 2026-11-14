#!/usr/bin/env python3
"""Audit *later GBIF re-indexing*, not butterfly establishment, following 2011-17 first indices.

Uses the exact original 225 x 51 species-island source, pre-2011 5 km port
plant inventories and source-conditioned risk join already preserved by PR #159.

This audit was added AFTER viewing the first-index interaction results. All
intervals and permutation/exact tests are exploratory, not confirmatory. A
plant without an indexed record is UNKNOWN, NOT demonstrably absent.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyze_port_host_ocean_interaction_20261114 import FILES, source_join


def source_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fisher_exact_2sided(a, b, c, d):
    """Two-sided fixed-margin hypergeometric exact p (unclustered diagnostic)."""
    n = a + b + c + d
    draws, wins = a + b, a + c

    def logchoose(x, y):
        if y < 0 or y > x:
            return -math.inf
        return math.lgamma(x + 1) - math.lgamma(y + 1) - math.lgamma(x - y + 1)

    def p_at(x):
        return math.exp(logchoose(wins, x) + logchoose(n - wins, draws - x)
                        - logchoose(n, draws))

    observed = p_at(a)
    lower = max(0, draws - (n - wins))
    upper = min(draws, wins)
    return min(1., sum(p_at(x) for x in range(lower, upper + 1)
                       if p_at(x) <= observed * (1 + 1e-9)))


def cells(rows, predicate):
    subset = [r for r in rows if predicate(r)]
    return {
        "n": len(subset),
        "later_reindexed": sum(r["reindexed"] for r in subset),
        "reindex_rate": (sum(r["reindexed"] for r in subset) / len(subset)
                         if subset else None),
        "islands": len({r["island"] for r in subset}),
        "butterfly_species": len({r["species"] for r in subset}),
        "mean_early_index_records": (
            sum(r["early_records"] for r in subset) / len(subset)
            if subset else None
        ),
    }


def island_bootstrap_rate_difference(rows, n_draws=4999, seed=20261114):
    """Resample complete destination-island clusters; differences in risk only."""
    by_island = defaultdict(list)
    for row in rows:
        by_island[row["island"]].append(row)
    islands = sorted(by_island)
    rng = np.random.default_rng(seed)
    delta = []
    for _ in range(n_draws):
        draw = rng.integers(0, len(islands), len(islands))
        counts = [0, 0]
        successes = [0, 0]
        for ix in draw:
            for row in by_island[islands[ix]]:
                k = row["introduced"]
                counts[k] += 1
                successes[k] += row["reindexed"]
        if counts[0] and counts[1]:
            delta.append(successes[1] / counts[1] -
                         successes[0] / counts[0])
    if not delta:
        return {"error": "No bootstrap draws with both groups"}
    return {
        "n_valid_draws": len(delta),
        "seed": seed,
        "percentile_95pct_interval_pp": [
            float(100 * np.quantile(delta, 0.025)),
            float(100 * np.quantile(delta, 0.975)),
        ],
        "method": "descriptive resample whole islands, not species or records",
    }


def evaluate(rows, *, bootstrap=True):
    intro, unindexed = (cells(rows, lambda r: r["introduced"] == z) for z in (1, 0))
    a, b = intro["later_reindexed"], intro["n"] - intro["later_reindexed"]
    c, d = unindexed["later_reindexed"], unindexed["n"] - unindexed["later_reindexed"]
    z = {
        "n_total": len(rows),
        "introduced_host_indexed": intro,
        "introduced_host_not_indexed_unknown_presence": unindexed,
        "risk_difference_pp": (
            100 * (intro["reindex_rate"] - unindexed["reindex_rate"])
            if intro["n"] and unindexed["n"] else None
        ),
        "unclustered_fisher_exact_p_2sided": (
            fisher_exact_2sided(a, b, c, d)
            if intro["n"] and unindexed["n"] else None
        ),
        "four_historical_record_categories": {
            "introduced_only_indexed": cells(rows,
                lambda r: r["introduced"] == 1 and r["native"] == 0),
            "native_and_introduced_indexed": cells(rows,
                lambda r: r["introduced"] == 1 and r["native"] == 1),
            "native_only_indexed": cells(rows,
                lambda r: r["introduced"] == 0 and r["native"] == 1),
            "no_compatible_host_indexed_NOT_absence": cells(rows,
                lambda r: r["introduced"] == 0 and r["native"] == 0),
        },
        "observed_two_group_table": {
            "introduced_reindexed": a, "introduced_not_reindexed": b,
            "not_introduced_indexed_reindexed": c,
            "not_introduced_indexed_not_reindexed": d,
        },
    }
    if bootstrap and intro["n"] and unindexed["n"]:
        z["island_cluster_bootstrap"] = island_bootstrap_rate_difference(rows)
    return z


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path,
                        default=Path("results/port_host_reindex_20261114"))
    args = parser.parse_args()
    root = args.repo_root
    butterfly_path = root / FILES["butterflies"]
    butterfly = source_csv(butterfly_path)
    hosts = source_csv(root / FILES["hosts"])
    effort_rows = source_csv(root / FILES["plant_effort"])
    effort = defaultdict(int)
    for row in effort_rows:
        effort[row["island"]] += int(row["records"])
    hosts_by_key = {(r["key"], r["island"]): r for r in hosts}
    risk, hashes, coefficients, support = source_join(root)
    earlier_conspecific = {
        (r["species"], r["island"]) for r in risk["2011_2017"]
    }

    candidates = []
    for r in butterfly:
        if int(r["eligible_2011_2017"]) != 1 or int(r["first_indexed_in_2011_2017"]) != 1:
            continue
        sp, isl = r["speciesKey"], r["landmass_id"]
        h = hosts_by_key.get((sp, isl))
        if h is None or effort[isl] <= 0:
            continue
        candidates.append({
            "island": isl,
            "species": sp,
            "name": r["species"],
            "introduced": int(int(h["introduced"]) > 0),
            "native": int(int(h["native"]) > 0),
            "unresolved_host_record": int(int(h["unresolved"]) > 0),
            "early_records": int(r["future2011_2017_indexed_records"]),
            "later_records": int(r["future2018_2022_indexed_records"]),
            "reindexed": int(int(r["future2018_2022_indexed_records"]) > 0),
            "pre2011_port_plant_records": effort[isl],
            "known_earlier_conspecific_on_other_original_island":
                int((sp, isl) in earlier_conspecific),
        })

    out = {
        "audit": "initial_index_to_later_reindex_20261114",
        "observation_windows": {
            "plant_port_records": "2001-2010",
            "first_butterfly_index": "2011-2017",
            "later_reindex": "2018-2022",
        },
        "input_sha256": hashes,
        "source_support": support,
        "audit_subsets": {},
        "not_estimated": [
            "new establishment", "adult population abundance", "observed larval feeding",
            "confirmed oviposition or breeding", "interisland dispersal",
            "actual ship carriage", "directed current or wind flight",
        ],
        "interpretation_limits": [
            "re-indexing by GBIF is not multi-year breeding confirmation",
            "no indexed introduced host does not imply absence of introduced hosts",
            "origin status is modern regional WCVP, not historic island introductions",
            "future records may share sampling effort / recorder and database bias",
            "early butterfly record counts differ between exposure groups",
            "15 historical introduced-host-indexed early-first cases in main 159",
            "multiple conditionally selected subgroups are sensitivity only",
            "Fisher exact p ignores island/species clustering; bootstrap resamples islands only",
            "no species-level historical true occupancy states are defined",
        ],
    }
    min_plants = 10
    main = [r for r in candidates if r["pre2011_port_plant_records"] >= min_plants]
    out["audit_subsets"]["main_pre2011_plant_records_ge10"] = evaluate(main)
    out["audit_subsets"]["earlier_conspecific_island_source_known"] = evaluate(
        [r for r in main
         if r["known_earlier_conspecific_on_other_original_island"] == 1]
    )
    out["audit_subsets"]["early_butterfly_records_ge2"] = evaluate(
        [r for r in main if r["early_records"] >= 2]
    )
    out["audit_subsets"]["historical_native_host_indexed"] = evaluate(
        [r for r in main if r["native"] == 1]
    )
    # Sensitivities below are not independent experiments or multiple confirmations.
    out["audit_subsets"]["all_pre2011_plant_record_positive"] = evaluate(candidates)
    out["audit_subsets"]["pre2011_plant_records_ge30"] = evaluate(
        [r for r in candidates if r["pre2011_port_plant_records"] >= 30]
    )
    out["introduced_case_index"] = [
        {k: r[k] for k in
         ("name", "species", "island", "native", "early_records", "later_records",
          "reindexed", "pre2011_port_plant_records",
          "known_earlier_conspecific_on_other_original_island")}
        for r in main if r["introduced"] == 1
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.out_dir / "audit.json"
    output.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    compact = {
        "output": str(output),
        "main": out["audit_subsets"]["main_pre2011_plant_records_ge10"],
        "source_conditioned": out["audit_subsets"]["earlier_conspecific_island_source_known"],
        "early_records_ge2": out["audit_subsets"]["early_butterfly_records_ge2"],
        "native_indexed_only_comparison": out["audit_subsets"]["historical_native_host_indexed"],
    }
    print(json.dumps(compact, ensure_ascii=False))


if __name__ == "__main__":
    main()
