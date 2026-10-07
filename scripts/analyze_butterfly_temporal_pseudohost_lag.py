#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


def quantile(values, p):
    xs = sorted(values)
    if not xs:
        return None
    pos = p * (len(xs) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(xs[lo])
    frac = pos - lo
    return float(xs[lo] * (1 - frac) + xs[hi] * frac)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queried-csv", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=99999)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--minimum-dated-pseudos-per-host", type=int, default=3)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-cell-csv", type=Path, required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.queried_csv.open(newline="", encoding="utf-8")))
    if not rows:
        raise RuntimeError("empty queried pseudo-host table")

    by_group = defaultdict(list)
    for row in rows:
        by_group[row["match_group"]].append(row)

    groups = {}
    for key, vals in by_group.items():
        actual = [r for r in vals if r["assignment_type"] == "actual"]
        pseudo = [r for r in vals if r["assignment_type"] == "pseudo"]
        if len(actual) != 1:
            raise RuntimeError(f"expected one actual row in {key}; got {len(actual)}")
        a = actual[0]
        ah = str(a["local_host_first_record_year"]).strip()
        pseudo_years = [
            int(float(r["local_host_first_record_year"]))
            for r in pseudo
            if str(r["local_host_first_record_year"]).strip()
        ]
        groups[key] = {
            "species": a["species"],
            "region": a["wgsrpd3_code"],
            "actual_host_name": a["candidate_name"],
            "actual_host_year": None if not ah else int(float(ah)),
            "butterfly_year": None if not str(a["butterfly_first_record_year_full"]).strip() else int(float(a["butterfly_first_record_year_full"])),
            "pseudo_years": pseudo_years,
            "pseudo_total": len(pseudo),
        }

    by_cell = defaultdict(list)
    for key, g in groups.items():
        by_cell[(g["species"], g["region"])].append((key, g))

    cell_rows = []
    eligible_cells = []
    for (species, region), gvals in sorted(by_cell.items()):
        butterfly_years = {g["butterfly_year"] for _, g in gvals if g["butterfly_year"] is not None}
        butterfly_year = None if not butterfly_years else min(butterfly_years)
        actual_years = [g["actual_host_year"] for _, g in gvals if g["actual_host_year"] is not None]
        full_group_support = all(
            g["actual_host_year"] is not None
            and len(g["pseudo_years"]) >= args.minimum_dated_pseudos_per_host
            and g["butterfly_year"] is not None
            for _, g in gvals
        )
        actual_earliest = None if not actual_years else min(actual_years)
        actual_lag = None if actual_earliest is None or butterfly_year is None else butterfly_year - actual_earliest
        row = {
            "species": species,
            "wgsrpd3_code": region,
            "actual_host_groups": len(gvals),
            "groups_with_minimum_pseudo_support": sum(
                len(g["pseudo_years"]) >= args.minimum_dated_pseudos_per_host for _, g in gvals
            ),
            "butterfly_first_record_year_full": butterfly_year,
            "actual_earliest_host_record_year": actual_earliest,
            "actual_record_lag_years": actual_lag,
            "eligible_for_identity_null": int(full_group_support),
        }
        cell_rows.append(row)
        if full_group_support:
            eligible_cells.append(((species, region), gvals, butterfly_year, actual_lag))

    if not eligible_cells:
        raise RuntimeError("no cells have complete actual and pseudo timing support")

    observed_lags = [int(x[3]) for x in eligible_cells]
    observed_median = float(statistics.median(observed_lags))
    observed_positive = sum(x > 0 for x in observed_lags)
    observed_nonnegative = sum(x >= 0 for x in observed_lags)

    rng = random.Random(args.seed)
    null_medians = []
    null_positive = []
    null_nonnegative = []
    cell_null_lags = defaultdict(list)

    for _ in range(args.permutations):
        lags = []
        for (species, region), gvals, butterfly_year, _actual_lag in eligible_cells:
            chosen = [rng.choice(g["pseudo_years"]) for _, g in gvals]
            pseudo_earliest = min(chosen)
            lag = int(butterfly_year) - int(pseudo_earliest)
            lags.append(lag)
            cell_null_lags[(species, region)].append(lag)
        null_medians.append(float(statistics.median(lags)))
        null_positive.append(sum(x > 0 for x in lags))
        null_nonnegative.append(sum(x >= 0 for x in lags))

    p_median = (1 + sum(x >= observed_median for x in null_medians)) / (args.permutations + 1)
    p_positive = (1 + sum(x >= observed_positive for x in null_positive)) / (args.permutations + 1)
    p_nonnegative = (1 + sum(x >= observed_nonnegative for x in null_nonnegative)) / (args.permutations + 1)

    cell_lookup = {(r["species"], r["wgsrpd3_code"]): r for r in cell_rows}
    for key, vals in cell_null_lags.items():
        r = cell_lookup[key]
        obs = int(r["actual_record_lag_years"])
        r["pseudo_lag_median"] = float(statistics.median(vals))
        r["pseudo_lag_q025"] = quantile(vals, 0.025)
        r["pseudo_lag_q975"] = quantile(vals, 0.975)
        r["actual_minus_pseudo_median_lag"] = obs - r["pseudo_lag_median"]
        r["cell_p_pseudo_lag_ge_actual"] = (1 + sum(v >= obs for v in vals)) / (len(vals) + 1)

    def summarize_subset(subset):
        if not subset:
            return {"cells": 0}
        keys = [(x[0][0], x[0][1]) for x in subset]
        obs_lags = [int(x[3]) for x in subset]
        obs_median = float(statistics.median(obs_lags))
        obs_positive = sum(x > 0 for x in obs_lags)
        obs_nonnegative = sum(x >= 0 for x in obs_lags)
        sim_medians, sim_positive, sim_nonnegative = [], [], []
        for b in range(args.permutations):
            vals = [cell_null_lags[key][b] for key in keys]
            sim_medians.append(float(statistics.median(vals)))
            sim_positive.append(sum(x > 0 for x in vals))
            sim_nonnegative.append(sum(x >= 0 for x in vals))
        return {
            "cells": len(subset),
            "species_regions": [f"{x[0][0]}|{x[0][1]}" for x in subset],
            "actual_median_record_lag_years": obs_median,
            "actual_positive_host_first_cells": obs_positive,
            "actual_nonnegative_host_first_cells": obs_nonnegative,
            "pseudo_median_lag_median": float(statistics.median(sim_medians)),
            "pseudo_median_lag_ci95": [quantile(sim_medians, 0.025), quantile(sim_medians, 0.975)],
            "pseudo_positive_cells_median": float(statistics.median(sim_positive)),
            "pseudo_positive_cells_ci95": [quantile(sim_positive, 0.025), quantile(sim_positive, 0.975)],
            "p_pseudo_median_lag_ge_actual": (1 + sum(x >= obs_median for x in sim_medians)) / (args.permutations + 1),
            "p_pseudo_positive_cells_ge_actual": (1 + sum(x >= obs_positive for x in sim_positive)) / (args.permutations + 1),
            "p_pseudo_nonnegative_cells_ge_actual": (1 + sum(x >= obs_nonnegative for x in sim_nonnegative)) / (args.permutations + 1),
        }

    post2017_cells = [x for x in eligible_cells if int(x[2]) >= 2018]
    post2017_summary = summarize_subset(post2017_cells)

    args.output_cell_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_cell_csv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(cell_rows[0].keys())
        extras = [
            "pseudo_lag_median","pseudo_lag_q025","pseudo_lag_q975",
            "actual_minus_pseudo_median_lag","cell_p_pseudo_lag_ge_actual"
        ]
        fields += [x for x in extras if x not in fields]
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(cell_rows)

    payload = {
        "schema": "chocho_butterfly_temporal_pseudohost_lag_null_v0.1",
        "status": "POSTHOC_IDENTITY_SPECIFIC_TEMPORAL_NULL",
        "seed": args.seed,
        "permutations": args.permutations,
        "minimum_dated_pseudos_per_host": args.minimum_dated_pseudos_per_host,
        "cells_total": len(cell_rows),
        "cells_eligible": len(eligible_cells),
        "actual": {
            "median_record_lag_years": observed_median,
            "positive_host_first_cells": observed_positive,
            "nonnegative_host_first_cells": observed_nonnegative,
        },
        "full_history_post2017_subset": post2017_summary,
        "matched_pseudohost_null": {
            "median_lag_median": float(statistics.median(null_medians)),
            "median_lag_ci95": [quantile(null_medians, 0.025), quantile(null_medians, 0.975)],
            "positive_cells_median": float(statistics.median(null_positive)),
            "positive_cells_ci95": [quantile(null_positive, 0.025), quantile(null_positive, 0.975)],
            "p_pseudo_median_lag_ge_actual": p_median,
            "p_pseudo_positive_cells_ge_actual": p_positive,
            "p_pseudo_nonnegative_cells_ge_actual": p_nonnegative,
        },
        "interpretation_rule": {
            "identity_specific_lag_supported": "Actual known hosts have a larger host-first lag than matched non-host alien plants under the conditional event-cell null.",
            "generic_alien_chronology": "Actual host timing falls within the matched pseudo-host timing distribution."
        },
        "claim_boundary": [
            "The nine focal cells were selected after apparent butterfly outcomes; this test is post-hoc and conditional on those cells.",
            "The null tests host identity/timing specificity, not whether introduced hosts caused butterfly arrival.",
            "First records are detection/digitization dates rather than establishment dates.",
            "A response-blind extension requires full-history candidate-cell outcomes before any temporal-lag claim enters the manuscript."
        ],
        "cells": cell_rows,
    }
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
