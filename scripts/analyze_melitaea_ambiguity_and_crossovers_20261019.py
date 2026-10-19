#!/usr/bin/env python3
"""Exploratory family-level host effect stability and measurement sensitivity.

Read MD5-verified 2020 Melitaea cinxia source. No synthetic outcomes
are ever reported as biological observation. Model simulations are
sensitivity-only and require unmeasured within-group dependence assumptions.

Run AFTER download_exact_cross_butterfly_sources_20261019.py.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

from analyze_cross_butterfly_lifestage_raw_20261019 import csv_records

SOURCE = "melitaea/Survival_data.csv"
DOI = "10.5061/dryad.kwh70rz2m"
TEMPERATURES = ("28", "30", "32", "34")
PLANTS = ("Pl", "Vs")


def source_rows():
    rows = csv_records(SOURCE)  # already requires exact source size and MD5
    assert len(rows) == 1680
    families = sorted({r["Family"] for r in rows}, key=int)
    assert len(families) == 14
    assert {r["Dead"] for r in rows} == {"0", "3"}
    group_n = defaultdict(int)
    for r in rows:
        group_n[(r["Family"], r["Temperature"], r["Plant"])] += 1
    assert len(group_n) == 112 and set(group_n.values()) == {15}
    ambiguous = [
        i for i, r in enumerate(rows)
        if (r["Dead"] == "0" and r["Time2"] != "NA")
        or (r["Dead"] == "3" and (r["Time2"] == "NA" or float(r["Time2"]) < float(r["Time"])))
    ]
    assert len(ambiguous) == 4, "The source's four contradictory event intervals changed"
    assert all(rows[i]["Dead"] == "0" for i in ambiguous)
    return rows, families, ambiguous


def effect_by_family(rows, families, uncertain_success: frozenset[int]):
    # As an uncertainty scenario, only the four anomalous rows may change
    # from event-coded success to failure. Denominators remain initial 15.
    counts = defaultdict(int)
    for i, r in enumerate(rows):
        if r["Dead"] == "0" and (r["Time2"] == "NA" or i in uncertain_success):
            counts[(r["Family"], r["Temperature"], r["Plant"])] += 1
    by_family = []
    by_temperature = {t: [] for t in TEMPERATURES}
    flip_n = 0
    for family in families:
        effects = []
        for temp in TEMPERATURES:
            pl = counts[(family, temp, "Pl")] / 15.0
            vs = counts[(family, temp, "Vs")] / 15.0
            diff = vs - pl
            effects.append(diff)
            by_temperature[temp].append(diff)
        flip_n += any(x > 0 for x in effects) and any(x < 0 for x in effects)
        by_family.append(statistics.mean(effects))
    return by_family, by_temperature, flip_n


def symmetric_signflip_diagnostic(effects):
    """Exact 2^14 numerical symmetry diagnostic; NOT an actual randomization p."""
    observed = abs(sum(effects))
    extreme = 0
    for signs in itertools.product((-1, +1), repeat=len(effects)):
        total = abs(sum(v * sign for v, sign in zip(effects, signs)))
        if total + 1e-10 >= observed:
            extreme += 1
    return extreme / 2**len(effects)


def bootstrap_family_mean(effects, seed=20261019, n=9999):
    rng = random.Random(seed)
    values = []
    count = len(effects)
    for _ in range(n):
        values.append(sum(effects[rng.randrange(count)] for _ in range(count)) / count)
    values.sort()
    return [round(values[int((n - 1) * q)], 8) for q in (0.025, 0.975)]


def all_binary_ambiguity_scenarios(rows, families, ambiguous):
    scenarios = []
    for bits in itertools.product((0, 1), repeat=len(ambiguous)):
        successes = frozenset(ambiguous[i] for i, v in enumerate(bits) if v)
        by_family, by_temperature, n_flips = effect_by_family(rows, families, successes)
        scenarios.append({
            "ambiguous_success_pattern": "".join(map(str, bits)),
            "estimate_Vs_minus_Pl": round(statistics.mean(by_family), 8),
            "family_crossovers": n_flips,
            "symmetric_signflip_two_sided": round(symmetric_signflip_diagnostic(by_family), 8),
        })
    assert len(scenarios) == 16
    return scenarios


def logodds_binomial_smooth(successes, n):
    return math.log((successes + 0.5) / (n - successes + 0.5))


def fit_illustrative_constant_family_host_effect(rows, families):
    """Plug-in null: per-family host odds constant across temps, with common
    all-family temperature host adjustments and flexible family×temp baseline.
    Fitted from this SAME dataset: a *posterior predictive-style diagnostic*,
    not a valid held-out discovery test or exact randomized inference.
    """
    counts = defaultdict(int)
    for r in rows:
        if r["Dead"] == "0":
            counts[(r["Family"], r["Temperature"], r["Plant"])] += 1

    across_temp = []
    for temp in TEMPERATURES:
        pl = sum(counts[(f, temp, "Pl")] for f in families)
        vs = sum(counts[(f, temp, "Vs")] for f in families)
        across_temp.append(logodds_binomial_smooth(vs, 210) -
                           logodds_binomial_smooth(pl, 210))
    temp_offsets = [beta - statistics.mean(across_temp) for beta in across_temp]

    cell_probs = []
    for family in families:
        pl_total = sum(counts[(family, t, "Pl")] for t in TEMPERATURES)
        vs_total = sum(counts[(family, t, "Vs")] for t in TEMPERATURES)
        family_beta = (
            logodds_binomial_smooth(vs_total, 60) -
            logodds_binomial_smooth(pl_total, 60)
        )
        by_temp = []
        for idx, t in enumerate(TEMPERATURES):
            m = (counts[(family, t, "Pl")] + counts[(family, t, "Vs")]) / 30
            odds_ratio = math.exp(family_beta + temp_offsets[idx])
            if m == 1:
                p_pl = p_vs = 1.0
            elif m == 0:
                p_pl = p_vs = 0.0
            else:
                lo, hi = 0.0, 1.0
                for _ in range(55):
                    p_pl = (lo + hi) / 2
                    p_vs = odds_ratio * p_pl / (1 - p_pl + odds_ratio * p_pl)
                    if (p_pl + p_vs) / 2 > m:
                        hi = p_pl
                    else:
                        lo = p_pl
                p_pl = (lo + hi) / 2
                p_vs = odds_ratio * p_pl / (1 - p_pl + odds_ratio * p_pl)
            by_temp.append((p_pl, p_vs))
        cell_probs.append(by_temp)
    return cell_probs


def simulate_reversal_count(null_probs, correlations=(0, 0.05, 0.15),
                            n=9999, seed=20261019):
    """Independent-binomial vs hypothetical beta-binomial group clustering.

    rho is assumed intragroup correlation among 15 larvae, NOT fitted.
    Conditional outcome probabilities are fitted to the same data. The
    tail fraction is a model-sensitivity calibration, NOT ecological p-value.
    """
    def binomial_15(p, rng):
        return sum(rng.random() < p for _ in range(15))

    def group_prob(p, rho, rng):
        if p in (0.0, 1.0) or rho == 0:
            return p
        a_plus_b = 1 / rho - 1
        return rng.betavariate(p * a_plus_b, (1 - p) * a_plus_b)

    summaries = []
    for i, rho in enumerate(correlations):
        rng = random.Random(seed + 101 * i)
        tallies = []
        for _ in range(n):
            family_crossovers = 0
            for family_cells in null_probs:
                plus = minus = False
                for p_pl, p_vs in family_cells:
                    pl = binomial_15(group_prob(p_pl, rho, rng), rng)
                    vs = binomial_15(group_prob(p_vs, rho, rng), rng)
                    if vs > pl:
                        plus = True
                    elif vs < pl:
                        minus = True
                if plus and minus:
                    family_crossovers += 1
            tallies.append(family_crossovers)
        tallies.sort()
        summaries.append({
            "assumed_larvae_within_group_correlation": rho,
            "simulated_mean_family_crossovers": round(statistics.mean(tallies), 5),
            "simulated_95pct_middle_range": [
                tallies[int((n - 1) * q)] for q in (0.025, 0.975)
            ],
            "simulated_fraction_ge_7": round(sum(z >= 7 for z in tallies) / n, 5),
            "simulated_fraction_ge_8": round(sum(z >= 8 for z in tallies) / n, 5),
            "simulations": n,
        })
    return summaries


def build(nsim=9999):
    rows, families, ambiguous = source_rows()
    primary, paired_temperature, primary_flip = effect_by_family(
        rows, families, frozenset(ambiguous))
    strict, _, strict_flip = effect_by_family(
        rows, families, frozenset())
    cases = all_binary_ambiguity_scenarios(rows, families, ambiguous)
    family_loocv = [
        statistics.mean([x for k, x in enumerate(primary) if k != i])
        for i in range(len(primary))
    ]
    model = fit_illustrative_constant_family_host_effect(rows, families)
    simulations = simulate_reversal_count(model, n=nsim)
    return {
        "schema": "chocho_melitaea_life_stage_uncertainty_crossover_v1",
        "source": DOI, "source_md5_checked_by": SOURCE,
        "status": "EXPLORATORY_SOURCE_CROSSCHECK_AND_MODEL_SENSITIVITY",
        "individual_larvae_n": len(rows),
        "mothers_n": len(families),
        "clutches_n": len({r["Clutch.ID"] for r in rows}),
        "end_point": "arrival alive at diapause, never adult recruitment",
        "ambiguous_records": [
            {k: rows[i][k] for k in ("Clutch.ID", "Family", "Temperature", "Plant", "Time", "Time2", "Dead")}
            for i in ambiguous
        ],
        "ambiguity_scenarios_n": len(cases),
        "all_scenarios": cases,
        "partial_identification_host_advantage": [
            min(x["estimate_Vs_minus_Pl"] for x in cases),
            max(x["estimate_Vs_minus_Pl"] for x in cases)
        ],
        "partial_identification_crossovers_range": [
            min(x["family_crossovers"] for x in cases),
            max(x["family_crossovers"] for x in cases)
        ],
        "partial_identification_signflip_diagnostic_range": [
            min(x["symmetric_signflip_two_sided"] for x in cases),
            max(x["symmetric_signflip_two_sided"] for x in cases)
        ],
        "primary_event_coded": {
            "mother_equal_weighted_Vs_minus_Pl": round(statistics.mean(primary), 8),
            "mother_cluster_bootstrap_95ci": bootstrap_family_mean(primary),
            "mother_LOO_point_range": [
                round(min(family_loocv), 8), round(max(family_loocv), 8)],
            "exact_symmetric_signflip_two_sided_DIAGNOSTIC_NOT_CAUSAL_P":
                round(symmetric_signflip_diagnostic(primary), 8),
            "families_with_opposite_host_effect_sign_across_four_temperatures": primary_flip,
            "family_mean_effects": [round(x, 8) for x in primary],
            "by_temperature_effects": {
                t: round(statistics.mean(paired_temperature[t]), 8)
                for t in TEMPERATURES
            }
        },
        "all_four_ambiguous_coded_not_successful_sensitivity": {
            "mother_equal_weighted_Vs_minus_Pl": round(statistics.mean(strict), 8),
            "mother_cluster_bootstrap_95ci": bootstrap_family_mean(strict, seed=20261020),
            "families_with_opposite_sign": strict_flip,
        },
        "independent_larva_null_and_group_dependence_sensitivity": {
            "description": "Fit family-specific constant host log odds, common host×temp offset and family×temp baseline from same data. Simulate observed direction switches; allow unmeasured within-group beta-binomial dependence rho. NOT an exact null test, not external validation.",
            "assumed_group_correlations_not_estimated": [0, 0.05, 0.15],
            "fits_same_data": True,
            "null_model_fitted_then_simulated": True,
            "results": simulations,
        },
        "interpretation_boundaries": [
            "16 ambiguity patterns are hypothetical recodings, NOT 16 studies.",
            "Source MD5 is checked by source_rows; no synthetic larva observation enters outcome summaries.",
            "A family sign change from two noisy binomial rates is not itself genetic reaction norm evidence.",
            "Fitted-null tail fractions depend on unmeasured shared group variation, cannot be promoted as a novel biological interaction or confirmatory significance.",
            "Host effects on diapause should not be pooled with field first-instar survival, lab survivor-only development, Vanessa death time, or global WCVP host availability.",
            "No oviposition choice, egg diversion, next-generation recruits, new host colonisation, or causal invasive-plant demographic effect was measured.",
        ]
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", help="JSON output filename")
    p.add_argument("--nsim", type=int, default=9999)
    args = p.parse_args()
    result = build(nsim=args.nsim)
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "mean_host_effect": result["primary_event_coded"]["mother_equal_weighted_Vs_minus_Pl"],
        "ambiguity_bounds": result["partial_identification_host_advantage"],
        "crossovers_primary": result["primary_event_coded"]["families_with_opposite_host_effect_sign_across_four_temperatures"],
        "crossovers_ambiguity_range": result["partial_identification_crossovers_range"],
        "hypothetical_null": result["independent_larva_null_and_group_dependence_sensitivity"]["results"]
    }, sort_keys=True))
