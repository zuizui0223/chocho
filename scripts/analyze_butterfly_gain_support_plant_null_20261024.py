#!/usr/bin/env python3
"""Does a negative butterfly gain vs single-host-support slope survive a plant-range null?

The source host and butterfly links are fixed. Accepted host plant introduced
regions are swapped within botanical family × WGSRPD Level1, holding the
plant-side degree structure and native host ranges intact. THIS IS NOT a
butterfly population, fitness, realized-resource or causal experiment.

Frozen protocol: provenance/reviewer_defenses/protocols/GAIN_SUPPORT_PLANT_RANGE_NULL_20261024.md
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from analyze_butterfly_shared_host_structural_null import (
    attempt_range_swaps, build_model, read_csv,
)

ROOT = Path(__file__).resolve().parents[1]
BLOBS = {
    "data/frozen/figure_sources/anthropogenic_species_metrics.csv":
        "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
    "data/frozen/figure_sources/host_contribution_metrics.csv":
        "4cf4851499ca2233ff2adadf1c4ee7f0e863f23d",
    "data/frozen/figure_sources/s1_resource_descriptors.csv":
        "09541b5894529cc8629c7515ea9022ded9edead2",
}
SEEDS = (20261024, 20261025)
DRAWS = 199
BURN_FACTOR = 5
SPACING_FACTOR = 2
OBSERVED_PRIMARY_SLOPE = -0.22071872904922496
OBSERVED_ADEQUATE_SLOPE = -0.2711287143300629
LEVEL1_SHA_COMMIT = "52da7828aba9d461dd133c27b3bd7a4407161f54"


def git_file_sha(raw: bytes):
    return hashlib.sha1(
        b"blob " + str(len(raw)).encode("ascii") + b"\x00" + raw
    ).hexdigest()


def frozen_csv(path: str):
    raw = (ROOT / path).read_bytes()
    observed_sha = git_file_sha(raw)
    if observed_sha != BLOBS[path]:
        raise ValueError(f"FROZEN_BLOB_CHANGED: {path}: {observed_sha}")
    return read_csv(ROOT / path)


def species_descriptors():
    mapping = {}
    for data in frozen_csv("data/frozen/figure_sources/anthropogenic_species_metrics.csv"):
        name = data["species"]
        if name in mapping:
            raise ValueError("Duplicate species in frozen main resource descriptors")
        mapping[name] = {
            "native": int(data["native_resource_units"]),
            "contemporary": int(data["contemporary_resource_units"]),
            "added": int(data["introduced_added_units"]),
            "host_count": int(data["resolved_host_species"]),
            "host_family_count": float(data["host_family_count"]),
            "adequate": data["host_taxonomy_lower_bound_adequate"] == "True",
        }
    hrows = frozen_csv("data/frozen/figure_sources/host_contribution_metrics.csv")
    h = {}
    for r in hrows:
        sp = r["species"]
        if sp in h:
            raise ValueError("Duplicate host-contribution row")
        h[sp] = r
    desc = frozen_csv("data/frozen/figure_sources/s1_resource_descriptors.csv")
    if len(mapping) != 239 or len(h) != 239 or len(desc) != 339 or set(h) != set(mapping):
        raise ValueError(f"Original species panel changed: {len(mapping)}/{len(h)}/{len(desc)}")
    return mapping, h, desc


def read_model(args):
    metrics, host_stats, descriptors = species_descriptors()
    lookup = json.loads(args.level1_json.read_text(encoding="utf-8"))
    if lookup.get("source_commit") != LEVEL1_SHA_COMMIT:
        raise ValueError("Level1 regional lookup source changed")
    model = build_model(
        descriptors, read_csv(args.insect_host_csv),
        read_csv(args.native_csv), read_csv(args.contemporary_csv),
        lookup["units"],
    )
    if len(model["species"]) != 239:
        raise ValueError("229/239 taxa source mismatch")
    if model["total_introduced_edges"] != 30373:
        raise ValueError(f"Original introduced plant edges changed: {model['total_introduced_edges']}")
    if not 28000 <= model["cumulative"][-1] <= 30373:
        raise ValueError("Insufficient source edge strata")
    for name, hostlist, nativebits in zip(
        model["species"], model["species_hosts"], model["native_species"]
    ):
        if name not in metrics:
            raise ValueError(f"Unmatched focal butterfly {name}")
        row = metrics[name]
        if (nativebits.bit_count() != row["native"] or
            len(hostlist) != row["host_count"]):
            raise ValueError(f"Frozen host-taxonomy/native mismatch: {name}")
    model["frozen_metrics"] = metrics
    model["frozen_host_stats"] = host_stats
    return model


def bits_for_introduced(host_introductions: dict[str, set[int]]):
    masks = {}
    for host, cols in host_introductions.items():
        b = 0
        for unit in cols:
            b |= 1 << unit
        masks[host] = b
    return masks


def individual_resource_states(model, host_introductions):
    masks = bits_for_introduced(host_introductions)
    states = []
    for species, hosts, nativebits in zip(
        model["species"], model["species_hosts"], model["native_species"]
    ):
        masked = 0
        multiples = 0
        for host in hosts:
            introduced = masks.get(host, 0) & ~nativebits
            multiples |= masked & introduced
            masked |= introduced
        gained = masked.bit_count()
        single = (masked & ~multiples).bit_count()
        if gained and single > gained:
            raise ValueError("Single known host support count exceeds total gain")
        states.append((species, gained, single, nativebits.bit_count()))
    return states


def fit(rows):
    """Matched original within-catalogue species-level adjusted OLS coefficient."""
    if len(rows) < 150:
        return None
    x = np.array([
        [1, r["gain"], r["log_native"], r["log_hosts"], r["families"]]
        for r in rows
    ], dtype=float)
    y = np.array([r["single_share"] for r in rows], dtype=float)
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    if rank != x.shape[1] or not np.isfinite(beta).all():
        return None
    return float(beta[1])


def regression_states(model, states, original_nonzero=None):
    rows = []
    total_gain = 0
    total_single = 0
    nonzero = set()
    for species, gained, single, native in states:
        meta = model["frozen_metrics"][species]
        if native != meta["native"]:
            raise ValueError(f"Native range unexpectedly moved: {species}")
        total_gain += gained
        total_single += single
        if gained <= 0:
            continue
        nonzero.add(species)
        rows.append({
            "species": species,
            "gain": math.log((native + gained) / native),
            "single_share": single / gained,
            "log_native": math.log(native),
            "log_hosts": math.log1p(meta["host_count"]),
            "families": meta["host_family_count"],
            "adequate": meta["adequate"],
        })
    full = fit(rows)
    adequate = fit([r for r in rows if r["adequate"]])
    if full is None or adequate is None:
        raise ValueError("Source panel loses full-rank 239-species adjusted regression")
    result = {
        "beta": full,
        "adequate_beta": adequate,
        "positive_species": len(nonzero),
        "adequate_positive_species": sum(r["adequate"] for r in rows),
        "all_added_units": total_gain,
        "single_host_supported_units": total_single,
    }
    if original_nonzero is not None:
        if original_nonzero <= nonzero:
            fixed = fit([r for r in rows if r["species"] in original_nonzero])
            result["fixed_original_positive_beta"] = fixed
        else:
            result["fixed_original_positive_beta"] = None
        result["original_positive_species_remaining"] = len(original_nonzero & nonzero)
    return result, nonzero


def source_checks(model, states):
    if sum(n for _, n, _, _ in states) != 14553:
        raise ValueError("Observed 14553 added butterfly-region units not reproduced")
    if sum(n for _, _, n, _ in states) != 8574:
        raise ValueError("Observed 8574 single-host-supported cells not reproduced")
    if sum(native for _, _, _, native in states) != 26530:
        raise ValueError("Observed 26530 native host cells not reproduced")
    for name, gain, single, _ in states:
        m = model["frozen_metrics"][name]
        h = model["frozen_host_stats"][name]
        if m["added"] != gain or m["contemporary"] != m["native"] + gain:
            raise ValueError(f"Individual geographic source count drift: {name}")
        if gain:
            expected = 1-float(h["fraction_added_units_with_multiple_contributing_hosts"])
            if abs(single/gain - expected) > 1e-9:
                raise ValueError(f"Individual host support share disagrees: {name}")
    out, nonzero = regression_states(model, states)
    if (len(nonzero) != 206 or
        abs(out["beta"] - OBSERVED_PRIMARY_SLOPE) > 1e-8 or
        abs(out["adequate_beta"] - OBSERVED_ADEQUATE_SLOPE) > 1e-8):
        raise ValueError(f"Observed source regression cannot reproduce reference: {out}")
    return out, nonzero


def family_region_margins(model, intro):
    fam = model["host_families"]
    lev = model["level1"]
    host_by_continent = Counter()
    region_family = Counter()
    introduced_edges = 0
    for host, units in intro.items():
        for u in units:
            if u in model["native_sets"].get(host, set()):
                raise ValueError("Introduced range intersects native plant footprint")
            host_by_continent[(host, lev[u])] += 1
            region_family[(u, fam[host])] += 1
            introduced_edges += 1
    return {
        "host_by_continent": host_by_continent,
        "region_family": region_family,
        "edges": introduced_edges,
    }


def margin_assert(model, introduced, reference):
    here = family_region_margins(model, introduced)
    if here != reference:
        raise ValueError("Botanical family / recipient or per-plant-continent margin broken")


def ess_estimate(history):
    """Positive sequence approximate ESS with the same procedure as existing nulls."""
    a = np.asarray(history, dtype=float)
    n = len(a)
    if n < 30 or not np.isfinite(a).all():
        return {"n": n, "ess": 0, "std": None, "median": None}
    x = a-a.mean()
    den = float(x@x)
    if den <= 1e-20:
        return {"n": n, "ess": 0, "std": 0, "median": float(np.median(a))}
    coefficients = []
    for lag in range(1, min(n//3, 30)):
        rho = float(x[:-lag]@x[lag:]/den)
        if rho <= 0:
            break
        coefficients.append(rho)
    return {
        "n": n,
        "ess": n/(1 + 2*sum(coefficients)),
        "lag1": float(x[:-1]@x[1:]/den),
        "median": float(np.median(a)),
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)),
        "early_late_drift": float(np.mean(a[-n//4:])-np.mean(a[:n//4])),
    }


def chain(model, seed, original_nonzero, reference_margins, draws=DRAWS):
    intro = {h: set(v) for h, v in model["introduced_sets"].items()}
    strata = copy.deepcopy(model["groups"])
    cumulative = model["cumulative"]
    E = model["total_introduced_edges"]
    rng = random.Random(seed)
    before = {h: frozenset(v) for h, v in intro.items()}
    accepted_burn = attempt_range_swaps(
        strata, cumulative, model["native_sets"], intro, BURN_FACTOR*E, rng
    )
    margin_assert(model, intro, reference_margins)
    moved = sum(len(before[h] ^ intro[h]) for h in intro)/(2*E)
    samples = []
    spacing_accept = []
    for j in range(draws):
        spacing_accept.append(attempt_range_swaps(
            strata, cumulative, model["native_sets"], intro, SPACING_FACTOR*E, rng
        ))
        if j in (0, draws//2, draws-1):
            margin_assert(model, intro, reference_margins)
        s, _ = regression_states(
            model, individual_resource_states(model, intro), original_nonzero
        )
        samples.append(s)
    stats = {
        key: ess_estimate([sample[key] for sample in samples])
        for key in ("beta", "adequate_beta")
    }
    all_values = [float(s["beta"]) for s in samples]
    all_adequate = [float(s["adequate_beta"]) for s in samples]
    fixed_values = [s["fixed_original_positive_beta"] for s in samples]
    fixed_eval = all(v is not None for v in fixed_values)
    return {
        "seed": seed,
        "draws": draws,
        "introduced_edges": E,
        "viable_edges": cumulative[-1],
        "accepted_burn": accepted_burn,
        "moved_fraction_after_burn": moved,
        "min_accepted_spacing": min(spacing_accept),
        "median_accepted_spacing": float(np.median(spacing_accept)),
        "sampling_statistics": stats,
        "full_panel_five_number": np.quantile(all_values,[0,.025,.5,.975,1]).tolist(),
        "adequate_panel_five_number": np.quantile(all_adequate,[0,.025,.5,.975,1]).tolist(),
        "fixed_positive_panel_evaluable_all_draws": fixed_eval,
        "fixed_positive_panel_five_number": (
            np.quantile(fixed_values,[0,.025,.5,.975,1]).tolist()
            if fixed_eval else None
        ),
        "original_positive_species_remaining_range": [
            min(s["original_positive_species_remaining"] for s in samples),
            max(s["original_positive_species_remaining"] for s in samples),
        ],
        "positive_species_range": [
            min(s["positive_species"] for s in samples),
            max(s["positive_species"] for s in samples),
        ],
        "adequate_positive_species_range": [
            min(s["adequate_positive_species"] for s in samples),
            max(s["adequate_positive_species"] for s in samples),
        ],
        "null_added_resource_units_range": [
            min(s["all_added_units"] for s in samples),
            max(s["all_added_units"] for s in samples),
        ],
        "null_single_supported_units_range": [
            min(s["single_host_supported_units"] for s in samples),
            max(s["single_host_supported_units"] for s in samples),
        ],
        "full_null_history": all_values,
        "adequate_null_history": all_adequate,
        "fixed_original_positive_history_if_complete": fixed_values if fixed_eval else None,
    }


def interpret_chains(model, observed, chains):
    quality = []
    for c in chains:
        primary = c["sampling_statistics"]["beta"]
        adequate = c["sampling_statistics"]["adequate_beta"]
        good = (
            c["accepted_burn"] >= 100 and
            c["moved_fraction_after_burn"] >= 0.10 and
            c["min_accepted_spacing"] >= 100 and
            primary["ess"] >= 30 and adequate["ess"] >= 30
        )
        quality.append({
            "seed": c["seed"],
            "qualified_local": bool(good),
            "source_movement_ge10pct": c["moved_fraction_after_burn"] >= 0.10,
            "min_swap_accept_ge100": c["min_accepted_spacing"] >= 100,
            "full_ess_ge30": primary["ess"] >= 30,
            "adequate_ess_ge30": adequate["ess"] >= 30,
        })
    median_full = [c["sampling_statistics"]["beta"]["median"] for c in chains]
    std_full = [c["sampling_statistics"]["beta"]["std"] for c in chains]
    median_adequate = [c["sampling_statistics"]["adequate_beta"]["median"] for c in chains]
    std_adequate = [c["sampling_statistics"]["adequate_beta"]["std"] for c in chains]
    concord_full = abs(median_full[0]-median_full[1]) <= 2*max(std_full)
    concord_adequate = abs(median_adequate[0]-median_adequate[1]) <= 2*max(std_adequate)
    qualified = all(c["qualified_local"] for c in quality) and concord_full and concord_adequate
    comparisons = {}
    for name, observed_slope, key in (
        ("primary", observed["beta"], "full_null_history"),
        ("adequate_catalogue", observed["adequate_beta"], "adequate_null_history"),
    ):
        per_chain = []
        for c in chains:
            values = np.asarray(c[key], dtype=float)
            per_chain.append({
                "observed_beta": observed_slope,
                "null_q025": float(np.quantile(values, .025)),
                "null_median": float(np.median(values)),
                "null_q975": float(np.quantile(values, .975)),
                "obs_minus_null_median": float(observed_slope-np.median(values)),
                "serial_left_tail_fraction_if_qualified": (
                    float((1+np.sum(values<=observed_slope))/(len(values)+1))
                    if qualified else None
                ),
                "serial_right_tail_fraction_if_qualified": (
                    float((1+np.sum(values>=observed_slope))/(len(values)+1))
                    if qualified else None
                ),
            })
        comparisons[name] = {
            "qualified": qualified,
            "source_conditioned_same_estimand": True,
            "chains": per_chain,
        }
    return {
        "numerically_qualified": qualified,
        "quality_per_chain": quality,
        "between_chain_full_median_agreement": concord_full,
        "between_chain_adequate_median_agreement": concord_adequate,
        "result": comparisons,
    }


def run(args):
    model = read_model(args)
    observed_states = individual_resource_states(model, model["introduced_sets"])
    observed, observed_nonzero = source_checks(model, observed_states)
    source_margins = family_region_margins(model, model["introduced_sets"])
    if source_margins["edges"] != 30373:
        raise ValueError("Original total added plant×region edges changed")
    chains = [
        chain(model, seed, observed_nonzero, source_margins, draws=args.draws)
        for seed in SEEDS
    ]
    interpretation = interpret_chains(model, observed, chains)
    return {
        "schema": "chocho_gain_vs_host_support_structural_plant_null_20261024",
        "status": (
            "QUALIFIED_POSTHOC_SOURCE_CONDITIONAL_PLANT_RANGE_NULL"
            if interpretation["numerically_qualified"]
            else "NOT_EVALUABLE_NUMERICAL_MIXING_GATE"
        ),
        "source": {
            "HOSTS_git_commit": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "WCVP_git_commit": "65bed76bae9d644ccb6ad200c05f9f5071d89e05",
            "level1_git_commit": LEVEL1_SHA_COMMIT,
            "local_input_git_blobs": BLOBS,
        },
        "focal_butterflies": 239,
        "source_botanical_introduced_host_area_edges": source_margins["edges"],
        "region_count": len(model["level1"]),
        "introduced_edges_in_viable_swap_strata": model["cumulative"][-1],
        "randomization": {
            "preserves": [
                "all accepted butterfly-to-host plant identities",
                "native plant×geographic-region support as structural exclusions",
                "each introduced plant's geographic count within a WGSRPD Level1 continent",
                "introduced plants per recipient region per botanical family",
                "global introduced plant×region edge budget",
            ],
            "does_not_preserve": [
                "per-butterfly added resource count; it changes with plant swaps",
                "plant niche climate or local distributional spatial contiguity",
                "actual butterfly occurrence, larval feeding, overwintering or demography",
            ],
            "seeds": list(SEEDS),
            "burn_factor_times_edges": BURN_FACTOR,
            "spacing_factor_times_edges": SPACING_FACTOR,
            "sample_draws_per_chain": args.draws,
            "posthoc_after_observed_slope": True,
        },
        "observed": observed,
        "numerical_evaluation": interpretation,
        "chains": chains,
        "interpretive_ceiling": (
            "A departure from this null would be plant-range/known-host geography only; "
            "not an identified butterfly biological preference or demographic stability. "
            "Nominal serial tails cannot be independent permutation p-values. "
            "The null is conditioned on botanical family and continent, not fine climate."
        ),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--insect-host-csv", type=Path, required=True)
    p.add_argument("--native-csv", type=Path, required=True)
    p.add_argument("--contemporary-csv", type=Path, required=True)
    p.add_argument("--level1-json", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--draws", type=int, default=DRAWS)
    a = p.parse_args()
    if a.draws != DRAWS:
        raise ValueError("Changed Monte Carlo sample size not permitted by frozen protocol")
    result = run(a)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    brief = {
        "status": result["status"],
        "observed_beta": result["observed"]["beta"],
        "observed_adequate_beta": result["observed"]["adequate_beta"],
        "null_comparison": result["numerical_evaluation"]["result"],
        "quality": result["numerical_evaluation"]["quality_per_chain"],
        "source_intro_edges": result["source_botanical_introduced_host_area_edges"],
    }
    print(json.dumps(brief, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
