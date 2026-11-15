#!/usr/bin/env python3
"""Audit 239-butterfly paired island/mainland native-host-gap filling heterogeneity.

Data are original frozen RESOURCE OPPORTUNITIES, not observed butterfly
occurrences, reproduction or island colonization. Do not use independent-species
p values as phylogenetically adjusted confirmatory evidence. Previous headline
+54.9%, island +118.9%, mainland +40.9% are preserved unchanged.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

SOURCE = {
    "by_species": "results/island_resource_gap_filling/by_species.csv",
    "traits": "data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "pooled": "results/island_resource_gap_filling/pooled.csv",
    "published": "results/island_resource_gap_filling/summary.json",
}
CLASS_ORDER = ("1_family", "2_families", "3_to_5_families", "6plus_families")
NUMERIC_COLS = ("cells", "native", "empty", "added", "unfilled")


def get_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def get_source(root, key):
    path = root / SOURCE[key]
    return (get_csv(path) if path.suffix == ".csv" else json.loads(path.read_text()),
            hashlib.sha256(path.read_bytes()).hexdigest())


def percentile(xs, q):
    xs = sorted(xs)
    position = (len(xs) - 1) * q
    i = int(position)
    fractional = position - i
    return xs[i] * (1 - fractional) + xs[min(i + 1, len(xs) - 1)] * fractional


def exact_sign_p(positive, negative):
    """Exploratory species-as-independent diagnostic; non-ties only."""
    n = positive + negative
    if n == 0:
        return None
    tail = sum(math.comb(n, k) for k in range(min(positive, negative) + 1))
    return min(1.0, 2.0 * tail / (2**n))


def summarize(rows, bootstrap=3999, seed=20261115):
    if not rows:
        return {"n": 0}
    diffs = [r["fill_delta_pp"] for r in rows]
    pos = sum(d > 0 for d in diffs)
    neg = sum(d < 0 for d in diffs)
    tied = sum(d == 0 for d in diffs)
    mean = sum(diffs) / len(diffs)
    rng = random.Random(seed)
    draws = []
    for _ in range(bootstrap):
        draws.append(sum(diffs[rng.randrange(len(diffs))]
                         for _ in range(len(diffs))) / len(diffs))
    return {
        "n_species": len(rows),
        "island_fill_advantaged": pos,
        "mainland_fill_advantaged": neg,
        "exact_ties": tied,
        "mean_island_fill_pct": sum(r["island_fill_pct"] for r in rows) / len(rows),
        "mean_mainland_fill_pct": sum(r["mainland_fill_pct"] for r in rows) / len(rows),
        "species_equal_island_minus_mainland_fill_pp": mean,
        "median_difference_pp": percentile(diffs, 0.5),
        "paired_species_bootstrap_95pct_pp": [percentile(draws, 0.025),
                                               percentile(draws, 0.975)],
        "bootstrap_draws": bootstrap,
        "source_species_sign_test_p_two_sided_exploratory_only":
            exact_sign_p(pos, neg),
        "source_species_sign_test_invalid_for_phylogenetically_independent_p":
            True,
        "native_island_zero_species": sum(r["island_native_cells"] == 0 for r in rows),
        "mean_intro_added_island_regions": sum(r["island_added_cells"] for r in rows) / len(rows),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--out-dir", type=Path,
                   default=Path("results/island_host_filling_heterogeneity_20261115"))
    a = p.parse_args()

    raw = {}
    hashes = {}
    for key in SOURCE:
        raw[key], hashes[key] = get_source(a.repo_root, key)
    by = raw["by_species"]
    trait = raw["traits"]
    original = raw["published"]
    pooled = raw["pooled"]

    if len(by) != 239 * 2 or len(trait) != 239 or len(pooled) != 2:
        raise ValueError("Original 239×(island,mainland) source or paired trait count changed")
    traits = {t["species"]: t for t in trait}
    if len(traits) != 239:
        raise ValueError("Duplicate trait species names")
    pairs = defaultdict(dict)
    sums = defaultdict(Counter)
    for b in by:
        name, group = b["species"], b["group"]
        if name not in traits or group not in {"island", "mainland"} or group in pairs[name]:
            raise ValueError("Changed butterfly name or geography class identity")
        data = {k: int(float(b[k])) for k in NUMERIC_COLS}
        if data["native"] + data["empty"] != data["cells"]:
            raise ValueError("Source native+empty != cells")
        if data["added"] + data["unfilled"] != data["empty"]:
            raise ValueError("Source added+unfilled != original gaps")
        if data["empty"] <= 0:
            raise ValueError("Cannot compare zero-denominator gap-filling species")
        expected_cells = 126 if group == "island" else 237
        if data["cells"] != expected_cells:
            raise ValueError("Regional denominator changed")
        data["fill"] = data["added"] / data["empty"]
        pairs[name][group] = data
        sums[group].update(data)
    if set(pairs) != set(traits):
        raise ValueError("Original 239 butterfly panel has changed")
    for name, vals in pairs.items():
        if len(vals) != 2:
            raise ValueError(f"Source paired island/mainland missing {name}")
    for published in pooled:
        g = published["geography"]
        for k in NUMERIC_COLS:
            if sums[g][k] != int(float(published[k])):
                raise ValueError(f"Published pooled {g}/{k} changed")
    if (sums["island"]["native"], sums["island"]["added"],
        sums["mainland"]["native"], sums["mainland"]["added"]) != (4829, 5741, 21215, 8686):
        raise ValueError("Frozen resource headline changed")

    result = []
    for species in sorted(pairs):
        i, m = pairs[species]["island"], pairs[species]["mainland"]
        t = traits[species]
        fam = t["host_breadth_stratum"]
        if fam not in CLASS_ORDER:
            raise ValueError(f"Unknown frozen host-family class: {fam}")
        result.append({
            "species": species,
            "host_breadth_stratum": fam,
            "host_family_count": t["host_family_count"],
            "resolved_host_species": t["resolved_host_species"],
            "taxonomy_adequate": t["host_taxonomy_lower_bound_adequate"],
            "island_native_cells": i["native"],
            "island_gap_cells": i["empty"],
            "island_added_cells": i["added"],
            "mainland_native_cells": m["native"],
            "mainland_gap_cells": m["empty"],
            "mainland_added_cells": m["added"],
            "island_fill_pct": 100 * i["fill"],
            "mainland_fill_pct": 100 * m["fill"],
            "fill_delta_pp": 100 * (i["fill"] - m["fill"]),
            "island_gain_pct_native_denominator": (
                100 * i["added"] / i["native"] if i["native"] > 0 else None
            ),
            "mainland_gain_pct_native_denominator": (
                100 * m["added"] / m["native"] if m["native"] > 0 else None
            ),
        })

    if (sum(r["fill_delta_pp"] > 0 for r in result),
        sum(r["fill_delta_pp"] == 0 for r in result),
        sum(r["fill_delta_pp"] < 0 for r in result)) != (
            original["island_fill_higher"], original["tied"],
            original["mainland_fill_higher"]
        ):
        raise ValueError("Previously published signs no longer reproduce")

    main = summarize(result)
    subset = [r for r in result if r["taxonomy_adequate"] == "True"]
    high = [r for r in result if r["island_native_cells"] > 0]
    breadth = {
        k: summarize([r for r in result if r["host_breadth_stratum"] == k],
                     seed=20261115 + n)
        for n, k in enumerate(CLASS_ORDER)
    }
    if len(subset) != 215 or main["island_fill_advantaged"] != 65:
        raise ValueError("Trait taxonomy sensitivity or primary support changed")

    # Counterfactual standardisation is a deterministic accounting identity,
    # NOT a causal intervention on butterfly ecology or islands.
    projected = sum(r["island_gap_cells"] * r["mainland_fill_pct"] / 100
                    for r in result)
    actual = sum(r["island_added_cells"] for r in result)
    gain = sum(r["island_native_cells"] for r in result)
    counter = {
        "description": "Apply each butterfly's mainland gap-fill fraction to the SAME "
                       "island's number of native-host-uncovered region slots",
        "actual_island_added_butterfly_region_cells": actual,
        "island_added_if_mainland_fill_fraction_reused": projected,
        "actual_minus_same_species_standardized": actual-projected,
        "native_only_island_butterfly_region_cells": gain,
        "actual_relative_gain_pct": 100*actual/gain,
        "standardized_relative_gain_pct": 100*projected/gain,
        "counterfactual_is_descriptive_identical_species_no_movement_or_mechanism": True,
    }
    top = sorted(result, key=lambda r:r["fill_delta_pp"], reverse=True)
    report = {
        "analysis": "239 paired species island-mainland larval host gap-filling heterogeneity",
        "source_sha256": hashes,
        "geographic_unit": "WGSRPD3 Level3 region (not physical single island, not butterfly occupancy)",
        "response": "Candidate native-host gap filled by introduced range of an already-known larval host",
        "main_equal_butterfly_weight": main,
        "taxonomy_adequate_only": summarize(subset, seed=20261116),
        "island_native_host_coverage_positive_only": summarize(high, seed=20261117),
        "host_family_breadth": breadth,
        "counterfactual_same_species_mainland_fill_applied_to_island_gap_pool": counter,
        "top_10_island_fill_advantage_species": [
            {"species":r["species"], "delta_pp":r["fill_delta_pp"],
             "island_added_regions":r["island_added_cells"],
             "island_native_regions":r["island_native_cells"]}
            for r in top[:10]
        ],
        "top_10_mainland_fill_advantage_species": [
            {"species":r["species"], "delta_pp":r["fill_delta_pp"],
             "island_added_regions":r["island_added_cells"],
             "mainland_added_regions":r["mainland_added_cells"]}
            for r in top[-10:]
        ],
        "original_pooled": {r["geography"]: r for r in pooled},
        "original_area_biogeography_matched": {
            "matched_blocks": original["matched_blocks"],
            "island_fill_pct": original["matched_island_fill_pct"],
            "mainland_fill_pct": original["matched_mainland_fill_pct"],
            "difference_pp": original["matched_diff_pp"],
            "descriptive_permutation_p": original["matched_descriptive_permutation_p"],
        },
        "limits": [
            "The results concern potential presence of a listed larval host plant in a botanical region, NOT realized feeding or plant fitness",
            "Species are non-independent due to phylogeny and diet-source sharing; species bootstrap and binomial p are descriptive ONLY",
            "The island/midcontinent regional classification and plant taxon catalogue have geographic source biases",
            "Large island relative gain is dominated by low native-host coverage; this is not preferential introduction or observed butterfly range expansion",
            "The same species list is used on island and mainland; the paired differences and breedth are not different populations",
            "Port-host richness and dated butterfly indices were measured on other supports; do not pool their probabilities",
            "All results are retrospective and exploratory; no pre-registered treatment effect or natural experiment",
        ]
    }
    a.out_dir.mkdir(parents=True, exist_ok=True)
    with (a.out_dir/"paired_species_contrasts.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(result[0]))
        writer.writeheader();writer.writerows(result)
    (a.out_dir/"summary.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"main":main,"adequate":report["taxonomy_adequate_only"],
                      "by_host_family":{k:{z:v[z] for z in (
                          "n_species","island_fill_advantaged","mainland_fill_advantaged",
                          "species_equal_island_minus_mainland_fill_pp")}
                          for k,v in breadth.items()},
                      "counterfactual":counter},ensure_ascii=False))


if __name__=="__main__":
    main()
