#!/usr/bin/env python3
"""Independent DOI-observation overlap with native/introduced larval-host geography.

Use the ORIGINAL 32 selected butterflies, all 126 strict insular WGSRPD3
regions and the independent 2026 GBIF DOI occurrence-index subset. The 159
positive island-species DOI indices do NOT yield 3,873 real absences.

Conditioning nulls preserve species & region observation margins (and optionally
species × broad biogeographic Level-1 observation totals). Thus the test asks
whether the OBSERVED positive records overlap host-resource geography more
than warranted by exactly the same record counts and regional indexing bias.

No estimate here is ecological presence/absence, true habitat suitability,
historic introduction dates, port transport, larval feeding, or colonization.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

SOURCES = {
    "resources": "results/island_resource_gap_filling/resource_cells.csv.gz",
    "regions": "results/island-mainland-resource-q1/region_metrics.csv",
    "butterflies": "data/frozen/figure_sources/occurrence_resource_validation_species.csv",
    "observed": "results/q3-observed-butterfly-island-opportunity-v02/observed_species_island_cells.csv",
    "summary": "results/q3-observed-butterfly-island-opportunity-v02/q3_archival_gbif_overlap_audit.json",
}
N_REPS = 2999
SWAP_THIN = 1600
SWAP_BURN = 80000
SEED = 20261115


def csv_data(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def source(root: Path):
    tables = {}
    hashes = {}
    for k, p in SOURCES.items():
        path = root / p
        if not path.is_file():
            raise FileNotFoundError(f"Frozen source file missing: {p}")
        hashes[k] = hashlib.sha256(path.read_bytes()).hexdigest()
        tables[k] = json.loads(path.read_text()) if path.suffix == ".json" else csv_data(path)
    return tables, hashes


def metrics(edges, values, n_species, n_regions):
    s = Counter(values[sp * n_regions + region] for sp, region in edges)
    return {"native_known": s[2], "introduced_only": s[1], "no_known_host": s[0],
            "indexed_cells": len(edges)}


def summarize_null(original, draws):
    if not draws:
        raise ValueError("No preserved-margin permutations")
    values = sorted(draws)
    mean = sum(values)/len(values)
    tail = (1 + sum(x >= original for x in values))/(1 + len(values))
    lower = (1 + sum(x <= original for x in values))/(1 + len(values))
    center = abs(original-mean)
    two = (1 + sum(abs(x-mean) >= center-1e-12 for x in values))/(1 + len(values))
    return {
        "observed": original,
        "null_mean": mean,
        "null_q025": values[int(.025*(len(values)-1))],
        "null_median": values[int(.5*(len(values)-1))],
        "null_q975": values[int(.975*(len(values)-1))],
        "observed_minus_null_mean": original-mean,
        "p_one_sided_enrichment": tail,
        "p_one_sided_depletion": lower,
        "p_two_sided_centered_on_sampled_mean": two,
        "null_sd": (sum((x-mean)**2 for x in draws)/len(draws))**.5,
    }


def conditional_chain(edges, statuses, realm_by_region, n_regions, restricted=False,
                      replications=N_REPS, burn=SWAP_BURN, thin=SWAP_THIN, seed=SEED):
    rng = random.Random(seed)
    mutable = list(edges)
    present = set(edges)
    n = len(edges)
    if len(present) != n or n < 2:
        raise ValueError("No unique DOI species-region edges")
    accepted = 0
    proposed = 0

    def propose(n_attempt):
        nonlocal accepted, proposed
        for _ in range(n_attempt):
            proposed += 1
            i = rng.randrange(n)
            j = rng.randrange(n)
            if i == j:
                continue
            sp1, region1 = mutable[i]
            sp2, region2 = mutable[j]
            if sp1 == sp2 or region1 == region2:
                continue
            if restricted and realm_by_region[region1] != realm_by_region[region2]:
                continue
            next1, next2 = (sp1, region2), (sp2, region1)
            if next1 in present or next2 in present:
                continue
            present.remove(mutable[i])
            present.remove(mutable[j])
            present.add(next1)
            present.add(next2)
            mutable[i] = next1
            mutable[j] = next2
            accepted += 1

    propose(burn)
    at_burn = accepted
    draws = {"introduced_only": [], "native_known": [], "no_known_host": []}
    for _ in range(replications):
        propose(thin)
        scored = metrics(mutable, statuses, None, n_regions)
        for field in draws:
            draws[field].append(scored[field])
    original = metrics(edges, statuses, None, n_regions)
    return {
        "preserved": (
            "Exact recorded-cell totals for each butterfly, each WGSRPD3 island region, "
            "and butterfly × Level1 realm" if restricted
            else "Exact recorded-cell totals for each butterfly and each WGSRPD3 island region"
        ),
        "n_samples": replications, "burnin_proposed_edge_swaps": burn,
        "between_samples_proposed_edge_swaps": thin,
        "accepted_burnin_swaps": at_burn,
        "accepted_total_swaps": accepted, "proposed_total_swaps": proposed,
        "acceptance_rate": accepted/proposed,
        "original": original,
        "counts": {k: summarize_null(original[k], v) for k, v in draws.items()},
        "caveat": "MCMC edge swaps sample conditional index relabelings; empirical p values "
                  "are exploratory and depend on chain mixing, source selection, and sample size.",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--out-dir", type=Path,
                   default=Path("results/island_observed_host_degree_null_20261115"))
    a = p.parse_args()

    x, hashes = source(a.repo_root)
    assert x["summary"]["validation_panel"]["gbif_validation_species"] == 32
    assert x["summary"]["validation_panel"]["island_butterfly_species_cells_with_any_gbif_record"] == 159

    species = [z["species"] for z in x["butterflies"]]
    assert len(species) == len(set(species)) == 32
    island_regions = sorted(r["area_code_l3"] for r in x["regions"] if r["class"] == "island")
    assert len(island_regions) == len(set(island_regions)) == 126
    sp_i = {z: n for n, z in enumerate(species)}
    reg_i = {z: n for n, z in enumerate(island_regions)}
    level1 = {z["area_code_l3"]: z["level1"] for z in x["regions"]
              if z["class"] == "island"}
    realms = [level1[r] for r in island_regions]
    n_region = len(island_regions)

    values = {}
    per_species = Counter()
    for r in x["resources"]:
        sp, region = r["species"], r["area_code_l3"]
        if sp not in sp_i or region not in reg_i:
            continue
        key = sp_i[sp] * n_region + reg_i[region]
        native = int(r["native_hosts"])
        contempor = int(r["contemporary"])
        if native < 0 or contempor not in (0,1) or (native > 0 and contempor == 0):
            raise ValueError("Source native/contemporary resource inconsistency")
        if key in values:
            raise ValueError("Duplicate original species-region resource cell")
        values[key] = 2 if native > 0 else (1 if contempor else 0)
        per_species[sp] += 1
    if len(values) != 32 * 126 or any(per_species[sp] != 126 for sp in species):
        raise ValueError(f"Original full 32×126 resource matrix incomplete, n={len(values)}")
    flattened = [values[k] for k in range(32*126)]

    original_observations = []
    original_observed_by_species = Counter()
    for o in x["observed"]:
        if o["species"] not in sp_i or o["area_code_l3"] not in reg_i:
            raise ValueError("External DOI observation not in source panel")
        if o["observed_butterfly"] != "1":
            raise ValueError("Original GBIF observation row not a positive")
        edge = (sp_i[o["species"]], reg_i[o["area_code_l3"]])
        if edge in original_observations:
            raise ValueError("Duplicated original DOI observed species-island")
        actual_added = flattened[edge[0]*n_region + edge[1]] == 1
        if int(o["introduced_host_resource_added"]) != int(actual_added):
            raise ValueError("Original dated DOI observation vs host resource mapping drift")
        original_observations.append(edge)
        original_observed_by_species[o["species"]] += 1
    if len(original_observations) != 159:
        raise ValueError("Original 159 positive GBIF DOI pairs changed")
    original = metrics(original_observations, flattened, None, n_region)
    if original["introduced_only"] != 13:
        raise ValueError("Original observed 13 introduced resource-overlap cells changed")
    population = Counter(flattened)
    if population[1] != x["summary"]["validation_panel"]["sampled_species_added_island_opportunity_cells"]:
        raise ValueError("Original 32-species added resource opportunity denominator changed")

    globalscope = conditional_chain(
        original_observations, flattened, realms, n_region, restricted=False, seed=SEED
    )
    within_realm = conditional_chain(
        original_observations, flattened, realms, n_region, restricted=True, seed=SEED+1
    )
    categorical = {
        "no_known_host": "No known larval-host species indexed in contemporary WCVP region",
        "introduced_only": "At least one known host is introduced; none native at WGSRPD3 region scale",
        "native_known": "At least one native known larval host in that WGSRPD3 region",
    }
    marginal = {
        k: {
            "opportunity_cells": population[code],
            "indexed_butterfly_cells": original[k],
            "GBIF_DOI_cell_index_fraction_not_occupancy_rate": original[k] / population[code],
            "semantic": categorical[k],
        }
        for k, code in (("no_known_host",0),("introduced_only",1),("native_known",2))
    }
    per_sp = []
    for sp in species:
        idx = sp_i[sp]
        group = Counter(flattened[idx*n_region:(idx+1)*n_region])
        observed = [edge for edge in original_observations if edge[0] == idx]
        found = Counter(flattened[idx*n_region + edge[1]] for edge in observed)
        per_sp.append({
            "species": sp,
            "GBIF_DOI_observed_island_regions": len(observed),
            "native_known_opportunity_regions": group[2],
            "introduced_only_opportunity_regions": group[1],
            "no_known_host_regions": group[0],
            "observed_native_known": found[2],
            "observed_introduced_only": found[1],
            "observed_without_known_host": found[0],
        })
    assert sum(z["GBIF_DOI_observed_island_regions"] for z in per_sp) == 159
    report = {
        "study": "Independent archived DOI butterfly records vs potential native/introduction host coverage on islands",
        "input_sha256": hashes,
        "GBIF_DOI": x["summary"]["gbif"]["doi"],
        "butterfly_panel": 32,
        "species_with_any_GBIF_DOI_island_records": sum(original_observed_by_species[s]>0 for s in species),
        "strict_island_botanical_regions": n_region,
        "all_species_region_cells": len(flattened),
        "source_observed_positive_butterfly_region_cells": len(original_observations),
        "opportunity_by_category": marginal,
        "null_preserving_butterfly_and_island_GBIF_DOI_margins": globalscope,
        "null_also_preserving_species_by_biogeographic_level1": within_realm,
        "scope_limitations": [
            "GBIF observation index 0 is NOT biological absence; result cannot be called occupancy probability",
            "The 2026 GBIF DOI archive returned 290 fewer records than original query and sampled 32 selected butterfly species",
            "The 239-butterfly global resource cohort is separate; 32 selected source taxa are not a representative replicate",
            "Contemporary WCVP origin/status may postdate butterfly occurrence; no true historical introduced-host timing is inferred",
            "Native-known, introduced-only and no-known-host are *source-catalogue statuses*, not observed plant availability or quality at larval locations",
            "WGSRPD3 regions can include multiple islands and distinct habitats; coastal-port-local host usage is not established",
            "Fixed species and island margins control total source-record opportunities but NOT all plant/GBIF source effort or ecological suitability",
            "Level1 restricted null preserves broad geography only; it does not independently control climate, host-plant density, or insect traits",
            "Bipartite swap null replicates are correlated Markov-chain samples and p values exploratory",
            "Existing 155-species known-butterfly forecast shows temporal/spatial gains but per-species loss rows are not contained in committed forecast outputs",
            "Do not claim that an association here independently validates colonization or a species-level host-prediction × introduced-range relationship",
        ],
    }
    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n")
    with (a.out_dir / "per_species_source_crosswalk.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(per_sp[0]))
        w.writeheader()
        w.writerows(per_sp)
    print(json.dumps({
        "scope": {
            "candidate_butterflies":len(species), "island_regions":n_region,
            "GBIF_index_positive_cells":len(original_observations),
            "by_host_source_status":marginal,
        },
        "fixed_margins": {
            "global": {k:v["observed_minus_null_mean"] for k,v in globalscope["counts"].items()},
            "within_realm": {k:v["observed_minus_null_mean"] for k,v in within_realm["counts"].items()},
            "two_sided": {
                "global":globalscope["counts"]["introduced_only"]["p_two_sided_centered_on_sampled_mean"],
                "within_realm":within_realm["counts"]["introduced_only"]["p_two_sided_centered_on_sampled_mean"]
            },
            "acceptance": [globalscope["acceptance_rate"],within_realm["acceptance_rate"]]
        },
        "output":str(a.out_dir),
    }, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
