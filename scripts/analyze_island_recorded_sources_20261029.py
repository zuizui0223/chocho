#!/usr/bin/env python3
"""Post-pilot conditional diagnostic of prior recorded-source coast proximity.

This uses the already observed 96-island panel. It does not estimate arrival,
colonization, establishment, causal port effects or a biological absence rate.
There are no network calls. Original input bytes and the new plan are pinned.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import sys
import zipfile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.stats import norm, rankdata

ROOT = Path(__file__).resolve().parents[1]
PANEL = "data/exploratory/island_expansion_20261029/"
TAXONOMY = "data/exploratory/island_taxonomy_20261029/"
PLAN = "provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_PLAN_20261029.md"
PINS = {
    PLAN: "67c25d1f28ff9f98def1122fc0fa0b4782062e08cf3dae5497cd3806795a82e6",
    PANEL + "expanded_island_panel_selection.json": "5ccc868b89973b8142a7fdb4b3014c4394f134596ae2ef58e1179eaab68f5243",
    PANEL + "expanded_island_panel_manifest.json": "612c4f5e981a045637fcd66cf1353241f5ef70ddbfd63da9719e92c504404040",
    PANEL + "expanded_island_panel_species_key_facets.csv": "061f9ce3cf21e4d64572456e0ffbbb88022fbccfe9aa924d5f02062027c97db6",
    TAXONOMY + "original239_to_colxr_bridge.csv": "cd8b2cb8e0279a5b940e46267c67f0bb9843d2e1c54f38bc405ec30a95fe6eb2",
    TAXONOMY + "taxonomy_bridge_manifest.json": "43cf640edaee99e12fa487947213c2451935e22f2396193492f1fff7ae5e85e9",
}
WINDOWS = ("2010_2017", "2018_2025")
SEEDS = (20261029, 20261030, 20261031, 20261032)
BURN_IN = 100_000
SAMPLES = 4_000
THIN = 50
LAZY_HOLD = .5
TIE_TOLERANCE = 1e-12


def sha(body):
    return hashlib.sha256(body).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def utc():
    return datetime.now(timezone.utc).isoformat()


def validate_occurrence_snapshot(roster, manifest):
    keys = {r["island_key"] for r in roster}
    if len(keys) != len(roster):
        raise ValueError("Duplicate roster identity")
    expected = {(k, w) for k in keys for w in WINDOWS}
    result = {}
    fields = ("SPECIES_KEY", "DATASET_KEY", "YEAR", "BASIS_OF_RECORD")
    for query in manifest["results"]:
        key = query["island_key"], query["window"]
        if key not in expected or key in result or query["status"] != "COMPLETE":
            raise ValueError("Missing, duplicated, out-of-roster or failed island-window query")
        if not query.get("first_page_recheck_passed"):
            raise ValueError("Independent first-page verification is missing")
        count = query["count"]
        if type(count) is not int or count < 0:
            raise ValueError("Invalid complete query count")
        facets = query["facets"]
        for field in fields:
            bins = facets[field]
            if any(not isinstance(k, str) or not k or type(v) is not int or v <= 0 for k, v in bins.items()):
                raise ValueError("Invalid positive facet bin")
            residual = query["unassigned_facet_records"][field]
            if type(residual) is not int or residual < 0 or sum(bins.values()) + residual != count:
                raise ValueError("Facet margin does not reconcile with complete count")
            if field in ("DATASET_KEY", "YEAR") and residual:
                raise ValueError("Required year/dataset margin is incomplete")
        lower, upper = map(int, query["window"].split("_"))
        if any(not lower <= int(y) <= upper for y in facets["YEAR"]):
            raise ValueError("Year facet outside its declared window")
        result[key] = dict(facets["SPECIES_KEY"])
    if set(result) != expected:
        raise ValueError("The complete frozen island-window Cartesian product is required")
    return result


def load_inputs(root):
    for relative, expected in PINS.items():
        if sha((root / relative).read_bytes()) != expected:
            raise ValueError(f"Pinned input changed: {relative}")
    selection = json.loads((root / PANEL / "expanded_island_panel_selection.json").read_bytes())
    manifest = json.loads((root / PANEL / "expanded_island_panel_manifest.json").read_bytes())
    if sha(canonical(selection["selection"])) != selection["selection_sha256"]:
        raise ValueError("Selection canonical content does not match its lock")
    if manifest["selection_sha256"] != selection["selection_sha256"] or manifest["status"] != "FINISHED":
        raise ValueError("Occurrence panel differs from finished frozen selection")
    if manifest["taxonomy"]["checklistKey"] != "7ddf754f-d193-4cc9-b351-99906754a03b" or manifest["taxonomy"]["taxonKey"] != "5G9":
        raise ValueError("COL XR Papilionoidea namespace changed")
    roster = selection["selection"]["selected"]
    if len(roster) != 96:
        raise ValueError("Expected the original fixed 96-island panel")
    counts = validate_occurrence_snapshot(roster, manifest)
    facet_counts = {key: {} for key in counts}
    with (root / PANEL / "expanded_island_panel_species_key_facets.csv").open() as handle:
        for row in csv.DictReader(handle):
            key = row["island_key"], row["window"]
            species = row["facet_key"]
            if key not in facet_counts or species in facet_counts[key]:
                raise ValueError("Unexpected or duplicated published species facet")
            facet_counts[key][species] = int(row["record_count"])
    if facet_counts != counts:
        raise ValueError("Published species facets disagree with the complete source manifest")
    with (root / TAXONOMY / "original239_to_colxr_bridge.csv").open() as handle:
        bridge = list(csv.DictReader(handle))
    receipt = json.loads((root / TAXONOMY / "taxonomy_bridge_manifest.json").read_bytes())
    if receipt["output_sha256"]["original239_to_colxr_bridge.csv"] != PINS[TAXONOMY + "original239_to_colxr_bridge.csv"]:
        raise ValueError("Bridge output does not match its official source receipt")
    eligible = [r for r in bridge if r["conservative_join_eligible"] == "True"]
    taxa = {r["current_accepted_key"]: r["original_name"] for r in eligible}
    if len(bridge) != 239 or len(eligible) != 230 or len(taxa) != 230:
        raise ValueError("Expected 230 unique conservative original-name/COL concepts")
    return roster, counts, taxa, selection


def read_roster_samples(source_cache, roster, source_pin):
    """Read existing validated source files without extraction or network writes."""
    import geopandas as gpd
    from shapely import orient_polygons, to_wkt
    from measure_island_port_geography_20261029 import (
        SHA, derived_key, dissolve_original, shoreline_xyz,
    )
    archive = source_cache / "gshhg-shp-2.3.7.zip"
    if source_pin["zip_sha256"] != SHA or sha(archive.read_bytes()) != SHA:
        raise ValueError("Original GSHHG archive differs from the frozen source")
    with zipfile.ZipFile(archive) as handle:
        members = [n for n in handle.namelist() if n.endswith("GSHHS_h_L1.shp")]
        if len(members) != 1:
            raise ValueError("Expected one original high-resolution L1 source")
        stem = members[0][:-4]
        for suffix in (".shp", ".shx", ".dbf", ".prj"):
            if sha(handle.read(stem + suffix)) != sha((source_cache / ("GSHHS_h_L1" + suffix)).read_bytes()):
                raise ValueError("Extracted GSHHG source does not match pinned archive bytes")
    source = gpd.read_file(source_cache / "GSHHS_h_L1.shp")
    land = dissolve_original(source)
    if len(source) != 144749 or len(land) != 144744:
        raise ValueError("Original source landmass census changed")
    needed = {r["island_key"]: r for r in roster}
    samples, geometry_receipts = {}, []
    for _, row in land.iterrows():
        key = derived_key(row.parent_source_ID)
        if key not in needed:
            continue
        geometry = orient_polygons(row.geometry, exterior_cw=False)
        if not geometry.is_valid or geometry.is_empty:
            raise ValueError("Invalid exact frozen geometry")
        geometry_sha = sha(to_wkt(geometry, rounding_precision=-1, trim=True).encode())
        if geometry_sha != needed[key]["geometry_sha256"]:
            raise ValueError("Exact frozen island geometry hash changed")
        samples[key], receipt = shoreline_xyz(geometry, 1.)
        geometry_receipts.append({"island_key": key, "geometry_sha256": geometry_sha, **receipt})
    if set(samples) != set(needed):
        raise ValueError("Original source is missing selected islands")
    return samples, sorted(geometry_receipts, key=lambda r: r["island_key"])


def measure_pair_distances(samples, step_km=1.):
    from measure_island_port_geography_20261029 import coast_gap_interval
    if not 0 < step_km <= 1:
        raise ValueError("Source-coast sampling bound must be at most 1 km")
    trees = {key: cKDTree(points) for key, points in samples.items()}
    result = {}
    for left, right in itertools.combinations(sorted(samples), 2):
        small, large = (left, right) if len(samples[left]) <= len(samples[right]) else (right, left)
        lower, upper = coast_gap_interval(samples[small], trees[large], step_km)
        if not 0 <= lower <= upper or upper - lower > step_km + 1e-8:
            raise ValueError("Invalid coast-to-coast sampling interval")
        result[left, right] = float(lower), float(upper)
    return result


def longitude_block(lon):
    lon = float(lon)
    if not math.isfinite(lon) or not -180 <= lon <= 180:
        raise ValueError("Invalid representative longitude")
    return f"longitude_{min(3, int((lon + 180) // 90))}"


def build_risk_set(roster, counts, eligible_taxa, distances, source_min_records=1, late_min_records=1):
    if source_min_records < 1 or late_min_records < 1:
        raise ValueError("Positive focal record thresholds required")
    keyed = {r["island_key"]: r for r in roster}
    if len(keyed) != len(roster):
        raise ValueError("Duplicate roster identities")
    expected = {(i, w) for i in keyed for w in WINDOWS}
    if set(counts) != expected:
        raise ValueError("Complete explicit zero/nonzero windows required for every island")
    targets = sorted(i for i in keyed if len(counts[i, WINDOWS[0]]) >= 3 and len(counts[i, WINDOWS[1]]) >= 4)
    blocks = {i: longitude_block(keyed[i]["representative_lon"]) for i in keyed}
    sources = {s: sorted(i for i in keyed if counts[i, WINDOWS[0]].get(s, 0) >= source_min_records)
               for s in sorted(eligible_taxa)}
    sources = {s: islands for s, islands in sources.items() if islands}
    row_labels = [(s, b) for s in sources for b in sorted({blocks[i] for i in targets})]
    row_index = {label: i for i, label in enumerate(row_labels)}
    shape = len(row_labels), len(targets)
    mask, state = np.zeros(shape, dtype=bool), np.zeros(shape, dtype=bool)
    values, lower_values = np.full(shape, np.nan), np.full(shape, np.nan)
    cells = []
    for ci, target in enumerate(targets):
        early, late = counts[target, WINDOWS[0]], counts[target, WINDOWS[1]]
        for species, islands in sources.items():
            if early.get(species, 0):
                continue
            if target in islands:
                raise ValueError("Baseline-absent target unexpectedly appears in source set")
            candidates = [(distances[tuple(sorted((target, source)))], source) for source in islands]
            upper, sampled_source = min((interval[1], source) for interval, source in candidates)
            lower = min(interval[0] for interval, _ in candidates)
            possible = sorted(source for interval, source in candidates if interval[0] <= upper)
            if not 0 <= lower <= upper or upper - lower > 1 + 1e-8:
                raise ValueError("Invalid minimum-over-sources sampling interval")
            ri = row_index[species, blocks[target]]
            mask[ri, ci] = True
            state[ri, ci] = late.get(species, 0) >= late_min_records
            values[ri, ci], lower_values[ri, ci] = math.log1p(upper), math.log1p(lower)
            cells.append({"row_index": ri, "column_index": ci, "species_key": species,
                          "original_name": eligible_taxa[species], "target_island": target,
                          "longitude_block": blocks[target], "target_early_total_species": len(early),
                          "target_late_total_species": len(late), "target_early_focal_records": 0,
                          "target_late_focal_records": late.get(species, 0), "late_positive": int(state[ri, ci]),
                          "source_min_records": source_min_records, "late_min_records": late_min_records,
                          "source_island_count": len(islands), "source_islands": "|".join(islands),
                          "sampled_nearest_source_island": sampled_source,
                          "sampled_nearest_source_baseline_records": counts[sampled_source, WINDOWS[0]][species],
                          "possible_nearest_source_islands": "|".join(possible),
                          "source_identity_ambiguous_at_sampling_resolution": int(len(possible) > 1),
                          "nearest_source_lower_km": lower, "nearest_source_upper_km": upper,
                          "log1p_source_lower_km": math.log1p(lower), "log1p_source_upper_km": math.log1p(upper)})
    validate_state(state, mask)
    if int(mask.sum()) != len(cells) or int(state.sum()) != sum(r["late_positive"] for r in cells):
        raise ValueError("Risk-mask/cell accounting discrepancy")
    return {"mask": mask, "state": state, "values": values, "lower_values": lower_values,
            "row_labels": row_labels, "column_labels": targets, "cells": cells,
            "summary": {"candidate_cells": len(cells), "positive_cells": int(state.sum()),
                        "primary_eligible_original_taxa": len(eligible_taxa),
                        "taxa_with_qualifying_baseline_source": len(sources),
                        "targets_with_fixed_recording_support": len(targets),
                        "targets_with_candidate_cells": int(np.count_nonzero(mask.sum(0))),
                        "targets_with_positive_cells": int(np.count_nonzero(state.sum(0))),
                        "species_with_candidate_cells": len({r["species_key"] for r in cells}),
                        "species_with_positive_cells": len({r["species_key"] for r in cells if r["late_positive"]}),
                        "species_block_rows": len(row_labels),
                        "positive_species_block_rows": int(np.count_nonzero(state.sum(1))),
                        "zero_margin_species_block_rows": int(np.count_nonzero(state.sum(1) == 0)),
                        "zero_margin_target_columns": int(np.count_nonzero(state.sum(0) == 0)),
                        "structural_zero_matrix_cells": int(mask.size - mask.sum())}}


def validate_state(state, mask, row_margins=None, column_margins=None):
    if state.ndim != 2 or state.shape != mask.shape or not np.isin(state, (0, 1)).all() or not np.isin(mask, (0, 1)).all():
        raise ValueError("Aligned binary state and risk mask required")
    if np.any(state & ~mask):
        raise ValueError("Positive state outside structural risk mask")
    if row_margins is not None and not np.array_equal(state.sum(1), row_margins):
        raise ValueError("Species/block row margins changed")
    if column_margins is not None and not np.array_equal(state.sum(0), column_margins):
        raise ValueError("Target-island column margins changed")


def cycle_length_probabilities(max_length):
    """Depends only on fixed positive margins/E, never current valid matchings."""
    if max_length < 2:
        return {}
    weights = {2: .8}
    if max_length >= 3:
        weights[3] = .15
    if max_length >= 4:
        weights.update({k: .05 / (max_length - 3) for k in range(4, max_length + 1)})
    total = sum(weights.values())
    return {k: v / total for k, v in weights.items()}


def apply_cycle(state, mask, ordered_edges):
    """One proposal; invalid selections are rejected without resampling."""
    rows = [edge[0] for edge in ordered_edges]
    columns = [edge[1] for edge in ordered_edges]
    k = len(rows)
    if k < 2 or len(set(rows)) != k or len(set(columns)) != k:
        return "repeated_row_or_column"
    if not all(state[r, c] for r, c in ordered_edges):
        raise ValueError("Proposal did not draw distinct existing positive edges")
    proposed = list(zip(rows, columns[1:] + columns[:1]))
    if not all(mask[r, c] for r, c in proposed):
        return "structural_zero"
    if any(state[r, c] for r, c in proposed):
        return "occupied_proposed_cell"
    for r, c in ordered_edges:
        state[r, c] = False
    for r, c in proposed:
        state[r, c] = True
    return "accepted"


def fiber_has_move(state, mask):
    """An alternating directed cycle is equivalent to another fiber state."""
    validate_state(state, mask)
    rows, columns = state.shape
    if not rows or not columns:
        return False
    addition_rows, addition_columns = np.where(mask & ~state)
    deletion_rows, deletion_columns = np.where(state)
    origin = np.r_[addition_rows, rows + deletion_columns]
    destination = np.r_[rows + addition_columns, deletion_rows]
    graph = csr_matrix((np.ones(len(origin)), (origin, destination)), shape=(rows + columns, rows + columns))
    _, labels = connected_components(graph, directed=True, connection="strong")
    # A cell has one orientation only, so no directed 2-cycle can exist.
    return bool(np.any(np.bincount(labels) >= 4))


def sample_chain(initial, mask, values, seed, burn_in=BURN_IN, samples=SAMPLES, thin=THIN):
    validate_state(initial, mask)
    if initial.sum() == 0 or burn_in < 0 or samples < 1 or thin < 1 or not np.isfinite(values[mask]).all():
        raise ValueError("Positive observed edges, finite allowed scores and valid run lengths required")
    state = initial.copy()
    row_margins, column_margins = state.sum(1), state.sum(0)
    edges = [tuple(map(int, edge)) for edge in zip(*np.where(state))]
    maximum = min(np.count_nonzero(row_margins), np.count_nonzero(column_margins), len(edges))
    probabilities = cycle_length_probabilities(maximum)
    lengths = list(probabilities)
    cumulative = np.cumsum(list(probabilities.values())).tolist()
    rng = random.Random(int(seed))
    outcomes, proposed_lengths, accepted_lengths = Counter(), Counter(), Counter()
    draws, states_seen = [], set()
    attempts = burn_in + samples * thin
    for attempt in range(1, attempts + 1):
        if rng.random() < LAZY_HOLD:
            outcome = "lazy_hold"
        elif not lengths:
            outcome = "no_possible_cycle_length"
        else:
            draw = rng.random()
            length = next((k for k, upper in zip(lengths, cumulative) if draw < upper), lengths[-1])
            proposed_lengths[length] += 1
            indices = rng.sample(range(len(edges)), length)
            chosen = [edges[index] for index in indices]
            outcome = apply_cycle(state, mask, chosen)
            if outcome == "accepted":
                new_columns = [edge[1] for edge in chosen[1:] + chosen[:1]]
                for index, edge, column in zip(indices, chosen, new_columns):
                    edges[index] = edge[0], column
                accepted_lengths[length] += 1
        outcomes[outcome] += 1
        if attempt > burn_in and (attempt - burn_in) % thin == 0:
            validate_state(state, mask, row_margins, column_margins)
            if len(set(edges)) != len(edges) or any(not state[r, c] for r, c in edges):
                raise ValueError("Mutable positive-edge list diverged from binary matrix")
            draws.append(float(np.mean(values[state])))
            states_seen.add(sha(state.tobytes()))
    if len(draws) != samples or sum(outcomes.values()) != attempts:
        raise ValueError("Attempted-step sample accounting failed")
    return np.asarray(draws), {"seed": int(seed), "attempted_steps": attempts, "burn_in_attempts": burn_in,
                              "sampling_every_attempts": thin, "stored_samples": samples,
                              "stored_margin_checks": samples, "accepted_moves": outcomes["accepted"],
                              "acceptance_fraction_all_attempts": outcomes["accepted"] / attempts,
                              "outcomes": dict(outcomes), "proposed_cycle_lengths": dict(proposed_lengths),
                              "accepted_cycle_lengths": dict(accepted_lengths),
                              "cycle_length_probabilities": probabilities, "distinct_stored_states": len(states_seen),
                              "sample_mean_T": float(np.mean(draws)), "sample_sd_T": float(np.std(draws, ddof=1))}


def split_chains(chains):
    chains = np.asarray(chains, dtype=float)
    if chains.ndim != 2 or chains.shape[0] < 2 or chains.shape[1] < 8 or not np.isfinite(chains).all():
        raise ValueError("At least two finite, equal-length chains with eight samples are required")
    half = chains.shape[1] // 2
    return np.concatenate((chains[:, :half], chains[:, -half:]), axis=0)


def rhat(chains):
    n = chains.shape[1]
    within = float(np.mean(np.var(chains, axis=1, ddof=1)))
    between = float(n * np.var(np.mean(chains, axis=1), ddof=1))
    if within <= 0:
        return None
    variance = (n - 1) / n * within + between / n
    return float(math.sqrt(variance / within))


def pooled_ess(chains):
    """Multi-chain Geyer initial-positive, initial-monotone paired sequence."""
    chains = np.asarray(chains, dtype=float)
    m, n = chains.shape
    within = float(np.mean(np.var(chains, axis=1, ddof=1)))
    variance = (n - 1) / n * within + float(np.var(chains.mean(axis=1), ddof=1))
    if within <= 0 or variance <= 0:
        return None
    centered = chains - chains.mean(axis=1, keepdims=True)
    fft_size = 1 << (2 * n - 1).bit_length()
    transformed = np.fft.rfft(centered, n=fft_size, axis=1)
    covariance = np.fft.irfft(transformed * transformed.conj(), n=fft_size, axis=1)[:, :n] / n
    rho = 1 - (within - covariance.mean(axis=0)) / variance
    rho[0] = 1.
    pair_sums = []
    for lag in range(0, n - 1, 2):
        pair = float(rho[lag] + rho[lag + 1])
        if pair < 0:
            break
        pair_sums.append(min(pair, pair_sums[-1]) if pair_sums else pair)
    tau = max(1., -1. + 2 * sum(pair_sums))
    return float(min(m * n, m * n / tau))


def assess_samples(chains, observed):
    chains = np.asarray(chains, dtype=float)
    split = split_chains(chains)
    count = int(chains.size)
    indicators = (chains <= observed + TIE_TOLERANCE).astype(float)
    tail_count = int(indicators.sum())
    result = {"stored_samples_total": count, "observed_T": float(observed),
              "null_mean_T": float(chains.mean()), "null_sd_T": float(chains.std(ddof=1)),
              "observed_to_null_geometric_mean_one_plus_gap_ratio": float(math.exp(observed - chains.mean())),
              "lower_tail_count_including_ties": tail_count, "other_tail_count": count - tail_count,
              "finite_mcmc_tail_fraction_plus_one": (tail_count + 1) / (count + 1),
              "qualified_tail_probability": None, "qualified_tail": False,
              "split_rhat": None, "pooled_ess": None, "tail_indicator_ess": None,
              "tail_probability_mcse": None, "gate_failures": []}
    if float(np.ptp(chains)) <= 1e-12:
        result["gate_failures"] = ["degenerate_statistic"]
        result["status"] = "DEGENERATE_STATISTIC_NO_QUALIFIED_TAIL"
        return result
    ranks = rankdata(split.ravel(), method="average").reshape(split.shape)
    rank_normal = norm.ppf((ranks - .375) / (split.size + .25))
    folded = np.abs(split - np.median(split))
    folded_ranks = rankdata(folded.ravel(), method="average").reshape(split.shape)
    folded_normal = norm.ppf((folded_ranks - .375) / (split.size + .25))
    raw_rhat, rank_rhat, folded_rhat = rhat(split), rhat(rank_normal), rhat(folded_normal)
    rhats = [x for x in (raw_rhat, rank_rhat, folded_rhat) if x is not None]
    result.update(raw_split_rhat=raw_rhat, rank_normalized_split_rhat=rank_rhat,
                  folded_rank_normalized_split_rhat=folded_rhat,
                  split_rhat=max(rhats) if rhats else None)
    raw_ess, rank_ess = pooled_ess(split), pooled_ess(rank_normal)
    result.update(raw_T_ess=raw_ess, rank_normalized_T_ess=rank_ess,
                  pooled_ess=min(raw_ess, rank_ess) if raw_ess is not None and rank_ess is not None else None)
    indicator_ess = pooled_ess(split_chains(indicators)) if 0 < tail_count < count else None
    result["tail_indicator_ess"] = indicator_ess
    if indicator_ess is not None:
        probability = tail_count / count
        result["tail_probability_mcse"] = math.sqrt(probability * (1 - probability) / indicator_ess)
    if result["split_rhat"] is None or result["split_rhat"] >= 1.05:
        result["gate_failures"].append("split_rhat_not_below_1_05")
    if result["pooled_ess"] is None or result["pooled_ess"] < 1000:
        result["gate_failures"].append("statistic_ess_below_1000")
    result["statistic_mixing_qualified"] = not result["gate_failures"]
    if min(tail_count, count - tail_count) < 20:
        result["gate_failures"].append("fewer_than_20_samples_in_one_tail")
    if indicator_ess is None or indicator_ess < 1000:
        result["gate_failures"].append("tail_indicator_ess_below_1000_or_undefined")
    if not result["gate_failures"]:
        result.update(qualified_tail=True, qualified_tail_probability=result["finite_mcmc_tail_fraction_plus_one"],
                      status="APPROXIMATE_CONDITIONAL_MCMC_TAIL_QUALIFIED")
    else:
        result["status"] = "MCMC_DIAGNOSTIC_OR_TAIL_PRECISION_GATE_FAILED"
    return result


def run_null(state, mask, values, label=""):
    validate_state(state, mask)
    if not state.any():
        return {"summary": {"status": "NO_POSITIVE_OUTCOMES", "qualified_tail": False,
                            "qualified_tail_probability": None, "observed_T": None}, "samples": [], "chains": []}
    observed = float(np.mean(values[state]))
    if not fiber_has_move(state, mask):
        return {"summary": {"status": "DETERMINISTIC_MARGIN_FIBER", "qualified_tail": False,
                            "qualified_tail_probability": None, "observed_T": observed,
                            "null_mean_T": observed, "gate_failures": ["no_alternating_cycle"]},
                "samples": [], "chains": []}
    draws, chains = [], []
    for number, seed in enumerate(SEEDS, 1):
        samples, report = sample_chain(state, mask, values, seed)
        draws.append(samples)
        chains.append(report)
        print(json.dumps({"stage": "recorded_source_chain_complete", "analysis": label,
                          "chain": number, "seed": seed, "accepted": report["accepted_moves"],
                          "attempted": report["attempted_steps"]}), flush=True)
    summary = assess_samples(np.array(draws), observed)
    return {"summary": summary, "samples": draws, "chains": chains}


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(source_cache, output, root=ROOT):
    started = utc()
    roster, counts, taxa, selection = load_inputs(root)
    output.mkdir(parents=True, exist_ok=True)
    pre_run = {"started_utc": started, "plan_sha256": PINS[PLAN], "input_sha256": PINS,
               "script_sha256": sha(Path(__file__).read_bytes()), "seeds": SEEDS,
               "burn_in_attempts_per_chain": BURN_IN, "samples_per_chain": SAMPLES,
               "attempts_between_samples": THIN, "lazy_hold_probability": LAZY_HOLD,
               "scenario_order": ["primary", "source_records_ge3", "late_records_ge2"]}
    (output / "recorded_source_run_lock.json").write_text(json.dumps(pre_run, indent=2, sort_keys=True) + "\n")
    samples, geometry_receipts = read_roster_samples(source_cache, roster, selection["selection"]["gshhg_source"])
    distances = measure_pair_distances(samples)
    if len(distances) != 96 * 95 // 2:
        raise ValueError("All fixed 96-island unordered coast pairs must be measured")
    distance_rows = [{"island_left": left, "island_right": right, "coast_gap_lower_km": lower,
                      "coast_gap_upper_km": upper} for (left, right), (lower, upper) in sorted(distances.items())]
    write_csv(output / "recorded_source_coast_pairs.csv", distance_rows, list(distance_rows[0]))
    support = [{"island_key": r["island_key"], "longitude_block": longitude_block(r["representative_lon"]),
                "early_recorded_species": len(counts[r["island_key"], WINDOWS[0]]),
                "late_recorded_species": len(counts[r["island_key"], WINDOWS[1]]),
                "fixed_target_recording_support": int(len(counts[r["island_key"], WINDOWS[0]]) >= 3
                    and len(counts[r["island_key"], WINDOWS[1]]) >= 4)} for r in roster]
    write_csv(output / "recorded_source_target_support.csv", support, list(support[0]))
    print(json.dumps({"stage": "exact_coast_pairs_complete", "pairs": len(distances),
                      "fixed_target_support": sum(r["fixed_target_recording_support"] for r in support)}), flush=True)
    scenarios, all_cells, all_margins, all_draws = {}, [], [], []
    for label, source_min, late_min in (("primary", 1, 1), ("source_records_ge3", 3, 1), ("late_records_ge2", 1, 2)):
        risk = build_risk_set(roster, counts, taxa, distances, source_min, late_min)
        state, mask = risk["state"], risk["mask"]
        outcome = run_null(state, mask, risk["values"], label)
        observed_lower = float(np.mean(risk["lower_values"][state])) if state.any() else None
        scenarios[label] = {"source_min_records": source_min, "late_min_records": late_min,
                            "support": risk["summary"], "diagnostic": outcome["summary"],
                            "observed_lower_bound_T_sampling_only": observed_lower,
                            "chains": outcome["chains"], "risk_mask_sha256": sha(mask.tobytes()),
                            "observed_binary_state_sha256": sha(state.tobytes()),
                            "row_labels": risk["row_labels"], "column_labels": risk["column_labels"],
                            "risk_mask_shape": list(mask.shape)}
        all_cells.extend({"analysis": label, **cell} for cell in risk["cells"])
        for index, (species, block) in enumerate(risk["row_labels"]):
            all_margins.append({"analysis": label, "dimension": "species_longitude_block", "index": index,
                                "unit_key": species + ":" + block, "positive_margin": int(state[index].sum()),
                                "allowed_cells": int(mask[index].sum())})
        for index, island in enumerate(risk["column_labels"]):
            all_margins.append({"analysis": label, "dimension": "target_island", "index": index,
                                "unit_key": island, "positive_margin": int(state[:, index].sum()),
                                "allowed_cells": int(mask[:, index].sum())})
        for chain, (seed, draws) in enumerate(zip(SEEDS, outcome["samples"]), 1):
            all_draws.extend({"analysis": label, "chain": chain, "seed": seed, "sample": number,
                              "attempted_step": BURN_IN + number * THIN, "mean_log1p_source_gap_km": float(value)}
                             for number, value in enumerate(draws, 1))
        print(json.dumps({"stage": "recorded_source_analysis_complete", "analysis": label,
                          "support": risk["summary"], "diagnostic": outcome["summary"]}), flush=True)
    if all_cells:
        write_csv(output / "recorded_source_risk_cells.csv", all_cells, list(all_cells[0]))
    write_csv(output / "recorded_source_fixed_margins.csv", all_margins,
              ["analysis", "dimension", "index", "unit_key", "positive_margin", "allowed_cells"])
    write_csv(output / "recorded_source_chain_samples.csv", all_draws,
              ["analysis", "chain", "seed", "sample", "attempted_step", "mean_log1p_source_gap_km"])
    report = {"schema": "chocho_postpilot_recorded_source_fixed_margin_diagnostic_20261029_v1",
              "status": "COMPLETE_POSTPILOT_EXPLORATORY_DIAGNOSTIC", "started_utc": started, "finished_utc": utc(),
              "input_sha256": PINS, "script_sha256": sha(Path(__file__).read_bytes()),
              "selection_content_sha256": selection["selection_sha256"],
              "gshhg_source": selection["selection"]["gshhg_source"], "geometry_receipts": geometry_receipts,
              "distance_method": "Every unordered fixed96 pair; exact unsimplified source exteriors sampled at <=1km arc spacing; artificial dateline cut edges excluded; spherical cKDTree minimum. Interval [max(0, sampled_upper-1km), sampled_upper] covers coast sampling only, not positional/spherical model error. Focal island excluded. Statistic uses sampled upper; observed lower-bound statistic retained separately.",
              "risk_set": "230 conservative original concepts; require qualifying early record on another selected island, no early focal target record, >=3 total early target taxa and >=4 total late target taxa. Target eligibility is fixed, independent of randomized focal labels. All 192 queries complete and species facets reconciled. Risk mask includes only the target's longitude-block row; all other cells are structural zeros.",
              "null_method": {"distribution": "Uniform binary matrices with fixed species-by-longitude-block row margins and fixed target-island column margins, restricted to the fixed risk mask.",
                              "proposal": "With probability0.5 hold. Otherwise draw k from a state-independent pmf: weights .80 for2, .15 for3, .05 equally across4..maximum; omit unavailable lengths and normalize. maximum=min(positive row count, positive column count, total positive edges), all fixed by margins. Draw ordered k edges uniformly without replacement from ALL E current positive edges; do NOT resample until valid. Reject repeated rows/columns, disallowed proposed cells or occupied proposed cells. Otherwise cyclically rotate columns and flip. All rejections count as attempted steps.",
                              "reversibility": "A valid ordered cycle and its reversed ordered cycle have equal k weight and identical 1/(E)_k draw probability. Invalid proposals add self-loop mass. Hence transition probabilities are symmetric and the uniform fiber law is stationary. Differences between any two feasible binary matrices decompose into alternating simple cycles; every such cycle length has positive proposal probability, giving connectivity. The explicit lazy hold supplies aperiodicity including a two-state fiber.",
                              "diagnostics": "Four fixed seeds;100000 attempted burn-in steps then4000 samples every50 attempts. Verify matrix/edge list, fixed mask and both margins at every stored sample. Split Rhat is max(raw, rank-normalized, folded rank-normalized); pooled ESS is min(raw and rank-normalized multichain initial-positive/monotone estimates). Tail indicator ESS also uses split chains. Gate nondegenerate statistic, Rhat<1.05, statistic ESS>=1000, both tail counts>=20 and indicator ESS>=1000. No run-length adaptation.",
                              "tail": "Directional closer-source tail includes T_null<=T_observed, with1e-12 absolute numerical tie tolerance. (count+1)/(draws+1) is an approximate MCMC tail estimate, not exact. Finite counts remain visible when qualification fails. MCSE uses the indicator ESS when defined.",
                              "refinements_before_sampling": "State-independent0.5lazy hold for aperiodicity; rank/folded Rhat checks; inclusive numerical tie handling; tail-indicator ESS and minimum20samples in each tail for precision. Frozen scientific plan unchanged.",
                              "seeds": SEEDS, "burn_in_attempts_per_chain": BURN_IN, "samples_per_chain": SAMPLES,
                              "attempts_between_samples": THIN, "lazy_hold_probability": LAZY_HOLD,
                              "total_attempts_per_nondegenerate_analysis": len(SEEDS) * (BURN_IN + SAMPLES * THIN)},
              "interpretation": ["All96outcomes were already known: post-pilot exploration, not temporal holdout validation or replacement for the frozen richness/port failure.",
                                 "A baseline nondetection is not absence; a later record is not first arrival, dispersal, establishment or abundance.",
                                 "At least3other recorded species is recording support, not visits, surveys or a detection model. Late recording support conditions the target population.",
                                 "Only the sampled96islands can be recorded source candidates; mainland and most island sources are absent. Occurrences do not verify breeding populations or direction of movement.",
                                 "Margins constrain target recording gains and species gains within coarse longitude blocks, not all habitat, climate, observer, taxonomy or source-sampling biases. Static port, mainland-distance and area main effects cannot be estimated under fixed island margins.",
                                 "The ratio is exp(observed mean log(1+gap) minus null mean of that statistic), a descriptive geometric-distance contrast, not odds, a causal effect or confidence interval.",
                                 "Source>=3records and late>=2records are separate declared sensitivities, with their own masks/margins; none is selected to replace the primary."],
              "analyses": scenarios,
              "output_sha256": {p.name: sha(p.read_bytes()) for p in sorted(output.glob("*.csv"))}}
    (output / "recorded_source_diagnostic.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    lines = ["# Prior recorded-source proximity: post-pilot diagnostic", "",
             "The fixed 96-island early/late outcomes were already observed. This is an exploratory conditional recording diagnostic, not a temporal holdout or a colonization test.", "",
             "The frozen plan is `provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_PLAN_20261029.md`. The original richness/port plan and failed gates are unchanged.", "",
             "| Analysis | Candidate cells | Later-record positives | Observed T | Null mean T | Geometric-distance ratio | Qualified approximate tail |", "|---|---:|---:|---:|---:|---:|---:|"]
    for label, analysis in scenarios.items():
        d, s = analysis["diagnostic"], analysis["support"]
        fmt = lambda key: "not qualified" if d.get(key) is None else f"{d[key]:.6g}"
        lines.append(f"| {label} | {s['candidate_cells']} | {s['positive_cells']} | {fmt('observed_T')} | {fmt('null_mean_T')} | {fmt('observed_to_null_geometric_mean_one_plus_gap_ratio')} | {fmt('qualified_tail_probability')} |")
    lines.extend(["", "T is mean log(1 + nearest sampled-source coast gap in km) among later-record positive cells. Source candidates had earlier focal records on other members of the sampled96. Upper coast-gap estimates supply T; the JSON retains sampling intervals and lower-bound observed T.", "",
                  "The null fixes species-by-longitude-block and target-island positive margins. Four fixed chains each use100,000 attempted burn-in steps and4,000 samples every50attempts. Rejections and0.5lazy holds count as attempts. Chains are never extended after observing a result. Full diagnostics, qualification failures and finite tail counts are in the JSON; samples, margins and the explicit allowed risk cells are supplied for independent review.", "",
                  "Three other recorded taxa do not establish visits or adequate detectability. Later-record outcomes are not arrival or establishment. The source pool omits mainland and many island populations. Fixed island margins prevent estimating standalone port, area or mainland-distance effects. The distance ratio is descriptive, not a causal effect or an odds ratio.", "",
                  "Replay from repository root:", "", "```bash", "/tmp/chocho-geo-venv/bin/python scripts/analyze_island_recorded_sources_20261029.py --source-cache /tmp/chocho-island-sources", "```", ""])
    (output / "README.md").write_text("\n".join(lines))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/exploratory/island_recorded_sources_20261029")
    args = parser.parse_args()
    run(args.source_cache, args.output_dir)
