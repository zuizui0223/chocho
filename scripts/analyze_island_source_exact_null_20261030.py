#!/usr/bin/env python3
"""Exact counting and independent sampling of the frozen conditional table fiber.

The ecological estimand, risk masks, margins and distance statistic are unchanged.
Integer completion counts give a rank/table bijection; only score arithmetic and
binomial quantiles use floating point. This remains a post-pilot recording study.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import itertools
import json
import math
import operator
from pathlib import Path
import platform
import random
import sys
import time

import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
PLAN = "provenance/archive/host_hindcast/ISLAND_EXACT_NULL_HOST_ROSTER_PLAN_20261030.md"
PLAN_SHA256 = "ce04b052a37ad4f07d2c81e87916dbf3addac3083a7117ebecea47545a4f99e9"
OLD = "data/exploratory/island_recorded_sources_20261029/"
PINS = {
    PLAN: PLAN_SHA256,
    OLD + "recorded_source_diagnostic.json": "372fec4a1035132c17eec67c5392c5eab0727d3c596edb2d4f2b45b88b6ecebb",
    OLD + "recorded_source_risk_cells.csv": "a8c7905f1ef6e187b2647d5ed39f1a5a67f901b94a5aab3a86b6358969f72fb4",
    OLD + "recorded_source_fixed_margins.csv": "e472a2dcb0977a02dda1f313bc66ab32d0c096a2b74d505c764791b1b15ecbd0",
    OLD + "recorded_source_coast_pairs.csv": "d83803dd3f976069533da40ff1c20f6602386b699de3dd744f37096e63f03545",
    OLD + "recorded_source_run_lock.json": "a0124d968516c1d9a5b5cffb1ba6d6dd7f5826cedccd81c188c9d2e295735619",
}
SCENARIOS = (("primary", 2026103001, 1, 1),
             ("source_records_ge3", 2026103002, 3, 1),
             ("late_records_ge2", 2026103003, 1, 2))
SAMPLES = 100_000
MAX_STATES = 100_000
MAX_TRANSITIONS = 2_000_000
TIE_TOLERANCE = 1e-12


def sha(body):
    return hashlib.sha256(body).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def verify_pinned_files(root, pins=PINS):
    for relative, expected in pins.items():
        if sha((Path(root) / relative).read_bytes()) != expected:
            raise ValueError(f"Pinned input changed: {relative}")


def integer_margins(values, length):
    try:
        result = tuple(operator.index(x) for x in values)
    except TypeError as error:
        raise ValueError("Margins must be integers") from error
    if len(result) != length or any(x < 0 for x in result):
        raise ValueError("Nonnegative margins of the required length are required")
    return result


class DPResourceLimitError(RuntimeError):
    """An incomplete count must never be used as a probability denominator."""


class ExactBlockFiber:
    """Labeled binary tables with structural zeros and both fixed margins.

    C(i, c) sums exact completion counts over the fixed row's allowed subsets.
    Column coordinates retain their labels. Zero-margin removal changes neither
    table multiplicity nor the reverse reconstruction. All row identities remain
    distinct even when their allowed patterns and margins coincide.
    """

    def __init__(self, mask, row_margins, column_margins,
                 max_states=MAX_STATES, max_transitions=MAX_TRANSITIONS):
        raw_mask = np.asarray(mask)
        if raw_mask.ndim != 2 or not np.isin(raw_mask, (0, 1)).all():
            raise ValueError("A two-dimensional binary mask is required")
        self.mask = raw_mask.astype(bool, copy=True)
        self.row_margins = integer_margins(row_margins, self.mask.shape[0])
        self.column_margins = integer_margins(column_margins, self.mask.shape[1])
        if sum(self.row_margins) != sum(self.column_margins):
            raise ValueError("Row and column margin totals differ")
        if max_states < 1 or max_transitions < 1:
            raise ValueError("Positive computation caps are required")
        self.max_states, self.max_transitions = max_states, max_transitions
        self.columns = tuple(j for j, margin in enumerate(self.column_margins) if margin)
        active_rows = [i for i, margin in enumerate(self.row_margins) if margin]
        allowed = {i: tuple(j for j, original in enumerate(self.columns) if self.mask[i, original])
                   for i in active_rows}
        self.rows = tuple(sorted(active_rows, key=lambda i: (
            math.comb(len(allowed[i]), self.row_margins[i])
            if self.row_margins[i] <= len(allowed[i]) else 0, len(allowed[i]), i)))
        self.patterns = tuple(tuple(itertools.combinations(allowed[i], self.row_margins[i]))
                              for i in self.rows)
        self.initial = tuple(self.column_margins[j] for j in self.columns)
        n, width = len(self.rows), len(self.columns)
        suffix_capacity = [[0] * width for _ in range(n + 1)]
        self.suffix_edges = [0] * (n + 1)
        for i in range(n - 1, -1, -1):
            suffix_capacity[i] = [suffix_capacity[i + 1][j] + (j in allowed[self.rows[i]])
                                  for j in range(width)]
            self.suffix_edges[i] = self.suffix_edges[i + 1] + self.row_margins[self.rows[i]]
        self.suffix_capacity = tuple(tuple(c) for c in suffix_capacity)
        self._counts, self._branch_cache = {}, {}
        self.memo_states = self.attempted_transitions = 0
        self.count = self._count(0, self.initial)

    def _count(self, i, residual):
        key = i, residual
        cached = self._counts.get(key)
        if cached is not None:
            return cached
        self.memo_states += 1
        if self.memo_states > self.max_states:
            raise DPResourceLimitError(f"Memo states exceed fixed cap {self.max_states}")
        if sum(residual) != self.suffix_edges[i] or any(
                c < 0 or c > cap for c, cap in zip(residual, self.suffix_capacity[i])):
            total = 0
        elif i == len(self.rows):
            total = int(not any(residual))
        else:
            total = 0
            for pattern in self.patterns[i]:
                self.attempted_transitions += 1
                if self.attempted_transitions > self.max_transitions:
                    raise DPResourceLimitError(f"Pattern transitions exceed fixed cap {self.max_transitions}")
                if any(residual[j] == 0 for j in pattern):
                    continue
                child = list(residual)
                for j in pattern:
                    child[j] -= 1
                total += self._count(i + 1, tuple(child))
        self._counts[key] = total
        return total

    def _branches(self, i, residual):
        """Positive-completion intervals; no branch uses floating weights."""
        key = i, residual
        if key not in self._branch_cache:
            cumulative, branches = 0, []
            for pattern in self.patterns[i]:
                if any(residual[j] == 0 for j in pattern):
                    continue
                child = list(residual)
                for j in pattern:
                    child[j] -= 1
                child = tuple(child)
                count = self._counts[i + 1, child]
                if count:
                    cumulative += count
                    branches.append((cumulative, pattern, child, count))
            if cumulative != self._counts[key]:
                raise AssertionError("Child completion counts do not sum to parent")
            self._branch_cache[key] = tuple(branches), tuple(b[0] for b in branches)
        return self._branch_cache[key]

    def unrank(self, rank):
        try:
            rank = operator.index(rank)
        except TypeError as error:
            raise ValueError("An integer rank is required") from error
        if not 0 <= rank < self.count:
            raise ValueError("rank is outside the complete table fiber")
        state = np.zeros(self.mask.shape, dtype=bool)
        residual = self.initial
        for i, original_row in enumerate(self.rows):
            branches, cutpoints = self._branches(i, residual)
            choice = bisect_right(cutpoints, rank)
            _, pattern, residual, _ = branches[choice]
            if choice:
                rank -= cutpoints[choice - 1]
            for j in pattern:
                state[original_row, self.columns[j]] = True
        if rank or any(residual):
            raise AssertionError("Unranking did not terminate at the unique empty suffix")
        return state

    def occupancy_counts(self):
        """Exact cell numerators under this block's uniform labeled-table law."""
        if not self.count:
            raise ValueError("An empty feasible fiber has no occupancy probabilities")
        occupancy = np.zeros(self.mask.shape, dtype=object)
        prefix = {self.initial: 1}
        for i, original_row in enumerate(self.rows):
            if sum(p * self._counts[i, c] for c, p in prefix.items()) != self.count:
                raise AssertionError("Prefix/suffix layer does not reconcile to the whole fiber")
            following = defaultdict(int)
            for residual, prefix_count in prefix.items():
                for _, pattern, child, count in self._branches(i, residual)[0]:
                    following[child] += prefix_count
                    through_count = prefix_count * count
                    for j in pattern:
                        occupancy[original_row, self.columns[j]] += through_count
            prefix = following
        if sum(prefix.values()) != self.count:
            raise AssertionError("Complete forward paths do not reconcile to fiber count")
        if list(occupancy.sum(1)) != [r * self.count for r in self.row_margins]:
            raise AssertionError("Exact occupancy row margins do not reconcile")
        if list(occupancy.sum(0)) != [c * self.count for c in self.column_margins]:
            raise AssertionError("Exact occupancy column margins do not reconcile")
        if any(not 0 <= x <= self.count for x in occupancy.ravel()):
            raise AssertionError("Cell occupancy numerator is outside its denominator")
        return occupancy

    def moments(self, values):
        """Mean/variance of the sum score, using centered floating mixtures."""
        values = np.asarray(values, dtype=float)
        if not self.count or values.shape != self.mask.shape or not np.isfinite(values[self.mask]).all():
            raise ValueError("A nonempty fiber and finite scores on the mask are required")
        cache = {}

        def visit(i, residual):
            key = i, residual
            if key in cache:
                return cache[key]
            if i == len(self.rows):
                return 0., 0.
            row = self.rows[i]
            pieces = []
            for _, pattern, child, count in self._branches(i, residual)[0]:
                child_mean, child_variance = visit(i + 1, child)
                shift = math.fsum(values[row, self.columns[j]] for j in pattern)
                pieces.append((count / self._counts[key], shift + child_mean, child_variance))
            mean = math.fsum(p * branch_mean for p, branch_mean, _ in pieces)
            variance = math.fsum(p * (branch_variance + (branch_mean - mean) ** 2)
                                 for p, branch_mean, branch_variance in pieces)
            if not math.isfinite(mean) or not math.isfinite(variance) or variance < 0:
                raise ValueError("Nonfinite or negative conditional score moment")
            cache[key] = mean, variance
            return cache[key]

        return visit(0, self.initial)


def unrank_product(fibers, rank):
    total = math.prod(f.count for f in fibers)
    if not isinstance(rank, int) or not 0 <= rank < total:
        raise ValueError("Global rank is outside the product fiber")
    result = []
    for fiber in fibers:
        rank, component = divmod(rank, fiber.count)
        result.append(fiber.unrank(component))
    if rank:
        raise AssertionError("Global rank failed to split into component ranks")
    return result


def draw_global_ranks(counts, seed, samples):
    if samples < 1 or not counts or any(type(c) is not int or c < 1 for c in counts):
        raise ValueError("Positive exact component counts and a positive draw count are required")
    total, rng = math.prod(counts), random.Random(int(seed))
    for _ in range(samples):
        # Random.randrange uses unbiased getrandbits/rejection for arbitrary ints.
        # Repeated ranks are retained: these are independent draws with replacement.
        yield rng.randrange(total)


def validate_table(state, mask, row_margins, column_margins):
    if state.shape != mask.shape or not np.isin(state, (0, 1)).all():
        raise ValueError("An aligned binary table is required")
    if np.any(state & ~mask):
        raise ValueError("Table has a positive cell outside the fixed mask")
    if not np.array_equal(state.sum(1), row_margins) or not np.array_equal(state.sum(0), column_margins):
        raise ValueError("Table changed a fixed row or column margin")


def in_closer_tail(value, observed):
    return value <= observed + TIE_TOLERANCE


def tail_summary(count, samples):
    if type(count) is not int or type(samples) is not int or samples < 1 or not 0 <= count <= samples:
        raise ValueError("Valid integer tail and independent-draw counts are required")

    def upper(alpha):
        if count == samples:
            return 1.
        if not count:
            return -math.expm1(math.log(alpha) / samples)
        return float(beta.ppf(1 - alpha, count + 1, samples - count))

    return {"lower_tail_count_including_ties": count, "independent_draws": samples,
            "qualified_bound": True, "one_sided_95_upper": upper(.05),
            "simultaneous_three_scenario_95_upper": upper(.05 / 3),
            "method": "One-sided exact binomial Clopper-Pearson; three-scenario family uses Bonferroni alpha=.05/3",
            "exact_enumerated_tail_probability": None}


def write_gzip_csv(path, rows, fields):
    """CSV gzip bytes omit filenames and timestamps, preserving decimal ranks."""
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def load_inputs(root=ROOT):
    """Reconstruct the exact old labeled masks, states, scores and both margins."""
    root = Path(root)
    verify_pinned_files(root)
    old = json.loads((root / OLD / "recorded_source_diagnostic.json").read_bytes())
    with (root / OLD / "recorded_source_risk_cells.csv").open(newline="", encoding="utf8") as handle:
        all_cells = list(csv.DictReader(handle))
    with (root / OLD / "recorded_source_fixed_margins.csv").open(newline="", encoding="utf8") as handle:
        all_margins = list(csv.DictReader(handle))
    expected_names = {s[0] for s in SCENARIOS}
    if set(old["analyses"]) != expected_names or {r["analysis"] for r in all_cells} != expected_names:
        raise ValueError("Original scenario identities changed")
    if {r["analysis"] for r in all_margins} != expected_names:
        raise ValueError("Original margin scenario identities changed")
    risks = {}
    for name, _, source_min, late_min in SCENARIOS:
        original = old["analyses"][name]
        shape = tuple(original["risk_mask_shape"])
        mask, state = np.zeros(shape, dtype=bool), np.zeros(shape, dtype=bool)
        values = np.full(shape, np.nan)
        cells = []
        for raw in (r for r in all_cells if r["analysis"] == name):
            cell = dict(raw)
            i, j = int(cell["row_index"]), int(cell["column_index"])
            positive = int(cell["late_positive"])
            if not 0 <= i < shape[0] or not 0 <= j < shape[1] or mask[i, j] or positive not in (0, 1):
                raise ValueError("Duplicated, invalid or out-of-range frozen risk cell")
            if tuple(original["row_labels"][i]) != (cell["species_key"], cell["longitude_block"]):
                raise ValueError("Frozen species/block row identity changed")
            if original["column_labels"][j] != cell["target_island"]:
                raise ValueError("Frozen target column identity changed")
            if (int(cell["source_min_records"]), int(cell["late_min_records"])) != (source_min, late_min):
                raise ValueError("Frozen source/late thresholds changed")
            value = float(cell["log1p_source_upper_km"])
            if not math.isfinite(value) or abs(value - math.log1p(float(cell["nearest_source_upper_km"]))) > 1e-12:
                raise ValueError("Frozen statistic does not match upper coast-gap score")
            mask[i, j], state[i, j], values[i, j] = True, positive, value
            cell.update(row_index=i, column_index=j, late_positive=positive)
            cells.append(cell)
        if sha(mask.tobytes()) != original["risk_mask_sha256"] or sha(state.tobytes()) != original["observed_binary_state_sha256"]:
            raise ValueError("Reconstructed frozen mask or observed table hash differs")
        margins = [r for r in all_margins if r["analysis"] == name]
        row_margins, column_margins = [], []
        for dimension, size, destination in (("species_longitude_block", shape[0], row_margins),
                                              ("target_island", shape[1], column_margins)):
            records = [r for r in margins if r["dimension"] == dimension]
            keyed = {int(r["index"]): r for r in records}
            if len(records) != size or set(keyed) != set(range(size)):
                raise ValueError("Frozen margin indices are incomplete or duplicated")
            for i in range(size):
                r = keyed[i]
                if dimension == "species_longitude_block":
                    key = ":".join(original["row_labels"][i])
                    allowed = int(mask[i].sum())
                else:
                    key = original["column_labels"][i]
                    allowed = int(mask[:, i].sum())
                if r["unit_key"] != key or int(r["allowed_cells"]) != allowed:
                    raise ValueError("Frozen margin label or allowed-cell total differs")
                destination.append(int(r["positive_margin"]))
        if len(margins) != sum(shape):
            raise ValueError("Unexpected frozen margin dimensions")
        validate_table(state, mask, row_margins, column_margins)
        if not state.any():
            raise ValueError("Frozen scenario unexpectedly has no positive outcomes")
        observed = math.fsum(values[state]) / int(state.sum())
        if abs(observed - original["diagnostic"]["observed_T"]) > TIE_TOLERANCE:
            raise ValueError("Observed statistic does not reproduce the original receipt")
        risks[name] = {"analysis": name, "mask": mask, "state": state, "values": values,
                       "row_margins": np.asarray(row_margins), "column_margins": np.asarray(column_margins),
                       "row_labels": original["row_labels"], "column_labels": original["column_labels"],
                       "cells": cells, "observed_T": original["diagnostic"]["observed_T"],
                       "observed_T_fsum_recomputed": observed,
                       "risk_mask_sha256": original["risk_mask_sha256"],
                       "observed_binary_state_sha256": original["observed_binary_state_sha256"]}
    return old, risks


def build_blocks(risk):
    column_blocks = {}
    for cell in risk["cells"]:
        column, block = int(cell["column_index"]), cell["longitude_block"]
        if column in column_blocks and column_blocks[column] != block:
            raise ValueError("A target column occurs in multiple longitude blocks")
        column_blocks[column] = block
    if set(column_blocks) != set(range(risk["mask"].shape[1])):
        raise ValueError("Every frozen target requires an explicit longitude block")
    blocks, covered = [], np.zeros_like(risk["mask"])
    for label in sorted({r[1] for r in risk["row_labels"]}):
        rows = [i for i, r in enumerate(risk["row_labels"]) if r[1] == label]
        columns = [j for j in range(risk["mask"].shape[1]) if column_blocks[j] == label]
        selector = np.ix_(rows, columns)
        mask = risk["mask"][selector]
        covered[selector] = mask
        fiber = ExactBlockFiber(mask, risk["row_margins"][rows], risk["column_margins"][columns],
                                max_states=MAX_STATES, max_transitions=MAX_TRANSITIONS)
        if not fiber.count:
            raise ValueError("Observed feasible table unexpectedly has a zero completion count")
        blocks.append({"longitude_block": label, "rows": rows, "columns": columns, "fiber": fiber})
    if not np.array_equal(covered, risk["mask"]):
        raise ValueError("Longitude blocks do not partition the frozen structural mask")
    return blocks


def table_from_global_rank(blocks, rank, shape):
    state = np.zeros(shape, dtype=bool)
    components = unrank_product([b["fiber"] for b in blocks], rank)
    for block, component in zip(blocks, components):
        state[np.ix_(block["rows"], block["columns"])] = component
    return state


def summarize_fiber(risk, blocks):
    total = math.prod(b["fiber"].count for b in blocks)
    edges = int(risk["state"].sum())
    full_occupancy = np.zeros(risk["mask"].shape, dtype=object)
    block_rows = []
    for block in blocks:
        fiber = block["fiber"]
        selector = np.ix_(block["rows"], block["columns"])
        values, observed = risk["values"][selector], risk["state"][selector]
        numerators = fiber.occupancy_counts()
        mean, variance = fiber.moments(values)
        occupancy_mean = math.fsum(int(numerators[i, j]) / fiber.count * values[i, j]
                                   for i, j in zip(*np.where(fiber.mask)))
        if abs(mean - occupancy_mean) > 1e-10:
            raise AssertionError("Independent occupancy and centered-mixture mean calculations differ")
        full_occupancy[selector] = numerators * (total // fiber.count)
        forced = int(np.count_nonzero((numerators == fiber.count) & fiber.mask))
        variable = int(np.count_nonzero((numerators > 0) & (numerators < fiber.count) & fiber.mask))
        block_edges = int(observed.sum())
        observed_sum = math.fsum(values[observed])
        block_rows.append({"longitude_block": block["longitude_block"],
                           "fiber_cardinality_decimal": str(fiber.count),
                           "allowed_cells": int(fiber.mask.sum()),
                           "active_allowed_cells": int(fiber.mask[np.ix_(fiber.rows, fiber.columns)].sum()),
                           "active_rows": len(fiber.rows), "active_columns": len(fiber.columns),
                           "positive_cells": block_edges, "forced_positive_cells": forced,
                           "variable_allowed_cells": variable, "variable_positive_cells": block_edges - forced,
                           "never_positive_allowed_cells": int(np.count_nonzero((numerators == 0) & fiber.mask)),
                           "memo_states": fiber.memo_states, "attempted_pattern_transitions": fiber.attempted_transitions,
                           "row_order_original_indices": [block["rows"][i] for i in fiber.rows],
                           "active_column_original_indices": [block["columns"][j] for j in fiber.columns],
                           "observed_sum_log1p_gap": observed_sum,
                           "conditional_mean_sum_log1p_gap": mean,
                           "conditional_variance_sum_log1p_gap": variance,
                           "occupancy_mean_sum_log1p_gap": occupancy_mean,
                           "observed_minus_null_contribution_to_full_T": (observed_sum - mean) / edges,
                           "exact_occupancy_row_and_column_reconciliation": True})
    if list(full_occupancy.sum(1)) != [int(r) * total for r in risk["row_margins"]]:
        raise AssertionError("Whole-fiber integer occupancy row totals changed")
    if list(full_occupancy.sum(0)) != [int(c) * total for c in risk["column_margins"]]:
        raise AssertionError("Whole-fiber integer occupancy column totals changed")
    mean = math.fsum(b["conditional_mean_sum_log1p_gap"] for b in block_rows) / edges
    variance = math.fsum(b["conditional_variance_sum_log1p_gap"] for b in block_rows) / edges ** 2
    occupancy_mean = math.fsum(int(full_occupancy[i, j]) / total * risk["values"][i, j]
                               for i, j in zip(*np.where(risk["mask"]))) / edges
    if abs(mean - occupancy_mean) > TIE_TOLERANCE:
        raise AssertionError("Whole-fiber independent mean checks differ")
    occupancy_rows = []
    for cell in risk["cells"]:
        i, j = int(cell["row_index"]), int(cell["column_index"])
        numerator = int(full_occupancy[i, j])
        status = "forced_positive" if numerator == total else "never_positive" if not numerator else "variable"
        occupancy_rows.append({"analysis": risk["analysis"], "row_index": i, "column_index": j,
                               "species_key": cell["species_key"], "original_name": cell["original_name"],
                               "target_island": cell["target_island"], "longitude_block": cell["longitude_block"],
                               "observed_late_positive": int(cell["late_positive"]),
                               "probability_numerator": str(numerator), "probability_denominator": str(total),
                               "occupancy_probability": numerator / total, "occupancy_status": status,
                               "log1p_source_upper_km": risk["values"][i, j]})
    forced = sum(b["forced_positive_cells"] for b in block_rows)
    summary = {"fiber_cardinality_decimal": str(total), "observed_T": risk["observed_T"],
               "null_mean_T": mean, "null_variance_T": variance, "null_sd_T": math.sqrt(variance),
               "occupancy_null_mean_T": occupancy_mean,
               "observed_to_null_geometric_mean_one_plus_gap_ratio": math.exp(risk["observed_T"] - mean),
               "support": {"candidate_cells": int(risk["mask"].sum()), "positive_cells": edges,
                           "forced_positive_cells": forced, "variable_positive_cells": edges - forced,
                           "variable_allowed_cells": sum(b["variable_allowed_cells"] for b in block_rows),
                           "never_positive_allowed_cells": sum(b["never_positive_allowed_cells"] for b in block_rows)},
               "moment_arithmetic": "Exact integer completion/occupancy counts; binary64 scores and centered weighted mean/variance arithmetic",
               "blocks": block_rows}
    return summary, occupancy_rows


def sample_scenario(risk, blocks, seed, samples, path):
    """Write all fixed-budget rank/T pairs and check every reconstructed table."""
    values = np.empty(samples, dtype=float)
    edges, checked, lower_tail_count = int(risk["state"].sum()), 0, 0
    started = time.perf_counter()

    def records():
        nonlocal checked, lower_tail_count
        for number, rank in enumerate(draw_global_ranks([b["fiber"].count for b in blocks], seed, samples), 1):
            state = table_from_global_rank(blocks, rank, risk["mask"].shape)
            validate_table(state, risk["mask"], risk["row_margins"], risk["column_margins"])
            value = math.fsum(risk["values"][state]) / edges
            if not math.isfinite(value):
                raise ValueError("Nonfinite sampled statistic")
            values[number - 1] = value
            checked += 1
            lower_tail_count += int(in_closer_tail(value, risk["observed_T"]))
            yield {"analysis": risk["analysis"], "sample": number, "seed": seed,
                   "global_rank": str(rank), "mean_log1p_source_gap_km": repr(value)}

    write_gzip_csv(path, records(), ["analysis", "sample", "seed", "global_rank", "mean_log1p_source_gap_km"])
    if checked != samples:
        raise RuntimeError("Fixed independent sample budget was not completed")
    mean = math.fsum(values) / samples
    return {"seed": seed, "stored_samples": samples, "stored_margin_checks": checked,
            "lower_tail_count_including_ties": lower_tail_count,
            "sample_mean_T": mean,
            "sample_variance_T": math.fsum((x - mean) ** 2 for x in values) / (samples - 1) if samples > 1 else None,
            "elapsed_seconds": time.perf_counter() - started,
            "sample_file": Path(path).name, "sample_file_sha256": sha(Path(path).read_bytes()),
            "qualified_bound": False}


def run(output, root=ROOT):
    root, output = Path(root), Path(output)
    verify_pinned_files(root)
    output.mkdir(parents=True, exist_ok=False)
    script_sha = sha(Path(__file__).read_bytes())
    lock = {"schema": "chocho_exact_conditional_source_null_run_lock_20261030_v1",
            "started_utc": utc(), "plan": PLAN, "plan_sha256": PLAN_SHA256,
            "input_sha256": PINS, "script_sha256": script_sha,
            "python_version": platform.python_version(), "numpy_version": np.__version__,
            "samples_per_scenario": SAMPLES, "sampling_with_replacement": True,
            "scenario_order_and_seeds": [{"analysis": name, "seed": seed} for name, seed, _, _ in SCENARIOS],
            "max_memo_states_per_block": MAX_STATES, "max_attempted_pattern_transitions_per_block": MAX_TRANSITIONS,
            "tail_tie_tolerance": TIE_TOLERANCE,
            "marginal_upper_alpha": .05, "simultaneous_three_scenario_upper_alpha": .05 / 3}
    lock_path = output / "exact_null_run_lock.json"
    lock_path.write_text(json.dumps(lock, indent=2, sort_keys=True, allow_nan=False) + "\n")
    lock_sha = sha(lock_path.read_bytes())
    analyses, occupancy_rows, block_csv_rows = {}, [], []
    stage = "load_pinned_inputs"
    try:
        old, risks = load_inputs(root)
        for name, seed, source_min, late_min in SCENARIOS:
            stage = name + ":exact_counts_and_moments"
            risk = risks[name]
            blocks = build_blocks(risk)
            summary, cells = summarize_fiber(risk, blocks)
            summary.update(source_min_records=source_min, late_min_records=late_min,
                           observed_T_fsum_recomputed=risk["observed_T_fsum_recomputed"],
                           risk_mask_sha256=risk["risk_mask_sha256"],
                           observed_binary_state_sha256=risk["observed_binary_state_sha256"],
                           old_mcmc_status=old["analyses"][name]["diagnostic"]["status"])
            analyses[name] = summary
            occupancy_rows.extend(cells)
            for block in summary["blocks"]:
                block_csv_rows.append({"analysis": name, **{k: v for k, v in block.items()
                                                            if not isinstance(v, list)}})
            print(json.dumps({"stage": "exact_counts_complete", "analysis": name,
                              "fiber_cardinality_decimal": summary["fiber_cardinality_decimal"],
                              "block_counts": [b["fiber_cardinality_decimal"] for b in summary["blocks"]]}), flush=True)
            stage = name + ":fixed_independent_draws"
            summary["sampling"] = sample_scenario(risk, blocks, seed, SAMPLES, output / f"{name}_rank_T_samples.csv.gz")
            print(json.dumps({"stage": "fixed_independent_draws_complete", "analysis": name,
                              "stored_samples": summary["sampling"]["stored_samples"],
                              "stored_margin_checks": summary["sampling"]["stored_margin_checks"]}), flush=True)
        stage = "final_source_and_completion_checks"
        verify_pinned_files(root)
        if sha(Path(__file__).read_bytes()) != script_sha or sha(lock_path.read_bytes()) != lock_sha:
            raise ValueError("Script or immutable pre-run lock changed during execution")
        if set(analyses) != {s[0] for s in SCENARIOS} or any(
                a["sampling"]["stored_samples"] != SAMPLES or a["sampling"]["stored_margin_checks"] != SAMPLES
                for a in analyses.values()):
            raise RuntimeError("All fixed scenario budgets and margin checks must complete")
        write_gzip_csv(output / "exact_cell_occupancy_counts.csv.gz", occupancy_rows, list(occupancy_rows[0]))
        write_gzip_csv(output / "exact_block_counts_and_moments.csv.gz", block_csv_rows, list(block_csv_rows[0]))
        for summary in analyses.values():
            summary["tail"] = tail_summary(summary["sampling"]["lower_tail_count_including_ties"], SAMPLES)
            summary["sampling"]["qualified_bound"] = True
            summary["sampling"]["sample_mean_minus_DP_mean_T"] = summary["sampling"]["sample_mean_T"] - summary["null_mean_T"]
            summary["sampling"]["sample_mean_MC_standard_error_from_DP"] = math.sqrt(summary["null_variance_T"] / SAMPLES)
        report = {"schema": "chocho_exact_conditional_source_null_20261030_v1",
                  "status": "COMPLETE_EXACT_COUNT_IID_MONTE_CARLO_BOUNDS", "qualified_bounds": True,
                  "started_utc": lock["started_utc"], "finished_utc": utc(),
                  "plan": PLAN, "plan_sha256": PLAN_SHA256, "run_lock_sha256": lock_sha,
                  "script_sha256": script_sha, "input_sha256": PINS,
                  "method": {"counting": "Exact Python-integer labeled-row suffix DP with structural zeros and fixed margins; longitude-block product",
                             "sampling": "One unbiased arbitrary-integer global rank per draw; mixed-radix block ranks and exact completion-count unranking; independent draws with replacement",
                             "uniformity_proof": "Each row-pattern branch occupies an integer interval of its child completion count. By induction every labeled block table has exactly one rank. Mixed-radix decomposition is a bijection between product ranks and block-rank tuples. Uniform independent ranks therefore give uniform independent whole tables.",
                             "count_cap_states_per_block": MAX_STATES, "count_cap_transitions_per_block": MAX_TRANSITIONS,
                             "score_arithmetic": "Stored binary64 log1p upper-distance scores; math.fsum statistic; centered mixture moments and independent exact-occupancy-weighted mean check",
                             "tail": "T_null <= original observed T + 1e-12; exact binomial one-sided95% upper bounds, plus simultaneous95% three-scenario Bonferroni bounds",
                             "not_enumerated_tail": True, "no_MCMC_ESS_or_tail_hit_gate": True,
                             "sample_validation": "Every one of300000 complete reconstructed binary tables checked against full structural mask and both original margins"},
                  "interpretation": ["Post-pilot numerical validation of the same conditional recording comparison; not an untouched ecological test.",
                                     "Original MCMC output, failed gates, rows, columns, distances, sources and margins remain unchanged.",
                                     "Exact refers to integer table counts and the uniform sampling law; score moments use floating arithmetic and tail probabilities retain finite Monte Carlo uncertainty.",
                                     "Zero observed tail hits never imply zero conditional tail probability. Bounds concern the fixed uniform-table probability, not an ecological effect-size interval.",
                                     "Forced positive cells contribute to every full statistic and the original positive-count denominator.",
                                     "Early nondetection is not absence; later recording is not arrival, dispersal or establishment. Mainland and most island sources, habitat/climate and observer processes remain omitted.",
                                     "The geometric-distance ratio is descriptive and conditional on the frozen roster/margins, not causal or a colonization probability."],
                  "analyses": analyses,
                  "output_sha256": {p.name: sha(p.read_bytes()) for p in sorted(output.glob("*.csv.gz"))}}
        (output / "exact_null_diagnostic.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
        return report
    except BaseException as error:
        failure = {"schema": "chocho_exact_conditional_source_null_20261030_v1",
                   "status": "INCOMPLETE_NO_QUALIFIED_BOUND", "qualified_bounds": False,
                   "failed_stage": stage, "error_type": type(error).__name__, "error": str(error),
                   "started_utc": lock["started_utc"], "finished_utc": utc(), "run_lock_sha256": lock_sha,
                   "script_sha256": script_sha, "plan_sha256": PLAN_SHA256, "input_sha256": PINS,
                   "partial_output_sha256": {p.name: sha(p.read_bytes()) for p in sorted(output.glob("*.csv.gz"))}}
        (output / "exact_null_diagnostic.json").write_text(json.dumps(failure, indent=2, sort_keys=True, allow_nan=False) + "\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/exploratory/island_source_exact_null_20261030")
    args = parser.parse_args()
    run(args.output_dir)
