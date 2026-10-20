#!/usr/bin/env python3
"""Post-protocol GBIF *record* turnover, never plant arrival or butterfly colonization.

Protocol: provenance/reviewer_defenses/protocols/GBIF_COUNTRY_RECORD_TURNOVER_20261020.md
Reads only frozen, public country aggregate counts and primary eligible cells.
Both windows are 8 years long; no 2026 partial year is used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

SCOPES = ("all_accepted_bases", "core")
WINDOWS = ("2010_2017", "2018_2025")
COLLISION = frozenset(("Oeneis bore", "Oeneis polixenes"))
EXPECTED_SOURCE_RECORDS = {"all_accepted_bases": 33_817_891, "core": 33_621_039}
EXPECTED_COUNTRY_ROWS = 23433
EXPECTED_FROZEN_MARGIN = (255, 3188, 1527, 48, 26)
SEED = 20261020


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_sources(csv_path: Path, manifest_path: Path, strata_path: Path):
    source = json.loads(manifest_path.read_text(encoding="utf-8"))
    if len(source["species_receipts"]) != 239:
        raise ValueError("Frozen original 239-name manifest missing")
    names = {}
    keys = defaultdict(set)
    for row in source["species_receipts"]:
        if row["status"] == "TAXON_UNRESOLVED":
            continue
        sp = row["species"]
        names[sp] = row
        keys[row["taxon"]["key"]].add(sp)
    colliding = {sp for ns in keys.values() if len(ns) > 1 for sp in ns}
    if colliding != COLLISION or len(names) != 231:
        raise ValueError("Accepted-key collisions or resolved species changed")
    species = frozenset(names) - COLLISION
    if len(species) != 229:
        raise ValueError("Expected 229 noncolliding taxa")
    counts = {}
    rows = 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or ()) != {"species", "window", "country", "basis_scope", "record_count"}:
            raise ValueError("Country facet CSV schema drift")
        for item in reader:
            rows += 1
            sp, ct, win, scope = item["species"], item["country"], item["window"], item["basis_scope"]
            if sp not in names or win not in WINDOWS or scope not in SCOPES or len(ct) != 2:
                raise ValueError("Invalid country/species/scope/window")
            key = (scope, win, sp, ct)
            if key in counts:
                raise ValueError("Duplicate facet entry")
            n = int(item["record_count"])
            if n <= 0:
                raise ValueError("Zero/negative count cannot be a positive country facet")
            if sp not in COLLISION:
                counts[key] = n
    if rows != EXPECTED_COUNTRY_ROWS:
        raise ValueError(f"Country facets truncated: {rows} vs {EXPECTED_COUNTRY_ROWS}")
    totals = {scope: sum(n for (s, _, _, _), n in counts.items() if s == scope)
              for scope in SCOPES}
    if totals != EXPECTED_SOURCE_RECORDS:
        raise ValueError(f"Frozen source-country totals changed: {totals}")
    strata = json.loads(strata_path.read_text(encoding="utf-8"))
    if len(strata) != EXPECTED_FROZEN_MARGIN[0]:
        raise ValueError("Frozen stratum count changed")
    margins = (len(strata), sum(x["N"] for x in strata), sum(x["A"] for x in strata),
               sum(x["K"] for x in strata), sum(x["O"] for x in strata))
    if margins != EXPECTED_FROZEN_MARGIN:
        raise ValueError(f"Frozen cell margins changed: {margins}")
    seen = set()
    for st in strata:
        if st["species"] not in species or st["genus"] != st["species"].split()[0]:
            raise ValueError("Frozen candidate species invalid")
        for c in st["cells"]:
            key = st["species"], c["country"]
            if key in seen:
                raise ValueError("Species-country candidate duplicated")
            seen.add(key)
            pooled = sum(counts.get(("all_accepted_bases", win, *key), 0)
                         for win in WINDOWS)
            if int(pooled > 0) != c["positive"]:
                raise ValueError("Frozen pooled presence differs from country facets")
    background = defaultdict(set)
    for (scope, win, sp, ct), n in counts.items():
        if n > 0:
            background[scope, win, ct].add(sp)
    return strata, counts, background, {
        "country_counts_sha256": source_hash(csv_path),
        "original_manifest_sha256": source_hash(manifest_path),
        "primary_strata_sha256": source_hash(strata_path),
        "source_taxa": len(species), "original_names": 239,
        "country_rows_including_excluded_aliases": rows, "source_record_totals": totals,
    }


def states_for(strata, counts, background, scope, effort_gate, focal_min):
    """Build a paired-effort-qualified dataset without changing host exposures."""
    by_st = []
    overall = {"added": Counter(), "control": Counter()}
    countries_with_paired_effort = set()
    for st in strata:
        sp = st["species"]
        cells = []
        for c in st["cells"]:
            ct = c["country"]
            other_early = len(background[scope, WINDOWS[0], ct] - {sp})
            other_late = len(background[scope, WINDOWS[1], ct] - {sp})
            if min(other_early, other_late) < effort_gate:
                continue
            countries_with_paired_effort.add(ct)
            e = counts.get((scope, WINDOWS[0], sp, ct), 0)
            l = counts.get((scope, WINDOWS[1], sp, ct), 0)
            early, late = int(e >= focal_min), int(l >= focal_min)
            state = f"{early}{late}"
            added = int(c["added"])
            overall["added" if added else "control"][state] += 1
            cells.append({"country": ct, "added": added, "early": early,
                          "late": late, "state": state, "early_count": e,
                          "late_count": l, "other_early": other_early,
                          "other_late": other_late})
        by_st.append({"species": sp, "genus": st["genus"],
                      "stratum": st["stratum"], "cells": cells})
    return by_st, {
        "paired_effort_eligible_countries": len(countries_with_paired_effort),
        "paired_effort_candidate_cells": sum(len(x["cells"]) for x in by_st),
        "states_added": dict(sorted(overall["added"].items())),
        "states_control": dict(sorted(overall["control"].items())),
    }


def build_direction_rows(by_st, direction):
    if direction not in ("late_only", "early_only"):
        raise ValueError("Unknown directional record-only outcome")
    eligibility = "early" if direction == "late_only" else "late"
    outcome = "late" if direction == "late_only" else "early"
    rows = []
    for st in by_st:
        cells = [c for c in st["cells"] if c[eligibility] == 0]
        N = len(cells)
        A = sum(c["added"] for c in cells)
        if not 0 < A < N:
            continue
        K = sum(c[outcome] for c in cells)
        O = sum(c[outcome] * c["added"] for c in cells)
        rows.append({"species": st["species"], "genus": st["genus"],
                     "stratum": st["stratum"], "N": N, "A": A, "K": K, "O": O,
                     "E": Fraction(A * K, N),
                     "positive_countries": [c["country"] for c in cells if c[outcome]],
                     "added_positive_countries": [c["country"] for c in cells if c[outcome] and c["added"]],
                     "positive_records": [{
                         "species": st["species"], "country": c["country"],
                         "botanical_level1": st["stratum"],
                         "added_host_country": c["added"],
                         "early_record_count": c["early_count"],
                         "late_record_count": c["late_count"],
                         "other_taxa_early": c["other_early"],
                         "other_taxa_late": c["other_late"]
                     } for c in cells if c[outcome]]})
    return rows


def exact_conditional_upper_tail(rows, observed):
    distribution = [Fraction(1)]
    for r in rows:
        N, A, K = (r[x] for x in ("N", "A", "K"))
        denom = math.comb(N, K)
        new = [Fraction(0)] * (len(distribution) + min(A, K))
        for z in range(max(0, K + A - N), min(A, K) + 1):
            p = Fraction(math.comb(A, z) * math.comb(N - A, K - z), denom)
            for idx, old in enumerate(distribution):
                new[idx + z] += old * p
        distribution = new
    if sum(distribution) != 1:
        raise RuntimeError("Conditional hypergeometric mass not normalized")
    return float(sum(distribution[observed:]))


def quantile(values, q):
    values = sorted(values)
    if not values:
        return None
    x = (len(values) - 1) * q
    lo = math.floor(x)
    hi = math.ceil(x)
    return values[lo] + (values[hi] - values[lo]) * (x - lo)


def cluster_ci(rows, label, draws):
    bygenus = defaultdict(lambda: [Fraction(0), 0])
    for row in rows:
        bygenus[row["genus"]][0] += Fraction(row["O"]) - row["E"]
        bygenus[row["genus"]][1] += row["K"]
    g = [value for _, value in sorted(bygenus.items())]
    if not g:
        return None, 0
    # Use common random numbers for overlapping source scopes: identical
    # stratum/genus data must produce identical confidence intervals.
    rng = random.Random(SEED)
    stats = []
    for _ in range(draws):
        choice = [g[rng.randrange(len(g))] for _ in range(len(g))]
        denominator = sum(x[1] for x in choice)
        if denominator:
            numerator = sum((x[0] for x in choice), Fraction(0))
            stats.append(100 * float(numerator / denominator))
    if len(stats) < draws * .9:
        return None, len(stats)
    return [round(quantile(stats, .025), 6),
            round(quantile(stats, .975), 6)], len(stats)


def result_for_direction(rows, label, draws=499):
    N, A = sum(r["N"] for r in rows), sum(r["A"] for r in rows)
    K, O = sum(r["K"] for r in rows), sum(r["O"] for r in rows)
    E = sum((r["E"] for r in rows), Fraction(0))
    informative = sum(0 < r["K"] < r["N"] for r in rows)
    positive_species = len({r["species"] for r in rows if r["K"] > 0})
    positive_genera = len({r["genus"] for r in rows if r["K"] > 0})
    enough = informative >= 50 and positive_species >= 25 and K >= 25
    ci, valid = cluster_ci(rows, label, draws) if K else (None, 0)
    pos = Counter(ct for r in rows for ct in r["positive_countries"])
    cell_details = sorted(
        (cell for r in rows for cell in r["positive_records"]),
        key=lambda x: (x["species"], x["country"], x["botanical_level1"]))
    if len(cell_details) != K or sum(x["added_host_country"] for x in cell_details) != O:
        raise ValueError("Transition record-level audit disagrees with conditional margins")
    return {
        "status": ("EXPLORATORY_INFORMATION_GATE_MET" if enough else
                   "BELOW_PROTOCOL_INFORMATION_THRESHOLD"),
        "candidate_species_region_strata": len(rows), "candidate_cells": N,
        "added_candidate_cells": A, "positive_record_cells": K,
        "observed_added_positive": O,
        "conditional_expected_added": round(float(E), 9),
        "excess_percentage_points": round(100 * float((O - E) / K), 6) if K else None,
        "informative_strata": informative, "positive_species": positive_species,
        "positive_genera": positive_genera,
        "genus_bootstrap_95_interval_pp": ci,
        "genus_bootstrap_valid_draws": valid,
        "exact_conditional_geographic_upper_tail": (round(exact_conditional_upper_tail(rows, O), 9)
                                                     if K else None),
        "positive_country_counts": dict(sorted(pos.items())),
        "positive_species_country_details": cell_details,
        "positive_window_single_record_cells": sum(
            x["late_record_count" if "late_only" in label else "early_record_count"] == 1
            for x in cell_details),
        "interpretation": "Indexed GBIF country records; NOT range expansion, colonization, breeding, host use, or a causal p-value",
    }


def run(strata, counts, background, source_info, *, bootstrap_draws=499):
    grid = []
    for scope in SCOPES:
        for effort_gate in (1, 3, 5):
            for focal_min in (1, 5):
                label = f"{scope}/other{effort_gate}/focal{focal_min}"
                by_st, states = states_for(strata, counts, background, scope, effort_gate, focal_min)
                first = build_direction_rows(by_st, "late_only")
                reverse = build_direction_rows(by_st, "early_only")
                grid.append({
                    "setting": label, "source_scope": scope,
                    "other_taxa_both_windows_min": effort_gate,
                    "focal_records_in_positive_window_min": focal_min,
                    "all_paired_states": states,
                    "late_only": result_for_direction(first, label+"/late_only", bootstrap_draws),
                    "early_only": result_for_direction(reverse, label+"/early_only", bootstrap_draws),
                })
    if len(grid) != 12:
        raise RuntimeError("Did not execute all prespecified 12 settings")
    p = next(x for x in grid if x["setting"] == "all_accepted_bases/other3/focal1")

    # Region deletions are explicitly POSTHOC probes of where the finite
    # indexed records come from, not stratified predictions of recruitment.
    base_states, _ = states_for(strata, counts, background,
                                "all_accepted_bases", 3, 1)
    caribbean = {"CU", "DO", "HT"}
    central = {"GT", "HN", "NI"}
    region_omissions = {
        "exclude_caribbean": caribbean,
        "exclude_central_america": central,
        "exclude_caribbean_and_central_america": caribbean | central,
    }
    regional_posthoc = {}
    for region_name, omitted in region_omissions.items():
        reduced = [{**st, "cells": [c for c in st["cells"]
                                   if c["country"] not in omitted]}
                   for st in base_states]
        regional_posthoc[region_name] = {
            direction: result_for_direction(
                build_direction_rows(reduced, direction),
                f"posthoc/{region_name}/{direction}", bootstrap_draws)
            for direction in ("late_only", "early_only")
        }
    reduced_level1 = [st for st in base_states if st["stratum"] != "8"]
    regional_posthoc["exclude_level1_8"] = {
        direction: result_for_direction(
            build_direction_rows(reduced_level1, direction),
            f"posthoc/exclude_level1_8/{direction}", bootstrap_draws)
        for direction in ("late_only", "early_only")
    }
    return {
        "schema": "gbif_country_two_eight_year_windows_record_turnover_exploratory_v1",
        "protocol": "provenance/reviewer_defenses/protocols/GBIF_COUNTRY_RECORD_TURNOVER_20261020.md",
        "protocol_status": "FROZEN_BEFORE_COUNTRY_TRANSITIONS_INSPECTED",
        "not_causal": True, "not_population_establishment": True,
        "observation_windows": WINDOWS,
        "source": source_info,
        "primary_setting": p["setting"], "primary": p,
        "sensitivity_grid": grid,
        "regional_posthoc_diagnostics": regional_posthoc,
        "information_gate": {"minimum_informative_strata": 50,
                             "minimum_positive_species": 25,
                             "minimum_positive_record_cells": 25},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--country-counts", required=True, type=Path)
    parser.add_argument("--original-manifest", required=True, type=Path)
    parser.add_argument("--frozen-strata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bootstrap-draws", type=int, default=499)
    args = parser.parse_args()
    if not 1 <= args.bootstrap_draws <= 2999:
        raise ValueError("Bootstrap draws must be in 1..2999")
    strata, counts, background, source_info = read_sources(
        args.country_counts, args.original_manifest, args.frozen_strata)
    result = run(strata, counts, background, source_info, bootstrap_draws=args.bootstrap_draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    compact = [{
        "setting": x["setting"], "states": x["all_paired_states"],
        "late_only": {k:v for k,v in x["late_only"].items() if k != "positive_country_counts"},
        "early_only": {k:v for k,v in x["early_only"].items() if k != "positive_country_counts"},
    } for x in result["sensitivity_grid"]]
    print("GBIF_COUNTRY_TEMPORAL_SUMMARY=" + json.dumps({
        "source": source_info, "primary_setting": result["primary_setting"],
        "settings": compact,
        "primary_positive_country_counts": {
            "late_only": result["primary"]["late_only"]["positive_country_counts"],
            "early_only": result["primary"]["early_only"]["positive_country_counts"]},
        "primary_positive_species_country_details": {
            "late_only": result["primary"]["late_only"]["positive_species_country_details"],
            "early_only": result["primary"]["early_only"]["positive_species_country_details"]},
        "regional_posthoc_diagnostics": {
            key: {direction: {k:v for k,v in result["regional_posthoc_diagnostics"][key][direction].items()
                              if k != "positive_country_counts"}
                  for direction in ("late_only","early_only")}
            for key in result["regional_posthoc_diagnostics"]}
        }, separators=(",", ":")))


if __name__ == "__main__":
    main()
