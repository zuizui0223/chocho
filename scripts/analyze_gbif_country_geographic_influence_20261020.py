#!/usr/bin/env python3
"""Audit posthoc geographic dependence of the completed country-grain analysis.

Uses ONLY frozen public species × Level-1 strata; no new host or occurrence
records, no new biological or causal hypothesis test. Country exclusions change
the conditional reference pool, so recompute AK/N after every exclusion.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import comb, isclose
from pathlib import Path

CARIBBEAN = frozenset(("CU", "DO", "HT"))
CENTRAL_AMERICA = frozenset(("GT", "HN", "NI"))
BASELINE = {"candidate_strata": 255, "candidate_cells": 3188,
            "added_cells": 1527, "positive_cells": 48,
            "observed_added_positive": 26, "informative_strata": 21}


def validate(frozen):
    if not isinstance(frozen, list) or len(frozen) != 255:
        raise ValueError("Frozen primary strata changed")
    seen = set()
    for row in frozen:
        key = (row["species"], row["stratum"])
        if key in seen or row["genus"] != row["species"].split()[0]:
            raise ValueError("Duplicate species/region or genus mismatch")
        seen.add(key)
        cells = row["cells"]
        if len({c["country"] for c in cells}) != len(cells):
            raise ValueError("Country repeated within a stratum")
        if any(c["positive"] not in (0, 1) or c["added"] not in (0, 1)
               or not 0 <= c["area_score"] <= 1 for c in cells):
            raise ValueError("Invalid frozen binary/area code")
        N = len(cells)
        A = sum(c["added"] for c in cells)
        K = sum(c["positive"] for c in cells)
        O = sum(c["added"] * c["positive"] for c in cells)
        if not (0 < A < N):
            raise ValueError("Frozen stratum has no exposure contrast")
        if (N, A, K, O) != (row["N"], row["A"], row["K"], row["O"]):
            raise ValueError("Stored primary margin mismatch")
        if not isclose(A * K / N, row["expected"], abs_tol=1e-12):
            raise ValueError("Stored primary expectation mismatch")


def exact_upper_tail(rows, observed):
    """Convolve conditional hypergeometrics exactly, using rational masses."""
    probabilities = [Fraction(1, 1)]
    for row in rows:
        N, A, K = row["N"], row["A"], row["K"]
        denom = comb(N, K)
        next_probs = [Fraction(0, 1)] * (len(probabilities) + min(A, K))
        for z in range(max(0, K - N + A), min(A, K) + 1):
            p = Fraction(comb(A, z) * comb(N - A, K - z), denom)
            for index, previous in enumerate(probabilities):
                next_probs[index + z] += previous * p
        probabilities = next_probs
    if sum(probabilities) != 1:
        raise ValueError("Conditional reference mass lost")
    return float(sum(probabilities[observed:]))


def round_nonnegative_fraction(value: Fraction, digits: int) -> float:
    """Half-even decimal rounding of a rational *before* floating conversion.

    Rounding float(Fraction) can round in the wrong direction near a tie;
    the exact Ireland-out expectation ends ...467499, not ...467500.
    """
    if value < 0 or digits < 0:
        raise ValueError("Expected nonnegative Fraction and decimal places")
    scale = 10 ** digits
    scaled = value * scale
    whole, remainder = divmod(scaled.numerator, scaled.denominator)
    if 2 * remainder > scaled.denominator or (
            2 * remainder == scaled.denominator and whole % 2):
        whole += 1
    return whole / scale


def summarize(frozen, *, omit_countries=frozenset(),
              omit_genera=frozenset(), included_level1=None,
              excluded_level1=frozenset()):
    rows = []
    for row in frozen:
        if row["genus"] in omit_genera or row["stratum"] in excluded_level1:
            continue
        if included_level1 is not None and row["stratum"] not in included_level1:
            continue
        cells = [c for c in row["cells"] if c["country"] not in omit_countries]
        N = len(cells)
        A = sum(c["added"] for c in cells)
        if N == 0 or not 0 < A < N:
            continue
        K = sum(c["positive"] for c in cells)
        O = sum(c["added"] * c["positive"] for c in cells)
        rows.append({"N": N, "A": A, "K": K, "O": O,
                     "genus": row["genus"], "E": Fraction(A * K, N)})
    N = sum(r["N"] for r in rows)
    A = sum(r["A"] for r in rows)
    K = sum(r["K"] for r in rows)
    O = sum(r["O"] for r in rows)
    E = sum((r["E"] for r in rows), Fraction(0))
    informative = sum(0 < r["K"] < r["N"] for r in rows)
    result = {
        "candidate_strata": len(rows), "candidate_cells": N,
        "added_cells": A, "positive_cells": K, "observed_added_positive": O,
        "expected_added_positive": round_nonnegative_fraction(E, 12),
        "excess_pp": round(100 * float((O - E) / K), 9) if K else None,
        "informative_strata": informative,
        "independent_positive_genera": len({r["genus"] for r in rows if r["K"] > 0}),
        "conditional_geographic_upper_tail": round(exact_upper_tail(rows, O), 12),
    }
    return result


def build(frozen):
    validate(frozen)
    cb, ca = CARIBBEAN, CENTRAL_AMERICA
    results = {
        "baseline": summarize(frozen),
        "exclude_caribbean": summarize(frozen, omit_countries=cb),
        "exclude_central_america": summarize(frozen, omit_countries=ca),
        "exclude_caribbean_and_central_america": summarize(frozen, omit_countries=cb | ca),
        "only_level1_8": summarize(frozen, included_level1={"8"}),
        "exclude_level1_8": summarize(frozen, excluded_level1={"8"}),
    }
    if any(results["baseline"][key] != value for key, value in BASELINE.items()):
        raise AssertionError("Not the locked primary country comparison")
    if not isclose(results["baseline"]["conditional_geographic_upper_tail"],
                   0.283249574739, abs_tol=1e-10):
        raise AssertionError("Original exact geographical null not reproduced")
    baseline = results["baseline"]["excess_pp"]
    positive = {}
    for row in frozen:
        for cell in row["cells"]:
            if cell["positive"]:
                pos = positive.setdefault(cell["country"], [0, 0])
                pos[0] += 1
                pos[1] += cell["added"]
    country_rows = []
    for country in sorted({c["country"] for row in frozen for c in row["cells"]}):
        r = {"country": country,
             "omitted_positive_cells": positive.get(country, [0, 0])[0],
             "omitted_added_positive_cells": positive.get(country, [0, 0])[1],
             **summarize(frozen, omit_countries={country})}
        r["delta_from_baseline_pp"] = round(r["excess_pp"] - baseline, 9) if r["excess_pp"] is not None else None
        country_rows.append(r)
    genus_rows = []
    for genus in sorted({row["genus"] for row in frozen}):
        r = {"genus": genus, **summarize(frozen, omit_genera={genus})}
        r["delta_from_baseline_pp"] = round(r["excess_pp"] - baseline, 9) if r["excess_pp"] is not None else None
        genus_rows.append(r)
    return {"schema": "gbif_country_geographic_dependence_posthoc_v1",
            "derived_date": "2026-10-20",
            "source_file": "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json",
            "source_sha": "53662700f6fade11c5c12ef7bacebe130d03cf30",
            "status": "POSTHOC_DIAGNOSTIC_BELOW_INFORMATION_GATE",
            "estimator": "Within the retained species × WGSRPD Level1 strata, discard specified country/genus cells, retain 0<A<N strata, then calculate 100*sum(O - AK/N)/sum(K). Exact conditional geographical upper tails are a reference, NOT biological or causal p-values.",
            "scope": "Not a prospective or independent test; source is the same 229-taxon country panel as primary. Zero GBIF records are not biological absences.",
            "information_minimum_required": 50, "results": results,
            "leave_one_country_out": country_rows, "leave_one_genus_out": genus_rows,
            "omission_groups": {"caribbean": ["CU", "DO", "HT"],
                                "central_america": ["GT", "HN", "NI"],
                                "level1_8_formal_name": "Southern America (includes Caribbean and Central America)"},
            "interpretation": "Geographic sign sensitivity is not evidence of location-specific biological effects. All retained data remain below information threshold."}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen-strata", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = build(json.loads(args.frozen_strata.read_text(encoding="utf-8")))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["results"], indent=2))


if __name__ == "__main__":
    main()
