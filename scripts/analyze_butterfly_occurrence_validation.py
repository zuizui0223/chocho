#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


TOTAL_WGSRPD3_UNITS = 369


def hypergeom_pmf(n_population: int, n_success: int, draws: int) -> list[float]:
    lo = max(0, draws - (n_population - n_success))
    hi = min(draws, n_success)
    denom = math.comb(n_population, draws)
    out = [0.0] * (draws + 1)
    for x in range(lo, hi + 1):
        out[x] = (
            math.comb(n_success, x)
            * math.comb(n_population - n_success, draws - x)
            / denom
        )
    return out


def convolve(a: list[float], b: list[float]) -> list[float]:
    out = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y:
                out[i + j] += x * y
    return out


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    outside = 0
    recovered = 0
    expected = 0.0
    distribution = [1.0]
    informative_species = 0
    for row in rows:
        n = int(row["observed_outside_native"])
        if n <= 0:
            continue
        native = int(row["native_host_units"])
        contemporary = int(row["contemporary_host_units"])
        N = TOTAL_WGSRPD3_UNITS - native
        K = contemporary - native
        x = int(row["outside_native_units_explained_by_introduced_host_ranges"])
        if not (0 <= K <= N and 0 <= x <= n <= N):
            raise RuntimeError(f"invalid hypergeometric counts for {row['species']}")
        outside += n
        recovered += x
        expected += n * K / N
        distribution = convolve(distribution, hypergeom_pmf(N, K, n))
        informative_species += 1
    upper_tail = sum(distribution[recovered:])
    observed_total = sum(int(r["butterfly_observed_units"]) for r in rows)
    outside_contemporary = sum(int(r["observed_outside_contemporary"]) for r in rows)
    return {
        "species": len(rows),
        "species_with_outside_native_occurrence": informative_species,
        "observed_species_x_unit_total": observed_total,
        "observed_outside_native_species_x_units": outside,
        "recovered_by_introduced_host_ranges": recovered,
        "fraction_outside_native_recovered": (
            None if outside == 0 else recovered / outside
        ),
        "random_geographic_overlap_expectation": expected,
        "observed_minus_random_expectation": recovered - expected,
        "conditional_combined_hypergeometric_upper_tail_p": upper_tail,
        "observed_outside_contemporary_species_x_units": outside_contemporary,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlap-json", type=Path, required=True)
    ap.add_argument("--quality-gate-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    overlap = json.loads(args.overlap_json.read_text(encoding="utf-8"))
    gate = json.loads(args.quality_gate_json.read_text(encoding="utf-8"))
    rows = {str(r["species"]): r for r in overlap["species"]}
    qualified = set(map(str, gate["qualified_species"]))

    payload = {
        "schema": "chocho_butterfly_occurrence_validation_v0.1",
        "status": "COMPLETE",
        "total_wgsrpd3_units": TOTAL_WGSRPD3_UNITS,
        "full_independent_panel": summarize([rows[n] for n in sorted(rows)]),
        "quality_qualified_panel": summarize(
            [rows[n] for n in sorted(qualified) if n in rows]
        ),
        "null_definition": (
            "Condition on each species' number of observed butterfly WGSRPD3 units "
            "outside its native host envelope. Under the diagnostic null, those units "
            "are sampled uniformly without replacement from all WGSRPD3 units outside "
            "that native envelope. The number falling in introduced-added host units "
            "is hypergeometric. Species-specific distributions are convolved exactly."
        ),
        "claim_boundary": {
            "diagnostic_geographic_null_not_dispersal_model": True,
            "occurrence_overlap_does_not_confirm_local_host_use": True,
            "gbif_nonobservation_is_not_true_absence": True,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
