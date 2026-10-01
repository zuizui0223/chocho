#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def load_level1(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for feature in payload.get("features", []):
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        region = str(props.get("LEVEL1_COD") or "").strip()
        if code and region:
            out[code] = region
    if len(out) < 300:
        raise RuntimeError("WGSRPD3 lookup drift")
    return out


def summarize_draws(draws, observed):
    draws = np.asarray(draws, dtype=float)
    return {
        "mean": float(np.mean(draws)),
        "q025": float(np.quantile(draws, 0.025)),
        "median": float(np.median(draws)),
        "q975": float(np.quantile(draws, 0.975)),
        "p_high": float((1 + np.sum(draws >= observed)) / (len(draws) + 1)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native-footprints-json", type=Path, required=True)
    ap.add_argument("--contemporary-footprints-json", type=Path, required=True)
    ap.add_argument("--unit-table-csv", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--null-replicates", type=int, default=99999)
    ap.add_argument("--bootstrap-replicates", type=int, default=50000)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    native_payload = json.loads(args.native_footprints_json.read_text(encoding="utf-8"))
    contemporary_payload = json.loads(
        args.contemporary_footprints_json.read_text(encoding="utf-8")
    )
    native = {sp: set(map(str, units)) for sp, units in native_payload["species"].items()}
    contemporary = {
        sp: set(map(str, units))
        for sp, units in contemporary_payload["species"].items()
    }
    level1 = load_level1(args.level3_geojson)
    universe = set(level1)

    observed_units = defaultdict(set)
    with args.unit_table_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["butterfly_observed"]):
                observed_units[row["species"]].add(row["wgsrpd3_code"])

    species_rows = []
    for species in sorted(contemporary):
        native_units = native.get(species, set()) & universe
        contemporary_units = contemporary[species] & universe
        added = contemporary_units - native_units
        outside = (observed_units.get(species, set()) & universe) - native_units
        if not outside:
            continue
        overlap = added & outside

        region_terms = []
        expected = 0.0
        for region in sorted(set(level1.values())):
            candidates = {
                code
                for code, value in level1.items()
                if value == region and code not in native_units
            }
            n = sum(level1[code] == region for code in added)
            K = sum(level1[code] == region for code in outside)
            N = len(candidates)
            if n:
                if n > N:
                    raise RuntimeError("added units exceed regional candidate pool")
                region_terms.append((N, K, n))
                expected += 0.0 if N == 0 else K * n / N

        species_rows.append(
            {
                "species": species,
                "outside_units": len(outside),
                "recovered_units": len(overlap),
                "fraction_recovered": len(overlap) / len(outside),
                "region_matched_expected_units": expected,
                "region_terms": region_terms,
            }
        )

    if len(species_rows) < 10:
        raise RuntimeError("too few informative species")

    rng = np.random.default_rng(20260928)
    B = args.null_replicates
    S = len(species_rows)
    simulated = np.zeros((B, S), dtype=np.int16)
    for j, row in enumerate(species_rows):
        for N, K, n in row["region_terms"]:
            simulated[:, j] += rng.hypergeometric(K, N - K, n, size=B).astype(np.int16)

    observed_recovered = np.asarray([r["recovered_units"] for r in species_rows], dtype=float)
    outside = np.asarray([r["outside_units"] for r in species_rows], dtype=float)
    observed_fraction = observed_recovered / outside
    simulated_fraction = simulated / outside[np.newaxis, :]

    observed_total = int(np.sum(observed_recovered))
    outside_total = int(np.sum(outside))
    sim_total = np.sum(simulated, axis=1)

    observed_mean_species_fraction = float(np.mean(observed_fraction))
    observed_median_species_fraction = float(np.median(observed_fraction))
    sim_mean_species_fraction = np.mean(simulated_fraction, axis=1)
    sim_median_species_fraction = np.median(simulated_fraction, axis=1)

    # Cluster bootstrap: resample species, retaining each species' numerator and denominator.
    Bboot = args.bootstrap_replicates
    boot_ratio = np.empty(Bboot, dtype=float)
    boot_mean_fraction = np.empty(Bboot, dtype=float)
    for i in range(Bboot):
        idx = rng.integers(0, S, size=S)
        boot_ratio[i] = np.sum(observed_recovered[idx]) / np.sum(outside[idx])
        boot_mean_fraction[i] = np.mean(observed_fraction[idx])

    leave_one_out = []
    for j, row in enumerate(species_rows):
        loo_observed = observed_total - int(observed_recovered[j])
        loo_outside = outside_total - int(outside[j])
        loo_draws = sim_total - simulated[:, j]
        leave_one_out.append(
            {
                "excluded_species": row["species"],
                "observed_recovered_units": loo_observed,
                "outside_native_units": loo_outside,
                "fraction_recovered": loo_observed / loo_outside,
                "region_matched_null": summarize_draws(loo_draws, loo_observed),
            }
        )

    by_name = {row["excluded_species"]: row for row in leave_one_out}
    complete = [
        {
            "species": row["species"],
            "recovered_units": row["recovered_units"],
            "outside_units": row["outside_units"],
        }
        for row in species_rows
        if row["recovered_units"] == row["outside_units"]
    ]

    payload = {
        "schema": "chocho_butterfly_occurrence_species_robustness_v0.1",
        "status": "POSTHOC_SPECIES_CLUSTER_ROBUSTNESS",
        "species": S,
        "observed": {
            "recovered_units": observed_total,
            "outside_native_units": outside_total,
            "pooled_fraction_recovered": observed_total / outside_total,
            "mean_species_fraction_recovered": observed_mean_species_fraction,
            "median_species_fraction_recovered": observed_median_species_fraction,
        },
        "species_level_null": {
            "mean_fraction_recovered": summarize_draws(
                sim_mean_species_fraction, observed_mean_species_fraction
            ),
            "median_fraction_recovered": summarize_draws(
                sim_median_species_fraction, observed_median_species_fraction
            ),
        },
        "species_cluster_bootstrap": {
            "replicates": Bboot,
            "pooled_fraction_recovered_q025": float(np.quantile(boot_ratio, 0.025)),
            "pooled_fraction_recovered_median": float(np.median(boot_ratio)),
            "pooled_fraction_recovered_q975": float(np.quantile(boot_ratio, 0.975)),
            "mean_species_fraction_q025": float(np.quantile(boot_mean_fraction, 0.025)),
            "mean_species_fraction_median": float(np.median(boot_mean_fraction)),
            "mean_species_fraction_q975": float(np.quantile(boot_mean_fraction, 0.975)),
        },
        "leave_one_species_out": leave_one_out,
        "leave_one_out_summary": {
            "minimum_fraction_recovered": min(r["fraction_recovered"] for r in leave_one_out),
            "maximum_fraction_recovered": max(r["fraction_recovered"] for r in leave_one_out),
            "maximum_region_matched_p_high": max(
                r["region_matched_null"]["p_high"] for r in leave_one_out
            ),
            "pyrgus_communis_excluded": by_name.get("Pyrgus communis"),
            "pieris_brassicae_excluded": by_name.get("Pieris brassicae"),
        },
        "complete_recovery_species": complete,
        "largest_single_species_contribution": max(
            species_rows, key=lambda r: r["recovered_units"]
        ),
        "claim_boundary": (
            "Species-cluster resampling treats butterfly species, not species x region units, "
            "as the replication unit. It still does not demonstrate local larval use or "
            "causal host-facilitated range expansion."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
