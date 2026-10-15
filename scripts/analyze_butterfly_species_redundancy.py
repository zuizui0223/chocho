#!/usr/bin/env python3
"""Descriptive butterfly-level redundancy audit from byte-frozen figure sources.

Post-hoc exploratory analysis. The records measure potential host-region
opportunity, not local larval use, demographic establishment or resilience.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def average_ranks(values: list[float] | np.ndarray) -> np.ndarray:
    """Deterministic midranks, including ties; no SciPy dependency."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all():
        raise ValueError("ranks require a finite 1-D array")
    order = np.argsort(x, kind="stable")
    ranked = np.empty(len(x), dtype=float)
    j = 0
    while j < len(x):
        k = j + 1
        while k < len(x) and x[order[k]] == x[order[j]]:
            k += 1
        ranked[order[j:k]] = 0.5 * (j + 1 + k)
        j = k
    return ranked


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) != len(y) or len(x) < 4:
        raise ValueError("inconsistent or too-small vectors")
    a = x - x.mean()
    b = y - y.mean()
    scale = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / scale) if scale > 1e-12 else float("nan")


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return correlation(average_ranks(x), average_ranks(y))


def partial_rank_correlation(
    x: np.ndarray, y: np.ndarray, controls: tuple[np.ndarray, ...]
) -> float:
    """Correlate rank residuals after linearly projecting out control ranks.

    This is an exploratory linear partial-rank association, NOT proof of
    conditional independence or a prospective confirmatory test.
    """
    if not controls:
        return spearman(x, y)
    z = np.column_stack(
        [np.ones(len(x))] + [average_ranks(v) for v in controls]
    )
    ranked_x = average_ranks(x)
    ranked_y = average_ranks(y)
    adjusted_x = ranked_x - z @ np.linalg.lstsq(z, ranked_x, rcond=None)[0]
    adjusted_y = ranked_y - z @ np.linalg.lstsq(z, ranked_y, rcond=None)[0]
    return correlation(adjusted_x, adjusted_y)


def quantile(x: list[float], p: float) -> float:
    return float(np.quantile(np.asarray(x, dtype=float), p))


def read_frozen_rows(path: Path) -> tuple[list[dict[str, object]], str]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open(newline="", encoding="utf-8") as f:
        raw_rows = list(csv.DictReader(f))
    if not raw_rows or len({r["species"] for r in raw_rows}) != len(raw_rows):
        raise ValueError("empty or duplicate species in input")
    result = []
    for row in raw_rows:
        added = int(row["introduced_added_units"])
        host_species = int(row["resolved_host_species"])
        host_families = float(row["host_family_count"])
        multi_raw = row["fraction_added_units_with_multiple_contributing_hosts"]
        if added < 0 or host_species < 0 or host_families < 0:
            raise ValueError("invalid counts")
        if added:
            if multi_raw == "":
                raise ValueError("added unit lacks contributor redundancy metric")
            multi = float(multi_raw)
            if not 0 <= multi <= 1:
                raise ValueError("redundancy proportion out of range")
        else:
            if multi_raw != "":
                raise ValueError("zero-added species must have no redundancy proportion")
            multi = None
        result.append(
            {
                "species": row["species"],
                "host_species": host_species,
                "host_families": host_families,
                "native": int(row["native_resource_units"]),
                "added": added,
                "multi": multi,
                "adequate": row["host_taxonomy_lower_bound_adequate"] == "True",
            }
        )
    return result, digest


def stratum(row: dict[str, object]) -> str:
    n = int(row["host_species"])
    if n == 1:
        return "1"
    if 2 <= n <= 4:
        return "2-4"
    if 5 <= n <= 9:
        return "5-9"
    return "10+"


def describe(rows: list[dict[str, object]]) -> dict[str, object]:
    expanded = [r for r in rows if int(r["added"]) > 0]
    if not expanded:
        raise ValueError("no expanded species")
    total_added = sum(int(r["added"]) for r in expanded)
    weighted_single = sum(
        float(1 - float(r["multi"])) * int(r["added"]) for r in expanded
    )
    group = {}
    for g in ("1", "2-4", "5-9", "10+"):
        xs = [r for r in expanded if stratum(r) == g]
        if not xs:
            raise ValueError("missing host-richness stratum")
        single = [1 - float(r["multi"]) for r in xs]
        group[g] = {
            "butterflies": len(xs),
            "median_single_supported_fraction": round(quantile(single, 0.5), 6),
            "majority_single_supported_butterflies": sum(f > 0.5 for f in single),
            "fully_single_supported_butterflies": sum(f == 1 for f in single),
        }

    def vector(key: str) -> np.ndarray:
        return np.array([float(r[key]) for r in expanded], dtype=float)

    fam = vector("host_families")
    host = vector("host_species")
    native = vector("native")
    added = vector("added")
    multi = vector("multi")
    associations = {
        "spearman_families_multi_support": spearman(fam, multi),
        "spearman_host_species_multi_support": spearman(host, multi),
        "partial_families_given_host_species": partial_rank_correlation(
            fam, multi, (host,)
        ),
        "partial_families_given_host_species_native_and_added_area":
            partial_rank_correlation(fam, multi, (host, native, added)),
        "partial_host_species_given_families": partial_rank_correlation(
            host, multi, (fam,)
        ),
    }

    richness = [
        r for r in expanded
        if int(r["host_species"]) >= 10 and int(r["added"]) >= 20
    ]
    rich_adequate = [r for r in richness if bool(r["adequate"])]
    examples = sorted(
        richness,
        key=lambda r: (
            -(1 - float(r["multi"])),
            -int(r["added"]),
            str(r["species"]),
        ),
    )[:8]
    return {
        "schema": "chocho_species_redundancy_exploratory_v0.1",
        "status": "POSTHOC_DESCRIPTIVE_AUDIT_NOT_CONFIRMATORY",
        "n_focal": len(rows),
        "n_expanded": len(expanded),
        "added_butterfly_region_units": total_added,
        "weighted_single_host_supported_units": int(round(weighted_single)),
        "weighted_single_host_supported_fraction": round(
            weighted_single / total_added, 6
        ),
        "species_median_single_host_supported_fraction": round(
            quantile([1 - float(r["multi"]) for r in expanded], 0.5), 6
        ),
        "strata": group,
        "rich_portfolio_probe": {
            "rule": ">=10 resolved host species and >=20 added resource regions",
            "butterflies": len(richness),
            "majority_single_supported": sum(
                1 - float(r["multi"]) > 0.5 for r in richness
            ),
            "adequate_butterflies": len(rich_adequate),
            "adequate_majority_single_supported": sum(
                1 - float(r["multi"]) > 0.5 for r in rich_adequate
            ),
            "examples": [
                {
                    "species": r["species"],
                    "recorded_hosts": r["host_species"],
                    "added_regions": r["added"],
                    "single_fraction": round(1 - float(r["multi"]), 6),
                }
                for r in examples
            ],
        },
        "rank_associations": {
            key: round(val, 6) for key, val in associations.items()
        },
        "interpretation": (
            "Many butterflies with rich recorded host portfolios still have "
            "little within-region redundancy in newly added host-resource "
            "opportunity. Host-family breadth contributes little descriptive "
            "rank association beyond resolved host-species richness. Both "
            "relationships are structurally constrained by the host portfolio."
        ),
        "boundaries": [
            "Percentages are structural regional resource opportunity, not actual plant use.",
            "This exploratory analysis was chosen after viewing resource results.",
            "Species are not phylogenetically independent.",
            "No null for per-region host-distribution architecture has yet been run.",
            "Unrecorded host links, geography and host phenology may change redundancy.",
            "Zero-added species have undefined redundancy, not zero redundancy.",
            "Partial rank residualization is descriptive, not an equivalence test.",
        ],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--species-csv", type=Path,
        default=Path("data/frozen/figure_sources/host_contribution_metrics.csv"),
    )
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()
    rows, digest = read_frozen_rows(args.species_csv)
    result = describe(rows)
    result["input_sha256"] = digest
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
