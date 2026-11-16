#!/usr/bin/env python3
"""Same-WGSRPD3 paired EuDiS-source composition: not a dispersal mechanism test."""
from __future__ import annotations
import argparse
import csv
import json
import random
import statistics
from collections import Counter
from pathlib import Path


def load(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def bootstrap_median_ci(values: list[float], n=4999) -> list[float] | None:
    if len(values) < 3:
        return None
    rng = random.Random(20261116)
    k = len(values)
    out = sorted(statistics.median(rng.choices(values, k=k)) for _ in range(n))
    return [out[int((n-1) * .025)], out[int((n-1) * .975)]]


def summarize(pairs: list[dict]) -> dict:
    vals = [p["native_minus_introduced_thalasso_pp"] for p in pairs]
    return {
        "same_regions": len(pairs),
        "median_native_minus_introduced_trait_pp": statistics.median(vals) if vals else None,
        "descriptive_region_bootstrap_95pp": bootstrap_median_ci(vals),
        "positive_deltas": sum(v > 0 for v in vals),
        "negative_deltas": sum(v < 0 for v in vals),
        "ties": sum(v == 0 for v in vals),
        "level1_counts": dict(sorted(Counter(str(p["level1"]) for p in pairs).items())),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rows = load(a.csv)
    by_code: dict[str, dict] = {}
    for r in rows:
        code = r["area_code_l3"]
        if r["introduced"] not in ("0", "1"):
            raise ValueError("Unexpected origin category")
        by_code.setdefault(code, {})[r["introduced"]] = r
    pairs = []
    for code, group in sorted(by_code.items()):
        if not {"0", "1"}.issubset(group):
            continue
        n, i = group["0"], group["1"]
        if n["eligible"].lower() != "true" or i["eligible"].lower() != "true":
            continue
        if n["class"] != i["class"] or n["level1"] != i["level1"]:
            raise ValueError("Source geographic mismatch within pair")
        n_fraction = float(n["thalasso_trait_fraction_among_EuDiS_matched"])
        i_fraction = float(i["thalasso_trait_fraction_among_EuDiS_matched"])
        pairs.append({
            "code": code, "class": n["class"], "level1": n["level1"],
            "native_matched": int(n["EuDiS_matched_species"]),
            "introduced_matched": int(i["EuDiS_matched_species"]),
            "native_trait_1": int(n["thalasso_trait_1"]),
            "introduced_trait_1": int(i["thalasso_trait_1"]),
            "native_fraction": n_fraction, "introduced_fraction": i_fraction,
            "native_minus_introduced_thalasso_pp": 100 * (n_fraction - i_fraction)
        })
    if not pairs:
        raise RuntimeError("No eligible same-region native vs introduced source comparisons")
    with (a.out_dir / "within_region_eudis_native_introduced_pairs.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(pairs[0]))
        writer.writeheader()
        writer.writerows(pairs)
    all_island = [p for p in pairs if p["class"] == "island"]
    mainland = [p for p in pairs if p["class"] == "mainland"]
    eu_island = [p for p in all_island if p["level1"] == "1"]
    outside_eu = [p for p in all_island if p["level1"] != "1"]
    result = {
        "status": "PAIRED_EUDIS_TRAIT_COMPOSITION_NOT_OCEAN_TRANSPORT",
        "threshold": "within identical region, both categories 20+ EuDiS matched species and >=5% WCVP coverage",
        "islands_all": summarize(all_island),
        "islands_TDWG_level1_Europe": summarize(eu_island),
        "islands_outside_Europe": summarize(outside_eu),
        "mainlands_all": summarize(mainland),
        "origin_records_are_region_specific": True,
        "intercomparison_notes": [
            "Native/introduced differences may reflect European-only EuDiS source ascertainment, not the sea-current supply mechanism.",
            "Inferred trait 0 means absent morphological syndrome, not proven absence of saltwater dispersal.",
            "Only same-region origin categories are paired; no coastal or historical ship transport observed.",
            "These geographic units and fractions are not the 51-island butterfly arrival-risk endpoints."
        ]
    }
    (a.out_dir / "paired_source_trait_composition.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("EUDIS_SAME_REGION_NATIVE_INTRODUCED_DIAGNOSTIC", json.dumps(
        result, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
