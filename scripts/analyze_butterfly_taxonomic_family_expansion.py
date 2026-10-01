#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import median


def load_taxonomy(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"Species", "Family"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("LepTraits taxonomy columns missing")
        for row in reader:
            species = str(row["Species"]).strip()
            family = str(row["Family"]).strip()
            if species and family:
                out[species] = family
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--anthropogenic-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    taxonomy = load_taxonomy(args.leptraits_csv)
    groups: dict[str, list[dict[str, float | int | str]]] = defaultdict(list)

    with args.anthropogenic_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            species = str(row["species"]).strip()
            family = taxonomy.get(species)
            if not family:
                raise RuntimeError(f"missing butterfly Family for {species}")
            groups[family].append({
                "species": species,
                "host_family_count": float(row["host_family_count"]),
                "native_resource_units": int(row["native_resource_units"]),
                "introduced_added_units": int(row["introduced_added_units"]),
                "log_resource_expansion": float(row["log_resource_expansion"]),
            })

    total = sum(len(v) for v in groups.values())
    if total != 239:
        raise RuntimeError(f"expected 239 butterflies, got {total}")

    families = []
    for family, rows in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        expanded = sum(int(r["introduced_added_units"]) > 0 for r in rows)
        families.append({
            "family": family,
            "species": len(rows),
            "expanded_species": expanded,
            "expanded_fraction": expanded / len(rows),
            "median_host_family_count": median(float(r["host_family_count"]) for r in rows),
            "median_native_resource_units": median(int(r["native_resource_units"]) for r in rows),
            "median_added_units": median(int(r["introduced_added_units"]) for r in rows),
            "median_log_resource_expansion": median(float(r["log_resource_expansion"]) for r in rows),
        })

    major = [x for x in families if x["species"] >= 10]
    payload = {
        "schema": "chocho_butterfly_taxonomic_family_expansion_v0.1",
        "status": "DESCRIPTIVE_BUTTERFLY_FAMILY_GENERALITY_CHECK",
        "resource_species": total,
        "butterfly_families": len(families),
        "major_family_threshold_n": 10,
        "major_families": major,
        "all_families": families,
        "interpretation": (
            "Descriptive taxonomic check of whether anthropogenic resource expansion "
            "is confined to one butterfly family. No family-level causal or comparative "
            "inference is made."
        ),
        "claim_boundary": (
            "Family summaries are descriptive and sample sizes differ strongly among "
            "families. Riodinidae is represented by one species and should not be "
            "interpreted as a family-level estimate."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
