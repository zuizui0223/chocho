#!/usr/bin/env python3
"""Source-only, outcome-blind HOSTS Location inventory (frozen 2026-10-27).

Read the exact historical HOSTS CSV without deleting location information.
No WCVP classification or new-ecological-result claim is made here.
Never upload/redistribute the original HOSTS source in this workflow.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

REQUIRED = {
    "Location", "Insect Genus", "Insect Species",
    "Hostplant Genus", "Hostplant Species",
}

def original_source_inventory(path: Path):
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    location = Counter()
    examples = Counter()
    stats = Counter()
    unmatched = Counter()
    header = None
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if not header or not REQUIRED <= set(header):
            raise ValueError("Pinned source HOSTS required fields are absent")
        for row in reader:
            stats["raw_records"] += 1
            val = (row.get("Location") or "").strip()
            insect = [(row.get(col) or "").strip() for col in
                      ("Insect Genus", "Insect Species")]
            plant = [(row.get(col) or "").strip() for col in
                     ("Hostplant Genus", "Hostplant Species")]
            if all(insect) and all(plant):
                stats["raw_species_level_insect_plant_rows"] += 1
                location[val or "<EMPTY>"] += 1
                if val:
                    stats["species_rows_with_location"] += 1
                    if ";" in val or "/" in val or "|" in val:
                        stats["location_with_explicit_separator"] += 1
                    if val.count(",") > 1:
                        stats["location_with_two_or_more_commas"] += 1
                else:
                    stats["species_rows_location_missing"] += 1
            else:
                stats["rows_missing_species_level_pair"] += 1
                if val:
                    unmatched[val] += 1
            if val and not all(insect):
                examples[val] += 1
    return {
        "schema": "chocho_source_hosts_location_inventory_20261027_v01",
        "source": {
            "original_upstream_repo": "globalbioticinteractions/HOSTS",
            "pinned_commit": "808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "raw_resource_sha256": sha,
            "raw_csv_bytes": path.stat().st_size,
            "columns": header,
        },
        "status": "SOURCE_INVENTORY_ONLY_NO_BIOLOGICAL_CLASSIFICATION",
        "record_flow": dict(stats),
        "distinct_nonempty_location_labels": len(location) - int("<EMPTY>" in location),
        "top_160_location_labels_species_rows": [
            {"label": lab, "rows": count}
            for lab, count in location.most_common(160)
        ],
        "top_20_locations_for_non_species_insect_rows": [
            {"label": lab, "rows": count}
            for lab, count in unmatched.most_common(20)
        ],
        "source_warnings": [
            "The original HOSTS location identifies where a source report states host use, not where the butterfly naturally originated.",
            "A Location string may denote country, region, broader biogeographic realm, province, or literature context, not a coordinate.",
            "Do not collapse multiple HOSTS records into a unique ecological observation before reviewing geography/reference duplication.",
        ],
        "retrieved_from_original_commit_source": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    obj = original_source_inventory(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "schema": obj["schema"],
        "row_flow": obj["record_flow"],
        "label_count": obj["distinct_nonempty_location_labels"],
        "first_40": obj["top_160_location_labels_species_rows"][:40],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
