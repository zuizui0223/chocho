#!/usr/bin/env python3
"""Describe observed-record support of the fixed 24-island pilot, without fits."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/exploratory/island_port_20261029"


def run():
    path = OUT / "island_butterfly_pilot_counts.csv"
    with path.open() as handle:
        records = list(csv.DictReader(handle))
    islands = defaultdict(list)
    for row in records:
        islands[(row["island_key"], row["stratum"])].append(row)
    if len(islands) != 24 or len(records) != 48:
        raise ValueError("Frozen pilot denominator changed")
    blocks = {}
    for (_, stratum), windows in islands.items():
        if len(windows) != 2 or {r["window"] for r in windows} != {"2010_2017", "2018_2025"}:
            raise ValueError("Duplicated or missing pilot island window")
        area, gap, port = stratum.split(":")
        block = blocks.setdefault((area, gap), {
            "area_class_km2": area, "continental_gap_class_km": gap,
            "near50_selected": 0, "far100_selected": 0,
            "near50_both_windows_ge3_species": 0, "far100_both_windows_ge3_species": 0,
            "near50_failed_or_missing_windows": 0, "far100_failed_or_missing_windows": 0})
        label = "near50" if port == "within50" else "far100"
        if port not in {"within50", "beyond100"}:
            raise ValueError("Unexpected source port class")
        block[label+"_selected"] += 1
        complete = all(r["status"] == "COMPLETE" for r in windows)
        block[label+"_failed_or_missing_windows"] += sum(r["status"] != "COMPLETE" for r in windows)
        block[label+"_both_windows_ge3_species"] += int(complete and all(
            int(r["recorded_species_bins"]) >= 3 for r in windows))
    rows = []
    for _, row in sorted(blocks.items()):
        row["both_port_classes_have_paired_record_support"] = int(
            row["near50_both_windows_ge3_species"] > 0 and row["far100_both_windows_ge3_species"] > 0)
        rows.append(row)
    with (OUT / "island_butterfly_pilot_design_support.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {"schema": "chocho_island_pilot_record_support_v1", "islands": len(islands),
        "source_count_table_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "geographic_blocks": len(rows),
        "blocks_with_both_port_classes_and_paired_ge3_species": sum(r["both_port_classes_have_paired_record_support"] for r in rows),
        "islands_with_paired_ge3_species": sum(r["near50_both_windows_ge3_species"]+r["far100_both_windows_ge3_species"] for r in rows),
        "definition": ">=3 recorded species bins in each period is a descriptive broad-target-group coverage screen. It is not survey effort, and is not the original focal-species >=3 OTHER taxa gate.",
        "interpretation": "Geographic balance at selection does not guarantee observed-record support. Do not impute zero butterfly abundance to empty record windows or claim port/host effects from this pilot.",
        "blocks": rows}
    (OUT / "island_butterfly_pilot_design_support.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "blocks"}, sort_keys=True))


if __name__ == "__main__":
    run()
