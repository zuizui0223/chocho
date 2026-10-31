#!/usr/bin/env python3
"""Species-only support overlap: dated harbor host links x original-coordinate remote-gap summaries.

This is an OBSERVABILITY audit, NOT an island-year joint test, resource effect,
ship pathway, butterfly colonization, establishment, or genuine relay analysis.
All inputs are tracked in PR #104 and pinned by Git blob SHA.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = "provenance/reviewer_defenses/results/"
LINKS = BASE + "island_dated_host_match_20261030/island_dated_host_match_links.csv"
EVENTS = BASE + "island_host_continuity_20261030/island_host_continuity_events.csv"
SPECIES = BASE + "original_island_continuity_20261030/global_original_continuity_species_counts.csv"
PINS = {
    LINKS: "88db773c74dea9cd2d58bde0b67a394d93226869",
    EVENTS: "efbc770469a720e2dea82478db2f6196ba03b353",
    SPECIES: "17bf86ad6bd3622d36cea738a363f3e76774cf04",
}

def locked_rows(root: Path, relpath: str):
    raw = (root / relpath).read_bytes()
    sha = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\x00" + raw).hexdigest()
    if sha != PINS[relpath]:
        raise ValueError(f"Source Git blob SHA changed: {relpath}: {sha}")
    return list(csv.DictReader(raw.decode("utf-8").splitlines()))

def audit(root: Path = ROOT):
    links = locked_rows(root, LINKS)
    if len(links) != 5:
        raise ValueError("Expected fixed five prior-host butterfly links")
    names = sorted({r["species"] for r in links})
    if len(names) != 5:
        raise ValueError("Expected five distinct butterfly taxa")
    keyset = {name: set() for name in names}
    for row in locked_rows(root, EVENTS):
        if row["species"] in keyset:
            keyset[row["species"]].add(row["speciesKey"])
    if any(len(v) != 1 for v in keyset.values()):
        raise ValueError("Unresolved or nonunique butterfly GBIF species keys")
    counts = locked_rows(root, SPECIES)
    if len(counts) != 2368:
        raise ValueError("Global original-coordinate species-count panel changed")
    remote = [r for r in counts if r["cohort"] == "both_period_opportunity"
              and r["distance_band"] == "500 km以上"]
    if len(remote) != 193:
        raise ValueError("Opportunity-qualified cohort size changed")
    bykey = {}
    for r in remote:
        key = r["speciesKey"]
        if key in bykey:
            raise ValueError("Duplicate species remote-distance summary")
        n = int(r["events"])
        if not 0 <= n <= int(r["total_events"]):
            raise ValueError("Invalid species remote-distance count")
        bykey[key] = n
    series = sorted((n for n in bykey.values() if n > 0), reverse=True)
    total = sum(series)
    top10 = sum(series[:10])
    if (total, len(series), top10) != (164, 86, 52):
        raise ValueError("Locked original remote-gap cohort changed")
    selected = [
        {"butterfly": name, "speciesKey": next(iter(keyset[name])),
         "remote_gap_500km_both_period_events_across_all_islands":
             bykey.get(next(iter(keyset[name])))}
        for name in names
    ]
    if any(row["remote_gap_500km_both_period_events_across_all_islands"] is None
           for row in selected):
        raise ValueError("Selected butterfly species absent from global distance cohort")
    overlap = sum(r["remote_gap_500km_both_period_events_across_all_islands"] > 0
                  for r in selected)
    nselected = sum(r["remote_gap_500km_both_period_events_across_all_islands"]
                    for r in selected)
    return {
        "schema": "chocho_remote_gap_prior_host_species_overlap_20261030_v1",
        "input_git_blob_pins": PINS,
        "remote_gap_cohort": {
            "name": "both_period_opportunity",
            "n_species": len(remote),
            "total_events_500km": total,
            "species_with_at_least_one_remote_gap": len(series),
            "maximum_remote_gaps_one_species": series[0],
            "top_ten_species_remote_events": top10,
            "top_ten_share": round(top10 / total, 6),
            "simpson_effective_species": round(total * total / sum(v*v for v in series), 6),
        },
        "selected_early_host_species": selected,
        "species_with_any_remote_gap": overlap,
        "remote_gap_events_in_these_selected_species_any_island": nselected,
        "joint_island_species_year_host_and_gap_event_status": "NOT_OBSERVABLE_FROM_TRACKED_AGGREGATES",
        "source_island_status_for_remote_species_counts": "NOT_CONTAINED_IN_SPECIES_COUNT_TABLE",
        "no_establishment_or_shipping_relay_inference": True,
    }

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = audit(args.root)
    content = json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content, encoding="utf-8")
        print(args.output)
    else:
        print(content, end="")

if __name__ == "__main__":
    main()
