#!/usr/bin/env python3
"""Source-coverage audit: author-selected migrant candidates and immature records.

No movement mechanism is inferred from review-derived author-selected taxa. The
record-status UNKNOWN group is NOT the "non-migratory" group. Stage observations
are a limited archived relay sample, NOT a global immature-record search.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

MIGRANT = "data/exploratory/public51_author_migrant_20261109/migration_gbif_author2025_crosswalk.csv"
STAGES = "results/island_breeding_relay/verified_stage_records.csv"
BUTTERFLIES = "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
HOSTS = "results/route_free_port_hosts/species_island_hosts.csv"
PLANT_EFFORT = "results/route_free_port_hosts/port_coverage.csv"
SOURCE_CODE = "REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"


def read_rows(root: Path, relative: str):
    path = root / relative
    b = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f)), hashlib.sha256(b).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--out-dir", type=Path, default=Path("results/port_host_reindex_20261114"))
    a = p.parse_args()
    tables = {}
    hashes = {}
    for k, path in (("butterflies", BUTTERFLIES), ("hosts", HOSTS),
                    ("plant_effort", PLANT_EFFORT), ("migrant", MIGRANT),
                    ("stages", STAGES)):
        tables[k], hashes[k] = read_rows(a.repo_root, path)
    mig = {x["speciesKey"]: x["author_2025_GBIF_requested_taxon_key_proxy"]
           for x in tables["migrant"]}
    host = {(x["key"], x["island"]): x for x in tables["hosts"]}
    effort = defaultdict(int)
    for x in tables["plant_effort"]:
        effort[x["island"]] += int(x["records"])
    stage = defaultdict(list)
    for x in tables["stages"]:
        stage[x["key"], x["island"]].append({
            "gbif_id": x["gbif_id"], "year": x["year"],
            "larva": x["larva"], "egg": x["egg"], "laying": x["laying"]
        })

    cases = []
    for x in tables["butterflies"]:
        if int(x["eligible_2011_2017"]) != 1 or int(x["first_indexed_in_2011_2017"]) != 1:
            continue
        k = x["speciesKey"], x["landmass_id"]
        if k not in host or effort[k[1]] < 10:
            continue
        status = mig.get(k[0], "NOT_LISTED_OR_TAXON_UNRESOLVED_UNKNOWN")
        cases.append({
            "species": x["species"], "speciesKey": k[0], "landmass_id": k[1],
            "host_indexed": int(host[k]["introduced"]) > 0,
            "candidate_migrant_author_selected": status == SOURCE_CODE,
            "other_mobility_status": "UNKNOWN" if status != SOURCE_CODE else "NOT_APPLICABLE",
            "later_reindexed": int(x["future2018_2022_indexed_records"]) > 0,
            "verified_immature_stage_in_selected_archived_relay_source": bool(stage.get(k)),
            "archived_stage_rows": len(stage[k]),
            "early_record_count": int(x["future2011_2017_indexed_records"]),
            "later_record_count": int(x["future2018_2022_indexed_records"]),
        })

    def group(items):
        return {
            "n": len(items),
            "author_selected_migrant_candidate": sum(
                x["candidate_migrant_author_selected"] for x in items
            ),
            "mobility_not_in_positive_list_unknown": sum(
                not x["candidate_migrant_author_selected"] for x in items
            ),
            "later_reindexed": sum(x["later_reindexed"] for x in items),
            "archived_immature_stage_match": sum(
                x["verified_immature_stage_in_selected_archived_relay_source"] for x in items
            ),
            "distinct_species": len({x["speciesKey"] for x in items}),
            "distinct_islands": len({x["landmass_id"] for x in items}),
        }

    present = [x for x in cases if x["host_indexed"]]
    unindexed = [x for x in cases if not x["host_indexed"]]
    assert len(cases) == 159 and len(present) == 15 and len(unindexed) == 144
    out = {
        "analysis": "port_host_and_author_selected_migrant_coverage",
        "source_sha256": hashes,
        "main": {"historic_introduced_host_indexed": group(present),
                 "historic_introduced_host_not_indexed_unknown_presence": group(unindexed)},
        "case_listing_introduced_host_indexed": present,
        "interpretation": [
            "Author-selected review-related positive list is indirect, NOT an independent measured migration phenotype",
            "Unlisted 225-list candidates have UNKNOWN dispersal status, not non-migratory",
            "Stage-record source was assembled for a distinct relay case study and is NOT systematic across 51 islands",
            "Stage nonoverlap does not imply larvae, eggs or breeding truly absent",
            "GBIF re-index is NOT established occupancy, population growth or annual persistence",
            "Host candidate statuses are present-day regional labels assigned to 2001-2010 port indices",
            "Do not use known-list share as a causal host x migratory ability interaction",
        ],
    }
    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / "candidate_migration_and_stage_source_audit.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"main": out["main"], "source_count": len(cases)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
