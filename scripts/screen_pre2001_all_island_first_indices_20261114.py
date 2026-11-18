#!/usr/bin/env python3
"""Historical pre-2001 GBIF check for ALL 2011-22 first-index events (225 x 51).

The 2001-2010 GBIF index baseline is *not* evidence of historical absence.
Every observed post-2001 first-index butterfly-island cell is screened using
the SAME archived GSHHG geometry and 1700-2000 GBIF year predicate.
Unknown pre-2001 indexing NEVER implies a genuinely new arrival.

Usage:
  python scripts/screen_pre2001_all_island_first_indices_20261114.py \
    --gshhg-zip ... --work ... --islands 9,19,... --batch-id 1 --out-dir ...
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from build_public51_butterfly_risk_20261109 import GBIF, GSHHG_SHA256, get_json, load_polygons

ROOT = Path(__file__).resolve().parents[1]
BUTTERFLIES = ROOT / "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
PERIODS = {
    "2011_2017": "first_indexed_in_2011_2017",
    "2018_2022": "first_indexed_in_2018_2022",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csvfile(path: Path):
    with path.open(newline="", encoding="utf8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        raise ValueError("Expected nonempty case list")
    with path.open("w", newline="", encoding="utf8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gshhg-zip", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--islands", required=True)
    parser.add_argument("--batch-id", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--max-queries", type=int, default=120)
    parser.add_argument("--wkt-max-chars", type=int, default=2300)
    return parser.parse_args()


def main():
    args = parse_args()
    original = csvfile(BUTTERFLIES)
    if len(original) != 225*51:
        raise ValueError("Original source 225 x 51 not found")
    if sum(int(x["first_indexed_in_2011_2017"]) for x in original) != 221 or \
       sum(int(x["first_indexed_in_2018_2022"]) for x in original) != 96:
        raise ValueError("Frozen 221 + 96 first-index source counts changed")

    ids = set(args.islands.split(","))
    original_islands = {x["landmass_id"] for x in original}
    if not ids or not ids <= original_islands:
        raise ValueError("Bad batch island list")
    selected = [
        x for x in original if x["landmass_id"] in ids
        for period, key in PERIODS.items() if int(x[key]) == 1
    ]
    if not selected:
        raise ValueError("Batch contains no prospective first-index cases")
    by_island = defaultdict(set)
    for x in selected:
        by_island[x["landmass_id"]].add(x["speciesKey"])

    geometry, area_audits, source_digest = load_polygons(
        args.gshhg_zip, args.work, set(by_island), args.wkt_max_chars
    )
    if source_digest != GSHHG_SHA256:
        raise ValueError("Source polygon checksum changed")
    total = sum(len(x) for x in geometry.values())
    if total > args.max_queries:
        raise ValueError(
            f"Batch would require {total} queries vs fixed maximum {args.max_queries}"
        )
    print("PRE2001_FIXED_ALL_FIRST_INDEX_BATCH", args.batch_id,
          "species_island_period_cases", len(selected),
          "polygon_fragment_queries", total, flush=True)

    found = defaultdict(int)
    fragment_failures = []
    receipts = []
    for island in sorted(by_island, key=int):
        keys = sorted(by_island[island], key=int)
        for frag_index, polygon in enumerate(geometry[island], 1):
            params = dict(
                geometry=polygon,
                year="1700,2000",
                taxonKey=keys,
                hasCoordinate="true",
                hasGeospatialIssue="false",
                occurrenceStatus="PRESENT",
                facet="speciesKey", facetLimit=1000, limit=0,
            )
            receipt = dict(
                island=island, fragment=frag_index,
                speciesKeys=len(keys),
                years="1700-2000",
                geometry_WKT_sha256=hashlib.sha256(polygon.encode()).hexdigest(),
                state="UNAVAILABLE",
            )
            try:
                response = get_json(GBIF + "/occurrence/search", params)
                n = int(response["count"])
                facets = {
                    str(f.get("field", "")).upper(): f.get("counts", [])
                    for f in response.get("facets", [])
                }
                if n > 0 and "SPECIES_KEY" not in facets:
                    raise ValueError("Nonzero count with missing GBIF SPECIES_KEY facet")
                counts = defaultdict(int)
                for item in facets.get("SPECIES_KEY", []):
                    k = str(item["name"])
                    if k not in by_island[island]:
                        # This can reveal an unsupported GBIF facet interpretation.
                        raise ValueError("GBIF returned speciesKey not in frozen requested set: "+k)
                    counts[k] += int(item["count"])
                if sum(counts.values()) > n:
                    raise ValueError("Species facet counts exceed total indexed records")
                for k, v in counts.items():
                    found[island, k] += v
                receipt.update({
                    "state":"COMPLETE",
                    "total_indexed":n,
                    "index_count_selected":int(sum(counts.values())),
                    "source_response_sha256":hashlib.sha256(json.dumps(
                        response, sort_keys=True, ensure_ascii=False
                    ).encode()).hexdigest(),
                })
            except Exception as exc:
                receipt.update(state="FAILED_UNKNOWN", error=repr(exc))
                fragment_failures.append(receipt)
            receipts.append(receipt)
            print("HIST_SHARD",args.batch_id,island,frag_index,
                  receipt["state"],receipt.get("index_count_selected"),flush=True)
            time.sleep(0.08)

    island_complete = {
        island:all(r["state"]=="COMPLETE" for r in receipts if r["island"]==island)
        for island in by_island
    }
    output = []
    for x in selected:
        island, speciesKey = x["landmass_id"], x["speciesKey"]
        this_period = next(
            period for period, key in PERIODS.items() if int(x[key]) == 1
        )
        complete = island_complete[island]
        count = found[island, speciesKey] if complete else None
        output.append({
            "species":x["species"],
            "speciesKey":speciesKey,
            "island":island,
            "apparent_first_index_period":this_period,
            "gbif_first_index_baseline":"2001-2010",
            "pre2001_indexed_records":count if complete else "",
            "historical_index_category":(
                "PRE2001_POSITIVE" if count and complete else
                "PRE2001_NONINDEX_UNKNOWN_OCCUPANCY" if complete else
                "INCOMPLETE_SOURCE_UNKNOWN"
            ),
            "not_actual_establishment":True,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / f"source_cases_batch_{args.batch_id}.csv", output)
    audit = {
        "status": "COMPLETE" if not fragment_failures else "PARTIAL_QUERY_FAILURES",
        "schema":"pre2001_earlier_butterfly_index_check_not_true_absence",
        "query_time_utc":datetime.now(timezone.utc).isoformat(),
        "batch":args.batch_id,
        "batch_island_ids":sorted(ids, key=int),
        "with_earlier_index_events":len(output),
        "pre2001_positive_count":sum(r["historical_index_category"]=="PRE2001_POSITIVE" for r in output),
        "source_butterflies_sha256":digest(BUTTERFLIES),
        "GSHHG_archive_sha256":source_digest,
        "geometry":area_audits,
        "query_count":len(receipts),
        "query_receipts":receipts,
        "query_failures":fragment_failures,
        "limits":[
            "GBIF 1700–2000 pre-index positive disproves inference of first presence during 2011–2022, not true historical population establishment",
            "GBIF zero indexes are not ecological absence or evidence of a newly colonizing butterfly",
            "The butterfly species-by-island first-index cells were already selected by 2001–2010 GBIF nonindex",
            "Original source taxonKey and GSHHG 2.3.7 landmass polygon geometry are unchanged",
            "Published historical inventories may extend beyond GBIF and should also be considered",
        ],
    }
    (args.out_dir / f"source_receipts_batch_{args.batch_id}.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf8"
    )
    print("COMPLETE_BATCH_SUMMARY",json.dumps({
        "batch":args.batch_id, "queried_fragments":len(receipts),
        "source_failures":len(fragment_failures), "original_future_events":len(output),
        "pre2001_indexed":audit["pre2001_positive_count"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
