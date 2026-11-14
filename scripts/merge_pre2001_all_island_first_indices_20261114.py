#!/usr/bin/env python3
"""Merge six original GSHHG pre-2001 historical-index shards, without false absence.

This is an observation-index audit of positive outcomes, not a matched
establishment probability study. Do not refit a 'corrected' colonization model:
all remaining source-negative risk cells would first need the same historical
screen, and 1700–2000 nonindex still does not prove no butterfly.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ORIGINAL = "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
HOST = "results/route_free_port_hosts/species_island_hosts.csv"
PORTS = "results/route_free_port_hosts/port_coverage.csv"


def read_csv(path):
    with path.open(newline="", encoding="utf8") as file:
        return list(csv.DictReader(file))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragments", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    raw = read_csv(REPO / ORIGINAL)
    if len(raw) != 11475:
        raise ValueError("Original fixed 225 x 51 cohort changed")
    expect = {(r["speciesKey"],r["landmass_id"],p)
        for r in raw for p,key in (
          ("2011_2017","first_indexed_in_2011_2017"),
          ("2018_2022","first_indexed_in_2018_2022"))
        if int(r[key]) == 1}
    if len(expect) != 317:
        raise ValueError(f"Unexpected original first-index set {len(expect)} vs 317")

    joined,shards=[],[]
    for k in range(1,7):
        cp = args.fragments / f"source_cases_batch_{k}.csv"
        jp = args.fragments / f"source_receipts_batch_{k}.json"
        if not cp.exists() or not jp.exists():
            raise FileNotFoundError(f"Missing source shard for batch {k}")
        sr=json.loads(jp.read_text(encoding="utf8"))
        if sr["status"] != "COMPLETE" or sr["query_failures"]:
            raise ValueError(f"Incomplete historical source shard {k}: {sr['status']}")
        if sr["source_butterflies_sha256"] != sha256(REPO / ORIGINAL):
            raise ValueError(f"Original 225 x 51 source drift in shard {k}")
        rows = read_csv(cp)
        if len(rows) != sr["with_earlier_index_events"]:
            raise ValueError(f"Shard file/receipt record count mismatch {k}")
        if sum(r["historical_index_category"]=="PRE2001_POSITIVE" for r in rows) != sr["pre2001_positive_count"]:
            raise ValueError(f"Source facet positive mismatch batch {k}")
        if sr["query_count"] != len(sr["query_receipts"]):
            raise ValueError("Source query receipt count mismatch")
        joined.extend(rows)
        shards.append({
            "id":k,
            "cases":len(rows),
            "positive_pre2001":sr["pre2001_positive_count"],
            "queries":sr["query_count"],
            "case_csv_sha256":sha256(cp),
            "source_receipt_sha256":sha256(jp),
        })

    actual = [(r["speciesKey"],r["island"],r["apparent_first_index_period"])
              for r in joined]
    if set(actual) != expect or len(actual)!=len(set(actual)):
        raise ValueError(f"6 shard total inconsistent: actual {len(actual)}, expected {len(expect)}")

    h = {(r["key"],r["island"]):r for r in read_csv(REPO / HOST)}
    effort=defaultdict(int)
    for row in read_csv(REPO/PORTS):
        effort[row["island"]]+=int(row["records"])
    introduced_focal = []
    for row in joined:
        key, isl = row["speciesKey"], row["island"]
        host=h.get((key,isl))
        row["historic_port_flora_indexed"] = bool(effort[isl] > 0)
        row["historic_port_plant_observation_count"] = effort[isl]
        row["candidate_intro_host_historic_port_indexed"] = bool(host and int(host["introduced"])>0)
        row["native_host_historic_port_indexed"] = bool(host and int(host["native"])>0)
        row["modern_origin_not_verified_2001_local"] = True
        row["pre2001_association_is_index_not_establishment"] = True
        if row["apparent_first_index_period"]=="2011_2017" and \
           effort[isl]>=10 and row["candidate_intro_host_historic_port_indexed"]:
            introduced_focal.append(row)
    if len(introduced_focal)!=15:
        raise ValueError(f"Focal 15 historical-host cases no longer matched: {len(introduced_focal)}")
    focal_positive = sum(x["historical_index_category"]=="PRE2001_POSITIVE" for x in introduced_focal)
    if focal_positive!=12:
        raise ValueError(f"Original 15 historical screen mismatch expected 12, got {focal_positive}")

    stats={}
    def summarize(rows):
        count=Counter(r["historical_index_category"] for r in rows)
        return {
            "total":len(rows),
            "pre2001_positive":count["PRE2001_POSITIVE"],
            "historically_unindexed_TRUE_PRESENCE_UNKNOWN":count["PRE2001_NONINDEX_UNKNOWN_OCCUPANCY"],
            "incomplete":count["INCOMPLETE_SOURCE_UNKNOWN"],
            "historical_record_fraction":count["PRE2001_POSITIVE"]/len(rows) if rows else None,
            "distinct_species":len({r["speciesKey"] for r in rows}),
            "islands":len({r["island"] for r in rows}),
        }
    stats["all_first_indexes"]=summarize(joined)
    for era in ["2011_2017","2018_2022"]:
        this=[r for r in joined if r["apparent_first_index_period"]==era]
        stats[era]=summarize(this)
        stats[f"{era}_candidate_intro_host_indexed"]=summarize([r for r in this if r["candidate_intro_host_historic_port_indexed"]])
    stats["focal_15_host_positive_early"]=summarize(introduced_focal)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    joined.sort(key=lambda z:(z["apparent_first_index_period"],int(z["island"]),z["species"]))
    out_csv=args.out_dir/"all_317_butterfly_island_historic_index_audit.csv"
    with out_csv.open("w", newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(joined[0]))
        w.writeheader()
        w.writerows(joined)
    report={
        "analysis":"225 butterfly x original51 island post-2001 event history audit",
        "source_6_shards":shards,
        "original_panel_sha256":sha256(REPO/ORIGINAL),
        "historic_plant_crosswalk_sha256":sha256(REPO/HOST),
        "port_effort_sha256":sha256(REPO/PORTS),
        "complete_source_status":True,
        "established_not_inferred":True,
        "statistics":stats,
        "historical_data_scope":"1700-2000 GBIF georeferenced PRESENT/no-geospatial-issue indexed records",
        "critical_finding":"Earlier source-positive GBIF records preclude interpreting 2011-2017 or 2018-2022 as true first island arrival",
        "no_earlier_record_is_not_absence":True,
        "why_not_corrected_establishment_model":[
            "Only the positive 2011-2022 first-index outcomes, not all unindexed risk cells, were retro-screened",
            "The 2001-2010 index baseline still selects which post-2001 outcomes count as 'first'",
            "Pre2001 GBIF nonindex does not prove true historical species absence",
            "All the original first-index outcome data are observation indices, not breeding or colonization",
            "The full comparison groups have non-identical data capture and some island origins are only contemporary regional status"
        ],
        "all_15_original_intro_cases":[
            {k:r[k] for k in ("species","speciesKey","island","apparent_first_index_period",
                              "pre2001_indexed_records","historical_index_category")}
            for r in introduced_focal
        ],
    }
    (args.out_dir/"analysis_summary.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("COMBINED_PRE2001_AUDIT",json.dumps({
        "source_status":True, "shards":6,"original_positive_cases":317,
        "stats":stats,"out_csv":str(out_csv)
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
