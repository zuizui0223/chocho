#!/usr/bin/env python3
"""Combine all 11,475 time-consistent butterfly island historical GBIF indices.

A later 'first GBIF record after no records in 1700–2010' is a source index
event, *not* verified colonization, breeding, or arrival. All 11,475 cells are
covered by the exact same pre-2001 extraction, without positive-outcome
conditioning. Source and geometry hashes are checked before merging.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
HISTORY_SHA256="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"


def read_csv(path):
    with path.open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sumcounts(rows,col):
    return sum(int(r[col]) for r in rows)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shards",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    source=read_csv(SOURCE)
    if len(source)!=11475:
        raise ValueError("Original 225 x 51 source not fixed")
    expected={(r["speciesKey"],r["landmass_id"]) for r in source}
    if len(expected)!=11475:raise ValueError("Original source contains duplicate key")

    out=[]
    shardreceipts=[]
    for n in range(1,7):
        cp=a.shards/f"full_history_batch_{n}.csv"
        jp=a.shards/f"full_history_receipt_batch_{n}.json"
        if not cp.exists() or not jp.exists():
            raise ValueError(f"Missing source batch {n}")
        data=json.loads(jp.read_text(encoding="utf8"))
        if data["source_status"]!="COMPLETE" or data["source_fragment_failures"]:
            raise ValueError(f"Incomplete source records in batch {n}")
        if data["GSHHG_sha256"]!=HISTORY_SHA256 or data["source_11_475_sha256"]!=sha(SOURCE):
            raise ValueError("Source geometry or 225 x 51 frozen bytes differ")
        rows=read_csv(cp)
        if len(rows)!=data["original_source_case_count"] or len(rows)%225:
            raise ValueError("Unexpected source shard cells")
        if sumcounts(rows,"original_2011_2017_first_index")!=data["original_early_first_count"] \
           or sumcounts(rows,"original_2018_2022_first_index")!=data["original_late_first_count"]:
            raise ValueError("Original butterfly denominators altered")
        if sumcounts(rows,"pre2001_any_index") != data["pre2001_positive_all_cells"]:
            raise ValueError("Pre2001 GBIF positive tally altered")
        out.extend(rows)
        shardreceipts.append({
            "batch":n, "cells":len(rows),
            "queries":data["source_fragment_queries"],
            "positive_pre2001":data["pre2001_positive_all_cells"],
            "csv_sha256":sha(cp),"source_receipt_sha256":sha(jp),
        })
    got={(r["speciesKey"],r["landmass_id"]) for r in out}
    if len(out)!=11475 or got!=expected:
        raise ValueError(f"All cells must be present once: observed {len(out)} unique {len(got)}")
    if any(r["pre2001_source_complete"]!="1" for r in out):
        raise ValueError("Incomplete historical source remained in full joined cells")
    if sumcounts(out,"original_2011_2017_first_index")!=221 \
       or sumcounts(out,"original_2018_2022_first_index")!=96:
        raise ValueError("Original 317 response window totals changed")

    stats={
        "all_11475":{
            "cells":len(out),
            "source_pre2001_indexed_cells":sumcounts(out,"pre2001_any_index"),
            "source_pre2001_no_index_TRUE_OCCUPANCY_UNKNOWN":
                len(out)-sumcounts(out,"pre2001_any_index"),
        }
    }
    for name,prefix in [("2011_2017","2011_2017"),("2018_2022","2018_2022")]:
        orig_key=f"original_{prefix}_first_index"
        new_key=f"indexed_first_{prefix}_longer_history"
        eligible_key=f"eligible_{prefix}_longer_history"
        subset=[x for x in out if x[orig_key]=="1"]
        old=sum(x["pre2001_any_index"]=="1" for x in subset)
        eligible=sumcounts(out,eligible_key)
        new=sumcounts(out,new_key)
        stats[name]={
            "original_2001_2010_baseline_first_indices":len(subset),
            "historical_pre2001_index_among_original_events":old,
            "historic_source_not_indexed_in_old_event_cells":len(subset)-old,
            "more_consistent_1700_2010_index_only_eligible_cells":eligible,
            "more_consistent_1700_2010_index_only_first_indices":new,
            "percent_original_first_indices_with_pre2001_index":100*old/len(subset),
            "confirmed_establishment_count":None,
            "source_first_index_not_colonization":True,
        }
        if len(subset)-old != new:
            raise ValueError("Long-history source eligibility and prior positive counts mismatch")

    # A contemporaneous source retrieval may drift. Consistency with the earlier
    # 317-event historical-only screen is required, or report drift explicitly.
    drift = {
        "prior_positive_only_screen_2011_2017":122,
        "new_full_screen_2011_2017":stats["2011_2017"]["historical_pre2001_index_among_original_events"],
        "prior_positive_only_screen_2018_2022":23,
        "new_full_screen_2018_2022":stats["2018_2022"]["historical_pre2001_index_among_original_events"],
    }
    stats["source_drift_comparison"]=drift
    stats["same_day_history_source_agrees_with_prior_317_screen"] = (
        drift["prior_positive_only_screen_2011_2017"]==drift["new_full_screen_2011_2017"]
        and drift["prior_positive_only_screen_2018_2022"]==drift["new_full_screen_2018_2022"]
    )

    out.sort(key=lambda r:(int(r["landmass_id"]),r["species"]))
    a.out_dir.mkdir(parents=True,exist_ok=True)
    cp=a.out_dir/"full_11475_same_source_pre2001_risk_cells.csv"
    with cp.open("w",newline="",encoding="utf8") as file:
        w=csv.DictWriter(file,fieldnames=list(out[0]))
        w.writeheader();w.writerows(out)
    report={
        "status":"SOURCE_COMPLETE",
        "original_panel_sha256":sha(SOURCE),
        "pre2001_years":"1700-2000",
        "source_window_2001_2010":"frozen original 2026-11-09 queried public GBIF",
        "GSHHG_original_sha256":HISTORY_SHA256,
        "six_batch_provenance":shardreceipts,
        "uniform_same_source_risk_screen":True,
        "outcomes_are_record_indices_not_actual_immigrations":True,
        "stats":stats,
        "limits":[
            "Historical 1700–2000 GBIF nonindex is not true absence",
            "GBIF indexing changes over time, six retrieval shards are not simultaneous",
            "The 51 landmasses are selected original ports with butterfly observations",
            "A source-consistent first index is not independent historical colonization or verified breeding",
            "A historical index does not document an exporting butterfly population for a later destination",
            "Regional contemporary WCVP introduced-origin labels are not local historic host introduction dates",
        ],
    }
    jp=a.out_dir/"full_risk_source_summary.json"
    jp.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("FULL_11475_HISTORY_REBUILD",json.dumps({
        "complete":True, "cases":len(out),
        "shards":len(shardreceipts),
        "query_fragments":sum(x["queries"] for x in shardreceipts),
        "statistics":stats,
        "output":str(cp),
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
