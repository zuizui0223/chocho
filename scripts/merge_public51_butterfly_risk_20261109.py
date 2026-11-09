#!/usr/bin/env python3
"""Merge six independently reindexed 225-butterfly x 51-island GBIF panels.

This is a new 2026 GBIF INDEX snapshot, not the author's private 83,573-row
original file, not biological absence, dispersal, reproduction, or migration.
No unknown original first-index identities are imputed.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_ISLANDS=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv"
ORIGINAL_TAXA=ROOT/"provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv"
SOURCE_ISLAND_GIT_BLOB="18cecde5bf3066f276380167f451f1fa5cc1caa1"
SOURCE_TAXA_GIT_BLOB="39fe84981cab2b0e47063399f0a5c0c9731e800d"
FIELDS=[
  "species","speciesKey","landmass_id","pre2001_2010_indexed_records",
  "future2011_2017_indexed_records","future2018_2022_indexed_records",
  "eligible_2011_2017","first_indexed_in_2011_2017",
  "eligible_2018_2022","first_indexed_in_2018_2022",
  "source_not_future_selected","index_not_true_absence_or_colonization"
]
OLD17=ROOT/"data/exploratory/ncomms_public17_access_20261108/full17_butterfly_island_225_atrisk_source.csv"
OLD17_BLOB="60e38e31d1d1b8896fb925aaa5ee0a405fbef2df"

def blob(p):
    data=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def csv_rows(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        fields=list(reader.fieldnames or [])
        return fields,list(reader)

def inspect_original():
    if blob(ORIGINAL_ISLANDS)!=SOURCE_ISLAND_GIT_BLOB or blob(ORIGINAL_TAXA)!=SOURCE_TAXA_GIT_BLOB:
        raise ValueError("Original island/taxon population source changed")
    _,islands=csv_rows(ORIGINAL_ISLANDS)
    _,species=csv_rows(ORIGINAL_TAXA)
    eligible={x["landmass_id"] for x in islands if x["eligible"]=="True"}
    taxa={x["taxon_key"]:x["species"] for x in species}
    if len(eligible)!=51 or len(taxa)!=225:
        raise ValueError("Original butterfly/island cohort changed")
    return eligible,taxa

def merge(part_dir,output,expect_parts=6):
    expected_islands,taxa=inspect_original()
    files=sorted(part_dir.rglob("batch_*_butterfly_island_225_atrisk.csv"))
    receipts=sorted(part_dir.rglob("batch_*_query_receipt.json"))
    if len(files)!=expect_parts or len(receipts)!=expect_parts:
        raise ValueError(f"Expected {expect_parts} complete source CSV/JSON shard pairs; got {len(files)}, {len(receipts)}")
    by_prefix={p.name.split("_butterfly_island_")[0]:p for p in files}
    receipt_prefix={p.name.split("_query_receipt")[0]:p for p in receipts}
    if set(by_prefix)!=set(receipt_prefix) or len(by_prefix)!=expect_parts:
        raise ValueError("Source batch filenames do not pair 1:1")
    seen=set()
    by_island=Counter()
    rows=[]
    audits=[]
    butterfly_filters=[]
    index_times=[]
    for prefix in sorted(by_prefix):
        p=by_prefix[prefix]
        src=json.loads(receipt_prefix[prefix].read_text(encoding="utf8"))
        fields,rr=csv_rows(p)
        if fields!=FIELDS:raise ValueError("Unexpected public original source schema")
        if src.get("schema")!="chocho_public_GBIF_unselected_225_butterfly_x_original51_island_window_v1":
            raise ValueError("Not the versioned 51-island public reindex source")
        if src.get("status")!="PUBLIC_REINDEXED_OCCURRENCE_HISTORY_NOT_TRUE_COLONIZATION_OR_HOST_MEDIATION":
            raise ValueError("GBIF source incompletely verified")
        if src.get("gshhg_shape_sha256")!="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41":
            raise ValueError("GSHHG source changed")
        pin=src.get("source_git_blob_pins",{})
        if pin.get("butterflies",{}).get("sha")!=SOURCE_TAXA_GIT_BLOB or pin.get("islands",{}).get("sha")!=SOURCE_ISLAND_GIT_BLOB:
            raise ValueError("Original source population hashes changed")
        ids={str(i) for i in src["islands_tested"]}
        if not ids<=expected_islands or not ids or len(rr)!=225*len(ids):
            raise ValueError("Invalid island shard original sample")
        if len(src.get("GBIF_index_query_receipts",[]))!=len(ids)*3 or len(src.get("geometry_audit",[]))!=len(ids):
            raise ValueError("Missing geometry/query receipt per island and period")
        if set(z["island_id"] for z in src["geometry_audit"])!=ids:
            raise ValueError("Geometry receipt does not cover source islands")
        butterfly_filters.append(src["butterfly_taxon_filter"].get("family_taxonkeys"))
        index_times.append(src.get("snapshot_timestamp_utc"))
        summaries=Counter()
        for x in rr:
            key=(x["speciesKey"],x["landmass_id"])
            if key in seen or x["landmass_id"] not in ids or x["speciesKey"] not in taxa or taxa[x["speciesKey"]]!=x["species"]:
                raise ValueError("Duplicate, foreign, or taxonomically shifted source cell")
            seen.add(key);by_island[x["landmass_id"]]+=1
            old,mid,new=[int(x[k]) for k in ("pre2001_2010_indexed_records",
                    "future2011_2017_indexed_records","future2018_2022_indexed_records")]
            if any(v<0 for v in (old,mid,new)):
                raise ValueError("Negative GBIF query result")
            binary={"eligible_2011_2017":int(old==0),
                "first_indexed_in_2011_2017":int(old==0 and mid>0),
                "eligible_2018_2022":int(old==0 and mid==0),
                "first_indexed_in_2018_2022":int(old==0 and mid==0 and new>0)}
            for k,v in binary.items():
                if int(x[k])!=v:raise ValueError("Invalid first-indexing risk rule")
                summaries[k]+=v
            if x["source_not_future_selected"]!="True" or x["index_not_true_absence_or_colonization"]!="True":
                raise ValueError("Future-selected or mislabeled historical response")
            rows.append(x)
        if len(rr)!=src["risk_cells"] or any(summaries[k]!=src[k] for k in binary):
            raise ValueError("Shard receipt and row counts disagree")
        audits.append({"prefix":prefix,"ids":sorted(ids,key=int),"risk_cells":len(rr),
            "csv_sha256":sha(p),"receipt_sha256":sha(receipt_prefix[prefix]),
            "index_timestamp_utc":src.get("snapshot_timestamp_utc"),
            "indexed_new_2011_2017":summaries["first_indexed_in_2011_2017"],
            "indexed_new_2018_2022":summaries["first_indexed_in_2018_2022"]})
    if not all(t==butterfly_filters[0] for t in butterfly_filters):
        raise ValueError("GBIF taxon filters changed across six source query jobs")
    if len(rows)!=11475 or len(seen)!=11475 or len(by_island)!=51 or any(v!=225 for v in by_island.values()):
        raise ValueError("Public 51-island source full population not reconstructed")
    rows.sort(key=lambda x:(int(x["landmass_id"]),int(x["speciesKey"])))
    totals={key:sum(int(x[key]) for x in rows) for key in FIELDS[6:10]}
    # Old independent 17-island result is a snapshot from October 31 and can
    # have drifted by November 9. Quantify disagreement; never equate sources.
    drift=None
    if OLD17.exists() and blob(OLD17)==OLD17_BLOB:
        _,old_rows=csv_rows(OLD17)
        old={(r["speciesKey"],r["landmass_id"]):r for r in old_rows}
        new={(r["speciesKey"],r["landmass_id"]):r for r in rows}
        if not set(old)<=set(new):raise ValueError("Old source's 17 islands not included in 51")
        cols=["pre2001_2010_indexed_records","future2011_2017_indexed_records",
              "future2018_2022_indexed_records"]
        drift={"old17_source_snapshot":"2026-10-31",
               "new51_source_snapshot":"2026-11-09",
               "comparable_original_species_island_cells":len(old),
               "changed_record_counts_by_period":{c:sum(x[c]!=new[k][c] for k,x in old.items()) for c in cols},
               "changed_binary_first_indexing":{c:sum(x[c]!=new[k][c] for k,x in old.items())
                   for c in ("first_indexed_in_2011_2017","first_indexed_in_2018_2022")},
               "note":"GBIF index snapshot differences are not biological trends"}
    output.mkdir(parents=True,exist_ok=True)
    table=output/"full51_225butterfly_period_index.csv"
    with table.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    receipt={
      "schema":"chocho_publicGBIF_225butterfly_original51_source_merged_20261109_v1",
      "scope":"NEW independent 2026 GBIF reindex, NOT the original private 83573-row parquet",
      "source_islands":51,"source_butterflies":225,"cells":11475,
      "source_island_git_blob":SOURCE_ISLAND_GIT_BLOB,
      "source_butterflies_git_blob":SOURCE_TAXA_GIT_BLOB,
      "source_GSHHG_sha256":"8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41",
      "island_counts":dict(sorted(by_island.items(),key=lambda z:int(z[0]))),
      "totals":totals,"shard_audits":audits,
      "query_snapshot_timestamps":index_times,
      "six_shards_not_simultaneous_snapshot":True,
      "original17_reindexed_snapshot_drift_audit":drift,
      "source_index_events_not_biological_arrival_or_colonization":True,
      "not_measured":["larval host plant introductions","shipping","actual butterfly movement",
                      "reproduction","true biological absence"],
      "merged_table_sha256":sha(table),
      "written_at_utc":datetime.now(timezone.utc).isoformat()
    }
    (output/"full51_source_validation.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"cells":11475,"totals":totals,"old17_drift":drift},indent=2,ensure_ascii=False))
    return receipt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--parts",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    ap.add_argument("--expected-shards",type=int,default=6)
    args=ap.parse_args()
    merge(args.parts,args.out_dir,args.expected_shards)

if __name__=="__main__":main()
