#!/usr/bin/env python3
"""Losslessly combine exactly four SHARDED transports of the prior partial 96.

No new butterfly records are inferred. Only fully retrieved species count as
observed, retaining the ORIGINAL GBIF ordinal window selection algorithm and
two temporal periods. Missing/partial taxa never become pseudo-absences.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SHARDS=((0,24),(24,48),(48,72),(72,96))
FIELDS=("species","window","wgsrpd3","basis","sampled_record_count")
SCHEMA="chocho_GBIF_96_partial_species_transportreread_20261019_v01"


def read_names(path):
    with Path(path).open(encoding="utf-8",newline="") as fh:
        rows=list(csv.DictReader(fh))
    assert len(rows)==96
    result={}
    for r in rows:
        sp=r["species"]
        k=int(r["accepted_gbif_species_key"])
        if sp in result or k<=0:
            raise ValueError("Wrong original taxon list or duplicate species")
        result[sp]=k
    return list(result.items())


def merge(shard_directory, frozen_taxa):
    expected=read_names(frozen_taxa)
    found={}
    totals=Counter()
    rows_by_shard=[]
    for shard,(start,stop) in enumerate(SHARDS):
        folder=Path(shard_directory)/f"gbif-96-shard-{shard}"
        receiptfile=folder/"gbif_96_recovery_receipts.json"
        countfile=folder/"gbif_96_recovered_region_basis_counts.csv"
        if not receiptfile.is_file() or not countfile.is_file():
            raise FileNotFoundError(f"Shard {shard} missing full metadata or count CSV")
        meta=json.loads(receiptfile.read_text(encoding="utf-8"))
        assert meta.get("schema")==SCHEMA,shard
        assert meta["species_attempted"]==(stop-start),shard
        assert (meta["shard_start_index"],meta["shard_end_exclusive"])==(start,stop)
        assert meta["original_ordinal_size"]==300
        assert meta["transport_only_subrequest_size"]==75
        assert meta["retained_original_source_sampling_protocol"] is True
        receipts=meta["species_receipts"]
        assert len(receipts)==stop-start
        scoped=dict(expected[start:stop])
        actual={r["species"]:r for r in receipts}
        if set(actual)!=set(scoped):
            raise ValueError(f"Shard {shard} contains non-frozen species or misses expected species")
        for sp,r in actual.items():
            if int(r["accepted_species_key"])!=scoped[sp]:
                raise ValueError(f"GBIF taxon key drift in {sp}")
            if sp in found:raise ValueError(f"Duplicate GBIF species across shards {sp}")
            if r["status"]=="COMPLETE":
                windows=r["windows"]
                if set(windows)!={"2010_2017","2018_2025"} or not all(z.get("complete") for z in windows.values()):
                    raise ValueError(f"False fully sampled GBIF species {sp}")
            elif r["status"]!="PARTIAL_TRANSPORT_NOT_ABSENCE":
                raise ValueError(f"Unknown source-missingness status for {sp}")
            found[sp]=r
        with countfile.open(encoding="utf-8",newline="") as fh:
            reader=csv.DictReader(fh)
            if tuple(reader.fieldnames or ())!=FIELDS:
                raise ValueError("GBIF source region count schema drift")
            local=Counter()
            seen=set()
            for r in reader:
                key=tuple(r[x] for x in FIELDS[:-1])
                n=int(r["sampled_record_count"])
                if n<=0 or key in seen or r["species"] not in actual or actual[r["species"]]["status"]!="COMPLETE":
                    raise ValueError("Duplicate/nonpositive/partial-species occurrence count")
                if r["window"] not in ("2010_2017","2018_2025"):
                    raise ValueError("GBIF period drift")
                local[key]+=n
                seen.add(key)
                if key in totals:raise ValueError("GBIF record overlap across distinct taxa")
                totals[key]=n
        if sum(local.values())!=meta["mapped_complete_records"]:
            raise ValueError(f"Shard mapped GBIF source count drift at {shard}")
        if len(local)!=meta["regional_basis_count_rows"]:
            raise ValueError(f"Shard regional count rows drift {shard}")
        if (sum(r["status"]=="COMPLETE" for r in actual.values())!=meta["complete_species"]
            or sum(r["status"]!="COMPLETE" for r in actual.values())!=meta["not_completed_species"]):
            raise ValueError(f"Shard status totals drift {shard}")
        rows_by_shard.append({
            "shard":shard,"start":start,"end":stop,
            "completed":meta["complete_species"],
            "incomplete":meta["not_completed_species"],
            "mapped_region_records":meta["mapped_complete_records"]
        })
    if len(found)!=96:
        raise ValueError("Not all 96 frozen source taxa were independently attempted")
    status=Counter(x["status"] for x in found.values())
    if status["COMPLETE"]+status["PARTIAL_TRANSPORT_NOT_ABSENCE"]!=96:
        raise ValueError("Missing original taxonomy status")
    result={
        "schema":SCHEMA,
        "source_original_239":39233336804,
        "source_transport":"4 fixed disjoint shards, 24 source taxa each, 75-record transport subrequests for original 300-record ordinal windows",
        "start_utc":None,
        "end_utc":datetime.now(timezone.utc).isoformat(),
        "species_attempted":96,
        "complete_species":status["COMPLETE"],
        "not_completed_species":status["PARTIAL_TRANSPORT_NOT_ABSENCE"],
        "regional_basis_count_rows":len(totals),
        "mapped_complete_records":sum(totals.values()),
        "original_ordinal_pages_per_window":2,
        "original_ordinal_size":300,
        "transport_only_subrequest_size":75,
        "retained_original_source_sampling_protocol":True,
        "species_receipts":[found[sp] for sp,_ in expected],
        "shard_source_audit":rows_by_shard,
        "partial_taxa_are_not_absent":True,
        "no_original_135_records_refetched":True,
        "sampling_index_live_GBIF_drift_possible":True,
    }
    return totals,result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shards-directory",type=Path,required=True)
    p.add_argument("--frozen-taxonomy",type=Path,required=True)
    p.add_argument("--output-counts",type=Path,required=True)
    p.add_argument("--output-manifest",type=Path,required=True)
    args=p.parse_args()
    rows,meta=merge(args.shards_directory,args.frozen_taxonomy)
    args.output_counts.parent.mkdir(parents=True,exist_ok=True)
    args.output_manifest.parent.mkdir(parents=True,exist_ok=True)
    with args.output_counts.open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f);writer.writerow(FIELDS)
        for key,n in sorted(rows.items()):writer.writerow([*key,n])
    args.output_manifest.write_text(json.dumps(meta,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "attempted":meta["species_attempted"],
        "completed":meta["complete_species"],
        "incomplete":meta["not_completed_species"],
        "source_mapped_records":meta["mapped_complete_records"],
        "shards":meta["shard_source_audit"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
