#!/usr/bin/env python3
"""Audit genuine 7737-row 17th-century Dutch intra-Asian cargo source.

Original DataverseNL 10.34894/5HGTCD file 569978, SHA1
d7cc3aff064fb63f58b8b549b1ca9808cb4b76d1; published 7,737 rows,
5 columns. Compare exact voyage ID keys to the separately SHA-checked
57,669 voyage legs only WHEN original cargo schema exposes the same key.
Do not treat transported commodity groups as living plants, larval hosts,
propagules or successful introductions without additional evidence.
Original ship/voyage cargo rows never leave GitHub runner.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
CARGO=Path("_private_history/voc_7737_cargo.tab")
VOYAGES=Path("_private_history/voc_57669_all_voyages.tab")
OUT=BASE/"voc_17c_7737_cargo_original_source_audit.json"
EXPECTED_SHA1="d7cc3aff064fb63f58b8b549b1ca9808cb4b76d1"
EXPECTED_BYTES=1389412
EXPECTED_ROWS=7737
EXPECTED_COLS=5
VOYAGE_SHA1="eed2479019a8570e392f37f3d781d5b0c67448b1"


def sha(path):
    h=hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path):
    with path.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        fields=list(reader.fieldnames or [])
        yield fields
        for row in reader:
            yield row


def analysis(cargo_path=CARGO,voyage_path=VOYAGES):
    if cargo_path.stat().st_size!=EXPECTED_BYTES or sha(cargo_path)!=EXPECTED_SHA1:
        raise RuntimeError("Original cargo source bytes or SHA1 mismatch")
    iterator=read_rows(cargo_path)
    fields=next(iterator)
    if len(fields)!=EXPECTED_COLS or len(set(fields))!=len(fields):
        raise RuntimeError("Original 5-field cargo schema not obtained")
    counts=Counter()
    per_column={k:Counter() for k in fields}
    nonempty=Counter()
    for row in iterator:
        counts["source_cargo_rows"]+=1
        for k in fields:
            v=str(row.get(k) or "").strip()
            if v:
                nonempty[k]+=1
                per_column[k][v]+=1
    if counts["source_cargo_rows"]!=EXPECTED_ROWS:
        raise RuntimeError("Original cargo source 7,737 records not fully read")
    idlike=[k for k in fields if "reisn" in k.casefold() or
            "voyage" in k.casefold() or "reis_nr" in k.casefold()]
    type_like=[k for k in fields if any(t in k.casefold() for t in
               ("cargo","lading","goods","commodity","product","goederen","waar"))]
    profiles={}
    for k in fields:
        unique=len(per_column[k])
        vals=per_column[k]
        # Only aggregated commodity/category names without person/voyage ID.
        is_safe_category=(k in type_like and unique<=5000)
        if not is_safe_category and not idlike and unique<=500:
            is_safe_category=True
        profile={"n_nonempty":nonempty[k],
                 "n_distinct":unique,
                 "top_count":vals.most_common(1)[0][1] if vals else 0}
        if is_safe_category:
            profile["frequent_category_values"]=[
                {"category":label,"rows":n}
                for label,n in vals.most_common(50) if n>=2 and len(label)<=90
            ]
        profiles[k]=profile
    join={
        "status":"NO_EXACT_VOYAGE_KEY_IN_CARGO_TABLE",
        "cargo_voyage_key_columns":idlike,
    }
    if idlike and voyage_path.exists() and sha(voyage_path)==VOYAGE_SHA1:
        with voyage_path.open("r",encoding="utf-8-sig",newline="") as f:
            vr=csv.DictReader(f,delimiter="\t")
            if "REISN,N,6,0" not in (vr.fieldnames or []):
                raise RuntimeError("Source voyage ID column not verified")
            sourceids={str(r["REISN,N,6,0"]).strip() for r in vr
                       if str(r.get("REISN,N,6,0") or "").strip()}
        if len(idlike)==1:
            cargoids=set(per_column[idlike[0]])
            both=sourceids&cargoids
            join={
                "status":"MATCHED_SOURCE_VOYAGE_ID_VOCABULARY" if both
                         else "CARGO_VOYAGE_IDS_NONOVERLAPPING",
                "source_voyage_ids":len(sourceids),
                "cargo_ids":len(cargoids),
                "intersecting_source_voyage_ids":len(both),
                "cargo_unique_ids_unmatched":len(cargoids)-len(both),
                "source_voyage_ids_without_cargo_id":len(sourceids)-len(both),
                "note":"Vocabulary overlap alone does not prove each cargo row refers to the same individual port-to-port leg",
            }
    return {
        "status":"SOURCE_VERIFIED_VOC_7737_CARGO_RECORDS_SCHEMA_NOT_LARVAL_HOST_TRANSPORT",
        "citation":"Leenstra and Parthesius 2022, Dutch Ships in Asian Waters in the 17th Century, DataverseNL DOI 10.34894/5HGTCD",
        "archive_file":569978,
        "original_cargo_source_sha1":EXPECTED_SHA1,
        "original_cargo_source_bytes":cargo_path.stat().st_size,
        "cargo_source_row_count":counts["source_cargo_rows"],
        "columns":fields,
        "source_record_column_profiles":profiles,
        "candidate_exact_voyage_id_join":join,
        "limitations":[
            "Original shipping commodities and merchandise records are not automatically living plant propagules or compatible butterfly larval host species.",
            "Cargo record taxonomy and voyage key relationship must be validated before cargo assignments to a particular island-to-island ship leg.",
            "Source 7,737 rows represent one historical dataset, not complete cargo manifests for all 57,669 voyages.",
            "No vessel, named ship, individual shipment date/quantity, cargo detail row or voyager is published here.",
            "No association with modern GloNAF/WCVP/GBIF butterfly outcomes is estimated from cargo field names.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cargo",type=Path,default=CARGO)
    p.add_argument("--voyages",type=Path,default=VOYAGES)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    j=analysis(args.cargo,args.voyages)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(j,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(j,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
