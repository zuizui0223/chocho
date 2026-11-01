#!/usr/bin/env python3
"""Audit ORIGINAL IISH / Dutch VOC intra-Asian voyage data: 57,669 records.

Leenstra & Parthesius (2022) Dutch Ships in Asian Waters in the 17th Century,
DOI 10.34894/5HGTCD / hdl:10622/YP2COP. DataverseNL version1 source
file id 569977, all voyages.tab, source JSON original SHA1
eed2479019a8570e392f37f3d781d5b0c67448b1, size 10,187,695.
It is distinct from Huygens DAS Europe–Asia 8,178 voyages and 2020 AIS rows.

No historic vessel/captain name or original individual journey is exported.
"""
from __future__ import annotations

from collections import Counter
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
SOURCE=Path("_private_history/voc_intra_asian_57669.tsv")
OUT=BASE/"voc_intra_asian_57669_original_source_schema.json"
EXPECTED_SHA1="eed2479019a8570e392f37f3d781d5b0c67448b1"
EXPECTED_BYTES=10187695
EXPECTED_RECORDS=57669
EXPECTED_COLUMNS=18


def sha1(path):
    return hashlib.sha1(Path(path).read_bytes()).hexdigest()


def read_source(path=SOURCE):
    raw=Path(path).read_bytes()
    try:
        decoded=raw.decode("utf-8-sig")
        encoding="utf-8-sig"
    except UnicodeDecodeError:
        decoded=raw.decode("latin1")
        encoding="latin1"
    # Published original file is tabular, but do not silently misparse HTML,
    # XML or Dataverse consent pages as observations.
    if decoded.lstrip().startswith(("<html","<?xml","<!DOCTYPE","{")):
        raise RuntimeError("Historical source URL returned metadata/access document instead of tabular voyages")
    reader=csv.DictReader(io.StringIO(decoded),delimiter="\t")
    columns=reader.fieldnames or []
    rows=list(reader)
    return rows,columns,encoding


def audit(path=SOURCE):
    checksum=sha1(path)
    size=path.stat().st_size
    rows,cols,encoding=read_source(path)
    if checksum!=EXPECTED_SHA1 or size!=EXPECTED_BYTES:
        raise RuntimeError(f"Published 2022 DataverseNL original file SHA/size mismatch: {checksum}, {size}")
    if len(rows)!=EXPECTED_RECORDS or len(cols)!=EXPECTED_COLUMNS:
        raise RuntimeError("Historic intra-Asian voyage 57669x18 source schema mismatched: "+
                           f"{len(rows)} records, {len(cols)} columns, {cols[:5]}")
    scope={}
    geohints=("place","port","from","to","depart","arrival","destination",
              "haven","herkomst","bestemming","loc","vertrek","aankomst","jaar","year","datum","date")
    for col in cols:
        values=[(r.get(col) or "").strip() for r in rows]
        observed=[x for x in values if x]
        hist=Counter(observed)
        entry={"rows_present":len(observed),"nonempty_distinct_values":len(hist)}
        # Port codes/labels are useful for finding source geography but no
        # individual voyage, ship, captain or crew identifier is published.
        low=col.casefold()
        if any(t in low for t in geohints) and not any(s in low for s in ("ship","vessel","captain","naam","person","schip","voyage_id")):
            entry["most_frequent_historical_place_or_date_values"]=[
                {"value":v,"count":n} for v,n in hist.most_common(25)
                if n>=2
            ]
        scope[col]=entry
    return {
        "schema":"chocho_IISG_VOC_intra_Asia_57669_original_schema_v1",
        "status":"SOURCE_VERIFIED_IISG_HISTORICAL_INTRA_ASIA_VOYAGE_ROWS_NOT_YET_PORT_MAPPED",
        "source_citation":"Leenstra, Menno; Parthesius, Robert (2022), Dutch Ships in Asian Waters in the 17th Century",
        "source_doi":"10.34894/5HGTCD",
        "source_IISG_hdl":"hdl:10622/YP2COP",
        "original_DataverseNL_file_id":569977,
        "expected_source_SHA1":EXPECTED_SHA1,
        "actual_source_SHA1":checksum,
        "original_file_size_bytes":size,
        "original_record_count":len(rows),
        "source_columns":cols,
        "columns_count":len(cols),
        "decoding":encoding,
        "column_profiles":scope,
        "license":"source metadata: CC BY-SA custom source terms; do not distribute original rows",
        "notices":[
            "Original 57,669 intra-Asia voyage records; NOT Huygens 8178 Europe–Asia voyage record set.",
            "Only source provenance/schema audited; no butterfly plant, colonial exposure or historical source–destination connection yet.",
            "Historical place names must be mapped to WPI 2019 original 51 islands with ambiguous names excluded.",
            "Voyage movement between islands, if independently recovered, does not prove host-plant or butterfly cargo.",
        ],
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,default=SOURCE)
    parser.add_argument("--output",type=Path,default=OUT)
    args=parser.parse_args()
    v=audit(args.source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(v,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:q for k,q in v.items() if k!="column_profiles"},indent=2))


if __name__=="__main__":
    main()
