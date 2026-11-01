#!/usr/bin/env python3
"""Inspect 57,669 original 17th-century VOC ship voyage source rows (metadata).

Leenstra and Parthesius (2022), IISG Dataverse
hdl:10622/YP2COP, public all voyages.tab (Dataverse fileId 10466).
The source table is NOT modern AIS and NOT proof of host-plant cargo.
No source individual vessel, named ship, individual dated row or original
historical voyage is published to the chocho repository.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import re

BASE=Path(__file__).resolve().parent
SOURCE=Path("_private_history/voc_all_voyages.tab")
OUT=BASE/"voc_17c_57669_source_schema.json"
EXPECTED_MD5="1a4d0867da62c1cb4ac4fad83f78d092"
EXPECTED_ROWS=57669
EXPECTED_COLS=18


def md5(path):
    h=hashlib.md5()
    with Path(path).open("rb") as f:
        for v in iter(lambda:f.read(1<<20),b""):
            h.update(v)
    return h.hexdigest()


def audit(path=SOURCE):
    if not path.exists():
        raise FileNotFoundError(path)
    original_checksum=md5(path)
    # Document exact original bytes, including source encoding, before any
    # downstream data reuse. The IISH Dataverse file browser reports an MD5.
    valid_checksum=(original_checksum==EXPECTED_MD5)
    raw=path.read_bytes()
    mode="utf-8-sig"
    try:
        decoded=raw.decode(mode)
    except UnicodeDecodeError:
        mode="latin-1"
        decoded=raw.decode(mode)
    snip=decoded[:5000]
    try:
        dialect=csv.Sniffer().sniff(snip,delimiters="\t,;|")
    except csv.Error:
        dialect=csv.excel_tab
    rows=list(csv.DictReader(io.StringIO(decoded),dialect=dialect))
    if not rows:
        raise RuntimeError("Original historical voyage file contains no parsed records")
    columns=list(rows[0])
    name_tokens=("from","to","port","loc","place","haven","departure","arrival",
                 "destination","herkomst","bestemming","afvaart","aankomst",
                 "city","harbor","harbour","vertrek","reis","date","year","tijd","voyage")
    candidate=[k for k in columns if any(t in k.casefold() for t in name_tokens)]
    coverage={}
    for key in columns:
        values=[(r.get(key) or "").strip() for r in rows]
        nonblank=sum(bool(x) for x in values)
        distinct=len(set(v for v in values if v))
        entry={"nonblank_rows":nonblank,"distinct_nonblank_values":distinct}
        # Limit published sample values to geography fields or numeric date/
        # year fields. Never publish raw vessel names or rare voyages.
        location_hint=any(t in key.casefold() for t in
             ("port","place","loc","haven","departure","arrival","destination",
              "afvaart","aankomst","herkomst","bestemming","from","to"))
        if location_hint and distinct>1:
            counts=Counter(v for v in values if v)
            # Avoid accidental publication of personal ship/voyage identifiers.
            if not any(s in key.casefold() for s in ("ship","vessel","naam","name","schip")):
                entry["most_frequent_source_place_names"]=[
                  {"name":place,"records":n} for place,n in counts.most_common(20)
                  if n>=2
                ]
        coverage[key]=entry
    result={
        "status":"SOURCE_ARCHIVE_SCHEMA_ONLY_UNTIL_VOYAGE_PLACE_KEYS_VERIFIED",
        "dataset":"Dutch Ships in Asian Waters in the 17th Century",
        "citation":"Leenstra and Parthesius (2022), hdl:10622/YP2COP",
        "download_file_dataverse_id":10466,
        "web_reported_rows":EXPECTED_ROWS,
        "web_reported_variables":EXPECTED_COLS,
        "web_reported_original_md5":EXPECTED_MD5,
        "observed_download_md5":original_checksum,
        "web_reported_original_md5_matches_download":valid_checksum,
        "source_content_bytes":len(raw),
        "decoded_encoding":mode,
        "detected_separator":dialect.delimiter,
        "source_rows":len(rows),
        "source_variables":len(columns),
        "columns":columns,
        "join_candidate_columns":candidate,
        "column_profiles":coverage,
        "limitations":[
            "Historical 17th-century VOC voyages are from one company/era, not global colonial shipping.",
            "At this stage no port-pair links have been matched to original WPI194 ports.",
            "Dataverse original-source MD5 must match before treating contents as the original archive.",
            "Voyages with matching place names do not confirm botanical consignments.",
            "No source expedition row, crew or historical ship identity is committed.",
        ]
    }
    # Do not claim a usable original until source checksum and grain verified.
    if len(rows)==EXPECTED_ROWS and len(columns)==EXPECTED_COLS and valid_checksum:
        result["status"]="SOURCE_VERIFIED_57669_HISTORICAL_VOC_VOYAGE_ROWS_SCHEMA"
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SOURCE)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    d=audit(args.source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in d.items() if k!="column_profiles"},
                     indent=2,ensure_ascii=False))
    if d["status"]!="SOURCE_VERIFIED_57669_HISTORICAL_VOC_VOYAGE_ROWS_SCHEMA":
        print("Source mismatch; archive data must not be entered as empirical historic voyage link evidence.")
        raise SystemExit(2)


if __name__=="__main__":
    main()
