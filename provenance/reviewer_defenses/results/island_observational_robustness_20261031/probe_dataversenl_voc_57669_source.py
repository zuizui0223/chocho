#!/usr/bin/env python3
"""Audit actual 57,669-row VOC intra-Asian voyages from DataverseNL mirror.

Leenstra and Parthesius (2022), hdl:10622/YP2COP,
source mirror doi:10.34894/5HGTCD; file id 569977 ("all voyages.tab").

The original IISG Dataverse API file 10466 returned XML metadata in a
previous source attempt. An XML response is never treated as a voyage.
Validate original published SHA-1, archived file size, exact 57,669 rows and
18 columns before claiming the historical from/to source is usable.

Only summary schema and common place strings are persisted. Individual
voyage data, ship names/IDs, dates, and crew data remain on the CI runner.
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
import zipfile

ROOT=Path(__file__).resolve().parent
LOCK=ROOT/"dataversenl_17c_original_57669_lock_20261101.json"
SOURCE=Path("_private_history/voc_57669_all_voyages.tab")
OUT=ROOT/"dataversenl_17c_57669_voyage_schema_audit.json"
EXPECTED_SIZE=10187695
EXPECTED_SHA1="eed2479019a8570e392f37f3d781d5b0c67448b1"
EXPECTED_ROWS=57669
EXPECTED_COLS=18


def checksum(data,kind):
    return hashlib.new(kind,data).hexdigest()


def parse_bytes(raw):
    if raw.lstrip().startswith((b"<?xml",b"<!DOCTYPE",b"<html",b"<HTML")):
        raise ValueError("Source access returned XML/HTML, not historic voyages")
    wrapped=False
    member=None
    payload=raw
    if raw.startswith(b"PK\x03\x04"):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            matches=[f for f in z.infolist()
                     if not f.is_dir() and
                     f.filename.lower().endswith((".tab",".csv",".tsv",".txt"))]
            if len(matches)!=1:
                raise ValueError("Original API ZIP has no single verifiable voyage table")
            member=matches[0].filename
            payload=z.read(matches[0])
            wrapped=True
    try:
        txt=payload.decode("utf-8-sig")
        encoding="utf-8-sig"
    except UnicodeDecodeError:
        txt=payload.decode("cp1252",errors="replace")
        encoding="cp1252-with-unmapped-byte-replacement"
    head=txt[:10000]
    try:
        dialect=csv.Sniffer().sniff(head,delimiters="\t,;|")
        delimiter=dialect.delimiter
    except csv.Error:
        delimiter="\t" if head.splitlines()[0].count("\t")>=5 else ","
    reader=csv.DictReader(io.StringIO(txt,newline=""),delimiter=delimiter)
    fields=[str(x or "").strip().lstrip("\ufeff") for x in (reader.fieldnames or [])]
    if not fields:
        raise ValueError("No record table column header")
    profiles={k:Counter() for k in fields}
    present=Counter()
    examples={}
    rows=0
    for row in reader:
        rows+=1
        for key in fields:
            val=(row.get(key) or "").strip() if isinstance(row.get(key),str) else ""
            if not val:
                continue
            present[key]+=1
            # The number of values is bounded; only precise geography
            # candidate metadata are released from the archival records.
            if len(profiles[key])<=350000:
                profiles[key][val]+=1
    columns={}
    place_tokens=("origin","destination","departure","arrival","from","to","port",
                  "harbour","harbor","haven","place","location","city","plaats",
                  "vertrek","aankomst","bestemming","afvaart","herkomst")
    private_tokens=("vessel","captain","crew","shipname","ship_name",
                    "person","mmsi","voyage_id","name_ship","schip","record_id")
    temporal_tokens=("date","year","datum","tijd","arrival","departure")
    for key in fields:
        k=key.casefold().replace(" ","_")
        vals=profiles[key]
        item={"nonblank_rows":present[key],
              "distinct_values":len(vals),
              "top_value_count":vals.most_common(1)[0][1] if vals else 0}
        if (any(t in k for t in place_tokens) and
            not any(t in k for t in private_tokens) and
            not any(t in k for t in ("date","time","datum","year","tijd")) and
            len(vals)<=5000):
            # Only repeated historical port names, not voyage identities.
            item["most_common_source_place_labels"]=[
                {"name":x,"source_rows":n}
                for x,n in vals.most_common(30) if n>=3 and len(x)<100
            ]
        if any(t in k for t in temporal_tokens):
            years=Counter()
            for val,n in vals.items():
                matches=re.findall(r"\b(1[56]\d\d|17\d\d|18\d\d)\b",val)
                if matches and len(matches)==1:
                    years[matches[0]]+=n
            item["source_years_identifiable"]=sum(years.values())
            item["source_year_min"]=min(years,default=None)
            item["source_year_max"]=max(years,default=None)
        columns[key]=item
    return {"source_rows":rows,"columns":fields,"column_count":len(fields),
            "separator":repr(delimiter),"encoding":encoding,
            "api_format_was_zip":wrapped,"archived_zip_member":member,
            "payload_bytes":len(payload),
            "per_column":columns}


def audit(source=SOURCE):
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock.get("schema")!="chocho_dataversenl_17c_asiatic_historic_voyages_source_gate_v1":
        raise RuntimeError("Original historical source lock missing")
    raw=source.read_bytes()
    receipt={
      "source_DataverseNL_DOI":"10.34894/5HGTCD",
      "historical_original_citation":"Leenstra and Parthesius (2022), hdl:10622/YP2COP",
      "source_file_ID":569977,
      "download_bytes":len(raw),
      "published_original_bytes":EXPECTED_SIZE,
      "download_SHA1":checksum(raw,"sha1"),
      "published_original_SHA1":EXPECTED_SHA1,
      "download_MD5":checksum(raw,"md5"),
      "original_exact_SHA1":checksum(raw,"sha1")==EXPECTED_SHA1,
      "original_exact_file_size":len(raw)==EXPECTED_SIZE,
      "source_metadata_lock_sha256":checksum(LOCK.read_bytes(),"sha256"),
      "status":"SOURCE_NOT_AUDITED",
    }
    try:
        result=parse_bytes(raw)
        receipt.update(result)
        correct_grain=(result["source_rows"]==EXPECTED_ROWS and
                       result["column_count"]==EXPECTED_COLS)
        receipt["source_reported_grain_matches"]=correct_grain
        receipt["status"]=(
           "PUBLIC_DATAVERSENL_ORIGINAL_57669_VOYAGES_SHA1_AND_GRAIN_VERIFIED"
           if receipt["original_exact_SHA1"] and correct_grain
           else "SOURCE_GRANULARITY_OR_ORIGINAL_SHA1_UNVERIFIED"
        )
    except (ValueError,csv.Error,UnicodeError,zipfile.BadZipFile) as exc:
        receipt["status"]="NOT_HISTORICAL_TABULAR_DATA"
        receipt["failure"]=f"{type(exc).__name__}: {str(exc)[:350]}"
    return receipt


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--output",type=Path,default=OUT)
    args=ap.parse_args()
    report=audit(args.source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",
                           encoding="utf8")
    print(json.dumps({k:v for k,v in report.items() if k!="per_column"},
                     indent=2,ensure_ascii=False))
    if report["status"]!="PUBLIC_DATAVERSENL_ORIGINAL_57669_VOYAGES_SHA1_AND_GRAIN_VERIFIED":
        raise SystemExit(2)


if __name__=="__main__":
    main()
