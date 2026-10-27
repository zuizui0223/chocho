#!/usr/bin/env python3
"""Source-only schema audit for original 2024 Zenodo nonnative insect dataset.

Turner et al. 2024 v0.2.1 International Non-native Insect Establishment Data,
Zenodo record 11519262. Exports only public file schemas/column names,
source hashes, row/grain counts and tiny non-identifying public previews.
No host-resource/insect-arrival ecological effects inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

SCHEMA="chocho_zenodo_11519262_international_nonnative_insect_species_region_source_schema_20261027_v01"
RECORD_ID=11519262
API=f"https://zenodo.org/api/records/{RECORD_ID}"
REQUIRED={"occurrence_table.csv","attribute_table.csv",
          "metatada_occurrence_table.csv","metatada_attribute_table.csv",
          "metatada_taxonomy_table.csv"}


def get(url,max_bytes=9_000_000):
    errors=[]
    for attempt in range(3):
        req=Request(url,headers={
            "User-Agent":"chocho-scientific-source-audit/1.0",
            "Accept":"application/json,text/csv,application/octet-stream,*/*",
        })
        try:
            with urlopen(req,timeout=65) as r:
                payload=r.read(max_bytes+1)
                if len(payload)>max_bytes:
                    raise ValueError("Unexpected Zenodo source exceeds size gate")
                return payload,r.url,r.headers.get("Content-Type")
        except (HTTPError,URLError,TimeoutError) as e:
            errors.append(f"{type(e).__name__}: {str(e)[:140]}")
            if attempt<2: time.sleep(2**attempt)
    raise RuntimeError("Zenodo original open source failed: "+str(errors))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def parse_csv(filename,raw):
    try:
        txt=raw.decode("utf-8-sig")
        encoding="utf-8-sig"
    except UnicodeDecodeError:
        txt=raw.decode("latin-1")
        encoding="latin-1 (byte-conserving fallback)"
    initial=txt[:3072]
    try:
        dialect=csv.Sniffer().sniff(initial,delimiters=",;\t")
        delimiter=dialect.delimiter
    except csv.Error:
        delimiter=","
    reader=csv.DictReader(io.StringIO(txt),delimiter=delimiter)
    if reader.fieldnames is None:
        raise ValueError("Original source lacks header "+filename)
    headers=reader.fieldnames
    if len(headers)!=len(set(headers)):
        raise ValueError("Duplicate column names "+filename)
    previews=[]
    missing=Counter()
    row_count=0
    for row in reader:
        row_count+=1
        if len(previews)<2:
            previews.append({k:(row.get(k) or "")[:115]
                             for k in headers[:18]})
        for k in headers:
            if row.get(k) is None or row.get(k)=="":
                missing[k]+=1
    return {
        "name":filename,
        "original_source_sha256":sha(raw),
        "original_source_bytes":len(raw),
        "original_text_encoding":encoding,
        "original_delimiter_repr":repr(delimiter),
        "original_header_columns":headers,
        "original_column_count":len(headers),
        "original_source_data_rows":row_count,
        "source_unlabeled_or_empty_value_counts":dict(missing),
        "first_two_public_source_rows_schema_preview":previews,
    }


def fetch():
    raw,url,kind=get(API)
    original=json.loads(raw)
    if original.get("id")!=RECORD_ID:
        raise ValueError("Not original Zenodo record ID")
    meta=original.get("metadata") or {}
    version=str(meta.get("version") or "")
    # Upstream publication calls this v0.2.1. Minor upstream
    # metadata differences should be reported, not silently changed.
    members=original.get("files") or []
    if isinstance(members,dict):
        members=list(members.values())
    names={r.get("key") or r.get("filename") for r in members}
    if not REQUIRED<=names:
        raise ValueError("2024 Zenodo original species/attribute metadata missing; available="+repr(sorted(names)))
    byname={}
    for item in members:
        name=item.get("key") or item.get("filename")
        if name in REQUIRED:
            link=(item.get("links") or {}).get("self") or (
                item.get("links") or {}).get("content")
            if not link:
                raise ValueError("No original public source content link for "+name)
            data,final,content_type=get(link)
            if int(item.get("size") or 0) and len(data)!=int(item["size"]):
                raise ValueError("Zenodo source declared member size mismatch "+name)
            checksum=item.get("checksum")
            if checksum and checksum.startswith("md5:"):
                if hashlib.md5(data).hexdigest()!=checksum.partition(":")[2]:
                    raise ValueError("Zenodo original md5 mismatch "+name)
            parsed=parse_csv(name,data)
            parsed["Zenodo_original_declared_checksum"]=checksum
            parsed["Zenodo_original_download_url"]=link
            parsed["source_final_url"]=final
            parsed["mime_type"]=content_type
            byname[name]=parsed
    return {
        "schema":SCHEMA,
        "status":"ORIGINAL_PUBLIC_ESTABLISHED_INSECT_OCCURRENCE_TABLES_SCHEMA_ONLY",
        "record_id":RECORD_ID,
        "original_record_metadata_sha256":sha(raw),
        "original_record_url":url,
        "original_dataset_DOI":meta.get("doi") or original.get("doi"),
        "original_dataset_title":meta.get("title"),
        "original_version":version,
        "original_dataset_created":original.get("created"),
        "original_dataset_published":meta.get("publication_date"),
        "source_original_declared_files":sorted(names),
        "required_source_files":byname,
        "scientific_inference_stop":[
            "The original occurrence table is an ALIEN ESTABLISHMENT POSITIVE LIST, not all unsuccessful insect introductions.",
            "A species unlisted in any recipient zone is NOT proven unestablished or biologically absent.",
            "2024 establishment evidence cannot be retroactively matched to a named 1994-2019 intercepted consignment unless source individual linkage exists.",
            "The exact column names and original source grain must be inspected before deriving species×zone outcomes.",
            "Original 239 butterfly HOSTS source is independent but WCVP accepted names may not match original Zenodo taxonomic synonym practices.",
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    x=fetch()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(x,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print(json.dumps({
        "status":x["status"],"record":x["record_id"],
        "original_doi":x["original_dataset_DOI"],
        "original_version":x["original_version"],
        "source_files":[{
           "name":k,"n":v["original_source_data_rows"],
           "columns":v["original_header_columns"],
           "example":v["first_two_public_source_rows_schema_preview"],
           "sha":v["original_source_sha256"]
        } for k,v in x["required_source_files"].items()],
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
