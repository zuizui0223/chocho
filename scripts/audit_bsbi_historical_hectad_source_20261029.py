#!/usr/bin/env python3
"""Source-only, read-only audit of BSBI Plant Atlas 2020 10km date-class archive.

No botanical occurrence rows, hectad locations, individual plant records,
original UKBMS transects or user identifiers are exported.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import zipfile

ZENODO_RECORD="https://zenodo.org/records/10998003"
SOURCE_NAME="PA2020_10by10kmSquareData.zip"
EXPECTED_MD5="66a8bb4305df484839587b5fe8acf1d2"
URLS=(
    f"{ZENODO_RECORD}/files/{SOURCE_NAME}?download=1",
    f"https://zenodo.org/api/records/10998003/files/{SOURCE_NAME}/content",
)
MAX_SOURCE=30_000_000
MAX_UNPACKED=650_000_000
MAX_CSV=250_000_000
SCHEMA="chocho_bsbi_plantatlas2020_hectad_dateclass_source_schema_20261029_v01"

def classify_header(m):
    columns=[str(x).strip()[:100] for x in m]
    signature=" ".join(columns).lower()
    return {"original_source_field_names":columns[:80],
            "plausibly_geographic":any(s in signature for s in ("hectad","grid","square","easting","northing","osgr")),
            "plausibly_taxonomic":any(s in signature for s in ("taxon","species","sci_name","scientific","taxon_name")),
            "plausibly_dated":any(s in signature for s in ("period","date","year","dc0","dc1","dc2","dc3","dc4"))}

def scan_member(z,info):
    name=info.filename
    suffix=Path(name).suffix.lower()
    item={"filename":name[:200], "suffix":suffix,
          "uncompressed_bytes":info.file_size}
    if suffix not in {".csv",".tsv",".txt"} or info.file_size>MAX_CSV:
        item["source_rows_scanned"]=False
        return item
    sep="\t" if suffix==".tsv" else ","
    with z.open(info) as stream:
        rdr=csv.reader(io.TextIOWrapper(stream,encoding="utf-8-sig",
                          errors="replace",newline=""),delimiter=sep)
        header=next(rdr,[])
        n=sum(1 for _ in rdr)
    item["source_rows_scanned"]=True
    item["source_tabular_rows"]=n
    item.update(classify_header(header))
    return item

def inspect_archive(data):
    if not data.startswith(b"PK\x03\x04"):
        raise ValueError("Original Zenodo archive must be a genuine ZIP, not an HTML/JSON response")
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        info=[x for x in z.infolist() if not x.is_dir()]
        if not info or sum(x.file_size for x in info)>MAX_UNPACKED:
            raise ValueError("ZIP archive empty or uncompressed source exceeds bound")
        for file in info:
            path=Path(file.filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("Unsafe original archive path")
        bad=z.testzip()
        if bad:
            raise ValueError("Original ZIP CRC invalid")
        members=[scan_member(z,x) for x in info]
    return {"original_archive_member_count":len(info),
            "file_type_counts":dict(sorted(Counter(x["suffix"] for x in members).items())),
            "original_source_members_schema":members[:60],
            "source_rows_scanned_total":sum(x.get("source_tabular_rows",0) for x in members),
            "source_tabular_files":sum(x.get("source_rows_scanned",False) for x in members),
            "no_plant_species_by_hectad_rows_exported":True}

def obtain():
    failures=[]
    for url in URLS:
        try:
            req=Request(url,headers={"User-Agent":"chocho-BSBI-historical-source-audit/1.0",
                                     "Accept":"application/zip,application/octet-stream,*/*"})
            with urlopen(req,timeout=110) as r:
                content=r.read(MAX_SOURCE+1)
                if r.status!=200:raise ValueError("Original Zenodo source not HTTP 200")
                domain=(r.url.split("/")[2]).lower()
                if "zenodo.org" not in domain and "cern.ch" not in domain:
                    raise ValueError("Original URL unexpected external redirect host")
            if len(content)>MAX_SOURCE:
                raise ValueError("Original source exceeds bounded 30MB")
            actual=hashlib.md5(content).hexdigest()
            if actual!=EXPECTED_MD5:
                raise ValueError(f"Unexpected source MD5: {actual}")
            return content,url,failures
        except (OSError,HTTPError,URLError,ValueError) as exc:
            failures.append({"source_url":url,"error_class":type(exc).__name__,
                             "reason":str(exc)[:160]})
    raise ValueError("Original Zenodo source unavailable/changed: "+str(failures)[:400])

def run(output):
    result={"schema":SCHEMA,"source_record":ZENODO_RECORD,
            "DOI":"10.5281/zenodo.10998003",
            "original_file":SOURCE_NAME,
            "original_publisher":"Botanical Society of Britain and Ireland / UK Centre for Ecology & Hydrology",
            "original_expected_md5":EXPECTED_MD5,
            "source_description":"Presence-only historical 10×10 km British/Irish botanical atlas by period, NOT plant absence",
            "host_taxa_accepted_or_human_management_evaluated":False,
            "source_rows_or_hectad_locations_exported":False}
    try:
        data,url,failures=obtain()
        result["retrieved_original_url"]=url
        result["original_archive_sha256"]=hashlib.sha256(data).hexdigest()
        result["original_archive_size_bytes"]=len(data)
        result["download_fallbacks"]=failures
        result.update(inspect_archive(data))
        result["status"]="SOURCE_VERIFIED_BSBI_HISTORICAL_GRID_SCHEMA"
    except (OSError,ValueError,zipfile.BadZipFile) as exc:
        result["status"]="SOURCE_BLOCKED_OR_NOT_VERIFIED"
        result["blocker_type"]=type(exc).__name__
        result["blocker"]=str(exc)[:400]
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:result.get(k) for k in ("status","original_archive_sha256",
        "original_archive_member_count","source_rows_scanned_total",
        "source_tabular_files","blocker")},sort_keys=True))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    r=run(a.output)
    if r["status"]!="SOURCE_VERIFIED_BSBI_HISTORICAL_GRID_SCHEMA":
        raise SystemExit(2)
