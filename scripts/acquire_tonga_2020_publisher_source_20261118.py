#!/usr/bin/env python3
"""Fail-closed source accession pilot for Fall & Drezner (2020), Ecology e02902.

Fetches only the named public publisher supplementary ZIP; does not fit models,
silently substitute a different dataset, read PDF metadata, or infer water vs
seawater traits. Output JSON separates retrieval success from model validity.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

NAME="ecy2902-sup-0001-DataS1.zip"
BASE="10.1002/ecy.2902"
URLS=(
    "https://esajournals.onlinelibrary.wiley.com/doi/suppl/"+BASE+"/supinfo/"+NAME,
    "https://onlinelibrary.wiley.com/doi/suppl/"+BASE+"/supinfo/"+NAME,
    "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2902&file="+NAME,
)
OUT=Path("tonga2020_source_gate")
OUT.mkdir(exist_ok=True)
MAX_BYTES=300_000
report={"source":"Fall & Drezner (2020), Ecology 101 e02902, DOI 10.1002/ecy.2902",
        "requested_archive":NAME,
        "source_status":"NOT_FETCHED",
        "attempts":[],
        "model_fit":False,
        "trait_ontology":"Not determined: requires original species-level columns and metadata, not assumed to be EuDiS thalassochory"}
for url in URLS:
    attempt={"url":url}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (research data archive audit; +https://github.com/zuizui0223/chocho)","Accept":"application/zip, application/octet-stream, */*"})
        with urllib.request.urlopen(req,timeout=12) as handle:
            attempt["http_status"]=handle.status
            attempt["content_type"]=handle.headers.get("Content-Type","")
            data=handle.read(MAX_BYTES+1)
            attempt["bytes_downloaded"]=len(data)
            attempt["final_url"]=handle.geturl()
        if len(data)>MAX_BYTES: raise ValueError("Downloaded bytes exceed 300KB source ZIP cap")
        if not data.startswith(b"PK"): raise ValueError("Response is not a ZIP signature")
        archive=zipfile.ZipFile(io.BytesIO(data))
        members=[]
        for info in archive.infolist():
            if info.file_size>25_000_000 or info.compress_size>MAX_BYTES:
                raise ValueError("Unexpected decompressed file size")
            if info.filename.startswith("/") or ".." in Path(info.filename).parts:
                raise ValueError("Unsafe ZIP member name")
            item={"filename":info.filename,"bytes":info.file_size}
            if info.filename.lower().endswith((".csv",".txt",".tsv")) and info.file_size<=3_000_000:
                raw=archive.read(info.filename)
                try:
                    lines=raw.decode("utf-8-sig").splitlines()
                except UnicodeDecodeError:
                    lines=raw.decode("latin-1").splitlines()
                item["line_count"]=len(lines)
                item["first_line"]=lines[0][:1000] if lines else ""
                item["delimiter_guess"]="tab" if item["first_line"].count("\t")>item["first_line"].count(",") else "comma"
            members.append(item)
        (OUT/NAME).write_bytes(data)
        report.update(source_status="ZIP_BYTES_VERIFIED",sha256=hashlib.sha256(data).hexdigest(),
                      zip_member_count=len(members),zip_members=members,
                      exact_source_url=url,received_bytes=len(data))
        attempt["result"]="ZIP_VALID"
        report["attempts"].append(attempt)
        break
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError, zipfile.BadZipFile, OSError) as err:
        attempt["result"]="NOT_FETCHED"
        attempt["error_class"]=type(err).__name__
        attempt["error_detail"]=str(err)[:320]
        report["attempts"].append(attempt)
        continue
(OUT/"source_gate.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print(json.dumps(report,indent=2,ensure_ascii=False))
if report["source_status"]!="ZIP_BYTES_VERIFIED":
    print("::warning::Tonga 2020 publisher source ZIP not acquired; do not assert any species-level result")
else:
    print("::notice::Tonga 2020 source ZIP acquired; schema/taxon/trait status requires separate source review")
