#!/usr/bin/env python3
"""Source-only Figshare audit: prior-published butterfly–parasitoid study.

Never treat a verified archive or pre-existing parasitism as new ecological discovery.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, re, zipfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ARTICLE_ID=14260211
META=f"https://api.figshare.com/v2/articles/{ARTICLE_ID}"
MAX_FILE=15_000_000
CSV_MAX=8_000_000

def request_bytes(url: str, limit=MAX_FILE, opener=urlopen):
    req=Request(url,headers={"User-Agent":"chocho-source-provenance-audit/1.0","Accept":"application/json, text/csv, */*"})
    with opener(req, timeout=30) as res:
        body=res.read(limit+1)
    if len(body)>limit or not body:
        raise ValueError("empty/oversized original source")
    return body

def parse_csv_preview(data:bytes, filename:str):
    txt=data.decode("utf-8-sig")
    if not txt.strip():
        raise ValueError("empty CSV")
    sniffer=csv.Sniffer()
    lines=txt.splitlines()
    delim="\t" if filename.lower().endswith(".tsv") else ","
    try:
        dialect=sniffer.sniff("\n".join(lines[:10]),delimiters=",;\t")
        delim=dialect.delimiter
    except csv.Error:
        pass
    reader=csv.reader(io.StringIO(txt),delimiter=delim)
    header=next(reader)
    if len(header)<2 or len(header)>300:
        raise ValueError("CSV lacks consistent multi-column header")
    count=0
    blank=0
    for row in reader:
        if not any(cell.strip() for cell in row):
            blank+=1
            continue
        if len(row)!=len(header):
            raise ValueError(f"row width mismatch {len(row)} expected {len(header)}")
        count+=1
    return {"columns":header,"records":count,"blank_rows":blank,"delimiter":"TAB" if delim=="\t" else delim}

def gate_inventory(meta:dict,fetch=request_bytes):
    if int(meta.get("id",-1)) != ARTICLE_ID:
        raise ValueError("wrong original Figshare article ID")
    fs=meta.get("files")
    if not isinstance(fs,list) or not fs:
        raise ValueError("missing deposited source files")
    receipt={"schema":"chocho_nettle_parasitoid_source_inventory_v01",
             "archive_id":ARTICLE_ID,"article_title":meta.get("title"),
             "doi":meta.get("doi"),"publication":meta.get("published_date"),
             "license":meta.get("license"),
             "study_paper":"Audusseau et al. 2021 Oikos 10.1111/oik.07953",
             "source_files":[],"unreviewed_files":[],"raw_butterfly_records_verified":False,
             "new_ecological_effect_estimated":False,
             "botanical_introduction_tested":False,"causal_heterospecific_transmission_tested":False,
             "GEB_PR38_unchanged":True}
    saw_raw=False
    for f in fs:
        name=str(f.get("name") or "")
        item={"name":name,"id":f.get("id"),"declared_size":f.get("size"),
              "declared_md5":f.get("supplied_md5") or f.get("computed_md5"),
              "download_url":f.get("download_url")}
        srcurl=f.get("download_url")
        if not name.lower().endswith((".csv",".tsv")):
            item["status"]="FILE_LISTED_NOT_PARSED"
            receipt["source_files"].append(item)
            continue
        if not srcurl or f.get("size",MAX_FILE+1)>CSV_MAX:
            item["status"]="MISSING_OR_OVERSIZED_SOURCE"
            receipt["source_files"].append(item)
            continue
        try:
            raw=fetch(srcurl,CSV_MAX)
            md5=hashlib.md5(raw).hexdigest()
            reported=item["declared_md5"]
            if reported and reported!=md5:
                raise ValueError("original deposited MD5 mismatch")
            item.update({"downloaded_sha256":hashlib.sha256(raw).hexdigest(),
                         "verified_md5":md5,"schema":parse_csv_preview(raw,name),
                         "status":"SOURCE_CSV_VALIDATED"})
            if "batch" in name.lower() or "monitor" in name.lower():
                saw_raw=True
        except (HTTPError,URLError,OSError,TimeoutError,ValueError,UnicodeError,csv.Error) as exc:
            item.update({"status":"SOURCE_CSV_ACCESS_OR_SCHEMA_FAILURE",
                         "error":f"{type(exc).__name__}: {str(exc)[:120]}"})
        receipt["source_files"].append(item)
    receipt["raw_butterfly_records_verified"]=bool(saw_raw)
    receipt["status"]="RAW_SOURCE_SCHEMA_VERIFIED_ONLY" if saw_raw else "SOURCE_INVENTORY_ONLY"
    return receipt

def audit():
    raw=request_bytes(META,1_000_000)
    meta=json.loads(raw)
    return gate_inventory(meta)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    try:
        result=audit()
    except (HTTPError,URLError,OSError,TimeoutError,ValueError,UnicodeError,json.JSONDecodeError) as e:
        result={"schema":"chocho_nettle_parasitoid_source_inventory_v01","status":"SOURCE_ACCESS_BLOCKED",
                "reason":f"{type(e).__name__}: {str(e)[:300]}",
                "new_ecological_effect_estimated":False,"GEB_PR38_unchanged":True}
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
