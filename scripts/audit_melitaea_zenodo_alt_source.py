#!/usr/bin/env python3
"""Audit alternative *independent* Melitaea source, not a DiLeo-v3 substitute.

No modelling; outputs only original ZIP provenance, safe file inventory and CSV
column names.  The Zenodo file belongs to Schulz et al. (2020), NOT DiLeo (2024).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

RECORD="4987060"
ARCHIVE="ECOG-04799.zip"
SOURCE_URLS=[
    f"https://zenodo.org/records/{RECORD}/files/{ARCHIVE}?download=1",
    f"https://zenodo.org/api/records/{RECORD}/files/{ARCHIVE}/content",
]
EXPECTED_MD5="69122a1d82b1fb970fb6638b02da3db4"
MAX_ARCHIVE=12_000_000
MAX_MEMBER=20_000_000
MAX_FILES=500
TEXT_EXTENSIONS={".csv",".tsv",".txt"}

def inventory(raw:bytes):
    if len(raw)>MAX_ARCHIVE or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("not a bounded ZIP source")
    output=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        info=z.infolist()
        if len(info)>MAX_FILES:
            raise ValueError("too many ZIP members")
        for member in info:
            if member.is_dir():
                continue
            path=Path(member.filename)
            item={"path":member.filename,"size":member.file_size}
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("unsafe ZIP path")
            if member.file_size>MAX_MEMBER:
                item["status"]="OVERSIZED_MEMBER_NO_READ"
            elif path.suffix.lower() in TEXT_EXTENSIONS and member.file_size<5_000_000:
                text=z.read(member).decode("utf-8-sig",errors="replace")
                header=text.splitlines()[:1]
                item["first_line"]=(header[0][:600] if header else "")
                if path.suffix.lower() in {".csv",".tsv"} and header:
                    delim="\t" if path.suffix.lower()==".tsv" else ","
                    item["columns"]=next(csv.reader(header,delimiter=delim))[:80]
                    item["approx_line_count"]=len(text.splitlines())
            output.append(item)
    return output

def fetch(opener=urlopen):
    attempts=[]
    for url in SOURCE_URLS:
        try:
            request=Request(url,headers={"User-Agent":"chocho-scientific-provenance-check/1.0"})
            with opener(request,timeout=35) as r:
                raw=r.read(MAX_ARCHIVE+1)
                status=getattr(r,"status",None)
            md5=hashlib.md5(raw).hexdigest()
            if md5!=EXPECTED_MD5:
                attempts.append({"url":url,"http_status":status,
                                 "problem":"MD5 mismatch","received_md5":md5})
                continue
            return raw,url,attempts
        except HTTPError as exc:
            attempts.append({"url":url,"http_status":exc.code,"problem":"HTTP denied"})
        except (OSError,URLError,TimeoutError) as exc:
            attempts.append({"url":url,"problem":type(exc).__name__})
    return None,None,attempts

def audit():
    raw,url,attempts=fetch()
    result={
        "schema":"chocho_alternative_melitaea_source_inventory_v01",
        "source":"Schulz, Vanhatalo & Saastamoinen 2020 Ecography",
        "source_doi":"10.5061/dryad.ksn02v707",
        "zenodo_record":RECORD,
        "separate_from_DiLeo_2024":True,
        "expected_archive_md5":EXPECTED_MD5,
        "attempts":attempts,
        "original_bytes_confirmed":raw is not None,
        "local_extinction_effect_estimated":False,
        "host_fallback_effect_estimated":False,
        "GEB_PR38_modified":False,
    }
    if raw is None:
        result["status"]="SOURCE_ACCESS_BLOCKED"
        return result
    result["source_url"]=url
    result["archive_sha256"]=hashlib.sha256(raw).hexdigest()
    try:
        result["inventory"]=inventory(raw)
        result["status"]="SOURCE_VERIFIED_INVENTORY_ONLY"
    except (ValueError,zipfile.BadZipFile) as exc:
        result["status"]="SOURCE_INVENTORY_INVALID"
        result["error"]=str(exc)
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    result=audit()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
