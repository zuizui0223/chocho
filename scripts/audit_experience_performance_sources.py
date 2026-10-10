#!/usr/bin/env python3
"""Byte-level audit of original 2019 monarch experiment files. No effect fit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

URLS = {
    "oviposition_experience.xlsx": "https://datadryad.org/downloads/file_stream/81836",
    "caterpillar_performance.xlsx": "https://datadryad.org/downloads/file_stream/81832",
}
MIN_SIZE = 7000
MAX_SIZE = 500000

def audit_file(name, url, opener=urlopen):
    result = {"name": name, "url": url, "status": "SOURCE_ACCESS_BLOCKED"}
    try:
        request = Request(url, headers={"User-Agent": "chocho-independent-source-audit/1.0"})
        with opener(request, timeout=20) as response:
            data = response.read(MAX_SIZE + 1)
            result["http_status"] = response.status
            result["reported_content_type"] = response.headers.get("Content-Type")
        if len(data) > MAX_SIZE or len(data) < MIN_SIZE or data[:4] != b"PK\x03\x04":
            result["reason"] = "Not a byte-verified .xlsx ZIP file within source-size bounds"
            return result
        result["byte_count"] = len(data)
        result["sha256"] = hashlib.sha256(data).hexdigest()
        result["status"] = "SOURCE_BYTES_VERIFIED_SCHEMA_UNCHECKED"
    except HTTPError as err:
        result["http_status"] = err.code
        result["reason"] = f"HTTP {err.code}"
    except (URLError, TimeoutError, OSError) as err:
        result["reason"] = "access error: " + type(err).__name__
    return result

def audit():
    files = [audit_file(name,url) for name,url in URLS.items()]
    return {"schema": "chocho_experience_performance_source_receipt_v01",
            "source_doi": "10.5061/dryad.8hd6764",
            "sources": files,
            "source_bytes_complete": all(f["status"]=="SOURCE_BYTES_VERIFIED_SCHEMA_UNCHECKED" for f in files),
            "row_level_schema_audited": False,
            "effect_estimated": False,
            "historical_host_removal_effect_estimated": False,
            "botanical_status": "Both maternal choice targets are native Asclepias incarnata subspecies; no introduced-host treatment in this source.",
            "main_GEB_manuscript_unchanged": True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    receipt=audit()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
