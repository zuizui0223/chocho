#!/usr/bin/env python3
"""Restore original public CC0 source bytes from Zenodo and fail if MD5 differs.

The browser-readable previews may remove BOM and normalize line breaks.
Only raw HTTPS file bytes are acceptable for final ecological reanalysis.
"""
from __future__ import annotations
import hashlib
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "data/exploratory/realized_host_20261019"
FILES = {
    "vanessa/larv.plants.csv": ("https://zenodo.org/records/10689711/files/larv.plants.csv?download=1",69772,"0af2cbb249494e11e8dc81a68dcf7bae"),
    "vanessa/larv.switched.csv": ("https://zenodo.org/records/10689711/files/larv.switched.csv?download=1",33479,"af98d718442e42d56a89d728478ff21e"),
    "danaus/field_performance.txt": ("https://zenodo.org/records/7465143/files/field_performance.txt?download=1",2764,"5437505ae948d0133927e44d18acc617"),
    "danaus/lab_performance.txt": ("https://zenodo.org/records/7465143/files/lab_performance.txt?download=1",1911,"4e0c3fd3866a9668bab975398a01d708"),
    "melitaea/Survival_data.csv": ("https://zenodo.org/records/4279915/files/Survival_data.csv?download=1",48291,"373acf78433faddd218bc9d9c427fc9a"),
}

def restore():
    for name,(url,expected_size,expected_md5) in FILES.items():
        p=DIR/name
        if p.exists():
            b=p.read_bytes()
            if len(b)==expected_size and hashlib.md5(b).hexdigest()==expected_md5:
                print(f"EXACT already cached: {name} ({expected_size})",flush=True)
                continue
        last=None
        for attempt in range(3):
            try:
                req=urllib.request.Request(url,headers={"User-Agent":"chocho-original-source-integrity-audit/0.1","Accept":"application/octet-stream"})
                with urllib.request.urlopen(req,timeout=40) as response:
                    b=response.read()
                checksum=hashlib.md5(b).hexdigest()
                if len(b)!=expected_size or checksum!=expected_md5:
                    raise ValueError(f"ORIGINAL_SOURCE_CHECKSUM_MISMATCH: {name} bytes={len(b)} md5={checksum}")
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_bytes(b)
                print(f"PASS {name}: bytes={len(b)}, MD5={checksum}",flush=True)
                break
            except (OSError,ValueError) as e:
                last=e
                if "MISMATCH" in str(e): raise
                if attempt<2: time.sleep(2*(attempt+1))
        else:
            raise RuntimeError(f"SOURCE_UNAVAILABLE {name}: {last}")
    print("ALL FIVE ORIGINAL SOURCE MD5 VERIFIED",flush=True)

if __name__=="__main__":
    restore()
