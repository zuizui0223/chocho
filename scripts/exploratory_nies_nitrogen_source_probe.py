#!/usr/bin/env python3
"""Inspect independent NIES GNINPO-DSM v2025 source without changing its meaning.

No made-up port nitrogen allocations; the source is COUNTRY x YEAR 2005–2022.
This is a source QA/provenance receipt, not a model or local N exposure.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

EXPECTED_MD5="67db41b6d892c5f1055b9403e194d046"
DOI="10.17595/20260227.001"


def inspect(raw: bytes) -> dict:
    md5=hashlib.md5(raw).hexdigest()
    if md5!=EXPECTED_MD5:
        raise ValueError(f"GNINPO official CSV MD5 mismatch: {md5}, expected {EXPECTED_MD5}")
    text=raw.decode("utf-8-sig")
    rows=list(csv.reader(io.StringIO(text)))
    if len(rows)<20:
        raise ValueError("Not enough rows: wrong download or HTML response")
    idx=next((i for i,r in enumerate(rows[:25]) if
              any("year" in s.lower() for s in r) and
              any("iso" in s.lower() or "country" in s.lower() for s in r)),None)
    if idx is None:
        idx=0  # Record actual schema, do not silently invent field mapping.
    h=rows[idx]
    data=[r for r in rows[idx+1:] if len(r)==len(h) and any(s.strip() for s in r)]
    return {"dataset":"GNINPO-DSM Version 2025","doi":DOI,"source_type":"COUNTRY_YEAR_MODELLED_NITROGEN_2005_2022",
            "sha256":hashlib.sha256(raw).hexdigest(),"md5":md5,"raw_bytes":len(raw),
            "source_rows":len(data),"header_position":idx,"columns":h,
            "sample_rows":data[:4],
            "columns_with_N_code":[c for c in h if any(k in c.upper() for k in ("UNM","PMN","NIS","TNE"))],
            "no_spatial_allocation_to_ports":True,
            "not_direct_N_load_to_terrestrial_vegetation":True,
            "needs_units_and_aggregation_check_before_regression":True}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    r=inspect(a.source.read_bytes())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("NIES_GNINPO_OFFICIAL_SOURCE_SCHEMA",json.dumps(r,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
