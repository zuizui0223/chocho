#!/usr/bin/env python3
"""Fail-closed independent audit of the public original GBIF224 RAW->annual reconstruction.

This verifies the exact historical DOI/zip checksum and aggregate source
receipt; it never commits GBIF raw records, individual locations or datasetKeys.
It does NOT refit the 155-species RF or the 224-species ecological models.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

DOI="10.15468/dl.umqmsf"
KEY="0013264-261023074918325"
ZIP_BYTES=1415705617
ZIP_SHA256="a1769481e91cb3002c25f69c6d8aa8e6d8d73a012c9495e8394032af02a9df18"
ROWS=14637152
TABLE_EXPECTED={
    "gbif_global224_annual_v1.species_region_year.csv.gz":91882,
    "gbif_global224_annual_v1.region_year_effort.csv":5390,
    "gbif_global224_annual_v1.species_region.csv.gz":9408,
}
COMPARE_KEYS=["records","year_counts_mapped","source_basis_of_record","coverage"]


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b=f.read(1024*1024*8)
            if not b:break
            h.update(b)
    return h.hexdigest()


def verify(original:Path,reconstructed:Path,source_zip:Path)->dict:
    if source_zip.stat().st_size!=ZIP_BYTES:
        raise ValueError("Original GBIF download size drift")
    raw_sha=sha256(source_zip)
    if raw_sha!=ZIP_SHA256:
        raise ValueError("Original GBIF archive SHA256 drift")
    old=json.loads(original.read_text(encoding="utf-8"))
    new=json.loads(reconstructed.read_text(encoding="utf-8"))
    for r in (old,new):
        if r.get("schema")!="chocho_global_GBIF_224_butterfly_369_region_year_effort_v1":
            raise ValueError("Source annual audit schema mismatch")
        if r.get("sampled_mode") or r["records"]["rows"]!=ROWS:
            raise ValueError("Not full source record count")
        src=r["source"]
        if (src["gbif_doi"]!=DOI or src["gbif_download_key"]!=KEY or
            src["GBIF_zip_size_bytes"]!=ZIP_BYTES or
            src["GBIF_archive_sha256"]!=ZIP_SHA256):
            raise ValueError("Original download identity changed")
    matched={k:old[k]==new[k] for k in COMPARE_KEYS}
    count_match={name:old["files"][name]["rows"]==new["files"][name]["rows"]==value
                 for name,value in TABLE_EXPECTED.items()}
    original_hashes={name:old["files"][name]["sha256"] for name in TABLE_EXPECTED}
    rebuild_hashes={name:sha256(reconstructed.parent/name) for name in TABLE_EXPECTED}
    # A .gz file may have timestamp/header differences even with equal row content.
    # Preserve byte-level result for both compressed and plain outputs.
    raw_hash_match={name:original_hashes[name]==rebuild_hashes[name]
                    for name in TABLE_EXPECTED}
    returned={
        "schema":"chocho_GEB_v05_original_GBIF224_raw_annual_independent_rebuild_v1",
        "status":("RAW_TO_ANNUAL_FULL_ORIGINAL_RECORDS_AND_COUNTS_MATCH" if
                  all(matched.values()) and all(count_match.values())
                  else "RAW_TO_ANNUAL_ORIGINAL_RECONSTRUCTION_DIFFERENCE"),
        "source_download":{"doi":DOI,"download_key":KEY,"rows":ROWS,
                           "zip_bytes":ZIP_BYTES,"zip_sha256":ZIP_SHA256},
        "annual_receipt_same_source":matched,
        "annual_table_rows_match_source":count_match,
        "table_original_hashes":original_hashes,
        "table_rebuilt_hashes":rebuild_hashes,
        "table_file_hash_match_source":raw_hash_match,
        "scope":"Exact original GBIF raw archive to 224-species regional annual/effort aggregates; NOT raw-to-155 RF, nor 224 observed ecology model refit; no source raw records shared.",
        "licensing":"GBIF original download and 1716 contributing datasets require source-specific citation and terms; raw archive not redistributed",
        "not_reviewer_deposit":True,
    }
    return returned


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-audit",type=Path)
    p.add_argument("--reconstructed-audit",type=Path)
    p.add_argument("--source-zip",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        assert ZIP_BYTES>10**9 and ROWS==14637152
        assert len(TABLE_EXPECTED)==3 and len(ZIP_SHA256)==64
        print("GEB_GBIF224_RAW_ANNUAL_AUDIT_SELFTEST_PASS")
        return
    if not (a.original_audit and a.reconstructed_audit and a.source_zip and a.out):
        p.error("--original-audit --reconstructed-audit --source-zip --out required")
    res=verify(a.original_audit,a.reconstructed_audit,a.source_zip)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":res["status"],
        "data_tables_rows_match":res["annual_table_rows_match_source"],
        "source_receipt_keys":res["annual_receipt_same_source"]}))
    if res["status"]!="RAW_TO_ANNUAL_FULL_ORIGINAL_RECORDS_AND_COUNTS_MATCH":
        raise SystemExit("GBIF224_RAW_ANNUAL_REPLAY_NOT_SOURCE_IDENTICAL")


if __name__=="__main__":
    main()
