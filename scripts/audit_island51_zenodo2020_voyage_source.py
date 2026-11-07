#!/usr/bin/env python3
"""Verify Zenodo Duvallet et al. 2020 AIS date-ship-port triplet archive and frozen WPI source.

Acquires public source only; does NOT assume the port IDs are original WPI OBJECTIDs
or that triplet adjacency equals a nonstop/biologically relevant vessel voyage.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import random
import re
import time
import urllib.request
import urllib.error
import zipfile

ROOT=Path(__file__).resolve().parents[1]
PORTS=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
PORTS_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
ZENODO="https://zenodo.org/records/10783698/files/triplet_2020.zip?download=1"
ZENODO_MD5="af27e1bef1f706e2e29ed53aee47ff10"
MAX_ARCHIVE_BYTES=85_000_000
MAX_UNCOMPRESSED_MEMBER_BYTES=2_000_000_000

def blob(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()

def original_port_roster():
    if blob(PORTS)!=PORTS_GIT_BLOB:
        raise ValueError("Original WPI194 geometry Git blob changed")
    with PORTS.open(encoding="utf-8",newline="") as f:
        data=list(csv.DictReader(f))
    if len(data)!=194 or len({int(r["port_id"]) for r in data})!=194:
        raise ValueError("Original WPI port identity cohort differs")
    if len({r["island_id"] for r in data})!=51:
        raise ValueError("Original 51 island groups missing")
    if sum(str(r["included"]).lower()=="true" for r in data)!=182:
        raise ValueError("Original 182 qualified marine ports differ")
    if any(r["landmass_id"]!=r["island_id"] for r in data):
        raise ValueError("Original WPI source island/landmass mismatch")
    for r in data:
        a=float(r["original_latitude"]);b=float(r["original_longitude"])
        if not(-90<=a<=90 and -180<=b<=180):
            raise ValueError("Invalid original WPI source coordinates")
    return data

def download(dest):
    if dest.exists():
        raw=dest.read_bytes()
        if hashlib.md5(raw).hexdigest()!=ZENODO_MD5:
            raise ValueError("Local Zenodo source checksum mismatch")
        return len(raw),"verified_local_cache"
    last=None
    for attempt in range(1,5):
        try:
            request=urllib.request.Request(
                ZENODO,headers={"User-Agent":"chocho-20261107-port-voyage-source-audit/1.0",
                                "Accept":"application/zip,application/octet-stream,*/*"})
            with urllib.request.urlopen(request,timeout=120) as response:
                length=response.headers.get("Content-Length")
                if length and int(length)>MAX_ARCHIVE_BYTES:
                    raise ValueError("External archive exceeds frozen size limit")
                with dest.open("wb") as output:
                    n=0
                    while True:
                        data=response.read(1024*1024)
                        if not data:break
                        n+=len(data)
                        if n>MAX_ARCHIVE_BYTES:
                            raise ValueError("External archive exceeds frozen download size limit")
                        output.write(data)
            check=hashlib.md5(dest.read_bytes()).hexdigest()
            if check!=ZENODO_MD5:
                raise ValueError("Publisher Zenodo MD5 mismatch: "+check)
            return n,"source_download_success"
        except (urllib.error.URLError,TimeoutError,ConnectionError) as exc:
            last=exc
            if dest.exists():dest.unlink()
            print("ZENODO_SOURCE_RETRY="+json.dumps({"attempt":attempt,"type":type(exc).__name__,
                "error":str(exc)[:180]},sort_keys=True),flush=True)
            if attempt<4:time.sleep(min(30,2**attempt))
    raise RuntimeError("Source Zenodo retrieval unsuccessful after four bounded attempts: "+str(last))

def masked_sample(text):
    # 9-digit Maritime Mobile Service Identity is not needed in the public audit preview.
    text=re.sub(r"(?<!\d)\d{9}(?!\d)","[REDACTED_9_DIGIT_VESSEL_ID]",text)
    return text[:450]

def check_zip(path):
    if not zipfile.is_zipfile(path):
        raise ValueError("Published Zenodo download is not a valid ZIP")
    items=[];total=0
    with zipfile.ZipFile(path,"r") as zipped:
        members=zipped.infolist()
        if not members or len(members)>10_000:
            raise ValueError("Unexpected no/too many ZIP source members")
        for n,member in enumerate(members):
            if member.is_dir():continue
            # No extraction to disk, no member path trust.
            if member.file_size>MAX_UNCOMPRESSED_MEMBER_BYTES:
                raise ValueError("Member too large for source audit: "+member.filename)
            total+=member.file_size
            lower=member.filename.lower()
            is_text=lower.endswith((".csv",".txt",".tsv",".jsonl",".json",".dat",".log"))
            preview=[]
            if is_text or member.file_size<80_000_000:
                with zipped.open(member) as f:
                    for j in range(4):
                        raw=f.readline(2000)
                        if not raw:break
                        line=raw.decode("utf-8-sig",errors="replace").rstrip("\r\n")
                        preview.append(masked_sample(line))
            items.append({
                "member_index":n,"member_name":member.filename,
                "member_compressed_bytes":member.compress_size,
                "member_uncompressed_bytes":member.file_size,
                "compression_method":member.compress_type,
                "original_archive_crc32":f"{member.CRC:08x}",
                "has_obvious_tabular_extension":is_text,
                "first_four_redacted_lines":preview,
                "first_four_tab_field_counts":[len(x.split("\t")) for x in preview],
                "first_four_comma_field_counts":[len(list(csv.reader([x]))[0]) for x in preview],
            })
    return items,total

def run(out_dir):
    out=Path(out_dir).resolve()
    if out.is_relative_to(ROOT/"results") or out.is_relative_to(ROOT/"data") or out.is_relative_to(ROOT/"provenance"):
        raise ValueError("Never overwrite original research data and results")
    ports=original_port_roster()
    out.mkdir(parents=True,exist_ok=True)
    zip_path=out/"triplet_2020_zenodo10783698_original.zip"
    n,source_status=download(zip_path)
    members,uncompressed=check_zip(zip_path)
    receipt={
        "schema":"chocho_original194_WPI_zenodo_2020_AIS_vessel_port_triplet_source_audit_v1",
        "status":"PUBLISHER_MD5_AND_ARCHIVE_MEMBER_SCHEMA_VERIFIED",
        "publisher":"Duvallet, Kerne, Morel and Pigne (2024), Zenodo",
        "source_doi":"10.5281/zenodo.10783698",
        "source_archive_url":ZENODO,
        "source_license":"CC BY 4.0",
        "retrieval_utc":datetime.now(timezone.utc).isoformat(),
        "source_archive_original_filename":"triplet_2020.zip",
        "original_md5_published":ZENODO_MD5,
        "original_md5_recomputed":hashlib.md5(zip_path.read_bytes()).hexdigest(),
        "source_sha256":hashlib.sha256(zip_path.read_bytes()).hexdigest(),
        "compressed_bytes":n,
        "uncompressed_bytes_across_members":uncompressed,
        "source_status":source_status,
        "member_count":len(members),
        "source_member_records":members,
        "source_WPI_original194_port_git_blob":PORTS_GIT_BLOB,
        "source_WPI_original_ports":len(ports),
        "source_WPI_original_islands":len({x["island_id"] for x in ports}),
        "source_WPI_marine_eligible":sum(str(x["included"]).lower()=="true" for x in ports),
        "port_identifier_correspondence_to_WPI_not_yet_verified":True,
        "source_is_2020_not_2001_2010_shipment_chronology":True,
        "no_ship_voyage_or_butterfly_community_model_fitted":True,
        "source_archive_is_referenced_by_DOI_not_copied_to_git":True,
    }
    (out/"original2020_AIS_triplet_schema_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ZENODO_VOYAGE_SOURCE_INSPECTION_JSON="+json.dumps({
       "status":receipt["status"],"sha256":receipt["source_sha256"],
       "archive_bytes":n,"member_count":len(members),
       "uncompressed_bytes":uncompressed,
       "members":[{"name":x["member_name"],"compressed":x["member_compressed_bytes"],
        "uncompressed":x["member_uncompressed_bytes"],
        "source_samples":x["first_four_redacted_lines"][:3]} for x in members],
       "original_WPI_ports":194,"original_WPI_islands":51},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    run(p.parse_args().out_dir)
