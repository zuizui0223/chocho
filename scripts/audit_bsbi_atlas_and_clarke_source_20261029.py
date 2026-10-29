#!/usr/bin/env python3
"""SHA-pinned original BSBI atlas download + published Dryad host table metadata.

Source-only: excludes botanical species observations, grid-cell rows, host use,
protected source rows, observer names and butterfly responses from output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

SOURCES = (
    ("bsbi_pre2019_10km_dateclasses", 10998003, "PA2020_10by10kmSquareData.zip", "66a8bb4305df484839587b5fe8acf1d2", 25_000_000),
    ("bsbi_uk_native_status", 11097371, "PA2020_NativeStatusData.zip", "a5f60929e6118b77094110d6ca77bb48", 2_000_000),
)
DRYAD_DOI = "doi:10.5061/dryad.1vhhmgr12"
MAX_UNCOMPRESSED = 350_000_000
MAX_ARCHIVE_MEMBERS = 30


def fetch(url: str, max_bytes: int) -> bytes:
    if not url.startswith(("https://zenodo.org/", "https://datadryad.org/")):
        raise ValueError("Unapproved original publisher domain")
    req = Request(url, headers={"User-Agent":"chocho-BSBI-UK-host-hindcast-source-audit/1.0", "Accept":"application/json,application/zip,*/*"})
    with urlopen(req, timeout=80) as response:
        if response.status != 200:
            raise ValueError("Original publisher returned non-200 HTTP")
        data = response.read(max_bytes + 1)
    if not data or len(data)>max_bytes:
        raise ValueError("Original publisher source unavailable/over bounded transfer")
    return data


def catalogue_members(archive: bytes):
    if archive[:4] != b"PK\x03\x04":
        raise ValueError("Original publisher ZIP is not a zip archive")
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        if z.testzip() is not None:
            raise ValueError("Original BSBI ZIP CRC mismatch")
        entries = [v for v in z.infolist() if not v.is_dir()]
        if not entries or len(entries)>MAX_ARCHIVE_MEMBERS or sum(v.file_size for v in entries)>MAX_UNCOMPRESSED:
            raise ValueError("Unexpected source ZIP member or unpacked size")
        metadata=[]
        for f in entries:
            p=Path(f.filename)
            if p.is_absolute() or ".." in p.parts:
                raise ValueError("Original source ZIP unsafe path")
            info={"basename": p.name[:180], "unpacked_bytes": f.file_size,
                  "original_suffix":p.suffix.lower()}
            if p.suffix.lower() in {".csv", ".tsv", ".txt"}:
                with z.open(f) as contents:
                    reader=csv.reader(io.TextIOWrapper(contents, encoding="utf-8-sig",errors="replace",newline=""),
                                      delimiter="\t" if p.suffix.lower()==".tsv" else ",")
                    info["header_fields"]=[str(s).strip()[:120] for s in next(reader, [])][:100]
                    # This stage does not open botanical species/cell outcome rows.
                    info["original_observation_rows_published"] = False
            metadata.append(info)
    return metadata


def inspect_zenodo(name: str, rid: int, filename: str, expected_md5: str, maxbytes: int):
    result={"key":name,"record_id":rid,"declared_file":filename,
            "reference":"https://doi.org/10.5281/zenodo."+str(rid),"expected_published_md5":expected_md5}
    try:
        meta=json.loads(fetch("https://zenodo.org/api/records/"+str(rid), 2_000_000))
        if str(meta.get("id"))!=str(rid):
            raise ValueError("Original Zenodo record ID drift")
        files=meta.get("files",[])
        candidates=[f for f in files if f.get("key","").split("/")[-1]==filename]
        if len(candidates)!=1:
            raise ValueError("Exactly one original published archive not found")
        f=candidates[0]
        result["metadata_file_size"]=f.get("size")
        result["metadata_checksum"]=f.get("checksum")
        uri=(f.get("links") or {}).get("content")
        if not isinstance(uri,str) or not uri.startswith("https://zenodo.org/"):
            uri="https://zenodo.org/records/"+str(rid)+"/files/"+filename+"?download=1"
        original=fetch(uri,maxbytes)
        md5=hashlib.md5(original).hexdigest()
        if md5 != expected_md5:
            raise ValueError("Original published source MD5 mismatch")
        result.update({"status":"ORIGINAL_MD5_VERIFIED_FORMAT_AND_HEADERS",
                       "source_zip_sha256":hashlib.sha256(original).hexdigest(),
                       "source_zip_bytes":len(original),
                       "member_source_schema":catalogue_members(original)})
    except (OSError,ValueError,TypeError,KeyError,zipfile.BadZipFile,HTTPError,URLError) as e:
        result["status"]="ORIGINAL_SOURCE_BLOCKED_OR_FAILED_QA"
        result["error_class"]=type(e).__name__
        result["error_reason"]=str(e)[:240]
    return result


def inspect_dryad():
    url="https://datadryad.org/api/v2/datasets/"+quote(DRYAD_DOI,safe="")
    out={"source_doi":DRYAD_DOI,"source_url":url}
    try:
        data=json.loads(fetch(url,2_000_000))
        out["metadata_title"]=str(data.get("title", ""))[:240]
        out["metadata_identifier"]=data.get("identifier")
        out["metadata_version_id"]=data.get("_links",{}).get("stash:version",{}).get("href")
        # Search dataset JSON for version href, but do not request restricted files.
        metadata=data.get("_links",{})
        out["public_link_relations"]=[str(k)[:90] for k in metadata][:40]
        out["status"]="ORIGINAL_DRYAD_PUBLISHED_METADATA_AVAILABLE"
    except (OSError,ValueError,TypeError,HTTPError,URLError) as e:
        out["status"]="ORIGINAL_DRYAD_METADATA_UNAVAILABLE"
        out["error_class"]=type(e).__name__
        out["error_reason"]=str(e)[:240]
    return out


def run(output: Path):
    zenodo=[inspect_zenodo(*args) for args in SOURCES]
    dryad=inspect_dryad()
    result={"schema":"chocho_UK_historical_plantatlas_and_hostlink_source_gate_20261029_v01",
            "source_stage":"BOTANICAL_ORIGINAL_SOURCE_SCHEMA_ONLY_NO_BUTTERFLY_OUTCOMES",
            "status":"BOTANICAL_SOURCE_VERIFIED" if zenodo[0]["status"]=="ORIGINAL_MD5_VERIFIED_FORMAT_AND_HEADERS" else "BOTANICAL_SOURCE_UNAVAILABLE",
            "zenodo":zenodo,"dryad":dryad,
            "source_EVIDENCE_facts":{"BSBI_10km_is_observed_plant_presences_not_absences":True,
             "BSBI_plants_recorded_after_1999_cannot_be_training_predictors":True,
             "native_status_does_not_show_human_planting_or_maintenance":True,
             "missing_plants_are_unknown_not_confirmed_absent":True,
             "no_ukbms_transect_record_or_botanical_hectad_row_exported":True}}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":result["status"],"zenodo":[{"key":x["key"],"status":x["status"],"members":x.get("member_source_schema"),"reason":x.get("error_reason")} for x in zenodo],
           "dryad_status":dryad["status"]},ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    outcome=run(p.parse_args().output)
    if outcome["status"]!="BOTANICAL_SOURCE_VERIFIED":raise SystemExit(2)
