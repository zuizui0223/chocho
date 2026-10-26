#!/usr/bin/env python3
"""Retry exact published Dryad CSV via official public website file_stream.

Dryad REST dataset ZIP yielded 401 without credentials; the user-facing public
download link is distinct and this is a SOURCE-TRANSPORT-ONLY v02 amendment.
Data rows do not leave the short-lived runner.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from audit_western_monarch_stage_host_phenology_source_20261026 import (
    DOI, ENC, BASE, MAIN, inspect_csv
)

VERSIONS=BASE+"/datasets/"+ENC+"/versions"
PUBLIC="https://datadryad.org/downloads/file_stream/"
MAX_BYTES=6_000_000
EXPECTED_VERSION_NUMBER=10


def get_json(url,cap=1_000_000):
    req=Request(url,headers={"Accept":"application/json",
                             "User-Agent":"chocho-monarch-original-published-metadata-gate/0.2"})
    with urlopen(req,timeout=65) as r:
        data=r.read(cap+1)
        if r.status!=200 or len(data)>cap:
            raise ValueError("Dryad metadata/versions exceed documented JSON bounds")
    return json.loads(data.decode("utf-8")), {
        "url":url,"sha256":hashlib.sha256(data).hexdigest(),
        "byte_size":len(data),"http_status":200,
    }


def recursive_objects(data):
    if isinstance(data,dict):
        yield data
        for v in data.values():
            yield from recursive_objects(v)
    elif isinstance(data,list):
        for v in data:
            yield from recursive_objects(v)


def object_id(d,resource="versions"):
    val=d.get("id")
    if isinstance(val,int):
        return val
    for v in (d.get("_links") or {}).values():
        if isinstance(v,dict):
            href=str(v.get("href") or "")
            m=re.search(rf"/{resource}/(\d+)",href)
            if m:
                return int(m.group(1))
    return None


def identify_current_version(json_obj):
    candidates=[]
    for obj in recursive_objects(json_obj):
        v=obj.get("versionNumber",obj.get("version"))
        i=object_id(obj,"versions")
        if i and str(v).isdigit():
            candidates.append({"versionNumber":int(v),"id":i,
                               "status":obj.get("versionStatus",obj.get("curationStatus")),
                               "publicationDate":obj.get("publicationDate")})
    unique={(x["id"],x["versionNumber"]):x for x in candidates}
    candidates=list(unique.values())
    matching=[x for x in candidates if x["versionNumber"]==EXPECTED_VERSION_NUMBER]
    if len(matching)!=1:
        raise ValueError(f"Version number 10 not uniquely verified: {len(matching)}")
    return matching[0],candidates


def identify_exact_file(obj):
    matches=[]
    for d in recursive_objects(obj):
        name=next((str(d.get(k)) for k in ("path","filename","fileName","name")
                   if isinstance(d.get(k),str)),None)
        if name and Path(name).name==MAIN:
            fid=object_id(d,"files")
            if fid:
                matches.append({"name":Path(name).name,"file_id":fid,
                                "size":d.get("size"),"digest":d.get("digest"),
                                "digestType":d.get("digestType"),
                                "status":d.get("status")})
    uniq={x["file_id"]:x for x in matches}
    if len(uniq)!=1:
        raise ValueError(f"Original Dryad published CSV file id not unique: {len(uniq)}")
    return next(iter(uniq.values()))


def download_web_file(file_id):
    url=PUBLIC+str(file_id)
    req=Request(url,headers={"Accept":"text/csv,application/octet-stream,*/*",
                             "User-Agent":"chocho-original-Dryad-public-website-file-stream/0.2"})
    with urlopen(req,timeout=90) as response:
        raw=response.read(MAX_BYTES+1)
        meta={
            "requested_public_website_url":url,
            "status_code":response.status,
            "content_type":response.headers.get("Content-Type"),
            "resolved_public_link":response.url.split("?")[0],
            "source_response_sha256":hashlib.sha256(raw).hexdigest(),
            "download_bytes":len(raw),
        }
    if len(raw)>MAX_BYTES or response.status!=200:
        raise ValueError("Source public file stream not a bounded public 200 response")
    if raw[:1500].lower().find(b"<html")>=0 or raw[:1500].lower().find(b"<doctype")>=0:
        raise ValueError("Official public Dryad download returned HTML rather than data file")
    if len(raw)<50_000:
        raise ValueError("File stream too small for original 2026 2.39MB scientific survey")
    return raw,meta


def run():
    out={
        "schema":"chocho_dryad_published_monarch_stage_source_v02_public_ui_20261026",
        "dataset_doi":DOI,
        "expected_original_survey_name":MAIN,
        "original_dataset_publication":"2026-10-08",
        "api_zip_previous_result":"HTTP 401 for documented API dataset ZIP; no data acquired",
        "source_transport_change_only":True,
    }
    try:
        meta,metainfo=get_json(BASE+"/datasets/"+ENC)
        out["metadata_source_sha256"]=metainfo["sha256"]
        if str(meta.get("identifier","")).removeprefix("doi:")!=DOI:
            raise ValueError("Original published Dryad DOI did not match")
        if int(meta.get("versionNumber",-1))!=EXPECTED_VERSION_NUMBER:
            raise ValueError("Original published Dryad version changed from 10")
        versions,info=get_json(VERSIONS)
        out["versions_response_sha256"]=info["sha256"]
        current,items=identify_current_version(versions)
        out["published_version"]=current
        out["known_versions"]=[{"v":x["versionNumber"],"id":x["id"]} for x in items]
        files,filesinfo=get_json(BASE+f"/versions/{current['id']}/files")
        out["filelist_response_sha256"]=filesinfo["sha256"]
        target=identify_exact_file(files)
        out["source_file"]=target
        raw,sourceinfo=download_web_file(target["file_id"])
        out["download"]=sourceinfo
        # no derived source numerical biological result before fully parsed exact rows
        out["source_quality"]=inspect_csv(raw)
        out["status"]=out["source_quality"]["status"]
    except (HTTPError,URLError,TimeoutError,OSError,ValueError,UnicodeError,
            json.JSONDecodeError) as e:
        out["status"]="ORIGINAL_PUBLISHED_FILE_UNAVAILABLE_OR_SCHEMA_UNIDENTIFIABLE"
        out["error_type"]=type(e).__name__
        out["http_status"]=getattr(e,"code",None)
        out["error_message"]=(str(e)[:180] if not isinstance(e,HTTPError) else None)
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    result=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "http":result.get("http_status"),
        "error":result.get("error_message"),
        "original_version":result.get("published_version"),
        "file":result.get("source_file"),
        "download":result.get("download"),
        "sampling_units":result.get("source_quality",{}).get("transect_round_sampling_units"),
        "egg_source_units":result.get("source_quality",{}).get("positive_egg_transect_round_plant_units"),
        "late_larvae_source_units":result.get("source_quality",{}).get("positive_late_transect_round_plant_units"),
        "dataset_gate":result.get("source_quality",{}).get("gate_counts"),
    },indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
