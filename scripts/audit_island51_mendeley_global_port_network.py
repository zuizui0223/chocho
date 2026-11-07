#!/usr/bin/env python3
"""Public Mendeley Data DOI 10.17632/vzzy3b9gg4.1 files, folders and small schema probe.

Only the advertised open-source research port-network/model-output data are
inspected. No download of protected files, no bearer tokens, no 51-island
biological outcome fit, and no modeled port network called vessel observations.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import time
from datetime import datetime,timezone
from urllib.parse import urlencode,quote,urlsplit
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import zipfile

ROOT=Path(__file__).resolve().parents[1]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
DATASET_ID="vzzy3b9gg4"
VERSION=1
HOST="https://api.data.mendeley.com"
DOI="10.17632/vzzy3b9gg4.1"
MAX_DOWNLOAD=120_000_000
MAX_TOTAL_DOWNLOAD=220_000_000
MAX_ZIP_MEMBER=350_000_000
MAX_ZIP_SAMPLE=8000

def git_blob(path):
    content=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()

def source_roster():
    if git_blob(WPI)!=WPI_GIT_BLOB:
        raise ValueError("Original 194 WPI port source Git blob changed")
    with WPI.open(encoding="utf-8",newline="") as f:
        data=list(csv.DictReader(f))
    if len(data)!=194 or len({r["island_id"] for r in data})!=51 or sum(r["included"].lower()=="true" for r in data)!=182:
        raise ValueError("Original 51-island/194-port/182-terminal cohort no longer fixed")
    return data

def request_bytes(url,allow_failure=True,maxbytes=5_000_000):
    last=None
    for attempt in range(3):
        try:
            req=Request(url,headers={"User-Agent":"chocho-original51-port-maritime-network-source-audit/1.0",
                                    "Accept":"application/json,application/octet-stream,*/*"})
            with urlopen(req,timeout=30) as response:
                status=getattr(response,"status",200)
                typ=response.headers.get("Content-Type","")
                length=response.headers.get("Content-Length","")
                if length and int(length)>maxbytes:
                    raise ValueError("Public resource is too large for source-only inspection")
                content=response.read(maxbytes+1)
                if len(content)>maxbytes:
                    raise ValueError("Read exceeded source-only inspection bytes")
                return {"status":status,"content_type":typ,"raw":content,
                        "redirected_host":urlsplit(response.geturl()).hostname,
                        "bytes":len(content)},None
        except HTTPError as exc:
            last={"http_status":exc.code,"kind":"HTTPError",
                  "reason":str(exc.reason)[:140],
                  "response_excerpt":exc.read(250).decode("utf-8",errors="replace")}
            if exc.code not in (408,429,500,502,503,504):
                break
        except (URLError,TimeoutError,ConnectionError,ValueError) as exc:
            last={"kind":type(exc).__name__,"reason":str(exc)[:180]}
            if isinstance(exc,ValueError):break
        if attempt<2:time.sleep(2*(attempt+1))
    return None,last

def get_json(path):
    url=HOST+path
    response,err=request_bytes(url)
    entry={"requested_path":path,"status":response["status"] if response else None,
           "downloaded_bytes":response["bytes"] if response else 0,
           "content_type":response["content_type"] if response else None,
           "failure":err}
    if response:
        try:
            data=json.loads(response["raw"])
            entry["valid_json"]=True
            entry["top_level_type"]=type(data).__name__
            entry["top_level_keys"]=list(data)[:30] if isinstance(data,dict) else []
            return data,entry
        except (UnicodeDecodeError,ValueError):
            entry["valid_json"]=False
    return None,entry

def normalize_lists(value):
    if isinstance(value,list):return value
    if isinstance(value,dict):
        for key in ("results","files","items","data","folders"):
            if isinstance(value.get(key),list):return value[key]
    return []

def records(files,source_label):
    out=[]
    for f in normalize_lists(files):
        if not isinstance(f,dict):continue
        info=f.get("content_details") or {}
        size=f.get("size",info.get("size",0))
        try:size=int(size or 0)
        except (TypeError,ValueError):size=0
        name=str(f.get("filename") or f.get("name") or "")
        out.append({
            "origin":source_label,"file_id":str(f.get("id") or ""),
            "name":name,"folder_id":str(f.get("folder_id") or ""),
            "bytes":size,
            "content_type":str(info.get("content_type") or f.get("content_type") or ""),
            "sha256_published":str(info.get("sha256_hash") or f.get("sha256_hash") or ""),
            "download_url_AVAILABLE":bool(info.get("download_url")),
            "public_download_url":info.get("download_url") or None
        })
    return out

def sample_binary(filename,raw):
    result={"sniffed_bytes":len(raw)}
    if zipfile.is_zipfile(io.BytesIO(raw)):
        result["format"]="ZIP"
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            result["members"]=[{"name":f.filename,"size":f.file_size,
                                "compressed_size":f.compress_size}
                               for f in z.infolist()[:60]]
            text=[f for f in z.infolist() if f.file_size<MAX_ZIP_MEMBER and
                  f.filename.lower().endswith((".csv",".tsv",".txt",".jsonl",".json"))]
            if text:
                with z.open(text[0]) as b:
                    result["first_member_sample"]=b.read(MAX_ZIP_SAMPLE).decode("utf-8-sig",errors="replace").splitlines()[:4]
        return result
    if raw[:4]==b"PAR1":
        return dict(result,format="PARQUET",magic="PAR1",schema_uninspected_without_optional_dependency=True)
    sample=raw[:MAX_ZIP_SAMPLE].decode("utf-8-sig",errors="replace")
    result["format"]="TEXT_CSV_OR_OTHER"
    result["first_lines"]=sample.splitlines()[:5]
    return result

def safe_file_probe(public,full):
    name=public["name"]
    return bool(re.search(r"(^|[^a-z])(port[_ -]?to[_ -]?port|maritime[_ -]?network|network|ports|port[_ -]?stat)",name.lower())) and
            full["bytes"]>0 and full["bytes"]<=MAX_DOWNLOAD)

def main(outdir):
    out=Path(outdir).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("data","results","provenance")):
        raise ValueError("Never overwrite original source data and scientific results")
    orig=source_roster()
    out.mkdir(parents=True,exist_ok=True)
    queries=[
      "/datasets/"+DATASET_ID+"?version=1",
      "/datasets/"+DATASET_ID+"/files?version=1&%24start=0&%24limit=100",
      "/datasets/publics/"+DATASET_ID+"/files?version=1&%24start=0&%24limit=100",
      "/datasets/"+DATASET_ID+"/folders?version=1",
      "/datasets/"+DATASET_ID+"/snapshot?version=1",
    ]
    manifests=[];meta={}
    file_list=[]
    for path in queries:
        data,receipt=get_json(path)
        manifests.append(receipt)
        meta[path]=data
        if data is not None:
            file_list.extend(records(data,path))
        print("MENDELEY_PUBLIC_ENDPOINT="+json.dumps({key:val for key,val in receipt.items() if key!="response_excerpt"},sort_keys=True),flush=True)
    first=meta[queries[0]]
    folders=[]
    for obj in (first,meta.get(queries[3])):
        if not isinstance(obj,dict):continue
        for fld in normalize_lists(obj.get("folders") or obj):
            if not isinstance(fld,dict):continue
            folders.append({"id":str(fld.get("id") or ""),"name":str(fld.get("name") or ""),
                            "parent_id":str(fld.get("parent_id") or "")})
    if isinstance(first,dict):
        file_list.extend(records(first.get("files"),"dataset.files"))
    # Folder index is returned only if endpoint allows an anonymous request.
    # Some Mendeley dataset files are in folder subtrees.
    folder_ids=sorted(set(x["id"] for x in folders if x["id"]))
    for folder_id in folder_ids[:30]:
        path="/datasets/"+DATASET_ID+"/files?version=1&folder_id="+quote(folder_id,safe="")
        data,receipt=get_json(path)
        manifests.append(receipt)
        if data is not None:
            file_list.extend(records(data,path))
    unique={}
    for file in file_list:
        key=file["file_id"] or file["name"]
        if key and key not in unique:
            unique[key]=file
    listed=sorted(unique.values(),key=lambda x:(x["name"].lower(),x["file_id"]))
    parent_names=[]
    if isinstance(first,dict):
        for f in ("name","description","version","id","licence","license"):
            if f in first:parent_names.append({f:str(first[f])[:2000]})
    samples=[]
    total=0
    candidates=[v for v in listed if safe_file_probe(v,v)]
    for item in candidates[:3]:
        if total+item["bytes"]>MAX_TOTAL_DOWNLOAD:break
        fid=item["file_id"]
        paths=[
          "/datasets/"+DATASET_ID+"/files/"+fid+"/file_downloaded?version=1",
          "/datasets/publics/"+DATASET_ID+"/files/"+fid+"?version=1",
        ] if fid else []
        response=None
        attempts=[]
        for path in paths:
            val,err=request_bytes(HOST+path,maxbytes=MAX_DOWNLOAD)
            attempts.append({"api_path":path,"status":val["status"] if val else None,"error":err})
            if val and val["raw"]:
                response=val
                break
        if response is None:
            samples.append({"name":item["name"],"file_id":fid,
                            "download_status":"NOT_ACCESSIBLE_PUBLICLY_IN_THIS_WORKFLOW",
                            "attempts":attempts})
            continue
        data=response["raw"]
        digest=hashlib.sha256(data).hexdigest()
        published=item["sha256_published"]
        if published and digest.lower()!=published.lower():
            raise ValueError("Original public downloaded model output sha256 does not match metadata for "+item["name"])
        total+=len(data)
        schema=sample_binary(item["name"],data)
        # No raw third-party file published in chocho; only a tiny schema receipt.
        samples.append({"name":item["name"],"file_id":fid,"download_status":"READABLE_VERIFIED_PUBLIC_FILE",
                        "bytes":len(data),"sha256":digest,
                        "publisher_sha256":published,"schema":schema,
                        "attempts":attempts})
    receipt={
       "schema":"chocho_mendeley_global1300_port_network_public_metadata_source_probe_v1",
       "status":"PUBLIC_SOURCE_AVAILABILITY_AND_SCHEMA_INSPECTION_COMPLETE",
       "reviewed_utc":datetime.now(timezone.utc).isoformat(),
       "article_doi":"10.1038/s41467-022-32070-0",
       "dataset_doi":DOI,
       "dataset_url":"https://data.mendeley.com/datasets/vzzy3b9gg4/1",
       "dataset_id":DATASET_ID,"dataset_version":VERSION,
       "source_license_from_public_page":"CC BY 4.0",
       "source_is_model_outputs_not_direct_observed_AIS_ship_voyages":True,
       "source_original_WPI_git_blob":WPI_GIT_BLOB,
       "original_WPI_ports":len(orig),"original_source_islands":51,
       "original_marine_eligible_ports":182,
       "public_api_endpoint_results":manifests,
       "dataset_metadata_fields":parent_names,
       "folder_index":folders,
       "public_file_metadata":[{k:v for k,v in x.items() if k!="public_download_url"} for x in listed],
       "source_file_count":len(listed),
       "selected_candidate_file_count":len(candidates),
       "selected_source_file_download_samples":samples,
       "all_public_file_urls_or_signed_download_tokens_omitted":True,
       "no_original_port_match_or_butterfly_model_performed":True,
       "model_network_is_not_confirmed_observed_ship_journeys":True,
    }
    (out/"source_metadata_and_schema_receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("MENDELEY_GLOBAL_PORT_NETWORK_SOURCE_JSON="+json.dumps({
      "status":receipt["status"],"dataset_id":DATASET_ID,
      "public_api_query_statuses":[{"path":x["requested_path"],
             "status":x["status"],"failure":x.get("failure")} for x in manifests[:5]],
      "folder_count":len(folders),"file_count":len(listed),
      "network_candidate_files":[x["name"] for x in candidates],
      "downloaded_verified_files":[{k:v for k,v in x.items()
          if k in ("name","download_status","bytes","sha256")}
          for x in samples],
      "original_WPI_ports":194,"original_islands":51},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",required=True,type=Path)
    main(ap.parse_args().out_dir)
