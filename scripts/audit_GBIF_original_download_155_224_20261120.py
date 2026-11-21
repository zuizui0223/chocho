#!/usr/bin/env python3
"""Read-only verify original GBIF download citation and its contributing datasets.

Source: immutable 2026-10-29 host-forecast README DOI 10.15468/dl.umqmsf.
The independent 155 forecast and 224 observed-community model can cite this
same original GBIF download, but this program must NOT claim the private 224
historical aggregate, the 225x51 facet requery, or 1596 port plants used it.

This script never creates or requests a new GBIF download or derived DOI.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import time
from pathlib import Path
from urllib.request import urlopen,Request
from urllib.error import HTTPError, URLError

DOI = "10.15468/dl.umqmsf"
API="https://api.gbif.org/v1/occurrence/download"
UUID=re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.I)

def get_json(url):
    error=None
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers={
                "Accept":"application/json",
                "User-Agent":"chocho-GEB-v05-readonly-citation-audit/2026-11-20"
            }),timeout=35) as response:
                if response.status!=200:
                    raise RuntimeError(f"Status {response.status} for official GBIF endpoint")
                raw=response.read()
            return json.loads(raw),hashlib.sha256(raw).hexdigest()
        except (HTTPError,URLError,TimeoutError,ValueError,OSError) as ex:
            error=ex
            if attempt<3 and not isinstance(ex,HTTPError) or (
                isinstance(ex,HTTPError) and ex.code in (408,429,500,502,503,504)):
                time.sleep(min(7,2**attempt))
                continue
            break
    raise RuntimeError(f"Official GBIF citation endpoint unavailable: {type(error).__name__} {str(error)[:200]}")

def providers_from(data):
    if isinstance(data,list):
        return data,True
    if isinstance(data,dict):
        for k in ("results","datasets","data"):
            if isinstance(data.get(k),list):
                complete=data.get("endOfRecords",not data.get("next"))
                return data[k],bool(complete)
    raise RuntimeError("Unexpected GBIF contributing-datasets API schema, fail-closed")

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    meta,meta_sha=get_json(API+"/"+DOI)
    if not isinstance(meta,dict):
        raise RuntimeError("Official GBIF download DOI metadata not JSON object")
    seen_doi=str(meta.get("doi","")).replace("https://doi.org/","")
    if seen_doi.lower()!=DOI:
        raise RuntimeError(f"GBIF returned unrelated download DOI: {seen_doi}")
    declared=meta.get("totalRecords") or meta.get("total_records") or meta.get("totalCount")
    if declared is not None and int(declared)!=14637152:
        raise RuntimeError(f"Original 14637152-record DOI and metadata count inconsistent: {declared}")
    # The official contributors endpoint is paginated; never mistake page 1
    # for the complete set of GBIF data publishers.
    items=[]
    provider_page_shas=[]
    offset=0
    pages=0
    for _ in range(150):
        dataset_page,page_sha=get_json(API+"/"+DOI+"/datasets?offset="+str(offset)+"&limit=300")
        page_rows,complete=providers_from(dataset_page)
        provider_page_shas.append(page_sha)
        items.extend(page_rows)
        pages+=1
        print("ORIGINAL_GBIF_PUBLISHER_PAGE",json.dumps({
            "page":pages,"offset":offset,"returned":len(page_rows),"complete":complete
        }),flush=True)
        if complete:
            break
        if not page_rows:
            raise RuntimeError("Original GBIF contributor paging returned zero rows before end")
        offset += len(page_rows)
    else:
        raise RuntimeError("Too many GBIF contribution pages without a complete response")
    provider_sha=hashlib.sha256("".join(provider_page_shas).encode()).hexdigest()
    attribution=[]
    for x in items:
        if not isinstance(x,dict):
            raise RuntimeError("Unexpected dataset contribution entry type")
        key=x.get("datasetKey") or x.get("key") or (
            x.get("dataset",{}).get("key") if isinstance(x.get("dataset"),dict) else None)
        count=x.get("recordCount",x.get("numberRecords",x.get("count")))
        if key is None or not UUID.fullmatch(str(key)):
            raise RuntimeError("Original download provider entry lacks valid datasetKey")
        if count is not None and (not str(count).isdigit() or int(count)<0):
            raise RuntimeError("Original GBIF dataset contribution has invalid count")
        attribution.append({"datasetKey":str(key),"recordCount":None if count is None else int(count)})
    if not attribution or len({x["datasetKey"] for x in attribution})!=len(attribution):
        raise RuntimeError("Empty or duplicate original GBIF download provider roster")
    totals=[x["recordCount"] for x in attribution if x["recordCount"] is not None]
    contribution_sum=sum(totals) if len(totals)==len(attribution) else None
    if contribution_sum is not None and declared is not None and contribution_sum !=int(declared):
        raise RuntimeError("Per-provider counts do not reconcile to source download; do not miscite counts")
    result={
      "schema":"chocho_original_GBIF155_224_download_publisher_attribution_v1",
      "status":"OFFICIAL_GBIF_DOWNLOAD_METADATA_AND_CONTRIBUTOR_LIST_VERIFIED_READONLY",
      "original_DOI":DOI,
      "original_record_count_verified":int(declared) if declared is not None else None,
      "same_original_download_used_by":[
          "155 source-eligible butterfly retrospective recording forecast",
          "224 indexed butterfly observed regional Jaccard model"],
      "original_source_years":"2001-2024",
      "dataset_contributors":attribution,
      "provider_count":len(attribution),
      "contributors_record_sum_if_available":contribution_sum,
      "official_metadata_response_sha256":meta_sha,
      "official_dataset_pages":pages,
      "official_dataset_page_sha256s":provider_page_shas,
      "official_dataset_pages_digest":provider_sha,
      "read_utc_note":"Retrieval time is execution time, not a new GBIF data download",
      "DOI_not_valid_for":[
          "historical private 224 species×landmass×year aggregated 14425536 records",
          "independent November 9 2026 225×51 speciesKey GBIF API facet reindex",
          "1596 WPI botanical port speciesKey GBIF API facets"],
      "no_GBIF_download_or_derived_DOI_created":True
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("GBIF_ORIGINAL155_224_SHARED_DOI_AND_PUBLISHER_DATASETS_VERIFIED",
          json.dumps({"DOI":DOI,"download_records":result["original_record_count_verified"],
                      "provider_count":len(attribution),"contributions_sum":contribution_sum}))
if __name__=="__main__":
    main()
