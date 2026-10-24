#!/usr/bin/env python3
"""Source metadata-only audit of independent US butterfly survey counts.

No site-level records downloaded at this stage, no ecological model run.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
API="https://api.figshare.com/v2/articles/27934602"
EXPECTED_ARTICLE_ID=27934602
EXPECTED_DOI="10.6084/m9.figshare.27934602"
ORIGINAL_FOCAL_GIT_BLOB="9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"


def get_meta():
    request=Request(API,headers={
        "Accept":"application/json",
        "User-Agent":"chocho-science-public-data-source-gate/1.0",
    })
    status={"status":"UNVERIFIED","url":API}
    try:
        with urlopen(request,timeout=40) as response:
            data=response.read(2_000_001)
            status["http_status"]=response.status
            status["content_type"]=response.headers.get("Content-Type","")
            if len(data)>2_000_000:
                status["status"]="METADATA_RESPONSE_TOO_LARGE"
                return status
        status["meta_sha256"]=hashlib.sha256(data).hexdigest()
        source=json.loads(data.decode("utf-8-sig"))
        if not isinstance(source,dict) or source.get("id")!=EXPECTED_ARTICLE_ID:
            status["status"]="SOURCE_ID_MISMATCH"
            return status
        fields=source.get("files",[])
        if not isinstance(fields,list):
            status["status"]="FILES_LIST_NOT_PRESENT"
            return status
        status.update({
            "status":"SOURCE_METADATA_FOUND",
            "id":source.get("id"),
            "title":source.get("title"),
            "doi":source.get("doi"),
            "version":source.get("version"),
            "published_date":source.get("published_date"),
            "modified_date":source.get("modified_date"),
            "license":source.get("license"),
            "is_public":source.get("is_public"),
            "files":[{
                "id":f.get("id"),
                "name":f.get("name"),
                "size":f.get("size"),
                "computed_md5":f.get("computed_md5"),
                "download_url":f.get("download_url"),
            }for f in fields],
            "total_download_bytes_if_all":sum(int(f.get("size") or 0) for f in fields),
            "file_count":len(fields),
        })
    except HTTPError as ex:
        status.update(status="HTTP_ERROR",http_status=ex.code)
    except (URLError,TimeoutError,OSError,UnicodeError,ValueError) as ex:
        status.update(status="SOURCE_ERROR",exception_class=type(ex).__name__)
    return status


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    source=ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    raw=source.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()
    if blob!=ORIGINAL_FOCAL_GIT_BLOB:
        raise ValueError("Focal 239 input is not source-locked")
    with source.open(encoding="utf-8",newline="") as h:
        names=[r["species"] for r in csv.DictReader(h)]
    if len(names)!=239 or len(set(names))!=239:
        raise ValueError("239 original butterflies differ")
    meta=get_meta()
    result={
        "schema":"chocho_sbus_figshare_count_source_gate_20261025_v01",
        "source":meta,
        "focal_species_n":len(names),
        "focal_source_sha":blob,
        "original_publication":"Henry et al (2025) Scientific Data DOI 10.1038/s41597-025-05513-8",
        "no_survey_rows_downloaded":True,
        "no_demographic_or_competition_effect_estimated":True,
        "status":(
            "READY_FOR_SIZE_BOUNDED_SOURCE_DOWNLOAD_AND_SCHEMA_QC"
            if meta.get("status")=="SOURCE_METADATA_FOUND"
            else "SOURCE_METADATA_UNAVAILABLE_DO_NOT_INFER"
        ),
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],"source_article":meta.get("id"),
        "original_doi":meta.get("doi"),
        "number_of_files":meta.get("file_count"),
        "file_sizes": [(z["name"],z["size"]) for z in meta.get("files",[])],
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
