#!/usr/bin/env python3
"""Check whether already cited global 205-butterfly GBIF download can be reused.

The published DOI https://doi.org/10.15468/dl.umqmsf underlies the previous
global-original butterfly index. A public download if still available is
preferable to burdening the GBIF API with 11,088 independent island queries.
No auth needed for a previously created, publicly retained download. This
metadata check makes no requests to create a new download or access private
occurrence files.
"""
from __future__ import annotations
import json,urllib.request,urllib.error
from datetime import datetime,timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
DOI="10.15468/dl.umqmsf"
QUERIES=[
 "https://api.gbif.org/v1/occurrence/download/"+DOI,
 "https://api.gbif.org/v1/occurrence/download/"+DOI.replace("/","%2F"),
 "https://www.gbif.org/occurrence/download/"+DOI,
]
def source(url):
    try:
        req=urllib.request.Request(url,headers={"Accept":"application/json",
            "User-Agent":"chocho-global-butterfly-original-GBIF-source-metadata-audit/20261031"})
        with urllib.request.urlopen(req,timeout=45) as r:
            head=r.read(350_000)
            response_url=r.url
            content_type=r.headers.get("Content-Type","")
            status=r.status
        if "json" in content_type:
            j=json.loads(head)
            return {"status":"RESPONDED_JSON","request_url":url,
                "response_url":response_url,"http":status,
                "metadata":{k:j.get(k) for k in (
                    "key","doi","status","size","totalRecords","downloadLink","format",
                    "created","modified","numberDatasets","source","eraseAfter")},
                "metadata_fields":sorted(j)[:40]}
        return {"status":"RESPONDED_NONJSON","request_url":url,
                "response_url":response_url,"http":status,"content_type":content_type,
                "content_excerpt":head[:220].decode("utf8",errors="replace")}
    except Exception as e:
        return {"status":"SOURCE_LOOKUP_FAILED","request_url":url,
                "error_type":type(e).__name__,"error":str(e)[:300]}
def main():
    result={"schema":"chocho_published_global_GBIF_butterfly_download_reavailability_v1",
            "citation_doi":"https://doi.org/"+DOI,
            "checked_at_utc":datetime.now(timezone.utc).isoformat(),
            "responses":[source(u) for u in QUERIES],
            "no_new_download_requested":True,
            "not_original_14million_source_content_verified":True,
            "not_equivalent_to_private_original_record_year_geography":True}
    (HERE/"original_global_GBIF_reuse_gate.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()
