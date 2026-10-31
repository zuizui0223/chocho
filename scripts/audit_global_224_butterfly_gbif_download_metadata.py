#!/usr/bin/env python3
"""Check whether the already-cited ~14.6M butterfly GBIF download is retrievable.

DOI 10.15468/dl.umqmsf is a historical global occurrence download cited by
chocho's island/host forecast work. This probe GETs *metadata only* and stores
only public provenance, sizes/status/availability. No records or coordinates
are accessed. It never guesses a download key, generates a fresh download,
or claims DOI contents can be equated with the GBIF 31-species pilot.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

DOI="10.15468/dl.umqmsf"
ENDPOINT="https://api.gbif.org/v1/occurrence/download/"+DOI
OUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/"
         "gbif_global_224_original_download_source_gate.json")

def main():
    req=Request(ENDPOINT,headers={"Accept":"application/json",
         "User-Agent":"chocho-source-locked-global-butterfly-macroecology/1.0"})
    try:
        with urlopen(req,timeout=30) as response:
            payload=json.loads(response.read(500000))
        if not isinstance(payload,dict):
            raise RuntimeError("Expected public GBIF occurrence download JSON")
        doi=str(payload.get("doi") or "")
        if DOI not in doi:
            raise RuntimeError("GBIF resolved to a different DOI: "+doi)
        status=str(payload.get("status") or "")
        key=str(payload.get("key") or "")
        size=payload.get("size")
        n=payload.get("totalRecords")
        link=str(payload.get("downloadLink") or "")
        if status!="SUCCEEDED":
            state="KNOWN_DOWNLOAD_NOT_SUCCEEDED"
        elif not key or size is None or n is None or not link:
            state="INCOMPLETE_PUBLIC_DOWNLOAD_METADATA"
        else:
            state="PUBLIC_GLOBAL_BUTTERFLY_SOURCE_METADATA_VERIFIED"
        receipt={
            "status":state,"source_doi":DOI,"source_metadata_endpoint":ENDPOINT,
            "gbif_public_download_key":key,"archival_status":status,
            "reported_total_records":n,"archive_size_bytes":size,
            "download_link":link,"archive_erase_after":payload.get("eraseAfter"),
            "original_created":payload.get("created"),
            "last_modified":payload.get("modified"),
            "format":payload.get("format"),
            "source_claim":"Original historical global butterfly/island source cited in chocho work, not a re-verified taxon list",
            "limits":[
                "Metadata record count is not a validated set of 224 taxonomically matched butterfly species.",
                "Archive availability requires direct 200 response/bytes check; link alone is not evidence it remains downloadable.",
                "No geographic, year or provider-date completeness inferred from metadata.",
                "This DOI is different from 10.15468/dl.pp5nc9, the 30 matched species pilot.",
                "Mass-download source cannot be automatically assumed to reproduce earlier local 14,637,152 coordinates.",
            ],
        }
    except (HTTPError,URLError,ValueError,TimeoutError) as e:
        receipt={
            "status":"GLOBAL_SOURCE_METADATA_UNAVAILABLE_NOT_NEGATIVE_BIOLOGICAL_RESULT",
            "source_doi":DOI,"source_metadata_endpoint":ENDPOINT,
            "error":type(e).__name__+": "+str(e)[:300],
        }
    receipt["checked_at_utc"]=datetime.datetime.now(
        datetime.timezone.utc).isoformat()
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(receipt,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
