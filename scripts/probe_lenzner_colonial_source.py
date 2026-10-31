#!/usr/bin/env python3
"""Inspect Lenzner et al. 2022 colonial-flora open data files and spatial keys.

Only fetch Zenodo *metadata*. Publication reports regional colonial occupancy
in alien-plant floras. No historical empire/plant observation is inferred
until region geometry and exact temporal exposure crosswalk are verified.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

SOURCE = "https://zenodo.org/api/records/6640264"
OUT = Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/"
           "lenzner_zenodo_colonial_metadata.json")


def inspect():
    req = Request(SOURCE, headers={"User-Agent":"chocho-research-reproducibility/1.0",
                                  "Accept":"application/json"})
    try:
        with urlopen(req, timeout=30) as reply:
            raw = reply.read(2_000_000)
            payload = json.loads(raw)
        files = payload.get("files") or []
        if isinstance(files, dict):
            files = files.get("entries") or []
        items=[]
        for file in files:
            items.append({
                "key":file.get("key") or file.get("filename"),
                "size":file.get("size"),
                "checksum":file.get("checksum"),
                "content_type":file.get("mimetype"),
                "public_link":(file.get("links") or {}).get("self") or
                              (file.get("links") or {}).get("content"),
            })
        metadata=payload.get("metadata") or {}
        result={
            "status":"SOURCE_METADATA_RETRIEVED_NOT_COLONIAL_CROSSWALK",
            "source_record":SOURCE,
            "source_doi":"10.5281/zenodo.6640264",
            "metadata_version":payload.get("id"),
            "title":metadata.get("title"),
            "date":metadata.get("publication_date"),
            "license":(metadata.get("license") or {}).get("id")
                      if isinstance(metadata.get("license"),dict)
                      else metadata.get("license"),
            "file_count":len(items),
            "files":items,
            "file_total_bytes":sum((x.get("size") or 0) for x in items),
            "limits":[
                "Metadata is not original region-by-empire exposure.",
                "GloNAF regional geometries are not identical to 369 WGSRPD level-3 codes.",
                "Colonial overlap must be historically date-resolved; modern country labels do not identify old island regions.",
                "Plant-flora similarities are already known from the cited study; new evidence must arise from consumer resource network and independent butterflies.",
            ],
        }
    except (HTTPError, URLError, ValueError, TimeoutError) as e:
        result={
            "status":"SOURCE_METADATA_UNAVAILABLE_NOT_NEGATIVE_EVIDENCE",
            "source_record":SOURCE,
            "source_doi":"10.5281/zenodo.6640264",
            "error":f"{type(e).__name__}: {str(e)[:350]}",
        }
    result["checked_at_utc"]=datetime.datetime.now(
        datetime.timezone.utc).isoformat()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                   encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__ == "__main__":
    inspect()
