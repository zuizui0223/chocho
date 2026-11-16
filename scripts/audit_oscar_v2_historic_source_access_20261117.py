#!/usr/bin/env python3
"""Anonymous source-access preflight for historical NASA OSCAR v2.0.

This reads PUBLIC CMR metadata, then attempts a limited byte-range probe of
the official data link. No Earthdata credentials, large downloads, velocities,
trajectories, or ecological inferences are involved.
"""
from __future__ import annotations

import argparse
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

COLLECTION = "C2098858642-POCLOUD"
SHORT_NAME = "OSCAR_L4_OC_FINAL_V2.0"
DATE = "2005-07-15"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"

def fetch_public_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={
        "User-Agent": "chocho-oscar-v2-source-access/1.0 (research)",
        "Accept": "application/json",
    })
    with urllib.request.urlopen(req, timeout=45) as response:
        raw = response.read(4_000_000)
        if len(raw) >= 4_000_000:
            raise RuntimeError("Unexpectedly large CMR metadata response")
        return json.loads(raw)

def safe_host(url: str) -> str:
    return urllib.parse.urlsplit(url).hostname or ""

def request_data_range(url: str) -> dict:
    """Only read first 16 bytes; never expose signed URLs or user credentials."""
    if not url.startswith("https://"):
        return {"status":"not_https","source_host":safe_host(url)}
    req = urllib.request.Request(url, headers={
        "User-Agent": "chocho-oscar-v2-source-access/1.0 (research)",
        "Range": "bytes=0-15",
        "Accept": "*/*"
    })
    try:
        with urllib.request.urlopen(req, timeout=45) as response:
            first = response.read(16)
            mime = response.headers.get("Content-Type", "")
            final_host = safe_host(response.geturl())
            magic = first.startswith(b"\x89HDF\r\n\x1a\n") or first.startswith(b"CDF")
            return {
                "status": "netcdf_magic" if magic else "not_netcdf_magic",
                "http_code": response.status, "content_type": mime,
                "source_host": safe_host(url), "final_host": final_host,
                "first_16_bytes_hex": first.hex(), "netcdf_magic_seen": magic,
                "range_header": response.headers.get("Content-Range"),
            }
    except urllib.error.HTTPError as e:
        return {
            "status": "http_error", "http_code": e.code,
            "source_host": safe_host(url), "final_host":safe_host(e.geturl()),
            "content_type":e.headers.get("Content-Type", ""),
            "netcdf_magic_seen":False
        }
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return {
            "status": "network_error", "error_type":type(e).__name__,
            "source_host":safe_host(url), "netcdf_magic_seen":False
        }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    query=urllib.parse.urlencode({
        "collection_concept_id": COLLECTION,
        "temporal": f"{DATE}T00:00:00Z,{DATE}T23:59:59Z",
        "page_size": 2,
    })
    meta=fetch_public_json(f"{CMR}?{query}")
    entries=meta.get("feed", {}).get("entry", [])
    data=[]
    for ent in entries:
        links=ent.get("links", [])
        real=[x for x in links if not x.get("inherited") and
              "data#" in str(x.get("rel","")).lower() and
              "s3#" not in str(x.get("rel","")).lower() and
              not x.get("title","").lower().startswith("metadata") and
              str(x.get("href","")).startswith("https://")]
        data.append({
           "granule_id":ent.get("id"), "time_start":ent.get("time_start"),
           "time_end":ent.get("time_end"), "data_link_count":len(real),
           "candidate_hosts":sorted({safe_host(x["href"]) for x in real}),
        })
    first_data_url=None
    for ent in entries:
        for x in ent.get("links", []):
            if not x.get("inherited") and "data#" in str(x.get("rel","")).lower() and
               "s3#" not in str(x.get("rel","")).lower() and
               str(x.get("href","")).startswith("https://"):
               first_data_url=x["href"]
               break
        if first_data_url:break
    probe=request_data_range(first_data_url) if first_data_url else {"status":"no_https_direct_data_link","netcdf_magic_seen":False}
    outcome={
      "status":"PUBLIC_SOURCE_ACCESS_PREFLIGHT_NOT_VELOCITY_ANALYSIS",
      "collection":SHORT_NAME,"collection_concept_id":COLLECTION,
      "source":"NASA PO.DAAC CMR official granule endpoint",
      "test_date":DATE,"interval_of_interest":"1993-01-01 through 2010-12-31",
      "granules_returned":len(entries),"granules":data,
      "anonymous_direct_data_probe":probe,
      "can_read_netcdf_anonymously":bool(probe.get("netcdf_magic_seen")),
      "is_directional_current_model":False,
      "notes":[
          "A CMR metadata match only proves that granule metadata is indexed, not that velocity arrays can be opened.",
          "If the protected PO.DAAC data URL redirects to Earthdata Login, do not bypass authentication.",
          "The CMR response may expose public S3 metadata links that are not the protected NetCDF file.",
          "No sea surface velocity, particle tracking, source-to-island drift times, ship OD, or butterfly outcomes were estimated."
      ],
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(outcome,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("OSCAR_SOURCE_ACCESS_PREFLIGHT",json.dumps({
        "granules":len(entries),"data_link_counts":[x["data_link_count"] for x in data],
        "probe":probe,"can_read_netcdf_anonymously":outcome["can_read_netcdf_anonymously"],
        "status":outcome["status"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
