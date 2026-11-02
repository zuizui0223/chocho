#!/usr/bin/env python3
"""Discover reusable historical empire source metadata; never infer from text guesses."""
from __future__ import annotations
import json
import urllib.request
from pathlib import Path

URL="https://zenodo.org/api/records/6640264"
def main():
    req=urllib.request.Request(URL,headers={"User-Agent":"chocho-colonial-resource-study/0.1 scientific-metadata-audit","Accept":"application/json"})
    try:
        with urllib.request.urlopen(req,timeout=90) as f:
            root=json.loads(f.read())
        files=root.get("files",[])
        data={
          "source_doi":"10.5281/zenodo.6640264",
          "source_id":root.get("id"),
          "title":root.get("metadata",{}).get("title"),
          "created":root.get("created"),
          "files": [{
               "key": f.get("key"),
               "size":f.get("size"),
               "checksum":f.get("checksum"),
               "download_url":f.get("links",{}).get("self"),
          } for f in files],
          "status":"SUCCESS",
        }
    except Exception as e:
        data={"status":"FETCH_FAILED","source_doi":"10.5281/zenodo.6640264",
              "error":str(e)}
    out=Path("provenance/reviewer_defenses/results/colonial_host_source_20261102")
    out.mkdir(parents=True,exist_ok=True)
    (out/"zenodo_source_inventory.json").write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(data,indent=2,ensure_ascii=False))
    # Successful source discovery only; if metadata endpoint is unavailable, no synthetic empire assignments.
    if data["status"]!="SUCCESS":raise RuntimeError(data["error"])
if __name__=="__main__":
    main()
