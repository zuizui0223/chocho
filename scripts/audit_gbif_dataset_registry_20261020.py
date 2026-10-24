#!/usr/bin/env python3
"""Posthoc GBIF dataset registry lookup for the six source-locked dataset keys.

Registry titles/publishers are provenance metadata, NOT deduplicated butterfly
observations, breeding evidence or independent biological replication.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

API="https://api.gbif.org/v1/dataset/"
KEY_RE = re.compile(r"^[a-f0-9]{8}(?:-[a-f0-9]{4}){3}-[a-f0-9]{12}$")


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def candidate_keys(facet_qc):
    doc=json.loads(Path(facet_qc).read_text(encoding="utf-8"))
    if doc.get("live_facets_complete")!=7 or len(doc.get("results",[]))!=7:
        raise ValueError("Seven complete preceding source facets required")
    keys=set()
    for row in doc["results"]:
        if row.get("status")!="COMPLETE_LIVE_FACET_QC" or row.get("source_drift_records")!=0:
            raise ValueError("Unsynchronized original/candidate GBIF source")
        keys.update(row["facet_counts"]["DATASETKEY"])
    if len(keys)!=6 or any(not KEY_RE.fullmatch(k) for k in keys):
        raise ValueError("Registry source-key cohort changed")
    return sorted(keys)


def parse_metadata(raw, requested_key):
    try: doc=json.loads(raw)
    except (ValueError, UnicodeDecodeError) as exc:
        raise ValueError("Truncated/invalid dataset registry response") from exc
    if doc.get("key")!=requested_key:
        raise ValueError("GBIF dataset key identity mismatch")
    title=doc.get("title")
    if not isinstance(title,str) or not title.strip():
        raise ValueError("GBIF registry missing dataset title")
    return {
        "dataset_key": requested_key,
        "title": title[:400],
        "publishing_organization_key": doc.get("publishingOrganizationKey"),
        "type": doc.get("type"),
        "doi": doc.get("doi"),
        "license": doc.get("license"),
        "homepage": doc.get("homepage"),
        "created": doc.get("created"),
        "modified": doc.get("modified"),
    }


def fetch_one(key,max_attempts=3):
    url=API+key
    errors=[]
    for attempt in range(max_attempts):
        try:
            req=Request(url,headers={"User-Agent":"chocho-country-dataset-provenance/1.0","Accept":"application/json"})
            with urlopen(req,timeout=40) as response:
                status=response.status
                raw=response.read(500000)
            if status!=200 or len(raw)>=500000:
                raise ValueError("Truncated or non-200 dataset metadata")
            item=parse_metadata(raw,key)
            return {"status":"COMPLETE_REGISTRY_METADATA","url":url,
                    "http_status":status,"retrieved_utc":utcnow(),
                    "raw_sha256":hashlib.sha256(raw).hexdigest(),
                    "attempts":attempt+1,**item}
        except (HTTPError,URLError,TimeoutError,OSError,ValueError) as exc:
            errors.append({"attempt":attempt+1,"kind":type(exc).__name__,"detail":str(exc)[:350]})
            if attempt<max_attempts-1:time.sleep(2*(attempt+1))
    return {"status":"REGISTRY_UNAVAILABLE_NO_INDEPENDENCE_INFERENCE",
            "dataset_key":key,"url":url,"retrieved_utc":utcnow(),"errors":errors}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--facet-audit",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    keys=candidate_keys(a.facet_audit)
    rows=[]
    for key in keys:
        rows.append(fetch_one(key))
        time.sleep(1)
    result={"schema":"gbif_country_dataset_registry_followup_posthoc_v1",
            "created_utc":utcnow(),"source_qc_sha256":hashlib.sha256(a.facet_audit.read_bytes()).hexdigest(),
            "n_selected_keys":len(keys),
            "n_complete":sum(r["status"]=="COMPLETE_REGISTRY_METADATA" for r in rows),
            "not_record_independence":True,
            "dataset_registry":rows}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("GBIF_DATASET_REGISTRY_SUMMARY="+json.dumps(result,separators=(",",":")),flush=True)
    if result["n_complete"]!=len(keys):
        raise SystemExit("Metadata unavailable for one or more datasets; do not infer missing")
if __name__=="__main__": main()
