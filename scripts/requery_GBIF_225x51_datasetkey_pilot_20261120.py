#!/usr/bin/env python3
"""Source-matched, NEW-timestamp GBIF datasetKey-facet probe for two original islands.

Uses original 2026-11-09 GSHHG 2.3.7 WKT fragmentation code and exact archived
source query filters, rather than guessing coastlines. New GBIF index provider
counts are NOT original 2026-11-09 publisher attribution and do not mint DOI.
Run on the pinned 2026 source checkout from an independent audit checkout.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import time
import uuid

ORIGINAL_COMMIT="365526cc70e57b8a6634ee29494bd678f95961ff"
ORIGINAL_BASE=Path("data/exploratory/public51_butterfly_index_20261109")
SOURCE_SCRIPT=Path("scripts/build_public51_butterfly_risk_20261109.py")
PERIODS={"pre2001_2010":(2001,2010),
         "first2011_2017":(2011,2017),
         "holdout2018_2022":(2018,2022)}
DEFAULT_ISLANDS=("9","276")
API="https://api.gbif.org/v1/occurrence/search"

def request_json(params):
    url=API+"?"+urlencode(params,doseq=True)
    last=None
    for attempt in range(4):
        try:
            req=Request(url,headers={"Accept":"application/json",
                     "User-Agent":"chocho-original51-current-index-publisher-probe-20261120"})
            with urlopen(req,timeout=90) as f:
                if f.status!=200:raise RuntimeError("GBIF HTTP "+str(f.status))
                obj=json.load(f)
            if not isinstance(obj,dict) or obj.get("error"):
                raise RuntimeError("Unusable GBIF facet response")
            return obj
        except (HTTPError,URLError,TimeoutError,ValueError,RuntimeError) as ex:
            last=ex
            if isinstance(ex,HTTPError) and ex.code not in (408,429,500,502,503,504):
                break
            if attempt<3:time.sleep(min(10,2**attempt))
    raise RuntimeError("GBIF indexed facet not available for query: "+str(type(last).__name__)+" "+str(last)[:180])

def dataset_facets(params):
    values={}
    count=None
    offset=0
    pages=0
    for _ in range(80):
        response=request_json({**params,"facetOffset":offset})
        total=int(response["count"])
        if count is not None and total!=count:
            raise ValueError("Index count drifted during provider-facet pagination")
        count=total
        facets={str(f.get("field","")).upper():f.get("counts",[])
                for f in response.get("facets",[])}
        if count>0 and "DATASET_KEY" not in facets:
            raise ValueError("GBIF did not return DATASET_KEY facets: "+str(sorted(facets)))
        part=facets.get("DATASET_KEY",[])
        for obj in part:
            key=str(obj.get("name",""))
            uuid.UUID(key)
            n=int(obj["count"])
            if n<0 or key in values:
                raise ValueError("Invalid or repeated datasetKey provider facet")
            values[key]=n
        pages+=1
        if len(part)<1000:
            break
        offset+=len(part)
    else:
        raise RuntimeError("More than 80000 GBIF providers in source-matched island-period")
    if sum(values.values())!=count:
        raise ValueError(f"Provider facets incomplete: sum {sum(values.values())}, indexed records {count}")
    return count,values,pages

def read_original_receipts(source_root):
    root=source_root/ORIGINAL_BASE/"source_shards"
    source={}
    butterfly_filter=None
    for p in sorted(root.glob("batch_*_query_receipt.json")):
        r=json.loads(p.read_text(encoding="utf8"))
        for island in r["islands_tested"]:
            island=str(island)
            if island in source:
                raise ValueError("Original source island appears in more than one shard")
            source[island]={
                "geometry":next(x for x in r["geometry_audit"] if str(x["island_id"])==island),
                "queries":{str(q["period"]):q for q in r["GBIF_index_query_receipts"] if str(q["island_id"])==island}
            }
        this_filter=r["butterfly_taxon_filter"]
        if butterfly_filter is None:
            butterfly_filter=this_filter
        elif butterfly_filter!=this_filter:
            raise ValueError("Inconsistent original GBIF butterfly taxon source filter")
    if len(source)!=51:
        raise ValueError("Original 51 island query receipts incomplete")
    return source,butterfly_filter

def probe(source_root,gshhg_zip,work_dir,islands):
    sys.path.insert(0,str((source_root/"scripts").resolve()))
    import build_public51_butterfly_risk_20261109 as original
    if set(PERIODS)!=set(original.PERIODS):
        raise ValueError("Original GBIF period definitions changed")
    src,filter_used=read_original_receipts(source_root)
    if not set(islands)<=set(src):
        raise ValueError("Requested island not present in original fixed source")
    shape,audits,digest=original.load_polygons(gshhg_zip,work_dir,set(islands),4300)
    if digest!=original.GSHHG_SHA256:
        raise ValueError("Original GSHHG archive source checksum mismatch")
    matched={}
    for audit in audits:
        island=str(audit["island_id"])
        old=src[island]["geometry"]
        if audit["WKT_sha256_fragments"]!=old["WKT_sha256_fragments"]:
            raise ValueError("Original 2026-11-09 island WKT geometry is not exactly reproduced")
        matched[island]=True
    if len(matched)!=len(islands):
        raise ValueError("Missing source-pinned geometry")
    taxon_keys=([filter_used["GBIF_taxonKey"]]
                if filter_used["method"]=="GBIF_SUPERFAMILY_TAXONKEY"
                else list(filter_used["family_taxonkeys"].values()))
    start=datetime.now(timezone.utc).isoformat()
    rows=[]
    for island in islands:
        old=src[island]
        if set(old["queries"])!=set(PERIODS):
            raise ValueError("Missing original island×3-period GBIF source receipt")
        for period,(y0,y1) in PERIODS.items():
            original_receipt=old["queries"][period]
            old_frags=original_receipt["fragment_receipts"]
            if len(old_frags)!=len(shape[island]):
                raise ValueError("Original WKT fragment count/GBIF source receipt mismatch")
            pcounts=Counter()
            new_total=old_total=0
            fragments=[]
            for i,wkt in enumerate(shape[island],start=1):
                f=old_frags[i-1]
                if hashlib.sha256(wkt.encode()).hexdigest()!=f["WKT_sha256"]:
                    raise ValueError("Source original WKT fragment SHA mismatch")
                params={"geometry":wkt,"year":f"{y0},{y1}","hasCoordinate":"true",
                        "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
                        "taxonKey":taxon_keys,"facet":"datasetKey","facetLimit":1000,"limit":0}
                total,counts,pages=dataset_facets(params)
                oldn=int(f["indexed_occurrences"])
                old_total+=oldn
                new_total+=total
                pcounts.update(counts)
                fragments.append({"fragment_number":i,"original_WKT_sha256":f["WKT_sha256"],
                    "original_2026_11_09_indexed_count":oldn,
                    "new_current_GBIF_indexed_count":total,
                    "new_current_GBIF_facets_pages":pages})
                print("GBIF_PROVIDER_FACET_PROBE",island,period,i,
                      "orig",oldn,"current",total,"datasetKeys",len(counts),flush=True)
            if sum(pcounts.values())!=new_total:
                raise ValueError("Provider source totals do not reconcile")
            rows.append({
                "original_landmass_id":island,"original_period":period,
                "original_source_GBIF_2026_11_09_index_count":old_total,
                "new_snapshot_index_count":new_total,
                "new_minus_original_index_count":new_total-old_total,
                "new_provider_count":len(pcounts),
                "original_source_2026_datasetKey_attribution":"NOT_CAPTURED",
                "new_snapshot_datasetkey_contributors":[{"datasetKey":k,"indexed_count":v}
                    for k,v in sorted(pcounts.items(),key=lambda kv:(-kv[1],kv[0]))],
                "source_fragments":fragments
            })
    return {
        "schema":"chocho_GBIF_20261120_new_snapshot_datasetKey_facet_after_20261109_source_replay_v1",
        "source_original_Git_commit":ORIGINAL_COMMIT,
        "original_2026_source":"GBIF occurrence/search source speciesKey-facet 2026-11-09; index only",
        "new_source":"Read-only current GBIF occurrence/search datasetKey facet on exact previous WKT geometries and source filters",
        "new_snapshot_start_utc":start,
        "new_snapshot_end_utc":datetime.now(timezone.utc).isoformat(),
        "original_GSHHG_237_sha256":digest,
        "source_taxon_filter":filter_used,
        "island_ids":list(islands),
        "n_new_index_queries":sum(len(z["source_fragments"]) for z in rows),
        "results":rows,
        "read_only_no_DOI_minted":True,
        "original_historical_dataset_publishers_recovered":False,
        "not_source_equivalent":"Current datasetKey facet is a newly retrieved GBIF index state; cannot be retroactively equated to original 2026-11-09 source dataset provider composition",
        "scientific_limit":"GBIF indexing is not true arrival, butterfly occupancy, larval performance, ship transport or plant-host colonization. Original simplified adjacent WKT borders could double-count boundary records."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--gshhg-zip",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--islands",default=",".join(DEFAULT_ISLANDS))
    a=p.parse_args()
    islands=tuple(sorted(set(a.islands.split(",")),key=int))
    j=probe(a.source_root,a.gshhg_zip,a.work,islands)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(j,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("NEW_GBIF_PROVIDER_FACET_DIAGNOSTIC_COMPLETED",
          json.dumps({"islands":len(islands),"island_periods":len(j["results"]),
                      "source_original_provider_unknown":True},ensure_ascii=False))
if __name__=="__main__":
    main()
