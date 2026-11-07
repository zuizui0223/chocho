#!/usr/bin/env python3
"""Acquire ALL original 194 WPI-port pre-2011 Plantae speciesKey facets from live public GBIF.

Not the missing original author's host/port parquet. Source dates are OBSERVATION years;
a GBIF index zero is not botanical absence. The 5km search differs from historical cache.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

ROOT=Path(__file__).resolve().parents[1]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
ARCHIVED=ROOT/"provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_plant_availability.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
ARCHIVED_BLOB="0a01831056c2cd59a4bd41bd1dbd4c84d77cf922"
API="https://api.gbif.org/v1/occurrence/search"
USER_AGENT="chocho-20261107-original194-early-Plantae-speciesKey-audit/1.0"
REQUEST_TIMEOUT=80
MAX_RETRIES=7
FREQUENCY_SECONDS=0.65
QUERY_CONST={
    "taxonKey":"6","year":"2001,2010",
    "hasCoordinate":"true","hasGeospatialIssue":"false",
    "occurrenceStatus":"PRESENT","facet":"speciesKey",
    "facetLimit":"1000","limit":"0",
}

def hash_blob(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def original_source_rows():
    if hash_blob(WPI)!=WPI_BLOB or hash_blob(ARCHIVED)!=ARCHIVED_BLOB:
        raise ValueError("The exact archived original WPI or early port-count source Git blob changed")
    g=read_csv(WPI)
    e=read_csv(ARCHIVED)
    if len(g)!=194 or len(e)!=194:
        raise ValueError("Original WPI194 or early GBIF194 inventory incomplete")
    geo={}
    for r in g:
        pid=int(r["port_id"]);iid=str(r["landmass_id"])
        if pid in geo or str(r["island_id"])!=iid:
            raise ValueError("Duplicated port source ID or different original island")
        a=float(r["original_latitude"]);b=float(r["original_longitude"])
        if not(-90<=a<=90 and -180<=b<=180) or not(math.isfinite(a) and math.isfinite(b)):
            raise ValueError("Invalid exact WPI2019 source port coordinates")
        geo[pid]={"island_id":iid,"port_id":pid,"latitude":a,"longitude":b,"country":r["country_code"],
                  "port_name":r["port_name"],"marine_eligible":r["included"].lower()=="true"}
    old={}
    for r in e:
        pid=int(r["port_id"])
        if pid in old or pid not in geo or geo[pid]["island_id"]!=str(r["landmass_id"]):
            raise ValueError("Early GBIF original194 port ID/island join disagrees with WPI")
        old[pid]=int(r["early_plant_records"])
        if old[pid]<0:raise ValueError("Negative original GBIF source record count")
    if len(geo)!=194 or len(old)!=194 or len({r["island_id"] for r in geo.values()})!=51:
        raise ValueError("Original 51 islands / 194 ports are not represented")
    if sum(x["marine_eligible"] for x in geo.values())!=182:
        raise ValueError("Original 182 marine-qualified port flags changed")
    if sum(old.values())!=155002:
        raise ValueError("Original saved 2001-2010 port-buffer record-index total changed")
    return geo,old

def query_params(port):
    params=dict(QUERY_CONST)
    params["geoDistance"]=f"{port['latitude']:.7f},{port['longitude']:.7f},5km"
    return params

def get_page(params):
    query=urlencode(params)
    url=API+"?"+query
    for i in range(MAX_RETRIES):
        try:
            req=Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
            with urlopen(req,timeout=REQUEST_TIMEOUT) as response:
                data=response.read()
                if response.status!=200:
                    raise ValueError("Unexpected status "+str(response.status))
            result=json.loads(data.decode("utf-8"))
            if not isinstance(result,dict) or "count" not in result or result.get("error"):
                raise ValueError("GBIF response missing count or contains error "+str(result.get("error"))[:140])
            return result,hashlib.sha256(data).hexdigest()
        except HTTPError as exc:
            if exc.code==400:
                snippet=exc.read(700).decode("utf-8",errors="replace")
                raise ValueError("Invalid GBIF query: "+snippet+" "+url[:320]) from exc
            if exc.code not in (408,425,429,500,502,503,504) or i+1==MAX_RETRIES:
                raise
            retry=exc.headers.get("Retry-After")
            try:delay=int(retry)
            except (TypeError,ValueError):delay=2**i
            time.sleep(min(120,max(2,delay)))
        except (URLError,TimeoutError,ConnectionError):
            if i+1==MAX_RETRIES:raise
            time.sleep(min(45,2**i+1))
    raise AssertionError("Unreachable GBIF request retries")

def fetch_port(port):
    const=query_params(port)
    count=None
    species={}
    pages=[]
    for offset in range(0,50001,1000):
        params=dict(const,facetOffset=str(offset))
        result,sha=get_page(params)
        n=int(result["count"])
        if n<0 or (count is not None and n!=count):
            raise ValueError("GBIF original port total changed within facet pagination")
        count=n
        facets=result.get("facets") or []
        arrays=[x.get("counts",[]) for x in facets
                if str(x.get("field","")).upper() in ("SPECIES_KEY","SPECIESKEY")]
        if len(arrays)>1:
            raise ValueError("Ambiguous duplicated GBIF SPECIES_KEY facet")
        values=arrays[0] if arrays else []
        if n>0 and not arrays:
            raise ValueError("GBIF returned plant records but no SPECIES_KEY facet")
        if len(values)>1000:
            raise ValueError("GBIF returned more than requested 1000 species facets")
        for r in values:
            key=int(r["name"]);records=int(r["count"])
            if key<=0 or records<=0 or key in species:
                raise ValueError("GBIF species key duplicate or invalid between pages")
            species[key]=records
        pages.append({"facet_offset":offset,"rows":len(values),"sha256":sha})
        if len(values)<1000:
            break
        if offset>=50000:
            raise ValueError("Exceeded 50,000 species facets at one WPI port")
        time.sleep(FREQUENCY_SECONDS)
    if sum(species.values())>count:
        raise ValueError("Species facets exceed source occurrence count")
    return {
       "port_id":port["port_id"],"island_id":port["island_id"],
       "source_country":port["country"],"query_parameters":const,
       "source_year_range":[2001,2010],
       "indexed_occurrence_count":count,"species_key_count":len(species),
       "records_without_species_key":count-sum(species.values()),
       "species_key_counts":[{"speciesKey":k,"indexed_records":v} for k,v in sorted(species.items())],
       "facet_pages":pages,"retrieved_at_utc":datetime.now(timezone.utc).isoformat(),
       "query_provides_indexed_presence_not_true_presence":True,
    }

def run(out_dir):
    destination=Path(out_dir).resolve()
    if any(destination.is_relative_to((ROOT/p).resolve()) for p in ("data","provenance","results")):
        raise ValueError("Never overwrite original result/source directories during acquisition")
    original,archive=original_source_rows()
    cache=destination/"port_json"
    cache.mkdir(parents=True,exist_ok=True)
    manifest=[]
    taxon_rows=[]
    for idx,pid in enumerate(sorted(original),1):
        port=original[pid]
        file=cache/(str(pid)+".json")
        if file.exists():
            item=json.loads(file.read_text(encoding="utf-8"))
            if (item.get("port_id")!=pid or item.get("island_id")!=port["island_id"] or
                item.get("query_parameters")!=query_params(port) or
                len(item.get("facet_pages",[]))<1):
                raise ValueError("Invalid incompatible acquisition checkpoint "+str(pid))
        else:
            time.sleep(FREQUENCY_SECONDS)
            item=fetch_port(port)
            file.write_text(json.dumps(item,sort_keys=True,separators=(",",":"))+"\n",encoding="utf-8")
        count=int(item["indexed_occurrence_count"])
        species=item["species_key_counts"]
        if len(species)!=item["species_key_count"] or sum(int(r["indexed_records"]) for r in species)>count:
            raise ValueError("Malformed saved GBIF original194 species facets")
        for r in species:
            taxon_rows.append((str(pid),port["island_id"],str(r["speciesKey"]),str(r["indexed_records"])))
        manifest.append({
            "port_id":pid,"island_id":port["island_id"],
            "original_early_port_buffer_records":archive[pid],
            "live_GBIF_2001_2010_plant_records":count,
            "live_faceted_species_keys":len(species),
            "live_non_species_key_records":int(item["records_without_species_key"]),
            "source_count_difference":count-archive[pid],
            "facet_pages":len(item["facet_pages"]),
            "query_5km_from_original_WPI_centroid":item["query_parameters"]["geoDistance"],
            "port_json_sha256":digest(file),
            "original_marine_eligible":bool(port["marine_eligible"]),
        })
        if idx%10==0 or idx in (1,194):
            print("GBIF194_PROGRESS="+json.dumps({"ports_complete":idx,"total_ports":194,
                   "last_port":pid,"speciesKey_rows":len(taxon_rows),
                   "counts_changed":sum(x["source_count_difference"]!=0 for x in manifest)},sort_keys=True),flush=True)
    if len(manifest)!=194 or len({x[0] for x in taxon_rows})<1:
        raise ValueError("Incomplete original194 botanical source")
    destination.mkdir(parents=True,exist_ok=True)
    with (destination/"GBIF_original194_port_2001_2010_coverage.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(manifest[0]))
        writer.writeheader();writer.writerows(manifest)
    with gzip.open(destination/"GBIF_original194_port_2001_2010_speciesKey.csv.gz","wt",newline="",encoding="utf-8",compresslevel=6) as f:
        writer=csv.writer(f)
        writer.writerow(["port_id","island_id","GBIF_speciesKey","indexed_records_2001_2010"])
        writer.writerows(sorted(taxon_rows,key=lambda r:(int(r[0]),int(r[2]))))
    islands=defaultdict(set)
    for pid,island,key,records in taxon_rows:
        islands[island].add(int(key))
    original_islands={x["island_id"] for x in original.values()}
    island_rows=[{
        "island_id":k,"original_port_count":sum(p["island_id"]==k for p in original.values()),
        "ports_with_plant_record":sum(x["island_id"]==k and x["live_GBIF_2001_2010_plant_records"]>0 for x in manifest),
        "unique_early_indexed_general_plant_speciesKeys":len(islands[k]),
        "sum_live_original_port_5km_records_NOT_unique_occurrences":sum(
            x["live_GBIF_2001_2010_plant_records"] for x in manifest if x["island_id"]==k),
        "sum_original_cached_early_records_NOT_unique_occurrences":sum(
            x["original_early_port_buffer_records"] for x in manifest if x["island_id"]==k),
    } for k in sorted(original_islands,key=int)]
    with (destination/"GBIF_original51_early_general_flora_islands.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(island_rows[0]))
        writer.writeheader();writer.writerows(island_rows)
    summary={
        "schema":"chocho_original194_direct_public_GBIF_Plant_speciesKey_2001_2010_v1",
        "status":"COMPLETE_ALL_194_SOURCE_PORTS_INDEPENDENT_GBIF_LIVE_ACQUISITION",
        "source_WPI_Git_blob":WPI_BLOB,"original_archive_Git_blob":ARCHIVED_BLOB,
        "GBIF_API":API,"query_constants":QUERY_CONST,"source_query_radius_km":5,
        "source_port_original_geo_count":194,"source_eligible_islands":51,
        "source_marine_eligible_ports":182,"source_period":[2001,2010],
        "source_retrieval_times_utc":[item["retrieved_at_utc"] for item in
              [json.loads((cache/(str(p)+".json")).read_text()) for p in (min(original),max(original))]],
        "live_total_indexed_port_buffer_records_can_overlap":sum(x["live_GBIF_2001_2010_plant_records"] for x in manifest),
        "original_cached_total_port_buffer_records":sum(archive.values()),
        "port_count_mismatches":sum(x["source_count_difference"]!=0 for x in manifest),
        "positive_index_ports":sum(x["live_GBIF_2001_2010_plant_records"]>0 for x in manifest),
        "zero_index_ports":sum(x["live_GBIF_2001_2010_plant_records"]==0 for x in manifest),
        "islands_with_speciesKey":sum(x["unique_early_indexed_general_plant_speciesKeys"]>0 for x in island_rows),
        "islands_with_no_speciesKey_not_true_absence":sum(x["unique_early_indexed_general_plant_speciesKeys"]==0 for x in island_rows),
        "total_port_speciesKey_rows":len(taxon_rows),
        "distinct_GBIF_speciesKeys_within_source_query":len({int(x[2]) for x in taxon_rows}),
        "query_live_is_not_original_cached_geometric_or_temporal_equivalence":True,
        "speciesKey_is_general_plant_index_not_butterfly_host":True,
        "overlapping_5km_buffers_do_not_provide_independent_record_counts":True,
        "port_identity_and_geography_uses_original_archived_WPI_not_current_mirror":True,
        "source_sha256":{
          p.name:digest(p) for p in destination.iterdir() if p.is_file() and p.name!="source_receipt.json"},
    }
    (destination/"source_receipt.json").write_text(json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GBIF194_SOURCE_RECEIPT_JSON="+json.dumps({
        k:v for k,v in summary.items() if k in (
           "status","source_port_original_geo_count","port_count_mismatches",
           "positive_index_ports","zero_index_ports","total_port_speciesKey_rows",
           "distinct_GBIF_speciesKeys_within_source_query","islands_with_speciesKey",
           "islands_with_no_speciesKey_not_true_absence",
           "live_total_indexed_port_buffer_records_can_overlap")},sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)
