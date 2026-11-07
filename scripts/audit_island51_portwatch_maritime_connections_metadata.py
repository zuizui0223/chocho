#!/usr/bin/env python3
"""Bounded public ArcGIS source *metadata* for PortWatch maritime connection tables.

Only service listing, port registry metadata and up to 4 potential connection
layer schemas, counts, at most 5 rows each. NOT a bulk IMF data download,
NO WPI194 joins or butterfly outcomes. Check IMF/third-party usage rights
before any systematic extraction.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone
from urllib.parse import urlencode,urlsplit
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import time

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
ORIGINAL_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
SERVICE_BASE="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services"
PORTS_LAYER=SERVICE_BASE+"/PortWatch_ports_database/FeatureServer/0"
SEARCH_BASE="https://www.arcgis.com/sharing/rest/search"
SERVICE_MAX_BYTES=15_000_000
META_MAX_BYTES=1_800_000
SAMPLE_MAX=5

def git_blob_sha(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def fixed_wpi():
    if git_blob_sha(ORIGINAL)!=ORIGINAL_GIT_BLOB:
        raise ValueError("Source WPI original194 Git blob modified")
    with ORIGINAL.open("rt",encoding="utf-8",newline="") as stream:
        rows=list(csv.DictReader(stream))
    if len(rows)!=194 or len({row["port_id"] for row in rows})!=194:
        raise ValueError("Original WPI194 ports no longer exact source")
    if len({row["island_id"] for row in rows})!=51:
        raise ValueError("Original WPI source landmasses not 51")
    if sum(row["included"].lower()=="true" for row in rows)!=182:
        raise ValueError("Original 182 marine-qualified ports missing")
    return rows

def metadata(url,params=None,maxbytes=META_MAX_BYTES):
    raw_url=url+("?"+urlencode(params) if params else "")
    if not (raw_url.startswith("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/") or
            raw_url.startswith("https://www.arcgis.com/sharing/rest/search")):
        raise ValueError("Only vetted IMF ArcGIS org and official Esri item catalog metadata permitted")
    receipt={"request_url":raw_url,"status":None,
             "bytes":0,"error":None,"returned_fields":[]}
    try:
        req=Request(raw_url,headers={"User-Agent":"chocho-original194-portwatch-marine-connection-metadata-audit/1.0",
                                    "Accept":"application/json"})
        with urlopen(req,timeout=30) as response:
            length=response.headers.get("Content-Length")
            if length and int(length)>maxbytes:
                raise ValueError("ArcGIS metadata response exceeds precommitted cap")
            data=response.read(maxbytes+1)
            if len(data)>maxbytes:
                raise ValueError("ArcGIS metadata response exceeds precommitted cap")
            receipt["status"]=getattr(response,"status",200)
            receipt["bytes"]=len(data)
        obj=json.loads(data.decode("utf-8"))
        if isinstance(obj,dict):
            receipt["returned_fields"]=sorted(obj)[:40]
            if obj.get("error"):
                receipt["error"]={"type":"ARCGIS_BODY_ERROR","data":str(obj["error"])[:500]}
                return None,receipt
        return obj,receipt
    except HTTPError as exc:
        receipt["status"]=exc.code
        receipt["error"]={"type":"HTTPError","detail":str(exc.reason)[:150]}
    except (URLError,TimeoutError,ConnectionError,ValueError) as exc:
        receipt["error"]={"type":type(exc).__name__,"detail":str(exc)[:200]}
    return None,receipt

def service_candidates(directory,search):
    entries=[]
    if isinstance(directory,dict):
        for s in directory.get("services") or []:
            name=str(s.get("name") or "")
            typ=str(s.get("type") or "")
            full=SERVICE_BASE+"/"+name+"/"+typ
            if typ not in ("FeatureServer","MapServer"):
                continue
            if any(k in name.lower() for k in ("maritime","connections","connection","shipping","network")):
                entries.append({"service_name":name,"type":typ,"source":"provider_directory","url":full})
    if isinstance(search,dict):
        for s in search.get("results") or []:
            title=str(s.get("title") or "")
            url=str(s.get("url") or "")
            if (any(k in title.lower() for k in ("maritime","connections","shipping","network")) and
                url.lower().startswith(SERVICE_BASE.lower())):
                entries.append({"service_name":title,"type":str(s.get("type") or ""),
                                "source":"official_esri_search","url":url})
    # Prefer specifically maritime connections.
    distinct={x["url"]:x for x in entries}
    return sorted(distinct.values(),key=lambda x:(
        0 if "maritime" in x["service_name"].lower() and
             "connection" in x["service_name"].lower() else
        1 if "connection" in x["service_name"].lower() else 2,
        x["service_name"].lower()))[:4]

def compact_fields(layer):
    return [{"name":str(f.get("name") or ""),"type":str(f.get("type") or ""),
             "alias":str(f.get("alias") or "")[:100]}
            for f in (layer.get("fields") or [])]

def safe_sample(result):
    rows=[]
    for feature in (result.get("features") or [])[:SAMPLE_MAX]:
        attributes=feature.get("attributes") or {}
        # Source metadata diagnostics ONLY, no personal vessel identifier.
        public={}
        for k,v in attributes.items():
            lower=k.lower()
            if any(secret in lower for secret in ("mmsi","imo","vessel_name","ship_name","callsign","ssvid")):
                continue
            if v is None or isinstance(v,(str,int,float,bool)):
                public[k]=str(v)[:100] if isinstance(v,str) else v
        rows.append(public)
    return rows

def source_probe(root):
    wpi=fixed_wpi()
    out=Path(root).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("results","data","provenance")):
        raise ValueError("Never overwrite original GBIF and maritime scientific data")
    out.mkdir(parents=True,exist_ok=True)
    fetched=[]
    directory,direct_req=metadata(SERVICE_BASE,{"f":"json"},SERVICE_MAX_BYTES)
    fetched.append(direct_req)
    search,search_req=metadata(SEARCH_BASE,{
        "f":"json","num":30,
        "q":"Maritime Connections PortWatch"
    },SERVICE_MAX_BYTES)
    fetched.append(search_req)
    registry,registry_req=metadata(PORTS_LAYER,{"f":"json"})
    fetched.append(registry_req)
    services=service_candidates(directory,search)
    network_entries=[]
    for cand in services:
        root_url=cand["url"]
        # First inspect service root: needed to choose tables/layers by provider ID.
        service,service_req=metadata(root_url,{"f":"json"})
        fetched.append(service_req)
        entry={"source_name":cand["service_name"],"service_url":root_url,
               "source_origin":cand["source"],"source_service_status":service_req["status"],
               "source_service_error":service_req["error"],"layers":[]}
        if not isinstance(service,dict):
            network_entries.append(entry)
            continue
        choices=(service.get("layers") or [])+(service.get("tables") or [])
        relevant=[v for v in choices if any(k in str(v.get("name") or "").lower()
                   for k in ("maritime","connection","network","shipping","port"))]
        for row in (relevant or choices)[:2]:
            layer_id=row.get("id")
            if not isinstance(layer_id,int):continue
            layer_url=root_url+"/"+str(layer_id)
            layer,layer_req=metadata(layer_url,{"f":"json"})
            fetched.append(layer_req)
            if not isinstance(layer,dict):
                entry["layers"].append({"id":layer_id,"error":layer_req["error"]})
                continue
            count,count_req=metadata(layer_url+"/query",{"where":"1=1","returnCountOnly":"true","f":"json"})
            fetched.append(count_req)
            example,sample_req=metadata(layer_url+"/query",{
                "where":"1=1","outFields":"*","returnGeometry":"false",
                "resultRecordCount":str(SAMPLE_MAX),"f":"json"
            })
            fetched.append(sample_req)
            example_rows=safe_sample(example) if isinstance(example,dict) else []
            entry["layers"].append({
                "id":layer_id,"name":str(layer.get("name") or ""),
                "type":str(layer.get("type") or ""),
                "objectIdField":str(layer.get("objectIdField") or layer.get("objectIdFieldName") or ""),
                "copyrightText":str(layer.get("copyrightText") or ""),
                "description":str(layer.get("description") or "")[:1000],
                "fields":compact_fields(layer),
                "record_count_from_provider":count.get("count") if isinstance(count,dict) else None,
                "provider_count_error":count_req["error"],
                "bounded_source_sample_rows":example_rows,
                "provider_example_error":sample_req["error"],
            })
        network_entries.append(entry)
    reg={}
    if isinstance(registry,dict):
        reg={"name":registry.get("name"),
             "type":registry.get("type"),
             "record_count_cap_per_page":registry.get("maxRecordCount"),
             "copyrightText":registry.get("copyrightText"),
             "fields":compact_fields(registry),
             "sourceLastEditDate":(registry.get("editingInfo") or {}).get("lastEditDate")}
    result={
        "schema":"chocho_original51_IMF_PortWatch_maritime_connections_metadata_only_v1",
        "status":"BOUNDED_PUBLIC_ARCGIS_METADATA_DISCOVERY_COMPLETE",
        "reviewed_utc":datetime.now(timezone.utc).isoformat(),
        "source":"International Monetary Fund IMF PortWatch / UN Global Platform",
        "source_portwatch_url":"https://portwatch.imf.org/",
        "source_official_registry":PORTS_LAYER,
        "source_arcgis_service_directory":SERVICE_BASE,
        "original_WPI2019_source_git_blob":ORIGINAL_GIT_BLOB,
        "original_WPI_ports":len(wpi),"original_WPI_islands":51,"marine_qualified_ports":182,
        "static_connections_155521_rows_other_paper_report_not_self_verified":True,
        "discovery_service_names":[x["service_name"] for x in services],
        "network_candidate_services":network_entries,
        "reference_ports_registry_schema":reg,
        "provider_request_receipts":fetched,
        "provider_service_total":len(directory.get("services") or []) if isinstance(directory,dict) else None,
        "original_service_search_total":search.get("total") if isinstance(search,dict) else None,
        "IMF_data_usage_policy":"https://www.imf.org/en/about/copyright-and-terms",
        "no_systematic_bulk_IMF_content_download_attempted":True,
        "no_large_port_network_source_rows_repushed_to_github":True,
        "at_most_5_example_rows_per_discovered_source_layer":True,
        "not_proof_of_cargo_or_butterfly_dispersal":True,
        "no_WPI_port_identity_mapping_or_biological_model_performed":True,
    }
    target=out/"portwatch_maritime_connections_metadata_receipt.json"
    target.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PORTWATCH_MARITIME_CONNECTIONS_METADATA_JSON="+json.dumps({
        "status":result["status"],
        "service_directory_status":direct_req["status"],
        "service_total":result["provider_service_total"],
        "search_status":search_req["status"],
        "search_results":result["original_service_search_total"],
        "port_registry_status":registry_req["status"],
        "port_registry_field_names":[f["name"] for f in reg.get("fields",[])],
        "network_candidate_services":[{"name":c["source_name"],
            "url":c["service_url"],"status":c["source_service_status"],
            "layers":[{"name":l.get("name"),"count":l.get("record_count_from_provider"),
                       "fields":[f["name"] for f in l.get("fields",[])],
                       "sample_rows":len(l.get("bounded_source_sample_rows",[]))}
                      for l in c["layers"]]} for c in network_entries],
        "no_bulk_port_connection_download":True},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,required=True)
    source_probe(p.parse_args().out_dir)
