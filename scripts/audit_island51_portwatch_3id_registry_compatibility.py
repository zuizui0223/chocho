#!/usr/bin/env python3
"""At most THREE IMF PortWatch sample maritime source IDs resolved to port registry IDs.

Only bounded source-metadata samples and 1 registry count + 1 <=5-row ID lookup.
No bulk IMF rows, no 194 WPI join, no butterfly modeling, no ship-voyage inference.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,urlopen
from urllib.parse import urlencode

ROOT=Path(__file__).resolve().parents[1]
METADATA=ROOT/"results/island51_portwatch_global_connections_20261107/portwatch_maritime_connections_metadata_receipt.json"
PORTS=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
PORTS_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
REGISTRY="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0"
CONNECTION="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/Spillover_Simulator_Maritime_Connections/FeatureServer/0"
EARTH=6371.0088
SAMPLE_ID_COUNT=3
MATCH_DISTANCE_KM=25.

def blob(path):
    b=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def fixed_scope():
    if blob(PORTS)!=PORTS_BLOB:raise ValueError("Original frozen 194 WPI port coordinates changed")
    with PORTS.open("r",newline="",encoding="utf-8") as f:
        x=list(csv.DictReader(f))
    if len(x)!=194 or len({r["island_id"] for r in x})!=51 or sum(r["included"].lower()=="true" for r in x)!=182:
        raise ValueError("Original 51 island/194 biological port/182 marine scope differs")
    return x

def haversine(a,b,c,d):
    p,q,r,s=[math.radians(float(z)) for z in (a,b,c,d)]
    v=math.sin((r-p)/2)**2+math.cos(p)*math.cos(r)*math.sin((s-q)/2)**2
    return 2*EARTH*math.asin(math.sqrt(max(0.,min(1.,v))))

def source_sample():
    snap=json.loads(METADATA.read_text(encoding="utf-8"))
    if (snap.get("status")!="BOUNDED_PUBLIC_ARCGIS_METADATA_DISCOVERY_COMPLETE" or
        snap.get("original_WPI_ports")!=194 or
        snap.get("original_WPI_islands")!=51):
        raise ValueError("Parent PortWatch source metadata receipt incomplete")
    target=next((s for s in snap["network_candidate_services"]
       if s["service_url"]+"/0"==CONNECTION),None)
    if not target:
        raise ValueError("Documented actual official PortWatch Maritime Connections source absent")
    layer=next((x for x in target["layers"] if int(x["id"])==0),None)
    if not layer or layer.get("record_count_from_provider")!=155521:
        raise ValueError("Public original 155521 source rows metadata changed")
    rows=layer.get("bounded_source_sample_rows") or []
    if len(rows)!=5:
        raise ValueError("Original bounded public source connection sample should have 5 rows")
    selected={}
    for r in rows:
        for prefix in ("from","to"):
            idkey=prefix+"_id"
            namekey=prefix+"_name"
            latkey=prefix+"_lat"
            lonkey=prefix+"_lon"
            key=str(r[idkey])
            if not key.startswith("port"):
                raise ValueError("Provider source connection ID does not use PortWatch port namespace")
            candidate={"portid":key,"source_connection_name":str(r[namekey]),
                "source_lat":float(r[latkey]),"source_lon":float(r[lonkey]),
                "source_iso3":str(r[prefix+"_iso3"])}
            if key in selected:
                previous=selected[key]
                if any(previous[f]!=candidate[f] for f in ("source_lat","source_lon","source_iso3")):
                    raise ValueError("The SAME original PortWatch portid has inconsistent example source coordinates")
            else:
                selected[key]=candidate
    first=list(selected.values())[:SAMPLE_ID_COUNT]
    if len(first)!=SAMPLE_ID_COUNT:
        raise ValueError("Fewer than three distinct source PortWatch identifiers in original five example links")
    return first

def fetch_public_metadata(params):
    url=REGISTRY+"/query?"+urlencode(params)
    if not url.startswith("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/"):
        raise ValueError("Never contact unverified source host for registry samples")
    request=Request(url,headers={"User-Agent":"chocho-original194-PortWatch-3-ID-source-audit/1.0","Accept":"application/json"})
    with urlopen(request,timeout=30) as resp:
        raw=resp.read(1000000+1)
        if len(raw)>1000000:
            raise ValueError("Registry source reply too large for THREE-record bounded sample")
    result=json.loads(raw.decode("utf-8"))
    if not isinstance(result,dict) or result.get("error"):
        raise ValueError("Provider port registry returned source error "+str(result.get("error"))[:200])
    return result

def main(out_dir):
    fixed_scope()
    original=source_sample()
    original_ids=[v["portid"] for v in original]
    count=fetch_public_metadata({"where":"1=1","returnCountOnly":"true","f":"json"})
    port_count=count.get("count")
    if not isinstance(port_count,int) or port_count<=0:
        raise ValueError("Official PortWatch registry count-only query unavailable")
    quoted=",".join("'"+v.replace("'","''")+"'" for v in original_ids)
    values=fetch_public_metadata({
        "where":"portid IN ("+quoted+")",
        "outFields":"portid,portname,country,ISO3,lat,lon,LOCODE",
        "resultRecordCount":"5","returnGeometry":"false","f":"json"
    })
    registry=[]
    for row in values.get("features") or []:
        item=(row or {}).get("attributes") or {}
        if item.get("portid") is not None:
            registry.append(item)
    by_id={}
    duplicates=[]
    for row in registry:
        key=str(row["portid"])
        if key in by_id:duplicates.append(key)
        by_id[key]=row
    checks=[]
    for source in original:
        key=source["portid"]
        registry_port=by_id.get(key)
        if registry_port:
            dist=haversine(source["source_lat"],source["source_lon"],
                registry_port.get("lat"),registry_port.get("lon"))
        else:
            dist=None
        checks.append({
            "connection_portid":key,
            "connection_source_port_name":source["source_connection_name"],
            "connection_source_iso3":source["source_iso3"],
            "connection_source_lon":source["source_lon"],
            "connection_source_lat":source["source_lat"],
            "registry_portid_match":bool(registry_port),
            "registry_port_name":registry_port.get("portname") if registry_port else None,
            "registry_country":registry_port.get("country") if registry_port else None,
            "registry_ISO3":registry_port.get("ISO3") if registry_port else None,
            "registry_LOCODE":registry_port.get("LOCODE") if registry_port else None,
            "registry_lat":registry_port.get("lat") if registry_port else None,
            "registry_lon":registry_port.get("lon") if registry_port else None,
            "original_provider_port_id_coordinate_distance_km":dist,
            "source_id_and_geography_illustratively_compatible":bool(
                registry_port and str(registry_port.get("ISO3"))==source["source_iso3"] and dist is not None and dist<=MATCH_DISTANCE_KM),
            "no_original_WPI194_port_identity_derived":True,
        })
    complete=all(x["source_id_and_geography_illustratively_compatible"] for x in checks) and not duplicates
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("results","data","provenance")):
        raise ValueError("Do not overwrite original provider or biological scientific source")
    out.mkdir(parents=True,exist_ok=True)
    report={
       "schema":"chocho_original51_IMF_PortWatch_3_sample_maritime_connection_to_registry_portid_v1",
       "status":"THREE_SAMPLE_PROVIDER_ID_JOIN_COMPLETED_NO_BULK_NETWORK_EXTRACTION",
       "source_connection_URL":CONNECTION,
       "source_registry_URL":REGISTRY,
       "source_provider":"International Monetary Fund PortWatch",
       "source_attribution":"Source: International Monetary Fund, IMF PortWatch",
       "source_utc":datetime.now(timezone.utc).isoformat(),
       "source_registry_record_count_from_count_only_query":port_count,
       "source_original_connection_edge_record_count_from_provider":155521,
       "source_3_portids_selected_before_query":original_ids,
       "number_registry_rows_returned":len(registry),
       "sample_identifier_duplicates":duplicates,
       "samples":checks,
       "all_three_sample_IDs_and_coordinates_compatible":complete,
       "source_original_WPI194_git_blob":PORTS_BLOB,
       "original_WPI_ports":194,"original_source_islands":51,
       "original_marine_eligible_ports":182,
       "provider_model_edges_not_vessel_level_2020_OD_trips":True,
       "source_IMF_rights_review_before_any_systematic_bulk_harvest":True,
       "no_original_WPI194_join_or_butterfly_model_performed":True,
       "no_bulk_provider_data_requested":True
    }
    (out/"source_bounded_connection_registry_3id_receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PORTWATCH_3_IDENTITY_SAMPLE_JSON="+json.dumps({
        "status":report["status"],"registry_count":port_count,"source_connection_rows":155521,
        "sample_provider_IDs":original_ids,"retrieved_registry_rows":len(registry),
        "all_3_source_ids_join_semantically":complete,
        "sample_checks":[{"portid":x["connection_portid"],
           "distance_km":x["original_provider_port_id_coordinate_distance_km"],
           "iso3_match":x["registry_ISO3"]==x["connection_source_iso3"],
           "compatible":x["source_id_and_geography_illustratively_compatible"]}
           for x in checks],
        "no_original_WPI_join":True},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    main(p.parse_args().out_dir)
