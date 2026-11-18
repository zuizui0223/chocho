#!/usr/bin/env python3
"""Independent standard-library source-only check of 3 PortWatch port-id reference joins."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
PARENT=OUT/"portwatch_maritime_connections_metadata_receipt.json"
TARGET=OUT/"source_bounded_connection_registry_3id_receipt.json"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"

def distance(a,b,c,d):
    a,b,c,d=[math.radians(float(v)) for v in (a,b,c,d)]
    v=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 2*6371.0088*math.asin(math.sqrt(max(0,min(1,v))))

def verify():
    p=json.loads(PARENT.read_text(encoding="utf-8"))
    r=json.loads(TARGET.read_text(encoding="utf-8"))
    raw=WPI.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()
    with WPI.open(encoding="utf-8",newline="") as f:
        original=list(csv.DictReader(f))
    if (blob!=WPI_BLOB or len(original)!=194 or len({v["island_id"] for v in original})!=51 or
        sum(v["included"].lower()=="true" for v in original)!=182):
        raise ValueError("Source original 194 WPI ports no longer the verified inputs")
    if (p["status"]!="BOUNDED_PUBLIC_ARCGIS_METADATA_DISCOVERY_COMPLETE" or
        r["status"]!="THREE_SAMPLE_PROVIDER_ID_JOIN_COMPLETED_NO_BULK_NETWORK_EXTRACTION" or
        p["original_WPI_ports"]!=194 or r["original_WPI_ports"]!=194 or
        r["source_original_connection_edge_record_count_from_provider"]!=155521):
        raise ValueError("PortWatch source-only discovery and 3-row ID namespace audit are inconsistent")
    service=next((x for x in p["network_candidate_services"]
          if x["service_url"]+"/0"==r["source_connection_URL"]),None)
    if service is None:
        raise ValueError("Connection source service differs from original authoritative metadata")
    layer=next((x for x in service["layers"] if x["id"]==0),None)
    if layer is None or layer["record_count_from_provider"]!=155521:
        raise ValueError("Original source connection count differs")
    sampled={}
    for connection in layer["bounded_source_sample_rows"][:5]:
        for kind in ("from","to"):
            key=str(connection[kind+"_id"])
            if key not in sampled:
                sampled[key]=connection
    sample_ids=list(sampled)[:3]
    if (len(r["samples"])!=3 or r["source_3_portids_selected_before_query"]!=sample_ids or
        r["number_registry_rows_returned"]>5):
        raise ValueError("Source sample port ID selection/post-query scope changed")
    if len({x["connection_portid"] for x in r["samples"]})!=3:
        raise ValueError("Repeated sampled PortWatch port ID")
    success=True
    for check in r["samples"]:
        pid=check["connection_portid"]
        c=sampled[pid]
        where="from" if c["from_id"]==pid else "to"
        if (check["connection_source_iso3"]!=c[where+"_iso3"] or
            not math.isclose(check["connection_source_lat"],c[where+"_lat"],abs_tol=1e-10) or
            not math.isclose(check["connection_source_lon"],c[where+"_lon"],abs_tol=1e-10)):
            raise ValueError("PortWatch connection example source identity differs")
        if check["registry_portid_match"]:
            d=distance(check["connection_source_lat"],check["connection_source_lon"],
                       check["registry_lat"],check["registry_lon"])
            if not math.isclose(d,check["original_provider_port_id_coordinate_distance_km"],
                                abs_tol=1e-7,rel_tol=1e-10):
                raise ValueError("Registry geographic join was not accurately calculated")
            compatible=d<=25 and check["connection_source_iso3"]==check["registry_ISO3"]
            if check["source_id_and_geography_illustratively_compatible"]!=compatible:
                raise ValueError("Source port ID+coordinate compatibility status differs")
        success &=bool(check["source_id_and_geography_illustratively_compatible"])
        if check["no_original_WPI194_port_identity_derived"] is not True:
            raise ValueError("Illustrative PortWatch portid was incorrectly equated to WPI")
    if r["all_three_sample_IDs_and_coordinates_compatible"]!=success:
        raise ValueError("Three PortWatch source ID summary does not match original records")
    report={
      "status":"ORIGINAL_THREE_PORTWATCH_PROVIDER_PORTIDS_SOURCE_REGISTRY_VERIFIED",
      "source_sample_IDs":sample_ids,
      "source_registry_row_count":r["source_registry_record_count_from_count_only_query"],
      "source_connections_table_row_count":155521,
      "source_sample_compatible_count":sum(x["source_id_and_geography_illustratively_compatible"] for x in r["samples"]),
      "no_194_WPI_port_name_match_claimed":True,
      "no_155k_arcgis_bulk_provider_download":True,
      "no_actual_ship_journey_count_or_ecological_prediction":True,
    }
    print("PORTWATCH_3ID_INDEPENDENT_SOURCE_AUDIT_JSON="+json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    verify()
