#!/usr/bin/env python3
"""Bounded four-WPI-source PortWatch directed modeled network edge diagnostic.

One count query for each of twelve ordered verified PortWatch port-ID pairs, and
one single-row example for each nonzero pair. Four incoming and four outgoing
degree count-only queries. Not real 2020 voyages, no 155k-row harvest and no
biological response models.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import time
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[1]
PORTS=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
SOURCE_SAMPLE=ROOT/"results/island51_portwatch_global_connections_20261107/WPI5_bounded_registry_source_pilot_receipt.json"
MODELED_LAYER="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/Spillover_Simulator_Maritime_Connections/FeatureServer/0"
ZENODO_VERIFIED=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/observed2020_verified_WPI_port_links.csv"
WPI_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
EXPECTED={
  "28":("port141","26","BELFAST","GB"),
  "428":("port1153","206","SANTA CRUZ DE TENERIFE","ES"),
  "817":("port647","85","LIMASSOL","CY"),
  "1416":("port254","31","COLOMBO","LK")
}
MAX_FEATURE_EXAMPLES=12
EARTH_RADIUS_KM=6371.0088
BUDGET_COUNTS=20
RETRYABLE={408,429,500,502,503,504}
OUT_FIELDS=("from_id","to_id","from_name","to_name","from_iso3","to_iso3",
            "from_lat","from_lon","to_lat","to_lon","distance","transit_days")

def git_blob(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def csv_rows(path):
    with Path(path).open("rt",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def geo_km(lat1,lon1,lat2,lon2):
    p,q,r,s=(math.radians(float(z)) for z in (lat1,lon1,lat2,lon2))
    v=math.sin((r-p)/2)**2+math.cos(p)*math.cos(r)*math.sin((s-q)/2)**2
    return 2*EARTH_RADIUS_KM*math.asin(math.sqrt(max(0,min(1,v))))

def fixed_sources():
    if git_blob(PORTS)!=WPI_GIT_BLOB:
        raise ValueError("Original 194-port WPI source SHA differs")
    rows=csv_rows(PORTS)
    if (len(rows)!=194 or len({r["port_id"] for r in rows})!=194 or
        len({r["island_id"] for r in rows})!=51 or
        sum(r["included"].lower()=="true" for r in rows)!=182):
        raise ValueError("Original 51-island/194-port WPI cohort invalid")
    byid={str(int(r["port_id"])):r for r in rows}
    parent=json.loads(SOURCE_SAMPLE.read_text(encoding="utf-8"))
    if parent.get("status")!="FIVE_ORIGINAL_WPI_SOURCE_PORT_PROVIDER_REGISTRY_SAMPLES_COMPLETE":
        raise ValueError("Prior exact-source WPI-to-PortWatch bounded pilot not verified")
    samples=parent.get("source_check_rows") or []
    if len(samples)!=5 or parent.get("original_WPI_sources_with_verified_provider_portid")!=4:
        raise ValueError("The frozen enriched WPI pilot has changed")
    source={}
    for row in samples:
        wpi=row["original_WPI_port_id"]
        if wpi in EXPECTED:
            pw,landmass,name,country=EXPECTED[wpi]
            original=byid[wpi]
            if (row.get("verified_provider_portid")!=pw or
                row.get("source_port_identity_verification_status")!="ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH" or
                row.get("provider_registry_sample_complete") is not True or
                original["island_id"]!=landmass or original["port_name"]!=name or
                original["country_code"]!=country):
                raise ValueError("One of the four original source-verified provider IDs no longer matches WPI")
            candidates=[r for r in row["all_provider_registry_sample_rows"] if str(r.get("portid"))==pw]
            if len(candidates)!=1 or candidates[0]["source_name_country_distance_all_match"] is not True:
                raise ValueError("Original WPI-to-provider source port identity audit missing exact registry entry")
            registry=candidates[0]
            source[wpi]={
               "wpi_id":wpi,"provider_id":pw,"island_id":landmass,
               "wpi_name":name,"wpi_country":country,
               "wpi_lat":float(original["original_latitude"]),
               "wpi_lon":float(original["original_longitude"]),
               "portwatch_registry_lat":float(registry["lat"]),
               "portwatch_registry_lon":float(registry["lon"]),
               "provider_ISO3":str(registry["ISO3"]),
            }
    if set(source)!=set(EXPECTED):
        raise ValueError("One of the 4 source-selected original port registry IDs is absent")
    edges=csv_rows(ZENODO_VERIFIED)
    if len(edges)!=94 or len({(r["origin_WPI_port_id"],r["dest_WPI_port_id"]) for r in edges})!=94:
        raise ValueError("Previous verified 2020 WPI 94 port edges source missing")
    2020_unused = 1
    old={(r["origin_WPI_port_id"],r["dest_WPI_port_id"]):int(r["source_observed_2020_triplet_count"])
         for r in edges}
    if any(n<=0 for n in old.values()):
        raise ValueError("Source confirmed 2020 matched OD count should always be positive")
    return source,old

def fetch_public(params):
    # All requests are exact-original-PortWatch-ID count-only or <=1 row samples.
    url=MODELED_LAYER+"/query?"+urlencode(params)
    attempts=[]
    for step in range(3):
        time.sleep(.16)
        req=Request(url,headers={"Accept":"application/json",
            "User-Agent":"chocho-20261107-four-verified-portwatch-WPI-source-edge-pilot/1.0"})
        try:
            with urlopen(req,timeout=35) as response:
                chunk=response.read(350001)
                if len(chunk)>350000:
                    raise ValueError("Source-only provider sample unexpectedly large")
                response_data=json.loads(chunk.decode("utf-8"))
                if not isinstance(response_data,dict) or response_data.get("error"):
                    raise ValueError("PortWatch layer returned error body "+str(response_data.get("error"))[:240])
                attempts.append({"status":response.status,"error":None,"bytes":len(chunk)})
                return response_data,{"success":True,"attempts":attempts}
        except HTTPError as e:
            attempts.append({"status":e.code,"error":str(e.reason)[:170]})
            if e.code not in RETRYABLE:
                break
            if step<2:time.sleep(1+step)
        except (URLError,TimeoutError,ConnectionError,ValueError) as e:
            attempts.append({"status":None,"error":str(e)[:190]})
            if isinstance(e,ValueError):
                break
            if step<2:time.sleep(1+step)
    return None,{"success":False,"attempts":attempts}

def main(outdir):
    orig,previous=fixed_sources()
    output=Path(outdir).resolve()
    if any(output.is_relative_to(ROOT/k) for k in ("results","data","provenance")):
        raise ValueError("Never overwrite original WPI source or earlier research")
    output.mkdir(parents=True,exist_ok=True)
    wpi_ids=tuple(EXPECTED)
    records=[]
    requests_used=0
    samples_used=0
    for a in wpi_ids:
        for b in wpi_ids:
            if a==b:continue
            query="from_id='"+orig[a]["provider_id"]+"' AND to_id='"+orig[b]["provider_id"]+"'"
            response,receipt=fetch_public({"where":query,"returnCountOnly":"true","f":"json"})
            requests_used+=1
            count=response.get("count") if isinstance(response,dict) else None
            complete=bool(receipt["success"] and type(count)==int and count>=0)
            examples=[]
            example_receipt=None
            if complete and count>0:
                sample,example_receipt=fetch_public({
                    "where":query,
                    "outFields":",".join(OUT_FIELDS),
                    "resultRecordCount":"1","returnGeometry":"false","f":"json"})
                if not (example_receipt["success"] and isinstance(sample,dict)
                        and isinstance(sample.get("features"),list) and len(sample["features"])==1):
                    complete=False
                else:
                    feature=sample["features"][0].get("attributes") or {}
                    if (feature.get("from_id")!=orig[a]["provider_id"] or
                        feature.get("to_id")!=orig[b]["provider_id"]):
                        complete=False
                    else:
                        try:
                            d1=geo_km(feature["from_lat"],feature["from_lon"],
                                      orig[a]["wpi_lat"],orig[a]["wpi_lon"])
                            d2=geo_km(feature["to_lat"],feature["to_lon"],
                                      orig[b]["wpi_lat"],orig[b]["wpi_lon"])
                        except (KeyError,ValueError,TypeError):
                            complete=False
                            d1=d2=None
                        if d1 is None or d1>25 or d2>25:
                            complete=False
                        examples=[{"from_id":feature.get("from_id"),
                          "to_id":feature.get("to_id"),
                          "from_name":feature.get("from_name"),
                          "to_name":feature.get("to_name"),
                          "from_iso3":feature.get("from_iso3"),
                          "to_iso3":feature.get("to_iso3"),
                          "distance":feature.get("distance"),
                          "transit_days":feature.get("transit_days"),
                          "from_source_WPI_geodesic_km":d1,
                          "to_source_WPI_geodesic_km":d2}]
                        samples_used+=1
            records.append({
               "original_WPI_origin":a,"original_WPI_destination":b,
               "original_from_island":orig[a]["island_id"],
               "original_to_island":orig[b]["island_id"],
               "portwatch_from_id":orig[a]["provider_id"],
               "portwatch_to_id":orig[b]["provider_id"],
               "portwatch_static_model_directed_record_count":count if complete else None,
               "PortWatch_count_query_complete":complete,
               "provider_count_request_receipt":receipt,
               "provider_sample_request_receipt":example_receipt,
               "provider_at_most_one_model_sample":examples,
               "verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED":
                     previous.get((a,b),0),
               "2020_missing_source_triplet_not_actual_no_ship":True,
               "portwatch_simulator_static_edge_not_actual_ship_journey":True,
            })
            print("PORTWATCH_FOUR_DIRECTED_MODEL_EDGE="+json.dumps({
                "origin":a,"destination":b,"original_Zenodo_matched_source_triplets":previous.get((a,b),0),
                "model_records":count if complete else None,"query_complete":complete
            },sort_keys=True),flush=True)
    degrees=[]
    for a in wpi_ids:
        for direction in ("from_id","to_id"):
            where=direction+"='"+orig[a]["provider_id"]+"'"
            response,receipt=fetch_public({"where":where,"returnCountOnly":"true","f":"json"})
            requests_used+=1
            count=response.get("count") if isinstance(response,dict) else None
            complete=receipt["success"] and type(count)==int and count>=0
            degrees.append({"wpi_id":a,"island_id":orig[a]["island_id"],
              "portwatch_id":orig[a]["provider_id"],"model_direction_field":direction,
              "provider_static_model_record_count":count if complete else None,
              "count_query_complete":bool(complete),"source_request_receipt":receipt})
    if requests_used!=BUDGET_COUNTS or len(records)!=12 or len(degrees)!=8 or samples_used>MAX_FEATURE_EXAMPLES:
        raise ValueError("Source-only 12 pairs plus 8 degrees query budget violated")
    if any(not x["PortWatch_count_query_complete"] for x in records):
        status="SOME_BOUNDED_PROVIDER_PAIR_QUERIES_INCOMPLETE_NO_ZERO_INFERRED"
    elif any(not x["count_query_complete"] for x in degrees):
        status="ALL_12_PAIR_QUERIES_BUT_SOME_DEGREE_QUERIES_INCOMPLETE"
    else:
        status="ALL_12_FIXED_DIRECTED_MODEL_PAIR_QUERIES_AND_8_DEGREES_COMPLETED"
    report={
      "schema":"chocho_original_WPI4_portwatch_12_directed_static_model_edge_pilot_v1",
      "status":status,"source_utc":datetime.now(timezone.utc).isoformat(),
      "source_provider":"International Monetary Fund IMF PortWatch",
      "source_network_layer":MODELED_LAYER,
      "source_attribution":"Source: International Monetary Fund, IMF PortWatch",
      "original_WPI_git_blob":WPI_GIT_BLOB,
      "source_enriched_WPI_port_count":4,"original_WPI_port_count":194,
      "original_island_count":51,"original_182_marine_qualified_port_count":182,
      "verified_4_WPI_source_registry_ids":[orig[k] for k in wpi_ids],
      "original_2020_Zenodo_94_WPI_directed_source_edges_verified":True,
      "paired_source_record_comparisons":records,
      "source_four_port_in_out_degree_count_queries":degrees,
      "n_pair_count_queries":12,"n_provider_degree_count_queries":8,
      "n_sample_feature_rows_requested":samples_used,
      "provider_model_edges_not_2020_individual_AIS_triplets":True,
      "no_155521_model_edge_bulk_download":True,
      "no_causal_island_butterfly_plant_or_host_dispersal_model_fitted":True,
      "original_2020_source_absence_not_no_shipping":True,
    }
    (output/"source_four_port_directed_model_edge_pilot.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("PORTWATCH_FOUR_PORT_LINK_AUDIT_JSON="+json.dumps({
       "status":status,
       "source_WPI_portids":[orig[k]["provider_id"] for k in wpi_ids],
       "n_directed_pair_conditions":len(records),
       "n_model_edges_positive":sum(x["PortWatch_count_query_complete"] and x["portwatch_static_model_directed_record_count"]>0 for x in records),
       "n_old_Zenodo_matched_2020_positive":sum(x["verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED"]>0 for x in records),
       "n_portwatch_model_positive_AND_zenodo_source_positive":sum(x["PortWatch_count_query_complete"] and x["portwatch_static_model_directed_record_count"]>0 and x["verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED"]>0 for x in records),
       "n_source_feature_rows":samples_used,"all_directed_queries_complete":all(x["PortWatch_count_query_complete"] for x in records)
    },sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",required=True,type=Path)
    main(p.parse_args().out_dir)
