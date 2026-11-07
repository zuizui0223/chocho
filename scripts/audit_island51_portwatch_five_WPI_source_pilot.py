#!/usr/bin/env python3
"""Five source-selected original WPI194 port vs provider PortWatch registry sample.

At most ten public ArcGIS queries and 75 source registry rows. No 155521-edge
download, no statistical 194-port coverage estimate or butterfly model.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import unicodedata
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode

ROOT=Path(__file__).resolve().parents[1]
PORTS=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
SOURCE_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
REGISTRY="https://services9.arcgis.com/weJ1QsnbMYJlCHdG/ArcGIS/rest/services/PortWatch_ports_database/FeatureServer/0"
IDS=("28","428","817","1416","1525")
EXPECTED={
  "28":("BELFAST","GB","26"),
  "428":("SANTA CRUZ DE TENERIFE","ES","206"),
  "817":("LIMASSOL","CY","85"),
  "1416":("COLOMBO","LK","31"),
  "1525":("KAHULUI","US","214"),
}
ISO3={"GB":"GBR","ES":"ESP","CY":"CYP","LK":"LKA","US":"USA"}
HALF_BOX_DEGREES=.35
MAX_ROWS_PER_BOX=15
EARTH_KM=6371.0088
MAX_CENTROID_KM=25.

def source():
    b=PORTS.read_bytes()
    if hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()!=SOURCE_BLOB:
        raise ValueError("Original WPI194 source hash modified")
    with PORTS.open(encoding="utf-8",newline="") as f:records=list(csv.DictReader(f))
    if len(records)!=194 or len({r["island_id"] for r in records})!=51 or sum(r["included"].lower()=="true" for r in records)!=182:
        raise ValueError("Original source port/island cohort differs")
    indexed={str(int(r["port_id"])):r for r in records}
    if len(indexed)!=194:
        raise ValueError("Original source port identifier duplication")
    for pid,expected in EXPECTED.items():
        row=indexed.get(pid)
        if row is None or (row["port_name"],row["country_code"],row["island_id"])!=expected:
            raise ValueError("Source-selected five original WPI port identities changed")
    return indexed

def normalize(text):
    normalized=unicodedata.normalize("NFKD",str(text))
    normalized=normalized.encode("ascii","ignore").decode().upper()
    return " ".join(re.findall("[A-Z0-9]+",normalized))

def distance(lat1,lon1,lat2,lon2):
    a,b,c,d=[math.radians(float(x)) for x in (lat1,lon1,lat2,lon2)]
    v=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 2*EARTH_KM*math.asin(math.sqrt(min(1.,max(0.,v))))

def api(params):
    url=REGISTRY+"/query?"+urlencode(params)
    if not url.startswith("https://services9.arcgis.com/weJ1QsnbMYJlCHdG/"):
        raise ValueError("Only original public IMF PortWatch registry may be queried")
    request=Request(url,headers={"User-Agent":"chocho-20261107-five-bounded-source-portwatch-WPI-examples/1.0",
                                "Accept":"application/json"})
    try:
        with urlopen(request,timeout=35) as res:
            b=res.read(500000+1)
            if len(b)>500000:raise ValueError("Sample provider response too large")
        data=json.loads(b)
        if not isinstance(data,dict) or data.get("error"):
            return None,{"kind":"ARCGIS_ERROR","detail":str(data.get("error"))[:250]}
        return data,None
    except (HTTPError,URLError,TimeoutError,ConnectionError,ValueError) as exc:
        return None,{"kind":type(exc).__name__,"detail":str(exc)[:250]}

def main(out_dir):
    indexed=source()
    out=Path(out_dir).resolve()
    if any(out.is_relative_to((ROOT/x).resolve()) for x in ("data","results","provenance")):
        raise ValueError("Do not overwrite original WPI and biological files")
    out.mkdir(parents=True,exist_ok=True)
    checks=[]
    for pid in IDS:
        w=indexed[pid]
        lat=float(w["original_latitude"]);lon=float(w["original_longitude"])
        where=(f"lat >= {lat-HALF_BOX_DEGREES:.10f} AND lat <= {lat+HALF_BOX_DEGREES:.10f} "
               f"AND lon >= {lon-HALF_BOX_DEGREES:.10f} AND lon <= {lon+HALF_BOX_DEGREES:.10f}")
        count,error=api({"where":where,"returnCountOnly":"true","f":"json"})
        count_n=count.get("count") if isinstance(count,dict) else None
        sample,sample_error=api({
           "where":where,
           "outFields":"portid,portname,country,ISO3,lat,lon,LOCODE",
           "resultRecordCount":str(MAX_ROWS_PER_BOX),
           "returnGeometry":"false","f":"json"
        })
        response=sample.get("features") if isinstance(sample,dict) else None
        items=[r.get("attributes") or {} for r in (response or [])]
        # Never treat a source truncated registry query as a complete identity search.
        complete=(isinstance(count_n,int) and count_n<=MAX_ROWS_PER_BOX and
                  count_n==len(items) and error is None and sample_error is None)
        match=[]
        diagnostics=[]
        for x in items:
            try:
                dist=distance(lat,lon,x.get("lat"),x.get("lon"))
            except (TypeError,ValueError):
                dist=None
            cname=normalize(x.get("portname") or "")
            source_name=normalize(w["port_name"])
            exact=cname==source_name and bool(cname)
            country_match=str(x.get("ISO3") or "")==ISO3[w["country_code"]]
            eligible=bool(exact and country_match and dist is not None and dist<=MAX_CENTROID_KM)
            if eligible:
                match.append((str(x["portid"]),dist))
            diagnostics.append({
              "portid":x.get("portid"),
              "portname":x.get("portname"),
              "ISO3":x.get("ISO3"),
              "LOCODE":x.get("LOCODE"),
              "lat":x.get("lat"),"lon":x.get("lon"),
              "original_WPI_distance_km":dist,
              "source_canonical_name_exact":exact,
              "original_WPI_countryISO3_match":country_match,
              "source_name_country_distance_all_match":eligible,
            })
        if not complete:status="INCOMPLETE_OR_FAILED_BOUNDED_PUBLIC_REGISTRY_QUERY"
        elif len(match)==1:status="ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH"
        elif len(match)>1:status="AMBIGUOUS_MULTIPLE_EXACT_PROVIDER_PORTS_IN_BOX"
        else:status="NO_EXACT_NAME_ISO3_WITHIN25KM_IN_BOUNDED_PROVIDER_REGISTRY"
        matched=match[0] if status=="ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH" else None
        checks.append({
          "original_WPI_port_id":pid,
          "original_WPI_port_name":w["port_name"],
          "original_WPI_country_ISO2":w["country_code"],
          "original_WPI_expected_provider_ISO3":ISO3[w["country_code"]],
          "original_WPI_island_id":w["island_id"],
          "original_WPI_lat":lat,"original_WPI_lon":lon,
          "requested_ArcGIS_where":where,
          "provider_box_returnCountOnly":count_n,
          "source_reg_count_query_error":error,
          "source_reg_sample_query_error":sample_error,
          "provider_sample_row_count":len(items),
          "provider_registry_sample_complete":complete,
          "all_provider_registry_sample_rows":diagnostics,
          "source_port_identity_verification_status":status,
          "verified_provider_portid":matched[0] if matched else None,
          "verified_provider_geodesic_km":matched[1] if matched else None,
          "retrospective_enriched_source_only_pilot_NOT_random_194_port_sample":True,
          "source_case_NOT_actual_ship_voyages_or_51_island_butterfly_prediction":True,
        })
        print("PORTWATCH_WPI5_SOURCE_PILOT_PROGRESS="+json.dumps({
           "original_WPI_port_id":pid,"original_name":w["port_name"],
           "box_source_count":count_n,"sample_rows":len(items),"status":status,
           "source_provider_portid":matched[0] if matched else None
        },sort_keys=True),flush=True)
    receipt={
      "schema":"chocho_portwatch_original_WPI194_five_preselected_source_name_country_geo_pilot_v1",
      "status":"FIVE_ORIGINAL_WPI_SOURCE_PORT_PROVIDER_REGISTRY_SAMPLES_COMPLETE",
      "source_provider":"International Monetary Fund IMF PortWatch",
      "source_registry_URL":REGISTRY,
      "source_checked_utc":datetime.now(timezone.utc).isoformat(),
      "source_WPI_original_git_blob":SOURCE_BLOB,
      "original_WPI_source_ports":194,"original_source_islands":51,"marine_source_ports":182,
      "source_preselected_WPI_ids":list(IDS),
      "source_WPI_sample_size":len(IDS),
      "provider_registry_queries_made":2*len(IDS),
      "maximum_source_registry_feature_rows_per_original_WPI_box":MAX_ROWS_PER_BOX,
      "source_record_example_rows_returned_total":sum(x["provider_sample_row_count"] for x in checks),
      "original_WPI_sources_with_verified_provider_portid":sum(bool(x["verified_provider_portid"]) for x in checks),
      "original_WPI_sources_with_truncated_or_failed_provider_query":sum(not x["provider_registry_sample_complete"] for x in checks),
      "source_check_rows":checks,
      "no_full_global_port_registry_or_155521_model_edges_downloaded":True,
      "source_pilot_is_selected_not_random_port_coverage_estimation":True,
      "provider_model_links_not_AIS_vessel_OD_records":True,
      "not_shipping_causal_inference_no_butterfly_model":True,
    }
    (out/"WPI5_bounded_registry_source_pilot_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PORTWATCH_WPI5_SOURCE_PILOT_JSON="+json.dumps({
      "source_original_WPI_IDs":list(IDS),
      "verified_provider_portids":{x["original_WPI_port_id"]:x["verified_provider_portid"] for x in checks},
      "source_boxes_complete":sum(x["provider_registry_sample_complete"] for x in checks),
      "provider_sample_rows":receipt["source_record_example_rows_returned_total"],
      "verified_count":receipt["original_WPI_sources_with_verified_provider_portid"],
      "no_bulk_provider_records":True},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,required=True)
    main(p.parse_args().out_dir)
