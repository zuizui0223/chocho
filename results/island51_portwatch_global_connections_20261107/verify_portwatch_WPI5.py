#!/usr/bin/env python3
"""Independently verify five original WPI → PortWatch port registry example identities.

No imports of original provider query/matching modules and no network requests.
Keep the intentionally enriched five WPI ports distinct from a random 194-port census.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import re
import unicodedata
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
RESULT=HERE/"WPI5_bounded_registry_source_pilot_receipt.json"
BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
FROZEN=("28","428","817","1416","1525")
COUNTRY={"GB":"GBR","ES":"ESP","CY":"CYP","LK":"LKA","US":"USA"}

def clean(value):
    base=unicodedata.normalize("NFKD",str(value)).encode("ascii","ignore").decode().upper()
    return " ".join(re.findall(r"[A-Z0-9]+",base))

def geo(a,b,c,d):
    p,q,r,s=[math.radians(float(x)) for x in (a,b,c,d)]
    z=math.sin((r-p)/2)**2+math.cos(p)*math.cos(r)*math.sin((s-q)/2)**2
    return 2*6371.0088*math.asin(math.sqrt(max(0.,min(1.,z))))

def verify():
    b=WPI.read_bytes()
    if hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()!=BLOB:
        raise ValueError("Original WPI source Git blob changed")
    with WPI.open(newline="",encoding="utf-8") as f:
        original=list(csv.DictReader(f))
    if len(original)!=194 or len({x["island_id"] for x in original})!=51 or sum(x["included"].lower()=="true" for x in original)!=182:
        raise ValueError("Original 194 WPI biological/182 marine/51 island source not intact")
    indexed={str(int(r["port_id"])):r for r in original}
    result=json.loads(RESULT.read_text())
    if result["status"]!="FIVE_ORIGINAL_WPI_SOURCE_PORT_PROVIDER_REGISTRY_SAMPLES_COMPLETE":
        raise ValueError("PortWatch original source sample not complete")
    if (tuple(result["source_preselected_WPI_ids"])!=FROZEN or
        len(result["source_check_rows"])!=5 or
        result["maximum_source_registry_feature_rows_per_original_WPI_box"]!=15 or
        result["provider_registry_queries_made"]!=10 or
        result["original_WPI_source_ports"]!=194 or
        result["original_source_islands"]!=51):
        raise ValueError("Source-only WPI5 original model-free protocol changed")
    count=0;verified=0;incomplete=0;positions=set()
    for row in result["source_check_rows"]:
        id=str(row["original_WPI_port_id"])
        if id not in FROZEN or id in positions:
            raise ValueError("Duplicate or non-source-selected original WPI port")
        positions.add(id)
        w=indexed[id]
        expected={
            "original_WPI_port_name":w["port_name"],
            "original_WPI_country_ISO2":w["country_code"],
            "original_WPI_expected_provider_ISO3":COUNTRY[w["country_code"]],
            "original_WPI_island_id":w["island_id"],
        }
        for key,value in expected.items():
            if row[key]!=value:
                raise ValueError("Original named WPI port or source country changed: "+key)
        a=float(w["original_latitude"]);b=float(w["original_longitude"])
        if not(math.isclose(a,row["original_WPI_lat"],abs_tol=1e-12) and
               math.isclose(b,row["original_WPI_lon"],abs_tol=1e-12)):
            raise ValueError("Original WPI centroids have been replaced")
        bounds=(float(a-.35),float(a+.35),float(b-.35),float(b+.35))
        where=(f"lat >= {bounds[0]:.10f} AND lat <= {bounds[1]:.10f} "
               f"AND lon >= {bounds[2]:.10f} AND lon <= {bounds[3]:.10f}")
        if row["requested_ArcGIS_where"]!=where:
            raise ValueError("Original source-only geographic box query changed")
        n=row["provider_box_returnCountOnly"]
        items=row["all_provider_registry_sample_rows"]
        if not isinstance(n,int) or n<0 or len(items)>15 or row["provider_sample_row_count"]!=len(items):
            raise ValueError("Source registry query counts not bounded")
        complete=n==len(items)
        if row["provider_registry_sample_complete"]!=complete:
            raise ValueError("Truncated PortWatch query was treated as complete")
        if complete and (row["source_reg_count_query_error"] or row["source_reg_sample_query_error"]):
            raise ValueError("Public PortWatch source returned error despite reported complete query")
        incomplete+=int(not complete)
        count+=len(items)
        matching=[]
        for source in items:
            lat=source["lat"];lon=source["lon"]
            if not(-90<=float(lat)<=90 and -180<=float(lon)<=180):
                raise ValueError("Provider port source outside valid Earth coordinates")
            d=geo(a,b,lat,lon)
            name_ok=clean(source["portname"])==clean(w["port_name"]) and bool(clean(source["portname"]))
            country_ok=source["ISO3"]==COUNTRY[w["country_code"]]
            okay=bool(d<=25 and name_ok and country_ok)
            if not math.isclose(d,source["original_WPI_distance_km"],abs_tol=1e-8,rel_tol=1e-10):
                raise ValueError("Source-to-provider port distance not independently reproducible")
            for field,exp in (("source_canonical_name_exact",name_ok),
                              ("original_WPI_countryISO3_match",country_ok),
                              ("source_name_country_distance_all_match",okay)):
                if source[field]!=exp:
                    raise ValueError("Provider port source example classification altered: "+field)
            if okay:matching.append((source["portid"],d))
        expected_status=("INCOMPLETE_OR_FAILED_BOUNDED_PUBLIC_REGISTRY_QUERY" if not complete
            else "ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH" if len(matching)==1
            else "AMBIGUOUS_MULTIPLE_EXACT_PROVIDER_PORTS_IN_BOX" if len(matching)>1
            else "NO_EXACT_NAME_ISO3_WITHIN25KM_IN_BOUNDED_PROVIDER_REGISTRY")
        if row["source_port_identity_verification_status"]!=expected_status:
            raise ValueError("WPI-to-provider source port identity status not supported")
        match=matching[0] if expected_status=="ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH" else None
        if row["verified_provider_portid"]!=(match[0] if match else None):
            raise ValueError("Unverified example was promoted to an IMF PortWatch ID")
        if match:
            if not math.isclose(row["verified_provider_geodesic_km"],match[1],abs_tol=1e-8):
                raise ValueError("WPI-to-provider saved distance differs")
            verified+=1
        elif row["verified_provider_geodesic_km"] is not None:
            raise ValueError("No matching provider ID may contain a verified geodesic")
        if row["retrospective_enriched_source_only_pilot_NOT_random_194_port_sample"] is not True:
            raise ValueError("Unrepresentative 5-port sample improperly described as random")
    if (positions!=set(FROZEN) or
        count!=result["source_record_example_rows_returned_total"] or
        verified!=result["original_WPI_sources_with_verified_provider_portid"] or
        incomplete!=result["original_WPI_sources_with_truncated_or_failed_provider_query"]):
        raise ValueError("Source pilot WPI matching totals differ")
    if result["no_full_global_port_registry_or_155521_model_edges_downloaded"] is not True:
        raise ValueError("Source-only metadata access permission restriction not kept")
    report={
      "schema":"chocho_portwatch_five_original_WPI_source_registry_independent_audit_v1",
      "status":"FIVE_ORIGINAL_WPI_AND_PROVIDER_PORT_ID_COORDINATES_VERIFIED_FROM_SAVED_SAMPLE",
      "source_original_ports":194,"original_islands":51,
      "source_selected_WPI_pilot_port_ids":list(FROZEN),
      "source_registry_example_rows_rechecked":count,
      "verified_WPI_to_provider_ID_pilot_examples":verified,
      "source_registration_not_found_or_ambiguous_examples":5-verified,
      "incomplete_geographic_box_queries":incomplete,
      "source_pilot_not_194_port_population_estimate":True,
      "no_full_155521_provider_network_rows_or_butterfly_models":True,
    }
    print("PORTWATCH_WPI5_INDEPENDENT_AUDIT_JSON="+json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    verify()
