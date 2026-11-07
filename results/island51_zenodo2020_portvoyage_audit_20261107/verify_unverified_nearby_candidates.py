#!/usr/bin/env python3
"""Independent source-only check of Zenodo2020-to-WPI194 matching and 27 nearby UNVERIFIED candidates.

Uses ONLY source-archived public port ledgers and original pinned 194 WPI coordinates.
Never promotes a geographical-neighbor candidate to an observed direct voyage.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
WPI_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
EARTH_R=6371.0088
MAX_MATCH_KM=25.

FIELDS=(
    "rank","external_port_id_NOT_WPI","external_port_name_UNVERIFIED",
    "external_normalized_port_name","source_lon","source_lat",
    "nearest_WPI_id_NOT_VERIFIED","nearest_WPI_name_NOT_VERIFIED",
    "nearest_WPI_country","nearest_original_island_id","geodesic_km",
    "within_5km","source_departure_events","source_arrival_events",
    "match_status","identity_verified"
)

def read(path):
    with Path(path).open("rt",encoding="utf-8",newline="") as stream:
        return list(csv.DictReader(stream))

def normalize(original):
    transliterated=unicodedata.normalize("NFKD",str(original))
    transliterated=transliterated.encode("ascii","ignore").decode("ascii").upper()
    words=re.findall("[A-Z0-9]+",transliterated)
    if not words or " ".join(words) in {"NAN","NONE","NULL","NA","UNKNOWN","?","-"}:
        return ""
    if len(words)>2 and words[:2]==["PORT","OF"]:
        words=words[2:]
    elif len(words)>1 and words[0]=="PORT":
        words=words[1:]
    while len(words)>1 and words[-1] in ("PORT","HARBOR","HARBOUR","TERMINAL","DOCK","QUAY"):
        words.pop()
    if not words or (len(words)==1 and words[0] in {"PORT","HARBOR","HARBOUR","TERMINAL","DOCK","QUAY"}):
        return ""
    return " ".join(words)

def haversine(lon_a,lat_a,lon_b,lat_b):
    a,b,c,d=[math.radians(float(v)) for v in (lat_a,lon_a,lat_b,lon_b)]
    h=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 2*EARTH_R*math.asin(math.sqrt(min(1.,max(0.,h))))

def assert_near(a,b,label):
    if not (math.isfinite(float(a)) and math.isfinite(float(b)) and
            math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-8)):
        raise ValueError(label+": "+str(a)+" is not "+str(b))

def verify(repo_results=HERE):
    output=Path(repo_results)
    source=WPI.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(source)).encode()+b"\0"+source).hexdigest()
    if blob!=WPI_GIT_BLOB:
        raise ValueError("Original WPI194 pinned source was modified")
    original=read(WPI)
    if len(original)!=194 or len({r["port_id"] for r in original})!=194 or len({r["island_id"] for r in original})!=51:
        raise ValueError("WPI original 194 ports / 51 islands incomplete")
    ports={str(int(r["port_id"])):r for r in original}
    by_name=defaultdict(list)
    for pid,r in ports.items():
        by_name[normalize(r["port_name"])].append(pid)
    original_receipt=json.loads((output/"source_matched_voyage_receipt.json").read_text(encoding="utf-8"))
    if original_receipt["status"]!="2020_ZENODO_SOURCE_WPI194_IDENTITIES_AUDITED_NO_BIOLOGICAL_MODELS":
        raise ValueError("Direct 2020 source producer receipt not verified")
    raw=read(output/"external2020_endpoints_to_WPI194_crosswalk.csv")
    if len(raw)!=1923:
        raise ValueError("Expected complete original external 1923 port endpoint list")
    status=Counter()
    closest=[]
    keys=set()
    for row in raw:
        uid=row["external_endpoint_hash"]
        if uid in keys:
            raise ValueError("Duplicate original external endpoint hash")
        keys.add(uid)
        name=normalize(row["source_port_name"])
        if row["normalized_port_name"]!=name:
            raise ValueError("External endpoint normalized port-name provenance differs")
        distances=[]
        for pid,w in ports.items():
            d=haversine(row["source_lon"],row["source_lat"],w["original_longitude"],w["original_latitude"])
            distances.append((d,int(pid),pid))
        distances.sort()
        distance_to_nearest,_,nearest=distances[0]
        if not name:
            # Source matcher intentionally refuses even geographic diagnostics
            # for missing/unknown port names, rather than issuing a nearest
            # WPI suggestion that a reviewer could mistake for verified ID.
            if (row["geography_only_nearest_WPI_id_NOT_VERIFIED"] or
                row["geography_only_nearest_distance_km"] or
                row["geo_only_5km_NOT_VERIFIED"].lower()!="false" or
                row["geo_only_10km_NOT_VERIFIED"].lower()!="false"):
                raise ValueError("Unnamed source endpoints must have NO geographic WPI suggestion")
        else:
            if nearest!=row["geography_only_nearest_WPI_id_NOT_VERIFIED"]:
                raise ValueError("Original nearest source WPI diagnostic ID wrong")
            assert_near(distance_to_nearest,row["geography_only_nearest_distance_km"],"nearest WPI haversine")
            if (row["geo_only_5km_NOT_VERIFIED"].lower()=="true")!=(distance_to_nearest<=5.):
                raise ValueError("Original 5-km geography diagnostic mismatch")
            if (row["geo_only_10km_NOT_VERIFIED"].lower()=="true")!=(distance_to_nearest<=10.):
                raise ValueError("Original 10-km geography diagnostic mismatch")
        candidates=sorted(
            ((d,int(pid),pid) for d,n,pid in distances if pid in by_name.get(name,[]) and d<=MAX_MATCH_KM),
            key=lambda z:(z[0],z[1])
        ) if name else []
        # The semantic port-identity rule MUST be checked independently of position.
        expected=("NO_RELIABLE_NAME" if not name else
                  "PRIMARY_NAME_AND_COORD_VERIFIED" if len(candidates)==1 else
                  "AMBIGUOUS_MULTIPLE_MATCHED_WPI_NAMES_NEARBY" if len(candidates)>1 else
                  "NO_EXACT_NAME_AND_COORD_MATCH")
        if row["primary_match_status"]!=expected:
            raise ValueError("Independent original exact-port-name + 25km identity rule failed")
        pid=row["verified_original_WPI_port_id"]
        if bool(pid)!=(expected=="PRIMARY_NAME_AND_COORD_VERIFIED"):
            raise ValueError("An unverified original node was promoted to WPI")
        if pid and pid!=candidates[0][2]:
            raise ValueError("An unrelated original WPI port was used as a verified endpoint")
        if pid:
            assert_near(row["original_port_identity_distance_km"],candidates[0][0],"verified geodesic")
        status[expected]+=1
        if expected=="NO_EXACT_NAME_AND_COORD_MATCH" and distance_to_nearest<=10.:
            closest.append((distance_to_nearest,float(row["external_port_id"]),row,ports[nearest],nearest))
    if status!=Counter(original_receipt["match_status_nodes"]):
        raise ValueError("Original source crosswalk node-status counts disagree")
    if len(closest)!=27 or sum(x[0]<=5 for x in closest)!=18:
        raise ValueError("Source match-failure neighborhood count not 27 within 10km / 18 within 5km")
    closest.sort(key=lambda x:(x[0],x[1]))
    archive=read(output/"source_only_10km_unverified_port_candidates.csv")
    if len(archive)!=27 or len({r["external_port_id_NOT_WPI"] for r in archive})!=27:
        raise ValueError("Unverified candidate audit rows missing/duplicated")
    for i,((distance,ext_number,ext,wpi,pid),saved) in enumerate(zip(closest,archive),1):
        expectations={
            "rank":i,"external_port_id_NOT_WPI":ext["external_port_id"],
            "external_port_name_UNVERIFIED":ext["source_port_name"],
            "external_normalized_port_name":ext["normalized_port_name"],
            "source_lon":ext["source_lon"],"source_lat":ext["source_lat"],
            "nearest_WPI_id_NOT_VERIFIED":pid,
            "nearest_WPI_name_NOT_VERIFIED":wpi["port_name"],
            "nearest_WPI_country":wpi["country_code"],
            "nearest_original_island_id":wpi["island_id"],
            "geodesic_km":distance,"within_5km":distance<=5,
            "source_departure_events":int(ext["raw_departures"]),
            "source_arrival_events":int(ext["raw_arrivals"]),
            "match_status":"UNVERIFIED_GEOGRAPHIC_PROXIMITY_ONLY",
            "identity_verified":False,
        }
        for field in FIELDS:
            actual=saved[field];expected=expectations[field]
            if field in ("source_lon","source_lat","geodesic_km"):
                assert_near(actual,expected,"geographic-neighbor diagnostic "+field)
            elif field in ("source_departure_events","source_arrival_events","rank"):
                if int(actual)!=expected:
                    raise ValueError("Unverified candidate source event/index count altered: "+field)
            elif field in ("within_5km","identity_verified"):
                if actual.lower()!=str(expected).lower():
                    raise ValueError("Unverified neighbor cannot be promoted to verified: "+field)
            elif actual!=str(expected):
                raise ValueError("Source-only unmatched geographic review row differs: "+field)
    report={
      "schema":"chocho_2020_AIS_194_WPI_rejected_geographic_match_independent_audit_v1",
      "status":"ALL_1923_SOURCE_NODES_AND_27_UNVERIFIED_NEARBY_WPI_CANDIDATES_VALIDATED",
      "source_nodes_verified":1923,
      "original_WPI_ports_verified":194,
      "primary_exact_name_plus25km_WPI_nodes":status["PRIMARY_NAME_AND_COORD_VERIFIED"],
      "names_not_matched_despite_valid_source_name":status["NO_EXACT_NAME_AND_COORD_MATCH"],
      "unnamed_external_source_nodes":status["NO_RELIABLE_NAME"],
      "unverified_source_names_within_5km_WPI":18,
      "unverified_source_names_within_10km_WPI":27,
      "unverified_nodes_promoted_to_WPI":0,
      "original2020_source_direct_voyage_counts_NOT_recomputed_or_changed":True,
      "no_fuzzy_name_or_nearest_neighbor_ship_route_assumptions":True,
    }
    print("ZENODO_WPI_UNVERIFIED_NEARBY_AUDIT_JSON="+json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--results-dir",type=Path,default=HERE)
    args=parser.parse_args()
    verify(args.results_dir)
