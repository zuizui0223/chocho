#!/usr/bin/env python3
"""Conservatively reconcile a public 2020 AIS origin-destination CSV with 194 original WPI ports.

Independent direct source: Duvallet et al. 2024, Zenodo 10783698.
Exact normalized port NAME plus <=25 km original WPI centroid, unique match only.
Never numerically equate external waiting-area IDs and WPI OBJECTIDs.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter,defaultdict
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import unicodedata
import zipfile

from audit_island51_zenodo2020_voyage_source import (
    ROOT,PORTS,PORTS_GIT_BLOB,ZENODO_MD5,original_port_roster
)

EARTH_RADIUS_KM=6371.0088
PRIMARY_MAX_DISTANCE_KM=25.
GEO_DIAGNOSTIC_CAP_KM=10.
DOI="10.5281/zenodo.10783698"
INPUT_MEMBER="triplet_2020.csv"
HEADER_REQUIRED=(
    "date-out","date-in","Port_ID_waiting_area_id_from","PortName_from",
    "lon_from","lat_from","Port_ID_waiting_area_id_to","PortName_to",
    "lon_to","lat_to",
)
UNKNOWN_NAMES={"","NAN","NA","NONE","NULL","UNKNOWN","?","-"}
NONINFORMATIVE={"PORT","HARBOR","HARBOUR","TERMINAL","DOCK","QUAY"}

def name_key(name):
    raw=unicodedata.normalize("NFKD",str(name))
    raw=raw.encode("ascii","ignore").decode("ascii").upper()
    tokens=re.findall(r"[A-Z0-9]+",raw)
    if not tokens or " ".join(tokens) in UNKNOWN_NAMES:
        return ""
    if len(tokens)>2 and tokens[:2]==["PORT","OF"]:
        tokens=tokens[2:]
    elif len(tokens)>1 and tokens[0]=="PORT":
        tokens=tokens[1:]
    while len(tokens)>1 and tokens[-1] in NONINFORMATIVE:
        tokens=tokens[:-1]
    if not tokens or (len(tokens)==1 and tokens[0] in NONINFORMATIVE):
        return ""
    return " ".join(tokens)

def km(lat1,lon1,lat2,lon2):
    a,b,c,d=[math.radians(float(x)) for x in (lat1,lon1,lat2,lon2)]
    v=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return EARTH_RADIUS_KM*2*math.asin(math.sqrt(max(0.,min(1.,v))))

def original():
    rows=original_port_roster()
    points={}
    by_name=defaultdict(list)
    for row in rows:
        pid=str(int(row["port_id"]))
        key=name_key(row["port_name"])
        if not key:
            raise ValueError("Original WPI has unmatchable port name: "+str(row["port_name"]))
        node={"port_id":pid,"island_id":str(row["island_id"]),
              "port_name":str(row["port_name"]),
              "country":str(row["country_code"]),
              "name_key":key,
              "latitude":float(row["original_latitude"]),
              "longitude":float(row["original_longitude"]),
              "marine_eligible":str(row["included"]).lower()=="true"}
        points[pid]=node
        by_name[key].append(node)
    if len(points)!=194 or len({x["island_id"] for x in points.values()})!=51:
        raise ValueError("Original frozen WPI source ports or island population differs")
    return points,by_name

def source_endpoint(raw_id,raw_name,raw_lon,raw_lat):
    try:
        lon=float(raw_lon);lat=float(raw_lat)
    except (ValueError,TypeError):
        return None
    if not(math.isfinite(lat) and math.isfinite(lon)) or not(-90<=lat<=90 and -180<=lon<=180):
        return None
    name=str(raw_name).strip()
    aid=str(raw_id).strip()
    norm=name_key(name)
    if not norm:
        # An external unknown waiting-area endpoint can still appear in raw source,
        # but cannot be assigned to an original WPI port from numeric ID alone.
        return (aid,name,round(lon,7),round(lat,7))
    return (aid,name,round(lon,7),round(lat,7))

def read_external(source_zip):
    raw=Path(source_zip).read_bytes()
    if hashlib.md5(raw).hexdigest()!=ZENODO_MD5:
        raise ValueError("Publisher original 2020 AIS archive MD5 does not match")
    states=Counter()
    origin=Counter()
    dest=Counter()
    raw_pairs=Counter()
    from_coverage=Counter()
    to_coverage=Counter()
    min_date=None
    max_date=None
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=[name for name in z.namelist() if not name.endswith("/")]
        if members!=[INPUT_MEMBER]:
            raise ValueError("Original Zenodo source member table structure changed")
        with z.open(INPUT_MEMBER) as stream:
            text=io.TextIOWrapper(stream,encoding="utf-8-sig",errors="replace",newline="")
            reader=csv.DictReader(text,delimiter=";")
            if not set(HEADER_REQUIRED).issubset(set(reader.fieldnames or [])):
                raise ValueError("Original public 2020 AIS origin/destination schema changed")
            for row in reader:
                states["source_rows"]+=1
                try:
                    departure=datetime.datetime.fromisoformat(str(row["date-out"]).strip())
                    arrival=datetime.datetime.fromisoformat(str(row["date-in"]).strip())
                except (TypeError,ValueError):
                    states["invalid_datetime"]+=1
                    continue
                if departure.year!=2020:
                    states["departure_not_2020"]+=1
                    continue
                if arrival<departure:
                    states["arrived_before_departure"]+=1
                    continue
                a=source_endpoint(row["Port_ID_waiting_area_id_from"],row["PortName_from"],
                                  row["lon_from"],row["lat_from"])
                b=source_endpoint(row["Port_ID_waiting_area_id_to"],row["PortName_to"],
                                  row["lon_to"],row["lat_to"])
                if a is None or b is None:
                    states["invalid_port_coordinate"]+=1
                    continue
                # The source table already claims a from/to journey; do not
                # reconstruct intermediate calls from vessel IDs.
                states["valid_2020_source_rows"]+=1
                origin[a]+=1
                dest[b]+=1
                raw_pairs[(a,b)]+=1
                min_date=min(min_date,departure) if min_date else departure
                max_date=max(max_date,departure) if max_date else departure
    if not states["source_rows"] or not states["valid_2020_source_rows"]:
        raise ValueError("No parsed 2020 original source voyages")
    return states,origin,dest,raw_pairs,min_date,max_date

def map_node(node,original_points,by_name):
    eid,source_name,lon,lat=node
    key=name_key(source_name)
    if not key:
        return {"status":"NO_RELIABLE_NAME","verified_wpi_id":None,"distance_km":None,
                "nearest_wpi_id":None,"nearest_km":None,"name_key":key}
    candidates=[]
    for wpi in by_name.get(key,[]):
        d=km(lat,lon,wpi["latitude"],wpi["longitude"])
        if d<=PRIMARY_MAX_DISTANCE_KM:
            candidates.append((d,wpi["port_id"]))
    candidates.sort(key=lambda z:(z[0],int(z[1])))
    nearest=sorted(((km(lat,lon,wpi["latitude"],wpi["longitude"]),wpi["port_id"])
                   for wpi in original_points.values()),key=lambda t:(t[0],int(t[1])))
    diagnostic=nearest[0] if nearest else (None,None)
    if len(candidates)==1:
        return {"status":"PRIMARY_NAME_AND_COORD_VERIFIED","verified_wpi_id":candidates[0][1],
                "distance_km":candidates[0][0],"nearest_wpi_id":diagnostic[1],
                "nearest_km":diagnostic[0],"name_key":key}
    if len(candidates)>1:
        return {"status":"AMBIGUOUS_MULTIPLE_MATCHED_WPI_NAMES_NEARBY","verified_wpi_id":None,
                "distance_km":None,"nearest_wpi_id":diagnostic[1],
                "nearest_km":diagnostic[0],"name_key":key}
    return {"status":"NO_EXACT_NAME_AND_COORD_MATCH","verified_wpi_id":None,
            "distance_km":None,"nearest_wpi_id":diagnostic[1],
            "nearest_km":diagnostic[0],"name_key":key}

def save_csv(path,rows,cols):
    with Path(path).open("wt",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def run(zip_path,out_dir):
    output=Path(out_dir).resolve()
    if any(output.is_relative_to((ROOT/s).resolve()) for s in ("data","provenance","results")):
        raise ValueError("Never overwrite previously frozen butterfly source/results")
    wp,by_name=original()
    states,origins,destinations,segments,dt_min,dt_max=read_external(zip_path)
    nodes=set(origins)|set(destinations)
    lookup={}
    rows=[]
    source_node_digest={}
    for node in sorted(nodes,key=lambda x:(x[0],x[1],x[2],x[3])):
        verdict=map_node(node,wp,by_name)
        stable=hashlib.sha256(json.dumps(node,separators=(",",":")).encode()).hexdigest()[:20]
        lookup[node]=verdict
        source_node_digest[node]=stable
        rows.append({
            "external_endpoint_hash":stable,
            "external_port_id":node[0],
            "source_port_name":node[1],
            "normalized_port_name":verdict["name_key"],
            "source_lon":node[2],"source_lat":node[3],
            "primary_match_status":verdict["status"],
            "verified_original_WPI_port_id":verdict["verified_wpi_id"] or "",
            "original_port_identity_distance_km":verdict["distance_km"] if verdict["distance_km"] is not None else "",
            "geography_only_nearest_WPI_id_NOT_VERIFIED":verdict["nearest_wpi_id"] or "",
            "geography_only_nearest_distance_km":verdict["nearest_km"] if verdict["nearest_km"] is not None else "",
            "geo_only_5km_NOT_VERIFIED":verdict["nearest_km"] is not None and verdict["nearest_km"]<=5,
            "geo_only_10km_NOT_VERIFIED":verdict["nearest_km"] is not None and verdict["nearest_km"]<=10,
            "raw_departures":origins[node],
            "raw_arrivals":destinations[node],
        })
    original_port_departs=Counter()
    original_port_arrivals=Counter()
    original_port_nodes=defaultdict(set)
    verified_node=0
    for node in nodes:
        match=lookup[node]
        if match["verified_wpi_id"] is not None:
            verified_node+=1
            pid=match["verified_wpi_id"]
            original_port_nodes[pid].add(source_node_digest[node])
            original_port_departs[pid]+=origins[node]
            original_port_arrivals[pid]+=destinations[node]
    wpi_edges=Counter()
    island_edges=Counter()
    source_legs=Counter()
    for (a,b),count in segments.items():
        wa=lookup[a]["verified_wpi_id"]
        wb=lookup[b]["verified_wpi_id"]
        if wa is None or wb is None:
            source_legs["one_or_both_endpoints_unverified"]+=count
            continue
        source_legs["both_endpoints_verified"]+=count
        if wa==wb:
            source_legs["same_WPI_port_both_endpoints"]+=count
            continue
        wpi_edges[(wa,wb)]+=count
        if wp[wa]["island_id"]==wp[wb]["island_id"]:
            source_legs["different_WPI_ports_same_island"]+=count
            continue
        island_edges[(wp[wa]["island_id"],wp[wb]["island_id"])]+=count
        source_legs["cross_original_island_voyages"]+=count
    output.mkdir(parents=True,exist_ok=True)
    save_csv(output/"external2020_endpoints_to_WPI194_crosswalk.csv",rows,list(rows[0]))
    portlist=[]
    for pid in sorted(wp,key=int):
        x=wp[pid]
        portlist.append({
           "original_WPI_port_id":pid,"island_id":x["island_id"],
           "country":x["country"],"source_port_name":x["port_name"],
           "original_marine_eligible":x["marine_eligible"],
           "verified_external_nodes":len(original_port_nodes[pid]),
           "observed_departures_at_verified_2020_endpoint":original_port_departs[pid],
           "observed_arrivals_at_verified_2020_endpoint":original_port_arrivals[pid],
           "verified_external_2020_endpoint":bool(original_port_nodes[pid]),
        })
    save_csv(output/"original194_source_matched_port_coverage.csv",portlist,list(portlist[0]))
    links=[]
    for (a,b),n in sorted(wpi_edges.items(),key=lambda x:(int(x[0][0]),int(x[0][1]))):
        links.append({
           "origin_WPI_port_id":a,"dest_WPI_port_id":b,
           "origin_island_id":wp[a]["island_id"],
           "dest_island_id":wp[b]["island_id"],
           "source_observed_2020_triplet_count":n,
           "same_original_island":wp[a]["island_id"]==wp[b]["island_id"]
        })
    save_csv(output/"observed2020_verified_WPI_port_links.csv",links,
        ["origin_WPI_port_id","dest_WPI_port_id","origin_island_id","dest_island_id",
         "source_observed_2020_triplet_count","same_original_island"])
    island_port_covered={r["island_id"] for r in portlist if r["verified_external_2020_endpoint"]}
    all_ids=sorted({x["island_id"] for x in wp.values()},key=int)
    unord=Counter()
    for (a,b),count in island_edges.items():
        t=tuple(sorted([a,b],key=int))
        unord[t]+=count
    pair_rows=[]
    for i,a in enumerate(all_ids):
        for b in all_ids[i+1:]:
            pair_rows.append({
                "island_a":a,"island_b":b,
                "both_original_islands_have_2020_source_matched_endpoint":a in island_port_covered and b in island_port_covered,
                "direct_observed_2020_source_voyages_in_either_direction":unord[(a,b)],
                "at_least_one_matched_2020_observation_between_source_islands":unord[(a,b)]>0,
                "zero_is_not_ship_journey_absence":True,
            })
    if len(pair_rows)!=1275:
        raise ValueError("Source 51-island unordered pair coverage not complete")
    save_csv(output/"original51_all1275_direct2020_source_pair_observation.csv",
        pair_rows,list(pair_rows[0]))
    receipt={
        "schema":"chocho_zenodo2020_actual_voyage_original_WPI194_conservative_name_geography_match_v1",
        "status":"2020_ZENODO_SOURCE_WPI194_IDENTITIES_AUDITED_NO_BIOLOGICAL_MODELS",
        "source_publisher":"Duvallet et al. 2024, Zenodo 10.5281/zenodo.10783698, CC BY 4.0",
        "published_archive_md5":ZENODO_MD5,
        "source_original194_WPI_git_blob":PORTS_GIT_BLOB,
        "source_original_194_WPI_ports":194,
        "source_original_182_marine_eligible":182,
        "source_original_51_landmasses":51,
        "linkage_name_rule":"NFKD ASCII uppercase canonical exact names; removable PORT/HARBOR/HARBOUR/TERMINAL edge tokens",
        "linkage_distance_max_km":PRIMARY_MAX_DISTANCE_KM,
        "external_port_ID_is_NOT_assumed_to_equal_WPI_OBJECTID":True,
        "source_raw_triplet_rows":states["source_rows"],
        "source_valid_2020_departure_arrival_geo_rows":states["valid_2020_source_rows"],
        "source_invalid_row_reasons":dict(states),
        "raw_unique_origin_destination_source_node_pairs":len(segments),
        "raw_unique_external_endpoint_nodes":len(nodes),
        "match_status_nodes":dict(Counter(x["status"] for x in lookup.values())),
        "verified_external_endpoint_nodes":verified_node,
        "verified_original_WPI_ports_with_one_or_more_external_endpoints":len(original_port_nodes),
        "original_islands_with_one_or_more_verified_2020_endpoints":len(island_port_covered),
        "WPI_port_departure_arrival_observation_events_at_verified_nodes":
             {"departures":sum(original_port_departs.values()),"arrivals":sum(original_port_arrivals.values())},
        "source_matched_triplet_status_counts":dict(source_legs),
        "distinct_WPI_origin_destination_port_pairs":len(wpi_edges),
        "distinct_original_landmass_directed_pairs":len(island_edges),
        "distinct_original_landmass_unordered_pairs_with_observed_2020_voyage":len(unord),
        "all_unordered_original_island_pairs_preserved":len(pair_rows),
        "source_2020_departure_range":[dt_min.isoformat() if dt_min else None,
                                      dt_max.isoformat() if dt_max else None],
        "coverage_gate_15_islands_and_80_direct_island_pairs":
           len(island_port_covered)>=15 and len(unord)>=80,
        "raw_source_missing_voyage_is_not_actual_zero_ship_traffic":True,
        "source_2020_postdates_early_2001_2010_butterflies_and_plants":True,
        "no_biological_model_or_colonization_inference":True,
    }
    (output/"source_matched_voyage_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ZENODO2020_WPI_MATCH_JSON="+json.dumps({k:receipt[k] for k in [
       "status","source_raw_triplet_rows","source_valid_2020_departure_arrival_geo_rows",
       "raw_unique_external_endpoint_nodes","match_status_nodes","verified_original_WPI_ports_with_one_or_more_external_endpoints",
       "original_islands_with_one_or_more_verified_2020_endpoints",
       "source_matched_triplet_status_counts","distinct_WPI_origin_destination_port_pairs",
       "distinct_original_landmass_directed_pairs",
       "distinct_original_landmass_unordered_pairs_with_observed_2020_voyage",
       "coverage_gate_15_islands_and_80_direct_island_pairs"
    ]},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--zip-path",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    x=p.parse_args()
    run(x.zip_path,x.out_dir)
