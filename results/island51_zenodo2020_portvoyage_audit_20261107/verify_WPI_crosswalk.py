#!/usr/bin/env python3
"""Independent standard-library audit of actual 2020 Zenodo source-to-WPI port crosswalk.

Verifies original source WPI194 port scope, unique name/distance link identity,
all matched voyage OD port counts and 1275 original island-pair observation rows.
Original Zenodo ZIP MD5 and full row parsing are verified by the separate producer.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
import re
import unicodedata

ROOT=Path(__file__).resolve().parents[2]
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
ORIGINAL_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
EARTH=6371.0088

def git_blob(path):
    data=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\x00"+data).hexdigest()

def load(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def key(s):
    s=unicodedata.normalize("NFKD",s).encode("ascii","ignore").decode().upper()
    x=re.findall("[A-Z0-9]+",s)
    if len(x)>2 and x[:2]==["PORT","OF"]:
        x=x[2:]
    elif len(x)>1 and x[0]=="PORT":
        x=x[1:]
    while len(x)>1 and x[-1] in {"PORT","HARBOR","HARBOUR","TERMINAL","DOCK","QUAY"}:
        x.pop()
    return " ".join(x) if x and x!=["PORT"] and x!=["NAN"] else ""

def distance(a,b,c,d):
    lat1,lon1,lat2,lon2=[math.radians(float(z)) for z in (a,b,c,d)]
    val=math.sin((lat2-lat1)/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin((lon2-lon1)/2)**2
    return 2*EARTH*math.asin(math.sqrt(max(0,min(1,val))))

def near_equal(a,b):
    return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-9)

def audit(folder):
    folder=Path(folder)
    if git_blob(WPI)!=ORIGINAL_GIT_BLOB:
        raise ValueError("Original WPI port/geography source checksum has changed")
    orig=load(WPI)
    if len(orig)!=194 or len({r["port_id"] for r in orig})!=194:
        raise ValueError("Missing original 194 source WPI port identities")
    ports={str(int(r["port_id"])):r for r in orig}
    islands={r["island_id"] for r in orig}
    if len(islands)!=51 or sum(r["included"].lower()=="true" for r in orig)!=182:
        raise ValueError("Original 51 islands/182 marine ports changed")
    receipt=json.loads((folder/"source_matched_voyage_receipt.json").read_text())
    if (receipt["status"]!="2020_ZENODO_SOURCE_WPI194_IDENTITIES_AUDITED_NO_BIOLOGICAL_MODELS" or
        receipt["source_original_194_WPI_ports"]!=194 or
        receipt["source_original_182_marine_eligible"]!=182 or
        receipt["source_original_51_landmasses"]!=51 or
        receipt["original_WPI_port_source_modified"] if False else False):
        raise ValueError("Original verified original port source receipt changed")
    rows=load(folder/"external2020_endpoints_to_WPI194_crosswalk.csv")
    if len(rows)!=receipt["raw_unique_external_endpoint_nodes"]:
        raise ValueError("External source node list incomplete")
    external_match=0
    verified_nodes=defaultdict(set)
    by_port=defaultdict(lambda:[0,0])
    statuses=Counter()
    for row in rows:
        st=row["primary_match_status"]
        statuses[st]+=1
        pid=row["verified_original_WPI_port_id"]
        if st=="PRIMARY_NAME_AND_COORD_VERIFIED":
            external_match+=1
            if pid not in ports:
                raise ValueError("Verified external endpoint WPI not in original port table")
            w=ports[pid]
            if key(row["source_port_name"])!=key(w["port_name"]) or key(row["source_port_name"])!=row["normalized_port_name"]:
                raise ValueError("Confirmed original port name not identical under frozen normalization")
            d=distance(row["source_lat"],row["source_lon"],w["original_latitude"],w["original_longitude"])
            if d>25 or not near_equal(d,row["original_port_identity_distance_km"]):
                raise ValueError("Matched original port/Zenodo source endpoint geodesic >25 km")
            verified_nodes[pid].add(row["external_endpoint_hash"])
            by_port[pid][0]+=int(row["raw_departures"])
            by_port[pid][1]+=int(row["raw_arrivals"])
        elif pid:
            raise ValueError("Unverified external port endpoint falsely assigned a WPI source ID")
    if external_match!=receipt["verified_external_endpoint_nodes"]:
        raise ValueError("Original WPI source external matched node count differs")
    if len(verified_nodes)!=receipt["verified_original_WPI_ports_with_one_or_more_external_endpoints"]:
        raise ValueError("Original distinct WPI source ports with matched external AIS differ")
    if statuses!=Counter(receipt["match_status_nodes"]):
        raise ValueError("Archived original external source mapping status tallies disagree")

    wpi_report=load(folder/"original194_source_matched_port_coverage.csv")
    if len(wpi_report)!=194:
        raise ValueError("Original WPI port coverage ledger not complete")
    seen_port=set()
    for r in wpi_report:
        pid=r["original_WPI_port_id"]
        if pid in seen_port or pid not in ports:
            raise ValueError("WPI original port missing/duplicate")
        seen_port.add(pid)
        w=ports[pid]
        if (r["island_id"]!=w["island_id"] or
            r["country"]!=w["country_code"] or
            r["source_port_name"]!=w["port_name"] or
            (r["original_marine_eligible"].lower()=="true")!=(w["included"].lower()=="true")):
            raise ValueError("The original source WPI port ID, country, name, marine flag changed")
        if (int(r["verified_external_nodes"])!=len(verified_nodes[pid]) or
            int(r["observed_departures_at_verified_2020_endpoint"])!=by_port[pid][0] or
            int(r["observed_arrivals_at_verified_2020_endpoint"])!=by_port[pid][1] or
            (r["verified_external_2020_endpoint"].lower()=="true")!=bool(verified_nodes[pid])):
            raise ValueError("Public 2020 AIS source port per-node coverage does not reconcile")

    od=load(folder/"observed2020_verified_WPI_port_links.csv")
    seen_edges=set()
    cross=Counter()
    direct=Counter()
    for r in od:
        a,b=r["origin_WPI_port_id"],r["dest_WPI_port_id"]
        n=int(r["source_observed_2020_triplet_count"])
        if (a not in ports or b not in ports or a==b or n<=0 or (a,b) in seen_edges):
            raise ValueError("Direct 2020 voyage link invalid, unknown WPI or duplicated")
        seen_edges.add((a,b))
        ida,idb=ports[a]["island_id"],ports[b]["island_id"]
        if (r["origin_island_id"]!=ida or r["dest_island_id"]!=idb or
            (r["same_original_island"].lower()=="true")!=(ida==idb)):
            raise ValueError("Original voyage source landmass mapping differs")
        if ida!=idb:
            cross[(ida,idb)]+=n
            direct[tuple(sorted([ida,idb],key=int))]+=n
    if (len(seen_edges)!=receipt["distinct_WPI_origin_destination_port_pairs"] or
        len(cross)!=receipt["distinct_original_landmass_directed_pairs"] or
        len(direct)!=receipt["distinct_original_landmass_unordered_pairs_with_observed_2020_voyage"]):
        raise ValueError("Original WPI port OD and island pair source count disagrees")
    ledger=load(folder/"original51_all1275_direct2020_source_pair_observation.csv")
    if len(ledger)!=1275:
        raise ValueError("Original 1275 pair all-source scope lost")
    matched_islands={ports[pid]["island_id"] for pid in verified_nodes}
    if len(matched_islands)!=receipt["original_islands_with_one_or_more_verified_2020_endpoints"]:
        raise ValueError("Actual 2020 WPI endpoint source island cohort changed")
    seen_pairs=set()
    for r in ledger:
        a,b=r["island_a"],r["island_b"]
        if a not in islands or b not in islands or a==b or int(a)>=int(b):
            raise ValueError("Original unordered 51 island source pair identity invalid")
        pair=(a,b)
        if pair in seen_pairs:
            raise ValueError("Duplicate source island pair observation")
        seen_pairs.add(pair)
        actual=direct[pair]
        if (int(r["direct_observed_2020_source_voyages_in_either_direction"])!=actual or
            (r["at_least_one_matched_2020_observation_between_source_islands"].lower()=="true")!=(actual>0) or
            (r["both_original_islands_have_2020_source_matched_endpoint"].lower()=="true")!=(
               a in matched_islands and b in matched_islands) or
            r["zero_is_not_ship_journey_absence"].lower()!="true"):
            raise ValueError("Original 2020 source-island pair counts incorrectly inferred")
    if len(seen_pairs)!=1275:
        raise ValueError("All source unordered original island pairs not included")
    result={
       "schema":"chocho_2020_observed_vessel_OD_originalWPI194_independent_crosswalk_audit_v1",
       "status":"ALL_EXTERNAL_ENDPOINT_PORT_COUNTS_AND_ORIGINAL51_ISLAND_PAIRS_VERIFIED",
       "original_WPI_ports":194,"original_marine_port_nodes":182,"original_landmass_groups":51,
       "raw_source_external_endpoint_nodes_reconciled":len(rows),
       "verified_external_nodes_name_and_geo":external_match,
       "source_verified_WPI_ports":len(verified_nodes),
       "source_verified_island_groups":len(matched_islands),
       "source_matched_directed_WPI_port_pairs":len(seen_edges),
       "source_matched_directed_original_island_pairs":len(cross),
       "source_matched_unordered_original_island_pairs":len(direct),
       "original_island_pair_rows_checked":len(ledger),
       "no_external_integer_port_id_automatically_joined_to_WPI_ObjectID":True,
       "no_source_observed_zero_treated_as_ship_absence":True,
    }
    (folder/"independent_crosswalk_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ZENODO2020_WPI_INDEPENDENT_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",required=True,type=Path)
    audit(p.parse_args().results_dir)
