#!/usr/bin/env python3
"""Independently reproduce verified 2020 WPI OD source and 27 UNVERIFIED what-if voyage coverage.

Reads publisher MD5 ZIP directly, with no import of the hypothetical runner.
A geography-only node candidate is NEVER treated as a confirmed source port.
"""
from __future__ import annotations
import argparse
import csv
import datetime
import hashlib
import json
import math
import zipfile
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
PARENT=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107"
WPI=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
PUB_MD5="af27e1bef1f706e2e29ed53aee47ff10"
SOURCE_MEMBER="triplet_2020.csv"

def read_csv(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def node(identity,name,lon,lat):
    return (str(identity).strip(),str(name).strip(),round(float(lon),7),round(float(lat),7))

def add_edge(pair_a,pair_b,ports,pair_counter):
    if not pair_a or not pair_b:
        return False
    if pair_a==pair_b:
        return True
    i,j=ports[pair_a],ports[pair_b]
    if i!=j:
        pair_counter[tuple(sorted([i,j],key=int))]+=1
    return True

def audit(raw_zip,out):
    path=Path(raw_zip);output=Path(out)
    contents=path.read_bytes()
    if hashlib.md5(contents).hexdigest()!=PUB_MD5:
        raise ValueError("Publisher Zenodo original 2020 AIS source MD5 incorrect")
    wpi=read_csv(WPI)
    if len(wpi)!=194 or len({r["port_id"] for r in wpi})!=194:
        raise ValueError("Original WPI source194 port IDs changed")
    ports={str(int(r["port_id"])):r["island_id"] for r in wpi}
    if len(ports)!=194 or len(set(ports.values()))!=51:
        raise ValueError("Original 51 island groups not retained")
    cross=read_csv(PARENT/"external2020_endpoints_to_WPI194_crosswalk.csv")
    near=read_csv(PARENT/"source_only_10km_unverified_port_candidates.csv")
    if len(cross)!=1923 or len(near)!=27:
        raise ValueError("Original source complete 1923 external nodes, 27 candidates not retained")
    real={};unmatched={}
    for r in cross:
        k=node(r["external_port_id"],r["source_port_name"],r["source_lon"],r["source_lat"])
        if k in real or k in unmatched:
            raise ValueError("Duplicated original external source endpoint")
        if r["primary_match_status"]=="PRIMARY_NAME_AND_COORD_VERIFIED":
            wpid=str(int(r["verified_original_WPI_port_id"]))
            if wpid not in ports:
                raise ValueError("Invalid original verified WPI port ID")
            real[k]=wpid
        else:
            if r["verified_original_WPI_port_id"]:
                raise ValueError("Unverified original source port was given a real WPI ID")
            unmatched[k]=r
    if len(real)!=46 or len(unmatched)!=1877:
        raise ValueError("Original source name and geography verification scope changed")
    hypothetical={}
    for r in near:
        k=node(r["external_port_id_NOT_WPI"],r["external_port_name_UNVERIFIED"],
               r["source_lon"],r["source_lat"])
        if k not in unmatched or k in hypothetical or r["identity_verified"].lower()!="false":
            raise ValueError("An unverified geographical node is outside the original source candidate cohort")
        old=unmatched[k]
        if old["primary_match_status"]!="NO_EXACT_NAME_AND_COORD_MATCH" or old["geo_only_10km_NOT_VERIFIED"].lower()!="true":
            raise ValueError("Geography candidate has source identity/10km mismatch")
        wpid=str(int(r["nearest_WPI_id_NOT_VERIFIED"]))
        if wpid!=str(int(old["geography_only_nearest_WPI_id_NOT_VERIFIED"])) or wpid not in ports:
            raise ValueError("A geographically nearest port was silently changed")
        hypothetical[k]=wpid
    if len(hypothetical)!=27 or set(real)&set(hypothetical):
        raise ValueError("Verified and UNVERIFIED source node sets overlap")
    assumed={**real,**hypothetical}
    r_events=Counter()
    h_events=Counter()
    r_pairs=Counter()
    h_pairs=Counter()
    r_port_edges=set()
    h_port_edges=set()
    trips=0
    with zipfile.ZipFile(path,"r") as z:
        if [n for n in z.namelist() if not n.endswith("/") ]!=[SOURCE_MEMBER]:
            raise ValueError("Original publisher ZIP member schema changed")
        with z.open(SOURCE_MEMBER) as raw:
            from io import TextIOWrapper
            reader=csv.DictReader(TextIOWrapper(raw,encoding="utf-8-sig",newline=""),delimiter=";")
            required={"date-out","date-in","Port_ID_waiting_area_id_from","PortName_from",
                "lon_from","lat_from","Port_ID_waiting_area_id_to","PortName_to","lon_to","lat_to"}
            if not required.issubset(set(reader.fieldnames or [])):
                raise ValueError("Original source port origin/destination schema changed")
            for row in reader:
                trips+=1
                dep=datetime.datetime.fromisoformat(row["date-out"].strip())
                arr=datetime.datetime.fromisoformat(row["date-in"].strip())
                if dep.year!=2020 or arr<dep:
                    raise ValueError("Original parsed 2020 source chronology inconsistent")
                a=node(row["Port_ID_waiting_area_id_from"],row["PortName_from"],row["lon_from"],row["lat_from"])
                b=node(row["Port_ID_waiting_area_id_to"],row["PortName_to"],row["lon_to"],row["lat_to"])
                for label,mapping,event,pairs,edges in (
                    ("verified",real,r_events,r_pairs,r_port_edges),
                    ("hypothetical",assumed,h_events,h_pairs,h_port_edges),
                ):
                    u=mapping.get(a);v=mapping.get(b)
                    if u is None or v is None:
                        event["one_or_both_unmapped"]+=1
                        continue
                    event["both_source_endpoints_mapped"]+=1
                    if u==v:
                        event["same_original_WPI_port"]+=1
                        continue
                    edges.add((u,v))
                    if ports[u]==ports[v]:
                        event["different_original_ports_same_island"]+=1
                    else:
                        event["different_original_islands"]+=1
                        pairs[tuple(sorted([ports[u],ports[v]],key=int))]+=1
                        if label=="hypothetical" and (a in hypothetical or b in hypothetical):
                            event["cross_island_triplet_requires_UNVERIFIED_nearby_node"]+=1
    if trips!=645581 or r_events["both_source_endpoints_mapped"]!=442 or r_events["different_original_islands"]!=133:
        raise ValueError("Original source actual-verified triplet event counts differ")
    if len(r_pairs)!=12 or len(r_port_edges)!=94:
        raise ValueError("Original 12 verified island pair types / 94 WPI links differ")
    prior=read_csv(PARENT/"original51_all1275_direct2020_source_pair_observation.csv")
    if len(prior)!=1275:
        raise ValueError("Original source all 1275 island pairs changed")
    saved=read_csv(output/"original51_verified_vs_UNVERIFIED_geonear_alias_2020_OD_pair_scenarios.csv")
    if len(saved)!=1275:
        raise ValueError("Hypothetical pair scenario did not retain complete source population")
    index={(x["island_a"],x["island_b"]):x for x in saved}
    if len(index)!=1275:
        raise ValueError("Hypothetical original island pairs duplicated")
    for p in prior:
        key=(p["island_a"],p["island_b"])
        if key not in index:
            raise ValueError("A source original island pair was removed")
        r=index[key]
        real_count=r_pairs[key]
        hypothetical_count=h_pairs[key]
        for label,actual in (
            ("previous_VERIFIED_source_events",real_count),
            ("hypothetical_27_geonear_alias_source_events_NOT_VERIFIED",hypothetical_count),
            ("additional_source_events_due_to_unverified_port_identity",hypothetical_count-real_count),
        ):
            if int(r[label])!=actual:
                raise ValueError("Source-only counterfactual pair event mismatch: "+label)
        if int(p["direct_observed_2020_source_voyages_in_either_direction"])!=real_count:
            raise ValueError("Original verified source paired voyages changed")
        if (r["previous_verified_direct_pair"].lower()=="true")!=(real_count>0) or (
            r["pair_exists_only_if_unverified_port_candidates_are_ASSUMED_same"].lower()=="true")!=(
                real_count==0 and hypothetical_count>0):
            raise ValueError("Hypothetical/verified source pair identity flag wrong")
        if r["zero_source_index_is_not_no_ship_traffic"].lower()!="true":
            raise ValueError("Missing source is incorrectly declared biologically no shipping")
    receipt=json.loads((output/"source_only_geonear27_whatif_receipt.json").read_text())
    expected={
        "complete_original_2020_source_triplet_rows":trips,
        "prior_verified_source_endpoint_nodes":len(real),
        "ADDITIONAL_UNVERIFIED_geography_only_endpoint_nodes":len(hypothetical),
        "hypothetically_mapped_source_nodes":len(assumed),
        "hypothetically_mapped_original_WPI_ports":len(set(assumed.values())),
        "hypothetically_mapped_original_island_groups":len({ports[p] for p in set(assumed.values())}),
        "previous_verified_both_endpoints_triplets":r_events["both_source_endpoints_mapped"],
        "hypothetical_both_endpoints_triplets":h_events["both_source_endpoints_mapped"],
        "previous_verified_cross_island_source_events":r_events["different_original_islands"],
        "hypothetical_cross_island_source_events":h_events["different_original_islands"],
        "previous_verified_unordered_island_pairs":len(r_pairs),
        "hypothetical_unordered_island_pairs":len(h_pairs),
        "previous_verified_WPI_directed_port_links":len(r_port_edges),
        "hypothetical_WPI_directed_port_links":len(h_port_edges),
        "additional_island_pair_types_requiring_unverified_node_assumptions":len(set(h_pairs)-set(r_pairs)),
    }
    for name,value in expected.items():
        if receipt[name]!=value:
            raise ValueError("Unverified 27-nearby port scenario receipt differs: "+name)
    if (not receipt["ACTUAL_SOURCE_GATE_REMAINS_FAILED_WITHOUT_VERIFIED_IDENTITIES"] or
        receipt["hypothetical_15islands_80pairs_threshold_reached_if_all_aliases_true"]!=(
           len({ports[p] for p in assumed.values()})>=15 and len(h_pairs)>=80)):
        raise ValueError("Source-only 80-pair coverage gate misinterpreted")
    result={
        "schema":"chocho_zenodo2020_hypothetical27_nearby_wpi_independent_raw_triplets_audit_v1",
        "status":"ORIGINAL_645581_SOURCE_ROWS_AND_1275_VERIFIED_VS_UNVERIFIED_PAIRS_RECONSTRUCTED",
        "publisher_source_rows_verified":trips,
        "verified_external_nodes":len(real),
        "geographic_candidate_nodes_still_UNVERIFIED":len(hypothetical),
        "original_verified_source_2020_cross_island_events":r_events["different_original_islands"],
        "hypothetical_source_cross_island_events_NOT_VERIFIED":h_events["different_original_islands"],
        "source_held_island_links_verified":len(r_pairs),
        "source_hypothetical_unverified_island_links":len(h_pairs),
        "original_source_island_pairs_individually_rechecked":len(saved),
        "source_27_geographic_candidates_not_promoted_to_actual_voyages":True,
        "no_butterfly_or_causal_ship_transport_claim":True,
    }
    (output/"independent_raw_source_audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ZENODO27_SOURCE_INDEPENDENT_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--zip-path",type=Path,required=True)
    parser.add_argument("--results-dir",type=Path,required=True)
    args=parser.parse_args()
    audit(args.zip_path,args.results_dir)
