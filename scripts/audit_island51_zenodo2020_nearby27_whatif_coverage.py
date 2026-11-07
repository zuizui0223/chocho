#!/usr/bin/env python3
"""ONLY a hypothetical coverage envelope from 27 UNVERIFIED geographically close AIS port nodes.

Preserves and reproduces original 46 WPI-verified nodes and 12 original source island-pair
links. The augmented 27 port assignments are knowingly NOT verified: no voyages or
butterfly model may be claimed from the added hypothetical links.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from reconcile_island51_zenodo2020_WPI_port_voyages import read_external,original
from audit_island51_zenodo2020_voyage_source import ZENODO_MD5

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107"
OUTPUT_REQUIRED_OLD={
   "source_event_rows":645581,
   "external_nodes":1923,
   "verified_node_count":46,
   "matched_2020_source_voyages_both":442,
   "cross_original_island_events":133,
   "unique_original_island_pair_links":12,
   "directed_WPI_links":94,
}
GEO_ONLY_RADIUS_KM=10.0
N_UNVERIFIED_NEARBY=27

def read_csv(path):
    with Path(path).open("rt",encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def source_node(row,which):
    if which=="near":
        return (row["external_port_id_NOT_WPI"],row["external_port_name_UNVERIFIED"],
                round(float(row["source_lon"]),7),round(float(row["source_lat"]),7))
    if which=="crosswalk":
        return (row["external_port_id"],row["source_port_name"],
                round(float(row["source_lon"]),7),round(float(row["source_lat"]),7))
    raise ValueError("Unknown original port source node role")

def fixed_port_assignment():
    ports,source_name_dict=original()
    cross=read_csv(INPUT/"external2020_endpoints_to_WPI194_crosswalk.csv")
    near=read_csv(INPUT/"source_only_10km_unverified_port_candidates.csv")
    initial=read_csv(INPUT/"original51_all1275_direct2020_source_pair_observation.csv")
    original_receipt=json.loads((INPUT/"source_matched_voyage_receipt.json").read_text())
    if (len(ports)!=194 or len(cross)!=1923 or len(near)!=N_UNVERIFIED_NEARBY or
        len(initial)!=1275 or
        original_receipt["verified_external_endpoint_nodes"]!=46):
        raise ValueError("Original source-selected 194/1923/27/1275 cohorts changed")
    verified={};seen=set()
    external_unmatched={}
    for row in cross:
        node=source_node(row,"crosswalk")
        if node in seen:
            raise ValueError("Duplicate source external endpoint")
        seen.add(node)
        if row["primary_match_status"]=="PRIMARY_NAME_AND_COORD_VERIFIED":
            pid=str(int(row["verified_original_WPI_port_id"]))
            if pid not in ports:
                raise ValueError("A source verified endpoint was outside original WPI194")
            verified[node]=pid
        else:
            if row["verified_original_WPI_port_id"]:
                raise ValueError("Unverified original source record has WPI assignment")
            external_unmatched[node]=row
    if len(verified)!=46 or len(seen)!=1923:
        raise ValueError("Original 46 verified WPI nodes or 1923 complete external records missing")
    provisional={}
    for row in near:
        node=source_node(row,"near")
        if node in provisional or node not in external_unmatched or node in verified:
            raise ValueError("The 27 geography-only source nodes were absent or duplicated")
        linked=external_unmatched[node]
        if (row["identity_verified"].lower()!="false" or
            linked["primary_match_status"]!="NO_EXACT_NAME_AND_COORD_MATCH" or
            linked["geo_only_10km_NOT_VERIFIED"].lower()!="true" or
            float(row["geodesic_km"])>GEO_ONLY_RADIUS_KM):
            raise ValueError("A 27 candidate was already verified or not strictly within 10 km")
        wpi=str(int(row["nearest_WPI_id_NOT_VERIFIED"]))
        if (wpi not in ports or
            str(int(linked["geography_only_nearest_WPI_id_NOT_VERIFIED"]))!=wpi or
            linked["source_port_name"]!=row["external_port_name_UNVERIFIED"]):
            raise ValueError("Source-only geography candidate does not match original WPI nearest ledger")
        provisional[node]=wpi
    if len(provisional)!=27 or set(provisional)&set(verified):
        raise ValueError("Original source-only hypothetical cohort changed")
    return ports,verified,provisional,original_receipt,initial

def aggregate(segments,ports,assignment,provisional):
    original_islands={p["island_id"] for p in ports.values()}
    event_status=Counter()
    wpi_edges=Counter()
    cross_island=Counter()
    directed_island=Counter()
    tentative_event_counts=Counter()
    # Publisher source endpoint is a pre-existing inferred voyage, NOT a confirmed
    # nonstop ship journey; source pair abundance can be incomplete/biased.
    for (a,b),n in segments.items():
        wa=assignment.get(a);wb=assignment.get(b)
        if wa is None or wb is None:
            event_status["unmapped_at_least_one_endpoint"]+=n
            continue
        event_status["both_endpoints_mapped"]+=n
        involved_unverified=a in provisional or b in provisional
        if involved_unverified:
            tentative_event_counts["source_events_involving_at_least_one_UNVERIFIED_node"]+=n
        if wa==wb:
            event_status["same_original_WPI_port"]+=n
            continue
        wpi_edges[(wa,wb)]+=n
        ia,ib=ports[wa]["island_id"],ports[wb]["island_id"]
        if ia==ib:
            event_status["different_WPI_ports_same_island"]+=n
            continue
        event_status["cross_original_island_events"]+=n
        directed_island[(ia,ib)]+=n
        pair=tuple(sorted((ia,ib),key=int))
        cross_island[pair]+=n
        if involved_unverified:
            tentative_event_counts["source_cross_island_events_involving_UNVERIFIED_node"]+=n
    if event_status["both_endpoints_mapped"]!=(
            event_status["same_original_WPI_port"]+
            event_status["different_WPI_ports_same_island"]+
            event_status["cross_original_island_events"]):
        raise ValueError("Source triplet aggregation outcomes not exhaustive")
    return event_status,wpi_edges,cross_island,directed_island,tentative_event_counts

def output_csv(path,rows):
    if not rows:
        raise ValueError("Source-only output has no row support")
    with Path(path).open("wt",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def run(zip_path,out_dir):
    out=Path(out_dir).resolve()
    if any(out.is_relative_to(ROOT/p) for p in ("data","results","provenance")):
        raise ValueError("Never overwrite source audited direct-port data")
    raw=Path(zip_path).read_bytes()
    if hashlib.md5(raw).hexdigest()!=ZENODO_MD5:
        raise ValueError("Source original 2020 AIS ZIP MD5 differs from publisher")
    original_ports,real,candidates,prior,source_pair_ledger=fixed_port_assignment()
    states,origin,destination,segments,first,last=read_external(zip_path)
    if states["source_rows"]!=OUTPUT_REQUIRED_OLD["source_event_rows"]:
        raise ValueError("Original publisher source 645581 rows changed")
    source_nodes=set(origin)|set(destination)
    if not set(real).issubset(source_nodes) or not set(candidates).issubset(source_nodes):
        raise ValueError("One of the 46 verified or 27 hypothetical source nodes is absent from the original 2020 AIS records")
    if len(set(origin)|set(destination))!=1923:
        raise ValueError("Original source endpoints differ from 1923 recorded nodes")
    base_status,base_ports,base_islands,base_direction,_=aggregate(segments,original_ports,real,{})
    if (base_status["both_endpoints_mapped"]!=442 or
        base_status["cross_original_island_events"]!=133 or
        len(base_islands)!=12 or len(base_ports)!=94 or len(base_direction)!=20):
        raise ValueError("Original verified WPI source voyages not exactly reproduced")
    optimistic={**real,**candidates}
    hypo_status,hypo_ports,hypo_islands,hypo_direction,provisional=aggregate(
        segments,original_ports,optimistic,candidates)
    if (hypo_status["both_endpoints_mapped"]<442 or
        hypo_status["cross_original_island_events"]<133 or
        len(hypo_islands)<12):
        raise ValueError("Adding hypothetical unverified source node labels reduced original observed support")
    unique_original_ports={str(int(x)) for x in optimistic.values()}
    counted_groups={original_ports[pid]["island_id"] for pid in unique_original_ports}
    expected_keys={(row["island_a"],row["island_b"]) for row in source_pair_ledger}
    if len(expected_keys)!=1275:
        raise ValueError("Original source 1275 island pairs not retained")
    pair_rows=[]
    new_pairs=0
    for row in source_pair_ledger:
        pair=(row["island_a"],row["island_b"])
        a=base_islands[pair]
        b=hypo_islands[pair]
        if int(row["direct_observed_2020_source_voyages_in_either_direction"])!=a or b<a:
            raise ValueError("Original validated port pair observation count changed")
        if a==0 and b>0:new_pairs+=1
        pair_rows.append({
           "island_a":pair[0],"island_b":pair[1],
           "previous_VERIFIED_source_events":a,
           "hypothetical_27_geonear_alias_source_events_NOT_VERIFIED":b,
           "additional_source_events_due_to_unverified_port_identity":b-a,
           "previous_verified_direct_pair":a>0,
           "pair_exists_only_if_unverified_port_candidates_are_ASSUMED_same":a==0 and b>0,
           "zero_source_index_is_not_no_ship_traffic":True,
           "historical_host_or_butterfly_colonization_not_inferred":True,
        })
    if len(pair_rows)!=1275 or new_pairs!=len(set(hypo_islands)-set(base_islands)):
        raise ValueError("Original island pair support or new unverified link types invalid")
    out.mkdir(parents=True,exist_ok=True)
    output_csv(out/"original51_verified_vs_UNVERIFIED_geonear_alias_2020_OD_pair_scenarios.csv",pair_rows)
    receipt={
        "schema":"chocho_2020_zenodo_originalWPI_27_unverified_name_mismatch_source_coverage_envelope_v1",
        "status":"SOURCE_ONLY_HYPOTHETICAL_27_NEARBY_NODE_COVERAGE_CALCULATED_NOT_VERIFIED_TRAFFIC",
        "source_original_Zenodo_DOI":"10.5281/zenodo.10783698",
        "source_original_Zenodo_MD5":ZENODO_MD5,
        "original_source_194_WPI_ports":194,"original_source_island_groups":51,
        "complete_original_2020_source_triplet_rows":states["source_rows"],
        "external_original_2020_source_nodes":1923,
        "prior_verified_source_endpoint_nodes":len(real),
        "ADDITIONAL_UNVERIFIED_geography_only_endpoint_nodes":len(candidates),
        "hypothetically_mapped_source_nodes":len(optimistic),
        "hypothetically_mapped_original_WPI_ports":len(unique_original_ports),
        "hypothetically_mapped_original_island_groups":len(counted_groups),
        "previous_verified_both_endpoints_triplets":base_status["both_endpoints_mapped"],
        "hypothetical_both_endpoints_triplets":hypo_status["both_endpoints_mapped"],
        "previous_verified_cross_island_source_events":base_status["cross_original_island_events"],
        "hypothetical_cross_island_source_events":hypo_status["cross_original_island_events"],
        "previous_verified_unordered_island_pairs":len(base_islands),
        "hypothetical_unordered_island_pairs":len(hypo_islands),
        "additional_island_pair_types_requiring_unverified_node_assumptions":new_pairs,
        "previous_verified_WPI_directed_port_links":len(base_ports),
        "hypothetical_WPI_directed_port_links":len(hypo_ports),
        "prior_source_15islands_80pairs_gate_passed":False,
        "hypothetical_15islands_80pairs_threshold_reached_if_all_aliases_true":(
            len(counted_groups)>=15 and len(hypo_islands)>=80
        ),
        "ACTUAL_SOURCE_GATE_REMAINS_FAILED_WITHOUT_VERIFIED_IDENTITIES":True,
        "hypothetical_event_status":dict(hypo_status),
        "hypothetical_unverified_node_event_influence":dict(provisional),
        "source_departure_range_2020":[first.isoformat(),last.isoformat()],
        "unverified_aliases_are_not_promoted_to_original_WPI_matches":True,
        "all_1275_original_source_island_pairs_retained":True,
        "zero_source_event_is_not_negative_ship_route_evidence":True,
        "no_butterfly_plant_or_causal_transport_model_fitted":True,
    }
    (out/"source_only_geonear27_whatif_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("ZENODO27_GEOONLY_ENVELOPE_RESULT_JSON="+json.dumps({
        k:receipt[k] for k in (
            "status","previous_verified_unordered_island_pairs",
            "hypothetical_unordered_island_pairs",
            "additional_island_pair_types_requiring_unverified_node_assumptions",
            "previous_verified_cross_island_source_events",
            "hypothetical_cross_island_source_events",
            "hypothetically_mapped_original_island_groups",
            "hypothetical_15islands_80pairs_threshold_reached_if_all_aliases_true",
            "ACTUAL_SOURCE_GATE_REMAINS_FAILED_WITHOUT_VERIFIED_IDENTITIES"
        )},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--zip-path",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    run(a.zip_path,a.out_dir)
