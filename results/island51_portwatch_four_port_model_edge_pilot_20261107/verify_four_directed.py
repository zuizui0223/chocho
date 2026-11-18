#!/usr/bin/env python3
"""Independent original-source four WPI ports / twelve provider modeled directed links auditor.

Read-only, no network and no fitting-module imports. Checks original WPI194 source
and direct Zenodo 2020 episode source edge table, prior PortWatch registry pilot,
all 12 model pair conditions, sample port IDs/geo, and 8 model degree counts.
A 0 provider model edge does not mean no actual vessel voyages.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
ORIGINAL=ROOT/"results/island51_host_circuit_atlas_20261106/original_port_geography.csv"
PILOT=ROOT/"results/island51_portwatch_global_connections_20261107/WPI5_bounded_registry_source_pilot_receipt.json"
REAL_ZENODO=ROOT/"results/island51_zenodo2020_portvoyage_audit_20261107/observed2020_verified_WPI_port_links.csv"
SRC_GIT_BLOB="a4681c3e3c1406d3a37f51464a298286e9d0ee4f"
SOURCE_WPI={
 "28":("port141","26","BELFAST","GB"),
 "428":("port1153","206","SANTA CRUZ DE TENERIFE","ES"),
 "817":("port647","85","LIMASSOL","CY"),
 "1416":("port254","31","COLOMBO","LK")
}
EARTH_R_KM=6371.0088

def read(path):
    with Path(path).open(encoding="utf-8",newline="") as stream:
        return list(csv.DictReader(stream))

def km(lat1,lon1,lat2,lon2):
    p,q,r,s=(math.radians(float(x)) for x in (lat1,lon1,lat2,lon2))
    a=math.sin((r-p)/2)**2+math.cos(p)*math.cos(r)*math.sin((s-q)/2)**2
    return 2*EARTH_R_KM*math.asin(math.sqrt(max(0.,min(1.,a))))

def source():
    raw=ORIGINAL.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if blob!=SRC_GIT_BLOB:
        raise ValueError("WPI original port cohort changed")
    rows=read(ORIGINAL)
    if (len(rows)!=194 or len({x["port_id"] for x in rows})!=194 or
        len({x["island_id"] for x in rows})!=51 or
        sum(x["included"].lower()=="true" for x in rows)!=182):
        raise ValueError("Original 51-island/194-WPI/182-marine source incomplete")
    wpi={str(int(row["port_id"])):row for row in rows}
    pilot=json.loads(PILOT.read_text(encoding="utf-8"))
    if pilot["original_WPI_sources_with_verified_provider_portid"]!=4:
        raise ValueError("Source verified four original WPI PortWatch identities no longer valid")
    matched={}
    for item in pilot["source_check_rows"]:
        sid=str(item["original_WPI_port_id"])
        if sid not in SOURCE_WPI:
            continue
        expected=SOURCE_WPI[sid]
        original=wpi[sid]
        if (item["verified_provider_portid"]!=expected[0] or
            original["island_id"]!=expected[1] or
            original["port_name"]!=expected[2] or
            original["country_code"]!=expected[3] or
            item["source_port_identity_verification_status"]!="ONE_UNIQUE_EXACT_NAME_ISO3_GEOGRAPHY_PROVIDER_MATCH" or
            item["provider_registry_sample_complete"] is not True):
            raise ValueError("One originally verified four WPI provider port ID incompatible")
        provider=[x for x in item["all_provider_registry_sample_rows"]
                  if str(x.get("portid"))==expected[0]]
        if len(provider)!=1:
            raise ValueError("Verified PortWatch provider registry sample no longer unique")
        r=provider[0]
        if km(original["original_latitude"],original["original_longitude"],r["lat"],r["lon"])>25:
            raise ValueError("Pilot PortWatch provider registry is too far from original WPI source")
        matched[sid]={"wpi":original,"portwatch":expected[0],"provider_registry":r}
    if set(matched)!=set(SOURCE_WPI):
        raise ValueError("Frozen 4 source provider identities incomplete")
    old=read(REAL_ZENODO)
    if len(old)!=94 or len({(x["origin_WPI_port_id"],x["dest_WPI_port_id"]) for x in old})!=94:
        raise ValueError("Original source matched2020 Zenodo directed WPI link table changed")
    old_links={(a["origin_WPI_port_id"],a["dest_WPI_port_id"]):int(a["source_observed_2020_triplet_count"])
               for a in old}
    return matched,old_links

def verify(out_dir):
    original,old=source()
    output=Path(out_dir)
    path=output/"source_four_port_directed_model_edge_pilot.json"
    current=json.loads(path.read_text(encoding="utf-8"))
    if (current["schema"]!="chocho_original_WPI4_portwatch_12_directed_static_model_edge_pilot_v1" or
        current["original_WPI_git_blob"]!=SRC_GIT_BLOB or
        current["source_enriched_WPI_port_count"]!=4 or
        current["original_WPI_port_count"]!=194 or
        current["original_island_count"]!=51 or
        current["original_182_marine_qualified_port_count"]!=182 or
        current["n_pair_count_queries"]!=12 or
        current["n_provider_degree_count_queries"]!=8 or
        current["n_sample_feature_rows_requested"]>12):
        raise ValueError("Bounded four-source-provider PortWatch result shape/source lock invalid")
    rows=current["paired_source_record_comparisons"]
    deg=current["source_four_port_in_out_degree_count_queries"]
    ids=list(SOURCE_WPI)
    if len(rows)!=12 or len(deg)!=8:
        raise ValueError("Expected all 12 ordered source comparisons and 8 network degrees")
    checked=set()
    valid=0;positive=0;previous_positive=0;intersect=0;examples=0
    for row in rows:
        a,b=str(row["original_WPI_origin"]),str(row["original_WPI_destination"])
        if a not in original or b not in original or a==b or (a,b) in checked:
            raise ValueError("Unknown, duplicated or self-linked originally selected WPI four-port pair")
        checked.add((a,b))
        o,d=original[a],original[b]
        if (row["original_from_island"]!=o["wpi"]["island_id"] or
            row["original_to_island"]!=d["wpi"]["island_id"] or
            row["portwatch_from_id"]!=o["portwatch"] or
            row["portwatch_to_id"]!=d["portwatch"] or
            int(row["verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED"])!=old.get((a,b),0)):
            raise ValueError("Source original WPI and verified 2020 Zenodo port link IDs differ")
        previous_positive+=old.get((a,b),0)>0
        completed=bool(row["PortWatch_count_query_complete"])
        q=row["provider_count_request_receipt"]
        if not isinstance(q,dict) or not isinstance(q.get("attempts"),list) or not q["attempts"]:
            raise ValueError("Source provider count query lacks API response ledger")
        if completed:
            n=row["portwatch_static_model_directed_record_count"]
            if type(n)!=int or n<0 or q["success"] is not True:
                raise ValueError("PortWatch static model edge count lacks verified provider response")
            valid+=1
            if n>0:
                positive+=1
                intersect+=old.get((a,b),0)>0
                sample=row["provider_at_most_one_model_sample"]
                if len(sample)!=1 or not isinstance(row["provider_sample_request_receipt"],dict) or not row["provider_sample_request_receipt"]["success"]:
                    raise ValueError("Positive provider model edge has no validated bounded example")
                x=sample[0]
                if x["from_id"]!=o["portwatch"] or x["to_id"]!=d["portwatch"]:
                    raise ValueError("Wrong provider port ID in bounded static model-edge example")
                if not (x["from_source_WPI_geodesic_km"] is not None and
                        x["to_source_WPI_geodesic_km"] is not None):
                    raise ValueError("Missing publisher model source endpoint geodesic")
                if float(x["from_source_WPI_geodesic_km"])>25 or float(x["to_source_WPI_geodesic_km"])>25:
                    raise ValueError("Provider edge source endpoint mismatches original WPI original port")
                if x["distance"] is not None and float(x["distance"])<0:
                    raise ValueError("Source provider model distance negative")
                if x["transit_days"] is not None and float(x["transit_days"])<0:
                    raise ValueError("Source provider model transit time negative")
                examples+=1
            elif row["provider_at_most_one_model_sample"]:
                raise ValueError("An empty provider static model-edge condition had nonempty example rows")
        elif row["portwatch_static_model_directed_record_count"] is not None:
            raise ValueError("Provider query failure was silently interpreted as source zero")
        if row["portwatch_simulator_static_edge_not_actual_ship_journey"] is not True or row["2020_missing_source_triplet_not_actual_no_ship"] is not True:
            raise ValueError("A source model edge was confused with a real vessel voyage/absence")
    if len(checked)!=12 or checked!={(a,b) for a in ids for b in ids if a!=b}:
        raise ValueError("Some ordered four original source pair conditions not evaluated")
    degseen=set()
    for row in deg:
        a=str(row["wpi_id"]);kind=row["model_direction_field"]
        if a not in original or kind not in ("from_id","to_id") or (a,kind) in degseen:
            raise ValueError("Unknown/duplicate count-only source model in/out-degree request")
        degseen.add((a,kind))
        if row["portwatch_id"]!=original[a]["portwatch"] or row["island_id"]!=original[a]["wpi"]["island_id"]:
            raise ValueError("Provider model count degree uses incorrect original WPI registry ID")
        if row["count_query_complete"] and (type(row["provider_static_model_record_count"])!=int or
                                            row["provider_static_model_record_count"]<0):
            raise ValueError("Provider degree count status or value inconsistent")
        if not row["count_query_complete"] and row["provider_static_model_record_count"] is not None:
            raise ValueError("Provider unavailable count silently treated as zero")
    if len(degseen)!=8:
        raise ValueError("Source provider degree count requests missing")
    report={
       "schema":"chocho_fourWPI_PortWatch_model_links_vs_original2020_Zenodo_source_audit_v1",
       "status":"ALL_12_FIXED_ORIGINAL_WPI_PORTWATCH_MODEL_CONDITIONS_AND_8_DEGREES_AUDITED",
       "original_WPI_source_ports_in_study":194,
       "original_source_verified_portid_sample":4,
       "source_directed_fixed_pairs":12,
       "source_count_queries_complete":valid,
       "PortWatch_static_model_edges_present_in_verified_four_source_pairs":positive,
       "old_Zenodo2020_source_positive_edges_in_same_four_WPI_ports":previous_positive,
       "cross_source_model_edge_and_Zenodo_source_both_positive":intersect,
       "provider_model_feature_rows_sampled":examples,
       "all_four_parked_source_port_IDs_preserved":True,
       "source_count_only_plus_at_most12_rows_not_bulk":True,
       "model_edge_not_direct_observed_ship_voyage":True,
       "no_model_fitted_to_island_butterfly_or_historical_plant_responses":True
    }
    (output/"independent_four_port_model_edge_audit.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("PORTWATCH_WPI4_DIRECTED_MODEL_INDEPENDENT_AUDIT_JSON="+json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",required=True,type=Path)
    verify(p.parse_args().results_dir)
