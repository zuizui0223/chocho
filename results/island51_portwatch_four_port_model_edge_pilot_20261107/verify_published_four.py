#!/usr/bin/env python3
"""Check 12 published PortWatch model-vs-2020 source comparisons and 4 degrees offline."""
import csv
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
NAMES={"28":"Belfast","428":"Santa Cruz de Tenerife","817":"Limassol","1416":"Colombo"}

def read(name):
    with (ROOT/name).open(encoding="utf-8",newline="") as stream:
        return list(csv.DictReader(stream))

def main():
    source=json.loads((ROOT/"source_four_port_directed_model_edge_pilot.json").read_text())
    rows=read("source_four_port_directed_model_edges_vs_2020_vendor_triplets.csv")
    degrees=read("source_four_WPI_model_in_out_edge_record_counts.csv")
    pair={(x["original_WPI_origin"],x["original_WPI_destination"]):x for x in source["paired_source_record_comparisons"]}
    deg={(x["wpi_id"],x["model_direction_field"]):x for x in source["source_four_port_in_out_degree_count_queries"]}
    idmap={x["wpi_id"]:x for x in source["verified_4_WPI_source_registry_ids"]}
    if source["status"]!="ALL_12_FIXED_DIRECTED_MODEL_PAIR_QUERIES_AND_8_DEGREES_COMPLETED" or len(rows)!=len(pair)!=12 or len(degrees)!=4 or len(deg)!=8:
        raise ValueError("Incomplete previously archived PortWatch pilot source")
    seen=set()
    for row in rows:
        a,b=row["from_wpi"],row["to_wpi"]
        if a not in NAMES or b not in NAMES or a==b or (a,b) in seen:
            raise ValueError("Unknown original ordered WPI IDs")
        seen.add((a,b))
        x=pair[(a,b)]
        if (row["from_port"]!=NAMES[a] or row["to_port"]!=NAMES[b] or
            row["from_island"]!=x["original_from_island"] or
            row["to_island"]!=x["original_to_island"] or
            int(row["PortWatch_model_directed_records"])!=x["portwatch_static_model_directed_record_count"] or
            int(row["Zenodo_2020_matched_source_triplet_rows"])!=x["verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED"] or
            row["source_pilot_query_success"].lower()!="true" or
            row["PortWatch_is_static_simulator_network_NOT_actual_AIS_voyage"].lower()!="true" or
            row["zero_Zenodo_source_matched_records_NOT_no_real_ships"].lower()!="true"):
            raise ValueError("Original pilot source-linked published table changed")
    if seen!={(a,b) for a in NAMES for b in NAMES if a!=b}:
        raise ValueError("12 original ordered source pairs do not match exactly")
    if len(idmap)!=4:
        raise ValueError("Original 4 provider identities missing")
    for row in degrees:
        i=row["WPI_id"]
        if i not in NAMES or row["source_port_name"]!=NAMES[i] or row["provider_portid"]!=idmap[i]["provider_id"] or row["original_island_id"]!=idmap[i]["island_id"]:
            raise ValueError("Original WPI provider port/in-out table mismatch")
        if (int(row["static_model_outgoing_record_count"])!=deg[(i,"from_id")]["provider_static_model_record_count"] or
            int(row["static_model_incoming_record_count"])!=deg[(i,"to_id")]["provider_static_model_record_count"]):
            raise ValueError("PortWatch outgoing/incoming modeled record count altered")
    result={
      "status":"ALL_12_WPI_DIRECTED_PAIR_AND_FOUR_PORT_PROVIDER_DEGREES_MATCH",
      "source_pairs":12,"original_WPI_ports":4,
      "PortWatch_modeled_positive_directions":sum(x["portwatch_static_model_directed_record_count"]>0 for x in pair.values()),
      "Zenodo2020_matched_source_positive_directions":sum(x["verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED"]>0 for x in pair.values()),
      "both_sources_positive_directions":sum(x["portwatch_static_model_directed_record_count"]>0 and x["verified_2020_Zenodo_WPI_source_triplet_count_OR_0_IF_NOT_MATCHED"]>0 for x in pair.values()),
      "not_actual_2020_traffic_or_butterfly_inference":True,
    }
    print("WPI4_PORTWATCH_PUBLISHED_OFFLINE_AUDIT_JSON="+json.dumps(result,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
