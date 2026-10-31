#!/usr/bin/env python3
"""A strict maritime *pairwise edge* admission gate for the 51-island cohort.

This is input QA, not an empirical shipping ecology result. Never fabricate
routes from port-level PLSCI scores or AIS density rasters. A missing edge is
treated as observed zero only when both island port sets have complete edge
coverage for the requested years. Ship movements and freight-model output
are kept in different variables, without unit conversion.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

PAIRS = Path(__file__).resolve().parents[1] / "island_host_network_20261030" / "island_host_network_pairs.csv"
OBSERVED = "AIS_DIRECT_PORT_PAIR"
MODELED = "MODEL_ALLOCATED_PORT_FREIGHT"
MODES = {OBSERVED, MODELED}


def rows(path):
    with Path(path).open(encoding="utf8", newline="") as f:
        return list(csv.DictReader(f))


def required(data, columns, name):
    if not data:
        raise ValueError(name + ": empty table is not evidence of zero links")
    absent = set(columns) - set(data[0])
    if absent:
        raise ValueError(name + " missing columns " + ",".join(sorted(absent)))


def scan(mapping, edge_rows, study_start, study_end, pair_rows):
    """Build island-pair direct link statistics with explicit source coverage.

    mapping columns: port_id,island_id,source_port_id,coverage_start_year,
    coverage_end_year,complete_observation. The source port ID may differ from
    WPI port ID, but must have been independently geocoded/cross-checked.

    edges columns: from_port_id,to_port_id,year,weight,measurement_type.
    Values are direct AIS vessel-link counts or model-allocated freight ONLY.
    """
    if study_start > study_end:
        raise ValueError("reversed time interval")
    required(mapping, ["port_id", "island_id", "source_port_id",
                       "coverage_start_year", "coverage_end_year",
                       "complete_observation"], "mapping")
    required(edge_rows, ["from_port_id", "to_port_id", "year", "weight",
                         "measurement_type"], "edges")
    required(pair_rows, ["island_a", "island_b", "host_jaccard",
                         "new_butterfly_jaccard", "earlier_butterfly_jaccard",
                         "port_distance_km"], "island pairs")

    source_to_island = {}
    per_island = {}
    for x in mapping:
        p, s, island = x["port_id"], x["source_port_id"], x["island_id"]
        if not p or not s or not island:
            raise ValueError("port IDs and island IDs are all mandatory")
        if s in source_to_island and source_to_island[s] != island:
            raise ValueError("source port is mapped to multiple islands")
        source_to_island[s] = island
        supported = (x["complete_observation"] == "True"
                     and int(x["coverage_start_year"]) <= study_start
                     and int(x["coverage_end_year"]) >= study_end)
        per_island.setdefault(island, []).append((s, supported))

    pair_link = {}
    drop = {"outside_period": 0, "unmatched_port": 0,
            "unknown_measurement": 0}
    for e in edge_rows:
        kind = e["measurement_type"]
        if kind not in MODES:
            drop["unknown_measurement"] += 1
            raise ValueError(
                "A port-node score, raster density or undocumented edge cannot be a route: "
                + kind)
        year = int(e["year"])
        if not (study_start <= year <= study_end):
            drop["outside_period"] += 1
            continue
        w = float(e["weight"])
        if not 0 <= w < float("inf"):
            raise ValueError("weights must be finite non-negative source units")
        a = source_to_island.get(e["from_port_id"])
        b = source_to_island.get(e["to_port_id"])
        if a is None or b is None:
            drop["unmatched_port"] += 1
            continue
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        q = pair_link.setdefault(key, {OBSERVED: 0.0, MODELED: 0.0,
                                      "route_edge_rows": 0})
        q[kind] += w
        q["route_edge_rows"] += 1

    out = []
    seen = set()
    for p in pair_rows:
        a, b = p["island_a"], p["island_b"]
        if a == b:
            raise ValueError("island pair cannot reference itself")
        key = (min(a, b), max(a, b))
        if key in seen:
            raise ValueError("duplicate island pair")
        seen.add(key)
        pa = per_island.get(a, [])
        pb = per_island.get(b, [])
        # A zero is observable only for fully covered island port sets.
        supported = bool(pa and pb and all(x[1] for x in pa + pb))
        links = pair_link.get(key, {})
        out.append({
            "island_a": a, "island_b": b,
            "study_start_year": study_start, "study_end_year": study_end,
            "n_matching_ports_a": len(pa), "n_matching_ports_b": len(pb),
            "island_pair_covered": supported,
            "observation_state": "COVERED" if supported else "UNKNOWN",
            "observed_direct_link_weight": (links.get(OBSERVED, 0.) if supported else None),
            "modeled_freight_link_weight": (links.get(MODELED, 0.) if supported else None),
            "matched_edge_rows": links.get("route_edge_rows", 0) if supported else None,
            "host_jaccard": float(p["host_jaccard"]),
            "new_butterfly_jaccard": float(p["new_butterfly_jaccard"]),
            "earlier_butterfly_jaccard": float(p["earlier_butterfly_jaccard"]),
            "port_distance_km": float(p["port_distance_km"]),
        })
    return {
        "status": "LINK_INPUT_QA_NOT_SHIPPING_EFFECT",
        "period": [study_start, study_end],
        "source_mapping_ports": len(source_to_island),
        "source_mapping_islands": len(per_island),
        "island_pairs": len(out),
        "covered_pairs": sum(x["island_pair_covered"] for x in out),
        "uncovered_pairs": sum(not x["island_pair_covered"] for x in out),
        "matched_link_rows": sum(x["matched_edge_rows"] or 0 for x in out),
        "source_rows_dropped": drop,
        "pair_rows": out,
        "limitations": [
            "AIS direct links and model-allocated freight are not interchangeable.",
            "A node-connectivity score or ship-density pixel is not a port-to-port edge.",
            "Unmatched island pairs are missing, not zero-travel connections.",
            "Port name alone is not adequate independent port ID geolocation.",
            "2019-2020 edges cannot be interpreted as pre-2011 host transports.",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--port-map", type=Path, required=True)
    p.add_argument("--edges", type=Path, required=True)
    p.add_argument("--start-year", type=int, required=True)
    p.add_argument("--end-year", type=int, required=True)
    p.add_argument("--output-json", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    args = p.parse_args()
    report = scan(rows(args.port_map), rows(args.edges),
                  args.start_year, args.end_year, rows(PAIRS))
    pair_rows = report.pop("pair_rows")
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
    with args.output_csv.open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
