#!/usr/bin/env python3
"""Read-only availability audit for maritime alien-host / island stepping-stone chain.

Descriptive cohort-support counts only. Does not estimate arrival, establishment,
historical ship-borne seeds or island-to-island movement. Original longitude/port
prediction outcomes are not modified. Source Git blobs are pinned.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = Path("provenance/reviewer_defenses/results")
PATHS = {
    "distance": BASE / "original_island_continuity_20261030/global_original_continuity_distance_summary.csv",
    "origin": BASE / "original_island_continuity_20261030/global_original_continuity_origin_summary.csv",
    "links": BASE / "island_dated_host_match_20261030/island_dated_host_match_links.csv",
    "plant_comparison": BASE / "island_dated_host_match_20261030/island_dated_host_match_comparison.csv",
    "events": BASE / "island_host_continuity_20261030/island_host_continuity_events.csv",
    "continuity": BASE / "island_host_continuity_20261030/island_host_continuity_summary.csv",
    "port_connectivity": BASE / "port_connectivity_archive_20261030/port_connectivity_join_audit.json",
}
PINS = {
    "distance": "90055903327e47c432252773bf966861c70a3e5b",
    "origin": "839ca8d9e9f443c08fdec8e08df9eec02a530189",
    "links": "88db773c74dea9cd2d58bde0b67a394d93226869",
    "plant_comparison": "15b2282b0cd255a2f3d518e495727353bee8b1e0",
    "events": "efbc770469a720e2dea82478db2f6196ba03b353",
    "continuity": "7405a41c8020022952ee58081c0ec35820564953",
    "port_connectivity": "fb9ab4b9da52411eca824518d76f99ad1b41b035",
}


def read_source(root: Path, name: str):
    path = root / PATHS[name]
    body = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(body)).encode() + b"\0" + body).hexdigest()
    if blob != PINS[name]:
        raise ValueError(f"Changed locked source: {name}")
    text = body.decode("utf-8")
    if name == "port_connectivity":
        return json.loads(text)
    return list(csv.DictReader(text.splitlines()))


def audit(root: Path = ROOT):
    dist = read_source(root, "distance")
    origin = read_source(root, "origin")
    links = read_source(root, "links")
    plant = read_source(root, "plant_comparison")
    events = read_source(root, "events")
    continuity = read_source(root, "continuity")
    port = read_source(root, "port_connectivity")

    distance = {}
    for cohort in ("all", "prior_opportunity", "both_period_opportunity"):
        rows = [r for r in dist if r["cohort"] == cohort]
        orig = [r for r in origin if r["cohort"] == cohort]
        if len(rows) != 4 or len(orig) != 2:
            raise ValueError(f"Missing original stage summary: {cohort}")
        total = int(rows[0]["total_events"])
        if sum(int(x["events"]) for x in rows) != total:
            raise ValueError("Distance bands do not reconcile")
        if sum(int(x["events"]) for x in orig) != total:
            raise ValueError("Prior-source landmass categories do not reconcile")
        # Prior source label refers to an earlier indexed occurrence,
        # not a documented breeding population on that island.
        by_bool = {r["prior_record_on_island"]: r for r in orig}
        far = next(r for r in rows if r["distance_band"] == "500 km以上")
        near = sum(int(r["events"]) for r in rows if r["distance_band"] in ("10 km未満", "10–100 km"))
        far_n = int(far["events"])
        near_island = int(by_bool["True"]["near100_events"])
        far_island = int(by_bool["True"]["far500_events"])
        if near_island + int(by_bool["False"]["near100_events"]) != near:
            raise ValueError("Nearer prior-landmass subtotals disagree")
        if far_island + int(by_bool["False"]["far500_events"]) != far_n:
            raise ValueError("Remote prior-landmass subtotals disagree")
        distance[cohort] = {
            "island_first_record_events": total,
            "apparent_nearest_prior_record_gap_ge500_km": far_n,
            "apparent_nearest_prior_record_gap_lt100_km": near,
            "far500_nearest_prior_indexed_record_on_island": far_island,
            "far500_nearest_prior_indexed_record_on_island_fraction": far_island / far_n,
            "near100_nearest_prior_indexed_record_on_island": near_island,
            "near100_nearest_prior_indexed_record_on_island_fraction": near_island / near,
            "far_minus_near_prior_island_fraction_points": 100 * (far_island/far_n - near_island/near),
        }

    if len(events) != 5665 or len(links) != 5:
        raise ValueError("Original butterfly and dated-host cohorts changed")
    event_by_key = {}
    for e in events:
        key = (e["species"], e["landmass_id"])
        if key in event_by_key:
            raise ValueError("Duplicate species×island five-year event")
        event_by_key[key] = e

    groups = Counter(e["host_group"] for e in events)
    if groups["NOT_EVALUABLE"] != 5333 or len(events) - groups["NOT_EVALUABLE"] != 332:
        raise ValueError("Retrospective port-host support changed")
    if len({l["species"] for l in links}) != 5:
        raise ValueError("Distinct dated-host butterfly count changed")
    edges = set()
    link_details = []
    for link in links:
        endpoints = []
        for label in ("island_a", "island_b"):
            key = (link["species"], link[label])
            edges.add(key)
            found = event_by_key.get(key)
            endpoints.append({
                "island": link[label],
                "has_2011_2017_five_year_window": found is not None,
                "first_indexed_year": int(found["first_record_year"]) if found else None,
                "any_followup_indexed_record": found["any_followup_record"] == "True" if found else None,
                "year5_indexed_record": found["record_in_year5"] == "True" if found else None,
            })
        link_details.append({
            "butterfly": link["species"],
            "early_known_host": link["host_name"],
            "both_island_host_record_window": "2001-2010",
            "paired_islands": endpoints,
        })
    dated_endpoints = [event_by_key[x] for x in edges if x in event_by_key]
    early = next(r for r in plant if r["plant_window"] == "Plants2001-2010")
    broad = next(r for r in plant if r["plant_window"] == "Plants2001-2024")
    if round(float(early["observed"]) * 1163) != 5 or round(float(broad["observed"])*1163) != 208:
        raise ValueError("Original early-vs-full port-host match comparison changed")
    if (len(edges), len(dated_endpoints)) != (10, 6):
        raise ValueError("Matched endpoint support changed")
    if (port["complete_plant_ports_with16of20_early_quarters"],
        port["complete_plant_islands_with16of20_early_quarters"],
        port["previous51_islands_with16of20_early_quarters"]) != (39, 27, 21):
        raise ValueError("Port chronology coverage changed")

    return {
        "schema": "chocho_maritime_alien_host_steppingstone_support_20261030_v1",
        "role": "observational_support_audit_not_a_mechanism_test",
        "all_source_git_blob_sha": PINS,
        "separate_cohorts_not_joint_risk_set": True,
        "original_coordinate_butterfly_history": distance,
        "early_port_host_link_support": {
            "common_new_butterfly_species_island_pair_denominator": 1163,
            "retrospective_plant_2001_2024_matches": 208,
            "genuinely_earlier_plant_2001_2010_matches": 5,
            "earlier_host_pairs": len(links),
            "earlier_host_butterfly_species": len({l["species"] for l in links}),
            "earlier_host_islands": len({l[k] for l in links for k in ("island_a", "island_b")}),
            "unique_species_island_endpoints_from_pairs": len(edges),
            "endpoints_observed_in_separate_2011_2017_followup_cohort": len(dated_endpoints),
            "endpoints_with_any_subsequent_record": sum(e["any_followup_record"] == "True" for e in dated_endpoints),
            "endpoints_with_record_in_year_five": sum(e["record_in_year5"] == "True" for e in dated_endpoints),
            "selected_endpoint_is_not_unbiased_success_comparison": True,
            "links": link_details,
        },
        "separate_five_year_followup_cohort": {
            "original_species_island_first_record_events": len(events),
            "port_host_group_matched_but_not_dated": len(events) - groups["NOT_EVALUABLE"],
            "not_evaluable": groups["NOT_EVALUABLE"],
            "host_groups": dict(sorted(groups.items())),
            "plant_window_for_host_groups": "2001-2024 (contains post-butterfly records)",
            "not_a_prior_host_colonization_effect": True,
        },
        "prior_shipping_connectivity_support": {
            "2006_2010_ports_with_historical_connectivity_and_complete_plant_inventory":
                port["complete_plant_ports_with16of20_early_quarters"],
            "such_islands": port["complete_plant_islands_with16of20_early_quarters"],
            "islands_overlap_prior_51_island_host_network": port["previous51_islands_with16of20_early_quarters"],
            "port_name_join_is_not_geo_verified": True,
            "liner_shipping_connectivity_is_not_vessel_route_or_seed_transport": True,
        },
        "mechanisms_identified_now": {
            "actual_seed_attachment_for_focal_host_species": False,
            "prior_confirmed_introduced_host_on_specific_island": False,
            "butterfly_local_breeding_or_establishment": False,
            "island_i_to_island_j_direct_movement": False,
            "joint_long_distance_x_earlier_host_x_recurrence_x_downstream_island": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.root)
    data = json.dumps(result, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        output = args.output if args.output.is_absolute() else args.root / args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(data, encoding="utf8")
        print(output)
    else:
        print(data)


if __name__ == "__main__":
    main()
