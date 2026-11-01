#!/usr/bin/env python3
"""Join observed *positive* 2020 AIS island links to historical empire sources.

Two historically different network objects, never equated:
(1) 2020 source AIS-derived from/to port voyage rows (incomplete source);
(2) Lenzner 2022 source-weighted regional history of four 1492-2010
    European empires, joined geographically via currently indexed WPI ports.

"Shared former imperial relation" is not the whole historical maritime
network: metropolitan countries can have zero *duration of colonization*
within their own mainland while connected to a former dependent territory.
Report TWO separate descriptive edge categories (shared duration and
metropole-country-to-historic-contact), never infer shipping era transfer.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
ORIGINAL=BASE.parent/"island_host_network_20261030/island_host_network_pairs.csv"
HISTORIC=BASE/"island_port_colonial_10km.islands.csv"
WPI=BASE/"wpi_51island_194ports_geocoded.csv"
COUNTRY_SHIPPING=BASE/"full_lsbci2010_island_pairs.csv"
LOCK=BASE/"ais_2020_island_voyage_lock_20261101.json"
METROPOLITAN={"GB":"years_GB","ES":"years_ESP","PT":"years_PRT","NL":"years_NL"}
EMPIRES=tuple(METROPOLITAN.values())


def read_csv(path):
    with Path(path).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def key(a,b):
    return tuple(sorted((a,b),key=int))


def hist(source):
    return (source is not None and source["eligible_for_history_test"]=="1")


def classify_pair(a,b,countries):
    if not (hist(a) and hist(b)):
        return None,None,[]
    shared=[name.removeprefix("years_") for name in EMPIRES
            if min(float(a[name]),float(b[name]))>=25]
    metropole=[]
    for first,second in ((a,b),(b,a)):
        codes=first["current_WPI_country_codes"].split(";")
        for code in codes:
            empire=METROPOLITAN.get(code)
            if empire is not None and float(second[empire])>=25:
                metropole.append(code)
    return bool(shared),bool(metropole),sorted(set(metropole))


def connected_components(edges):
    adjacency=defaultdict(set)
    for a,b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    visited=set()
    components=[]
    for node in sorted(adjacency,key=int):
        if node in visited:
            continue
        group=[]
        stack=[node]
        while stack:
            item=stack.pop()
            if item in visited:
                continue
            visited.add(item)
            group.append(item)
            stack+=list(adjacency[item]-visited)
        components.append(sorted(group,key=int))
    return sorted(components,key=lambda x:(-len(x),int(x[0])))


def analyse(km=5):
    if km not in (2,5,10):
        raise ValueError("Must use predeclared AIS WPI match radii")
    if not LOCK.exists():
        raise RuntimeError("Predeclared 2020 positive AIS source lock missing")
    ais=read_csv(BASE/f"ais_2020_wpi51_island_{km}km.island_pairs.csv")
    orig={key(r["island_a"],r["island_b"]):r for r in read_csv(ORIGINAL)}
    history={r["landmass_id"]:r for r in read_csv(HISTORIC)}
    ship={key(r["island_a"],r["island_b"]):r for r in read_csv(COUNTRY_SHIPPING)}
    countries=defaultdict(set)
    for w in read_csv(WPI):
        countries[w["landmass_id"]].add(w["country"])
    if len(ais)!=1275 or len(orig)!=1275 or len(ship)!=1275:
        raise RuntimeError("Source island-pair or shipping cohort changed")
    matched=[]
    for r in ais:
        if r["has_positive_2020_AIS_observed_from_to_pair"]!="1":
            if r["ais_2020_direct_voyages_observed_in_sample"]:
                raise RuntimeError("Unobserved AIS source pair was given a route weight")
            continue
        a,b=r["island_a"],r["island_b"]
        original=orig[key(a,b)]
        x,y=history.get(a),history.get(b)
        common,metro,metro_country=classify_pair(x,y,countries)
        shipping=ship[key(a,b)]
        if not (float(r["ais_2020_direct_voyages_observed_in_sample"])>0):
            raise RuntimeError("Nonpositive direct harbor voyage")
        matched.append({
            "island_a":a,"current_WPI_country_a":";".join(sorted(countries[a])),
            "island_b":b,"current_WPI_country_b":";".join(sorted(countries[b])),
            "observed_2020_AIS_voyage_source_rows":int(r["ais_2020_direct_voyages_observed_in_sample"]),
            "port_distance_km":float(original["port_distance_km"]),
            "candidate_introduced_host_jaccard":float(original["host_jaccard"]),
            "newly_GBIF_indexed_butterfly_jaccard":float(original["new_butterfly_jaccard"]),
            "old_GBIF_indexed_butterfly_jaccard":float(original["earlier_butterfly_jaccard"]),
            "both_islands_have_Lenzner_source_history":int(hist(x) and hist(y)),
            "shared_four_empire_duration_25years":"" if common is None else int(common),
            "metropolitan_source_country_to_hist_contact_25years":"" if metro is None else int(metro),
            "metropolitan_source_country_if_connected":";".join(metro_country),
            "2010_UNCTAD_country_bilateral_score":(
                shipping["lsbci_2010"]
                if shipping["status"]=="INDEX_2010_OBSERVED" else ""),
            "2010_UNCTAD_coverage_status":shipping["status"],
        })
    if not matched:
        raise RuntimeError("No positive direct AIS island links in source, skip ecology")
    links=[key(z["island_a"],z["island_b"]) for z in matched]
    components=connected_components(links)
    top=max(matched,key=lambda z:z["observed_2020_AIS_voyage_source_rows"])
    excl=[r for r in matched if r is not top]
    audit={
        "schema":"chocho_island_historic_colonial_and_2020_AIS_positive_links_v1",
        "status":"SOURCE_BACKED_POSITIVE_ROUTE_PAIRS_X_HOST_BUTTERFLY_HISTORY_NOT_CAUSAL",
        "data_sources":{
            "historical_Lenzner":"10.5281/zenodo.6640264",
            "AIS_2020_Duvallet":"10.5281/zenodo.10783698",
            "source_original_51island_pairs_sha256":sha(ORIGINAL),
            "geocoded_51_island_history_sha256":sha(HISTORIC),
            "observed_2020_AIS_positive_pairs_sha256":
                sha(BASE/f"ais_2020_wpi51_island_{km}km.island_pairs.csv"),
            "country_shipping_2010_sha256":sha(COUNTRY_SHIPPING),
            "AIS_before_result_design_lock_sha256":sha(LOCK),
        },
        "fixed_AIS_unique_port_match_radius_km":km,
        "original_island_pairs":1275,
        "positive_direct_2020_AIS_island_pair_count":len(matched),
        "distinct_linked_islands":len({i for edge in links for i in edge}),
        "recorded_AIS_route_row_total":sum(r["observed_2020_AIS_voyage_source_rows"] for r in matched),
        "links_with_positive_host_candidate_similarity":sum(r["candidate_introduced_host_jaccard"]>0 for r in matched),
        "links_with_positive_new_GBIF_butterfly_similarity":sum(r["newly_GBIF_indexed_butterfly_jaccard"]>0 for r in matched),
        "both_islands_historic_source_count":sum(r["both_islands_have_Lenzner_source_history"]==1 for r in matched),
        "colonial_shared_25years_positive_link_pairs":sum(r["shared_four_empire_duration_25years"]==1 for r in matched),
        "metropolitan_country_to_contact_positive_link_pairs":sum(r["metropolitan_source_country_to_hist_contact_25years"]==1 for r in matched),
        "2020_link_graph_components_on_positive_records":components,
        "largest_source_route_edge":{
            "islands":key(top["island_a"],top["island_b"]),
            "raw_source_voyage_rows":top["observed_2020_AIS_voyage_source_rows"],
            "fraction_of_all_matched_source_interisland_route_rows":
                top["observed_2020_AIS_voyage_source_rows"]/
                sum(r["observed_2020_AIS_voyage_source_rows"] for r in matched),
        },
        "after_largest_link_removed":{
            "positive_link_pairs":len(excl),
            "positive_route_rows":sum(r["observed_2020_AIS_voyage_source_rows"] for r in excl),
            "pairs_with_host_similarity_gt_zero":sum(r["candidate_introduced_host_jaccard"]>0 for r in excl),
            "pairs_with_new_butterfly_similarity_gt_zero":sum(r["newly_GBIF_indexed_butterfly_jaccard"]>0 for r in excl),
        },
        "unverified_false_routes_or_absences":[
            "Do not compare observed AIS route links against 1250+ unobserved island pairs as if route-negative.",
            "A positive direct 2020 AIS row does not identify which goods, plant propagules or butterfly life stages were carried.",
            "AIS source coverage concentrated in some seas and ship types; observed links are not a global connection census.",
            "2019 WPI harbor and 2020 AIS voyages are centuries later than much of the colonial contact source data.",
            "Source historical empire metrics cover only GB, ESP, PRT, NL; French/other empires absent.",
            "Metropolitan-country link is defined using a *current* WPI country code and a Lenzner historical-area proxy, not documented 1800s port operation.",
            "Candidate plants use 2001-2024 indexing, butterfly response uses 2011-2022 first GBIF indexing; not guaranteed time precedence.",
            "12 link observations reuse only 9 island nodes and offer no independent geospatial causality proof.",
        ],
    }
    return matched,audit


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--wpi-radius-km",type=int,default=5)
    ap.add_argument("--output-prefix",type=Path,
                    default=BASE/"ais_2020_island_colonial_overlap_5km")
    args=ap.parse_args()
    pairs,result=analyse(args.wpi_radius_km)
    args.output_prefix.parent.mkdir(parents=True,exist_ok=True)
    with Path(str(args.output_prefix)+".positive_links.csv").open(
        "w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(pairs[0]))
        w.writeheader()
        w.writerows(pairs)
    Path(str(args.output_prefix)+".audit.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
