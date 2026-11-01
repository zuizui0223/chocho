#!/usr/bin/env python3
"""Source-locked VOC cargo commodity GENUS compatibility with 239 HOSTS butterflies.

Original 17c VOC historical cargo descriptions are commodity group mentions
(e.g. pepper or cinnamon), not taxonomically confirmed viable plant species.
The exact genus lookup is a screening crosswalk to modern recorded larval host
genus names for 239 previously frozen butterfly species. Never infer
shipment of living larval host plants, 1600s botanical introduction,
butterfly arrival, or consumption at an island.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter, defaultdict
from pathlib import Path

BASE=Path(__file__).resolve().parent
HOSTS=Path("_external/HOSTS/resource.csv")
PANEL=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
LOCK=BASE/"voc_17c_commodity_to_HOSTS_genus_lock_20261101.json"
CARGO=BASE/"voc_17c_cargo_on_verified_island_routes_v1.json"
OUT=BASE/"voc_17c_cargo_HOSTS_genus_239_screen_v1.json"


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def csv_rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        yield from csv.DictReader(f)


def screen(hosts=HOSTS,panel=PANEL,lock=LOCK,cargo=CARGO):
    pre=json.loads(lock.read_text(encoding="utf8"))
    if pre["schema"]!="chocho_17c_VOC_commodity_name_to_focal239_HOSTS_genus_screen_v1":
        raise RuntimeError("Original pre-response commodity genus mapping lock changed")
    named=[r["species"].strip() for r in csv_rows(panel)]
    if len(named)!=239 or len(set(named))!=239:
        raise RuntimeError("Original GEB 239 butterfly panel changed")
    selection={name.casefold() for name in named}
    cats=pre["genus_mapping"]
    setgenus={g.casefold() for values in cats.values() for g in values}
    b_by_genus=defaultdict(set)
    hosts_by_genus=defaultdict(set)
    focal_hosts=defaultdict(set)
    n=Counter()
    for row in csv_rows(hosts):
        n["source_HOSTS_rows"]+=1
        butterfly=(str(row.get("Insect Genus") or "").strip()+" "+
                   str(row.get("Insect Species") or "").strip()).strip()
        if butterfly.casefold() not in selection:
            continue
        n["focal239_source_host_rows"]+=1
        genus=str(row.get("Hostplant Genus") or "").strip()
        host_ep=str(row.get("Hostplant Species") or "").strip()
        if not genus or not host_ep:
            n["focal_host_genus_or_epithet_missing"]+=1
            continue
        focal_hosts[butterfly.casefold()].add((genus.casefold(),host_ep.casefold()))
        if genus.casefold() in setgenus:
            b_by_genus[genus.casefold()].add(butterfly)
            hosts_by_genus[genus.casefold()].add((genus.casefold(),host_ep.casefold()))
            n["source_focal_candidate_genus_host_rows"]+=1
    if len(focal_hosts)!=239:
        raise RuntimeError("Original 239 HOSTS butterfly coverage changed: "+str(len(focal_hosts)))
    voyage=json.loads(cargo.read_text(encoding="utf8"))
    if voyage["verified_historic_island_edges"]!=9:
        raise RuntimeError("Exact 9 original historic direct-link results drifted")
    out_categories={}
    for cat,genera in cats.items():
        g0=[g.casefold() for g in genera]
        setbutterflies=set().union(*(b_by_genus.get(g,set()) for g in g0))
        sethosts=set().union(*(hosts_by_genus.get(g,set()) for g in g0))
        out_categories[cat]={
            "historical_source_linked_cargo_keyword_rows":
                voyage["on_historic_island_routes_cargo_category_keyword_row_counts"].get(cat,0),
            "predefined_possible_modern_host_genera":genera,
            "actual_frozen239_HOSTS_butterflies_with_at_least_one_matching_host_genus":
                len(setbutterflies),
            "focal_HOSTS_recorded_distinct_host_genus_species_names":
                len(sethosts),
            "focal_HOSTS_butterfly_names_preview":
                sorted(setbutterflies)[:12],
            "source_has_any_239_host_genus_overlap":bool(setbutterflies),
        }
    route=[]
    for p in voyage["pair_specific_aggregates"]:
        categories=p["source_cargo_category_counts"]
        hits=[x for x in categories
              if x in out_categories and
              out_categories[x]["source_has_any_239_host_genus_overlap"]]
        route.append({
            "island_a":p["island_a"],"island_b":p["island_b"],
            "source_voyage_leg_rows":p["observed_original_voyage_leg_rows"],
            "historical_cargo_rows_joined_exact_ID":p["source_cargo_rows_joined_by_unrepeated_REISN"],
            "commodity_categories_with_ANY_globally_recorded_focal239_HOSTS_genus":
                sorted(hits),
            "number_such_categories":len(hits),
            "number_of_commodity_keyword_mentions_in_such_categories":
                sum(categories[x] for x in hits),
            "is_biological_propagule_or_local_larval_host_proven":False,
        })
    return {
        "schema":"chocho_VOC_17c_cargo_commodity_genus_vs_GEB239_HOSTS_screen_v1",
        "status":"SOURCE_LOCKED_BROAD_GENUS_COMPATIBILITY_SCREEN_NOT_HISTORIC_BUTTERFLY_RESOURCE_INTRODUCTION",
        "sources":{
            "HOSTS_pinned_commit":pre["host_network_commit"],
            "HOSTS_original_resource_file_sha256":digest(hosts),
            "GEB239_frozen_panel_sha256":digest(panel),
            "historical_cargo_screen_receipt_sha256":digest(cargo),
            "mapping_lock_sha256":digest(lock),
        },
        "focal_source_coverage":dict(n),
        "focal239_with_host_genus_records":len(focal_hosts),
        "matched_historical_commodity_categories":sum(
            int(x["historical_source_linked_cargo_keyword_rows"]>0)
            for x in out_categories.values()),
        "categories_with_both_historical_cargo_mentions_and_source_HOSTS_genus":
            [name for name,x in out_categories.items()
             if x["historical_source_linked_cargo_keyword_rows"]>0
             and x["source_has_any_239_host_genus_overlap"]],
        "candidate_genus_mapping_by_commodity_category":out_categories,
        "source_nine_historical_pair_candidate_genus_overlap":route,
        "limitations":pre["limits"]+[
            "The original HOSTS genus taxonomy is not guaranteed to be taxonomically reconciled to current WCVP accepted species.",
            "The matching species universe is the globally pooled fixed239, not necessarily the actual butterfly species recorded on an island.",
            "Cargo commodity genus associations are source-informed hypotheses and cannot be treated as records of exact species cargo.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--hosts",type=Path,default=HOSTS)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    result=screen(hosts=args.hosts)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k not in
          ("limitations","candidate_genus_mapping_by_commodity_category",
           "source_nine_historical_pair_candidate_genus_overlap")},
          ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
