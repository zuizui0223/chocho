#!/usr/bin/env python3
"""Route-specific 17c VOC cargo genus × 2011–2017 GBIF-indexed island butterflies.

SOURCE LEVELS:
- Original VOC 1595–1698 island routes: 9 positive links, with unambiguous
  shipment cargo recorded for 7. Freight commodity words are NOT live host
  plants, the exact species being transported, or goods offloaded on the leg.
- Verified VOC commodity -> contemporary HOSTS genus compatibility:
  original fixed 239 species, exact selected genus names (screen receipt).
- The GBIF independent outcome is a first-indexed butterfly SPECIES on a
  particular island in 2011–2017; no HOSTS information enters that first-
  index outcome. It is NOT confirmed arrival, local host-use or colonization.
- Diagnostic: permute seven complete commodity category PACKETS among the
  7 original cargo-supported historic route edges (7!=5040 assignments).
  This is a source-conditional label test, NOT a colonial trade effect p value.

Never infer negatives for historic links with no indexed cargo data or
19th/20th/21st century absence of a biological record. No raw historic
voyage row, cargo description, or GBIF point is published.
"""
from __future__ import annotations

from collections import Counter,defaultdict
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

BASE=Path(__file__).resolve().parent
SOURCE_CARGO=BASE/"voc_17c_cargo_on_verified_island_routes_v1.json"
SOURCE_GENUS=BASE/"voc_17c_cargo_HOSTS_genus_239_screen_v1.json"
SOURCE_EVENTS=BASE.parent/"island_host_continuity_20261030/island_host_continuity_events.csv"
SOURCE_LOCK=BASE/"voc_17c_cargo_island_specific_butterfly_match_lock_20261102.json"
OUTPUT=BASE/"voc_17c_cargo_actual_island_butterfly_first_records_v1.json"
PAIR_CSV=BASE/"voc_17c_cargo_actual_island_butterfly_first_records_v1.csv"
TARGET_ISLANDS=("12","17","19","31","43","44")
PERIOD=range(2011,2018)
EXPECTED_GBIF_SOURCE_EVENTS=5665
EXPECTED_CARGO=7737


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def source_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def load_source_receipts(cargo=SOURCE_CARGO,screen=SOURCE_GENUS,lock=SOURCE_LOCK):
    study=json.loads(lock.read_text(encoding="utf8"))
    if study.get("schema")!="chocho_VOC_cargo_local_butterflies_exact_route_null_lock_v1":
        raise RuntimeError("Original pre-response source lock mismatch")
    c=json.loads(cargo.read_text(encoding="utf8"))
    g=json.loads(screen.read_text(encoding="utf8"))
    if c.get("verified_historic_island_edges")!=9 or c.get("historical_source_original_cargo_rows")!=EXPECTED_CARGO:
        raise RuntimeError("Nine authenticated original cargo/historic route cohort changed")
    if g.get("focal239_with_host_genus_records")!=239:
        raise RuntimeError("Frozen 239-host genus screen changed")
    cargo_rows=c.get("pair_specific_aggregates") or []
    genus_rows=g.get("candidate_genus_mapping_by_commodity_category") or {}
    if len(cargo_rows)!=9:
        raise RuntimeError("Nine archived island-route pair aggregates absent")
    butterfly_by_cat={}
    for cat,data in genus_rows.items():
        count=data["actual_frozen239_HOSTS_butterflies_with_at_least_one_matching_host_genus"]
        names=data["focal_HOSTS_butterfly_names_preview"]
        if count!=len(names):
            raise RuntimeError("Source names preview is truncated for "+cat+"; must rescan original HOSTS")
        if len(set(x.casefold() for x in names))!=len(names):
            raise RuntimeError("Duplicate names under original HOSTS genus screen for "+cat)
        butterfly_by_cat[cat]={name.casefold() for name in names}
    edges={}
    for r in cargo_rows:
        a,b=str(r["island_a"]),str(r["island_b"])
        k=tuple(sorted((a,b),key=int))
        if k in edges or a==b or not set(k)<=set(TARGET_ISLANDS):
            raise RuntimeError("Source original 9 island-pair labels or IDs changed")
        cats=r["source_cargo_category_counts"]
        if sum(cats.values())!=r["island_pair_cargo_mentions_in_predeclared_plant_commodity_categories"]:
            raise RuntimeError("Historical original cargo categories no longer reconcile")
        if not set(cats)<=set(butterfly_by_cat):
            raise RuntimeError("New cargo plant category lacks predeclared genus mapping")
        supported={cat for cat,count in cats.items() if count>0 and butterfly_by_cat[cat]}
        expected=g["source_nine_historical_pair_candidate_genus_overlap"]
        line=next((x for x in expected if
            tuple(sorted((x["island_a"],x["island_b"]),key=int))==k),None)
        if line is None or set(line["commodity_categories_with_ANY_globally_recorded_focal239_HOSTS_genus"])!=supported:
            raise RuntimeError("Precomputed global-genus cargo source mismatch at route "+str(k))
        edges[k]={
            "historic_voyage_leg_rows":r["observed_original_voyage_leg_rows"],
            "historic_cargo_indexed_rows":r["source_cargo_rows_joined_by_unrepeated_REISN"],
            "categories":cats,
            "original_category_mentions":sum(cats.values()),
            "globally_host_compatible_categories":supported,
        }
    if sum(x["historic_voyage_leg_rows"] for x in edges.values())!=1193:
        raise RuntimeError("Original authentic VOC direct island route totals changed")
    if sum(x["historic_cargo_indexed_rows"] for x in edges.values())!=504:
        raise RuntimeError("Original cargo source ID join totals changed")
    if sum(x["historic_cargo_indexed_rows"]>0 for x in edges.values())!=7:
        raise RuntimeError("Exact seven cargo-indexed historic route count changed")
    return edges,butterfly_by_cat,study


def independently_indexed_butterflies(path=SOURCE_EVENTS):
    events=source_rows(path)
    if len(events)!=EXPECTED_GBIF_SOURCE_EVENTS:
        raise RuntimeError("Original independent first-indexed butterfly event source changed")
    new=defaultdict(set)
    seen=set()
    year_counts=Counter()
    dropped=Counter()
    for row in events:
        island=str(row["landmass_id"])
        sp=str(row["species"]).strip().casefold()
        try:
            year=int(row["first_record_year"])
        except ValueError:
            dropped["missing_first_record_year"]+=1
            continue
        if (sp,island) in seen:
            raise RuntimeError("Duplicated source butterfly first-index species×island event")
        seen.add((sp,island))
        if not year in PERIOD:
            raise RuntimeError("Original 2011–2017 butterfly first-record cohort drifted")
        year_counts[year]+=1
        if island in TARGET_ISLANDS:
            new[island].add(sp)
    for island in TARGET_ISLANDS:
        if not new.get(island):
            raise RuntimeError("Historic island has no independently first indexed butterfly records: "+island)
    if sum(len(new[island]) for island in TARGET_ISLANDS)!=49:
        raise RuntimeError("Original six historical source-island butterfly-first-indexed events changed")
    return dict(new),{"global_original_new_events":len(events),
                      "all_original_year_counts":dict(sorted(year_counts.items())),
                      "source_target_six_islands":{k:len(new[k]) for k in TARGET_ISLANDS},
                      "source_other_unmapped":dict(dropped)}


def source_route_stats(pair,packet,first,butterfly_by_cat):
    """Measure source category compatibility with independently indexed butterflies."""
    a,b=pair
    aa,bb=first[a],first[b]
    compat=set().union(*(butterfly_by_cat.get(cat,set()) for cat in packet))
    a_match=aa&compat
    b_match=bb&compat
    both=aa&bb&compat
    available=sorted(cat for cat in packet if
        (aa|bb)&butterfly_by_cat.get(cat,set()))
    both_cat=sorted(cat for cat in packet if
        aa&butterfly_by_cat.get(cat,set()) and
        bb&butterfly_by_cat.get(cat,set()))
    return {
        "has_either_endpoint_with_compatible_record":int(bool(a_match or b_match)),
        "has_both_endpoints_each_with_compatible_record":int(bool(a_match) and bool(b_match)),
        "has_same_compatible_recorded_butterfly_in_both_islands":int(bool(both)),
        "number_compatible_source_commodity_categories":len(available),
        "number_both_endpoints_compatible_commodity_categories":len(both_cat),
        "compatible_recorded_butterfly_species_a":sorted(a_match),
        "compatible_recorded_butterfly_species_b":sorted(b_match),
        "same_recorded_butterfly_species_both":sorted(both),
        "source_commodity_categories_with_any_local_match":available,
        "source_commodity_categories_compatible_at_both_islands":both_cat,
    }


def exact_permutation_diagnostic(edges,first,butterfly_by_cat):
    selected=[(pair,v) for pair,v in sorted(edges.items(),key=lambda x:x[0])
              if v["historic_cargo_indexed_rows"]>0]
    if len(selected)!=7:
        raise RuntimeError("Exactly seven original cargo-supported link packets required")
    endpoints=[pair for pair,_ in selected]
    packet=[set(v["globally_host_compatible_categories"]) for _,v in selected]
    # Important: the whole original cargo packet is permuted across fixed
    # historic route endpoint pairs, no reallocation of individual commodity
    # counts or cargo-containing archival voyage rows is performed.
    stats=("has_either_endpoint_with_compatible_record",
           "has_both_endpoints_each_with_compatible_record",
           "has_same_compatible_recorded_butterfly_in_both_islands")
    def score(perm):
        aggregate=Counter()
        for idx,pair in enumerate(endpoints):
            r=source_route_stats(pair,packet[perm[idx]],first,butterfly_by_cat)
            for key in stats:
                aggregate[key]+=r[key]
            aggregate["number_compatible_source_commodity_categories"]+=r[
                "number_compatible_source_commodity_categories"]
        return dict(aggregate)
    obs=score(tuple(range(7)))
    samples={key:Counter() for key in (*stats,"number_compatible_source_commodity_categories")}
    total=0
    for perm in itertools.permutations(range(7)):
        result=score(perm)
        total+=1
        for key in samples:
            samples[key][result.get(key,0)]+=1
    if total!=5040:
        raise RuntimeError("Full exact 7! labeled commodity packet permutations not evaluated")
    summary={}
    for key,hist in samples.items():
        number=obs.get(key,0)
        values=[n for n,count in sorted(hist.items()) for _ in range(count)]
        summary[key]={
            "observed":number,
            "null_mean":sum(k*c for k,c in hist.items())/total,
            "null_min":min(hist),
            "null_max":max(hist),
            "null_q025":values[int(.025*(total-1))],
            "null_q975":values[int(.975*(total-1))],
            "conditional_finite_permutation_right_tail_fraction":sum(c for k,c in hist.items() if k>=number)/total,
            "null_histogram":dict(sorted(hist.items())),
        }
    return selected,summary,total


def analyze():
    edges,by_cat,lock=load_source_receipts()
    first,source_audit=independently_indexed_butterflies()
    selected,permutations,perm_n=exact_permutation_diagnostic(edges,first,by_cat)
    outputs=[]
    for (a,b),route in sorted(edges.items(),key=lambda x:x[0]):
        r=source_route_stats((a,b),route["globally_host_compatible_categories"],first,by_cat)
        outputs.append({
            "island_a":a,"island_b":b,
            "historic_17c_island_voyage_legs":route["historic_voyage_leg_rows"],
            "original_cargo_records_linked":route["historic_cargo_indexed_rows"],
            "island_a_first_GBIF_indexed_butterflies_2011_2017":len(first[a]),
            "island_b_first_GBIF_indexed_butterflies_2011_2017":len(first[b]),
            "globally_candidate_larval_host_genera_commodity_categories":";".join(
                sorted(route["globally_host_compatible_categories"])),
            "historic_cargo_commodity_mentions":route["original_category_mentions"],
            "any_local_butterfly_with_same_modern_known_host_genus":r[
                "has_either_endpoint_with_compatible_record"],
            "both_islands_have_some_compatible_butterfly":r[
                "has_both_endpoints_each_with_compatible_record"],
            "same_compatible_butterfly_GBIF_first_indexed_in_both_islands":r[
                "has_same_compatible_recorded_butterfly_in_both_islands"],
            "compatible_local_species_first_index_island_a":";".join(
                r["compatible_recorded_butterfly_species_a"]),
            "compatible_local_species_first_index_island_b":";".join(
                r["compatible_recorded_butterfly_species_b"]),
            "same_compatible_species_both_islands":";".join(
                r["same_recorded_butterfly_species_both"]),
            "commodity_categories_matching_locally_indexed_species":";".join(
                r["source_commodity_categories_with_any_local_match"]),
            "source_row_status":("ORIGINAL_CARGO_INDEXED"
                if route["historic_cargo_indexed_rows"] else "NO_MATCHED_CARGO_ROWS_NOT_NO_PLANTS"),
        })
    supported=[v for v in outputs if v["original_cargo_records_linked"]>0]
    audit={
        "schema":"chocho_VOC_cargo_HostGenus_actual_GBIF_island_butterflies_v1",
        "status":"SOURCE_CONDITIONAL_LOCAL_GENUS_OVERLAP_SCREEN_NOT_BIOLOGICAL_TRANSPORT",
        "sources":{
            "original_nine_history_cargo_sha256":sha256(SOURCE_CARGO),
            "original_239_HOSTS_genus_screen_sha256":sha256(SOURCE_GENUS),
            "original_GBIF_island_first_2011_2017_sha256":sha256(SOURCE_EVENTS),
            "pretest_lock_sha256":sha256(SOURCE_LOCK),
        },
        "original_authentic_VOC_island_links":len(edges),
        "historic_links_with_cargo_rows":len(supported),
        "islands_with_indexed_butterfly_first_record_set":len(first),
        "butterflies_source_event_quality":source_audit,
        "primary_observed":{
            "cargo_supported_links_with_at_least_one_local_compatible_butterfly":sum(
                r["any_local_butterfly_with_same_modern_known_host_genus"] for r in supported),
            "cargo_supported_links_with_compatible_species_at_each_endpoint":sum(
                r["both_islands_have_some_compatible_butterfly"] for r in supported),
            "cargo_supported_links_with_same_compatible_butterfly_indexed_at_both_endpoints":sum(
                r["same_compatible_butterfly_GBIF_first_indexed_in_both_islands"] for r in supported),
        },
        "exact_seven_factorial_commodity_packet_permutation":permutations,
        "number_labeled_packet_assignments":perm_n,
        "derived_island_pair_rows":outputs,
        "limitations":[
            "A VOC commodity category (e.g. rice/indigo) is not taxonomic proof that a live source HOSTS plant species was carried on the particular island leg.",
            "The matched contemporary known larval HOSTS genus does not establish that this butterfly larva used it on the listed island, now or historically.",
            "GBIF 2011–2017 first-indexed species events are not actual colonization or population change, and do not cover 2018–2022 events.",
            "Seven cargo-supported historic links share six island nodes; seven factorial exchangeability is untested across asymmetric shipping routes and geographies.",
            "The exact label-permutation right tail is NOT a valid causal colonial shipping p-value, and no negative route controls were created.",
            "No cargo indexed at 2 historic links is UNKNOWN, not no cargo.",
            "All nine historic positive links and the original 51-island host × butterfly association remain intact.",
        ],
    }
    return audit,outputs


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUTPUT)
    p.add_argument("--pairs",type=Path,default=PAIR_CSV)
    args=p.parse_args()
    receipt,rows=analyze()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",
                           encoding="utf8")
    with args.pairs.open("w",encoding="utf8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "source":receipt["sources"],
        "observed":receipt["primary_observed"],
        "permutation":receipt["exact_seven_factorial_commodity_packet_permutation"],
        "source_status":receipt["status"]},indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
