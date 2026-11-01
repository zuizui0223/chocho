#!/usr/bin/env python3
"""Evidence-constrained historic VOC cargo-to-original island voyage join.

Uses DOI 10.34894/5HGTCD two original files (57,669 original voyage legs
and 7,737 cargo rows), both source-checksum locked. Cargo descriptions are
free text; matched words are COMMODITY MENTIONS, not identified viable plant
propagules, larval host species or demonstrated colonial introductions.

Source-level IDs link to island pairs only if exactly one archival voyage
row has that REISN key. No source vessel names, individual cargo text, dates
or voyage rows are committed.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter,defaultdict
from pathlib import Path
from derive_dataversenl_17c_VOC_island_routes import (
    load_aliases, source_year, normalize_port, PORT_FROM, PORT_TO,
    DATE_FROM, DATE_TO, SOURCE_KEY, EXPECTED_SHA1, EXPECTED_ROWS)

BASE=Path(__file__).resolve().parent
CARGO=Path("_private_history/voc_7737_cargo.tab")
VOYAGES=Path("_private_history/voc_57669_all_voyages.tab")
LOCK=BASE/"voc_17c_cargo_commodity_lock_20261101.json"
ROUTES=BASE/"voc_57669_17c_interisland_routes_v1.audit.json"
OUT=BASE/"voc_17c_cargo_on_verified_island_routes_v1.json"
EXPECTED_CARGO_SHA1="d7cc3aff064fb63f58b8b549b1ca9808cb4b76d1"


def sha(path,kind="sha1"):
    h=hashlib.new(kind)
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def source_rows(path):
    with path.open("r",encoding="utf-8-sig",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")


def classify(text,categories):
    t=str(text or "").casefold()
    return [k for k,words in categories.items() if any(w in t for w in words)]


def aggregate(voyages=VOYAGES,cargo=CARGO,lock=LOCK,reference=ROUTES):
    if sha(voyages)!=EXPECTED_SHA1 or sha(cargo)!=EXPECTED_CARGO_SHA1:
        raise RuntimeError("Original verified VOC voyage/cargo source checksum mismatch")
    locked=json.loads(lock.read_text(encoding="utf8"))
    if locked["schema"]!="chocho_17c_VOC_archived_cargo_commodity_category_lock_v1":
        raise RuntimeError("Original commodity classification lock missing")
    types=locked["categories"]
    ref=json.loads(reference.read_text(encoding="utf8"))
    if (ref["historical_17c_positive_source_island_pairs"]!=9 or
        ref["historic_direct_island_pair_voyage_leg_rows_total"]!=1193):
        raise RuntimeError("Previously audited 9 historic island links drifted")
    alias,_,_=load_aliases(include_broad=False)
    route_by_key=defaultdict(set)
    voyage_id_row_count=Counter()
    observed_edges=Counter()
    source_chronology=Counter()
    n=0
    for r in source_rows(voyages):
        n+=1
        key=str(r.get(SOURCE_KEY) or "").strip()
        if key: voyage_id_row_count[key]+=1
        a=alias.get(normalize_port(r.get(PORT_FROM)))
        b=alias.get(normalize_port(r.get(PORT_TO)))
        if not a or not b or a==b:
            continue
        start=source_year(r.get(DATE_FROM))
        finish=source_year(r.get(DATE_TO))
        if start is None:
            source_chronology["mapped_two_island_departure_year_unknown"]+=1
            continue
        if finish is None:
            source_chronology["mapped_two_island_arrival_year_unknown"]+=1
        elif finish<start:
            source_chronology["mapped_two_island_arrival_year_before_departure"]+=1
        else:
            source_chronology["mapped_two_island_year_order_nonnegative"]+=1
        edge=tuple(sorted((a,b),key=int))
        observed_edges[edge]+=1
        if key: route_by_key[key].add(edge)
    if n!=EXPECTED_ROWS or len(observed_edges)!=9 or sum(observed_edges.values())!=1193:
        raise RuntimeError("Historical complete voyage source mismatches earlier 9/1193 receipt")
    # Full source IDs are unique for most legs. Ambiguous IDs are completely
    # excluded rather than re-used across different legs.
    safely_mapped={
        k:next(iter(edges)) for k,edges in route_by_key.items()
        if len(edges)==1 and voyage_id_row_count[k]==1
    }
    mapping_other=Counter()
    for k,edges in route_by_key.items():
        if len(edges)>1: mapping_other["mapped_id_multiple_island_pairs"]+=1
        elif voyage_id_row_count[k]>1: mapping_other["mapped_id_reused_in_source"]+=1
    seen_id=set()
    raw_match=Counter()
    source_categories=Counter()
    matched_categories=Counter()
    island_categories=defaultdict(Counter)
    island_rows=Counter()
    match_exact_ids=set()
    plant_cargo_ids=set()
    joined_rows_with_text=0
    total_cargo=0
    cargo_column_nonblank=0
    ambiguous_cargo_id=0
    for r in source_rows(cargo):
        total_cargo+=1
        k=str(r.get(SOURCE_KEY) or "").strip()
        text=str(r.get(locked["source_cargo_text_column"]) or "").strip()
        if text: cargo_column_nonblank+=1
        tags=classify(text,types)
        for cat in tags: source_categories[cat]+=1
        if k in route_by_key and k not in safely_mapped:
            ambiguous_cargo_id+=1
        if k not in safely_mapped:
            continue
        match_exact_ids.add(k)
        edge=safely_mapped[k]
        island_rows[edge]+=1
        if text: joined_rows_with_text+=1
        for cat in tags:
            matched_categories[cat]+=1
            island_categories[edge][cat]+=1
            plant_cargo_ids.add(k)
    if total_cargo!=7737:
        raise RuntimeError("Historical cargo rows not fully scanned")
    rows=[]
    for (a,b),nlegs in sorted(observed_edges.items(),key=lambda x:(int(x[0][0]),int(x[0][1]))):
        row={
            "island_a":a,"island_b":b,
            "observed_original_voyage_leg_rows":nlegs,
            "source_cargo_rows_joined_by_unrepeated_REISN":island_rows[(a,b)],
            "island_pair_cargo_mentions_in_predeclared_plant_commodity_categories":
                sum(island_categories[(a,b)].values()),
            "source_cargo_category_counts":dict(sorted(island_categories[(a,b)].items())),
        }
        rows.append(row)
    return {
        "schema":"chocho_17c_VOC_exact_cargo_7737_to_9_historic_island_routes_v1",
        "status":"SOURCE_VERIFIED_HISTORIC_VOYAGE_CARGO_KEYWORDS_NOT_LIVE_PLANT_OR_BUTTERFLY_TRANSPORT",
        "source":{
            "DOI":"10.34894/5HGTCD",
            "voyage_original_SHA1":sha(voyages),
            "cargo_original_SHA1":sha(cargo),
            "prelocked_commodity_lexicon_SHA256":sha(lock,"sha256"),
            "source_archived_historic_island_pair_receipt_SHA256":sha(reference,"sha256"),
        },
        "historical_source_original_voyage_rows":n,
        "historical_source_original_cargo_rows":total_cargo,
        "historical_original_cargo_descriptions_present":cargo_column_nonblank,
        "verified_historic_island_edges":len(observed_edges),
        "verified_historic_island_voyage_leg_rows":sum(observed_edges.values()),
        "unique_original_2island_voyage_ids":len(route_by_key),
        "source_route_id_quality":dict(mapping_other),
        "source_cargo_rows_linked_to_unambiguous_two_island_leg":sum(island_rows.values()),
        "unique_cargo_voyage_ids_linked_to_unambiguous_two_island_leg":len(match_exact_ids),
        "linked_cargo_rows_with_descriptions":joined_rows_with_text,
        "linked_source_cargo_ids_containing_any_declared_plant_commodity":len(plant_cargo_ids),
        "cargo_rows_from_ambiguous_or_repeated_island_route_voyage_ids":ambiguous_cargo_id,
        "all_original_cargo_category_keyword_row_counts":dict(sorted(source_categories.items())),
        "on_historic_island_routes_cargo_category_keyword_row_counts":dict(sorted(matched_categories.items())),
        "historic_leg_original_year_order_quality":dict(source_chronology),
        "pair_specific_aggregates":rows,
        "limitations":locked["classification_downgrade"]+[
            "Source REISN is a voyage-level archival key, but commodity details may apply to a ship's whole itinerary, not the specific island-to-island movement.",
            "Voyage IDs reused across rows or ambiguously mapped are excluded from the strict joined cargo evidence.",
            "The fixed plant commodity keyword list is an incomplete lexical screen, not a validated species-level botanical dataset.",
            "No positive cargo keyword result can itself prove viable plant propagules, larval host availability or actual butterfly transport.",
            "Other historic 17th-century shipping company and cargo records are not covered.",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--voyages",type=Path,default=VOYAGES)
    ap.add_argument("--cargo",type=Path,default=CARGO)
    ap.add_argument("--output",type=Path,default=OUT)
    args=ap.parse_args()
    r=aggregate(args.voyages,args.cargo)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in r.items()
         if k not in ("limitations","pair_specific_aggregates")},indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
