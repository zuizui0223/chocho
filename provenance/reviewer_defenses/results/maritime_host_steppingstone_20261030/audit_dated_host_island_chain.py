#!/usr/bin/env python3
"""Triage previously dated shared-host island pairs for possible time-ordering.

The archived butterfly follow-up table is limited to *2011-2017 index events*.
Absence from it does NOT prove later first arrival or biological nondetection.
Port-distance is a port-representative geodesic, not true dispersal or ship route.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=Path("provenance/reviewer_defenses/results")
FILES={
    "links": ("island_dated_host_match_20261030/island_dated_host_match_links.csv",
              "88db773c74dea9cd2d58bde0b67a394d93226869"),
    "pairs": ("island_host_network_20261030/island_host_network_pairs.csv",
              "5663d2ffe42419f5c4496b92c9f3bbf76175907b"),
    "events": ("island_host_continuity_20261030/island_host_continuity_events.csv",
               "efbc770469a720e2dea82478db2f6196ba03b353")
}


def inputs(root: Path, key: str):
    p,blob=FILES[key]
    data=(root/BASE/p).read_bytes()
    git_blob=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if git_blob!=blob:
        raise ValueError(f"Input not pinned: {key}")
    return list(csv.DictReader(io.StringIO(data.decode("utf-8"))))


def canonical_pair(a: str,b: str):
    if not a.isdecimal() or not b.isdecimal() or a==b:
        raise ValueError("Invalid distinct island identifiers")
    return tuple(sorted((a,b),key=int))


def year_observation(row: dict|None, island: str):
    if row is None:
        return {
            "island":island,
            "first_year_in_2011_2017_followup_table":None,
            "followup_years_in_first_five_years":None,
            "record_in_followup_year5":None,
            "first_year_not_in_followup_table_is_not_evidence_of_absence":True,
        }
    first=int(row["first_record_year"])
    if first<2011 or first>2017:
        raise ValueError("Unexpected initial window")
    years=([int(x) for x in row["followup_year_list"].split(";")]
           if row["followup_year_list"] else [])
    if any(y<=first or y>first+5 for y in years) or len(set(years))!=len(years):
        raise ValueError("Invalid followup years")
    if bool(years)!=(row["any_followup_record"]=="True"):
        raise ValueError("Followup flag mismatch")
    return {
        "island":island,
        "first_year_in_2011_2017_followup_table":first,
        "followup_years_in_first_five_years":years,
        "record_in_followup_year5":row["record_in_year5"]=="True",
        "first_year_not_in_followup_table_is_not_evidence_of_absence":False,
    }


def audit(root: Path = ROOT):
    links=inputs(root,"links");pairs=inputs(root,"pairs");events=inputs(root,"events")
    if (len(links),len(pairs),len(events))!=(5,1275,5665):
        raise ValueError("Source cohort has drifted")
    pm={}
    for p in pairs:
        k=canonical_pair(p["island_a"],p["island_b"])
        if k in pm:
            raise ValueError("Duplicate island pair")
        pm[k]=p
    em={}
    for e in events:
        k=(e["species"],e["landmass_id"])
        if k in em:
            raise ValueError("Duplicate butterfly island-year record")
        em[k]=e
    rows=[]
    for link in links:
        a,b=canonical_pair(link["island_a"],link["island_b"])
        pa=pm.get((a,b))
        if not pa:
            raise ValueError("Early host evidence pair missing from fixed island network")
        if int(pa["shared_new_species"])<1:
            raise ValueError("Dated host link outside shared-new butterfly risk set")
        end_a=year_observation(em.get((link["species"],a)),a)
        end_b=year_observation(em.get((link["species"],b)),b)
        have_both=(end_a["first_year_in_2011_2017_followup_table"] is not None and
                   end_b["first_year_in_2011_2017_followup_table"] is not None)
        ordered=False
        upstream=None
        follows_before_second=None
        if have_both:
            fa=end_a["first_year_in_2011_2017_followup_table"]
            fb=end_b["first_year_in_2011_2017_followup_table"]
            ordered=fa!=fb
            if ordered:
                u,v=(end_a,end_b) if fa<fb else (end_b,end_a)
                upstream=u["island"]
                follows_before_second=any(y<v["first_year_in_2011_2017_followup_table"]
                                          for y in u["followup_years_in_first_five_years"])
        d=float(pa["port_distance_km"])
        if not math.isfinite(d) or d<0:
            raise ValueError("Nonfinite port representative distance")
        rows.append({
            "butterfly":link["species"],
            "catalogue_compatible_prior_host":link["host_name"],
            "both_islands_have_2001_2010_port_host_records":True,
            "island_a":a, "island_b":b,
            "port_representative_great_circle_km":d,
            "not_a_measured_sea_crossing_or_actual_route":True,
            "endpoints":[end_a,end_b],
            "both_initial_record_years_available_in_followup_table":have_both,
            "ordered_first_record_years_in_available_2011_2017_table":ordered,
            "earlier_record_island_if_both_verified":upstream,
            "upstream_has_re_record_before_other_island_if_both_verified":follows_before_second,
            "unobserved_or_unknown_endpoint_year_not_inferred":True,
        })
    if len({r["butterfly"] for r in rows})!=5 or len({canonical_pair(r["island_a"],r["island_b"]) for r in rows})!=4:
        raise ValueError("Focal link identity changed")
    if sum(r["both_initial_record_years_available_in_followup_table"] for r in rows)!=1:
        raise ValueError("Temporal linkage support changed")
    if sum(r["port_representative_great_circle_km"]<500 for r in rows)!=3:
        raise ValueError("Small-distance link support changed")
    return {
        "schema":"chocho_prior_host_island_record_chain_screen_20261030_v1",
        "source_git_blob_pins":{k:blob for k,(_,blob) in FILES.items()},
        "status":"EXPLORATORY_CHRONOLOGY_TRIAGE_NOT_RELAY_EVIDENCE",
        "links":rows,
        "butterfly_species":5,
        "unique_island_pairs":4,
        "shorter_than_500km_port_representative_pairs_counting_butterfly":3,
        "both_record_years_2011_2017_available":1,
        "ordered_and_upstream_re_recorded_before_second":sum(
            r["upstream_has_re_record_before_other_island_if_both_verified"] is True for r in rows),
        "none_of_these_is_evidence_of_actual_ship_transport_or_island_to_island_dispersal":True,
        "null_event_not_proof_of_date_2018_2022":True,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    s=json.dumps(audit(args.root),indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        dst=args.output if args.output.is_absolute() else args.root/args.output
        dst.parent.mkdir(parents=True,exist_ok=True)
        dst.write_text(s,encoding="utf-8")
        print(dst)
    else:
        print(s)


if __name__=="__main__":
    main()
