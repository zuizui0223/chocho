#!/usr/bin/env python3
"""Resolve first indexed years of five prior-host island links from local raw year parquet.

Requires the original work/barriers/global_original_species_island_year.parquet;
the published PR #104 provides only aggregates and 2011-2017 follow-up rows.
This does not infer biological arrival, breeding, ship routing, or dispersal.
"""
from __future__ import annotations
from collections import defaultdict
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=Path("provenance/reviewer_defenses/results")
SCREEN=BASE/"maritime_host_steppingstone_20261030/dated_host_island_chain_screen.json"
FOLLOWUP=BASE/"island_host_continuity_20261030/island_host_continuity_events.csv"
FOLLOWUP_GIT_BLOB="efbc770469a720e2dea82478db2f6196ba03b353"


def species_keys_from_followup(root: Path=ROOT):
    f=root/FOLLOWUP
    data=f.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if sha != FOLLOWUP_GIT_BLOB:
        raise ValueError("Source butterfly key mapping changed")
    mapping=defaultdict(set)
    for row in csv.DictReader(data.decode("utf8").splitlines()):
        mapping[row["species"]].add(str(row["speciesKey"]))
    return {k:next(iter(v)) for k,v in mapping.items() if len(v)==1}


def full_year_resolution(screen: dict, names_to_keys: dict, occurrence_year_rows):
    if screen["schema"]!="chocho_prior_host_island_record_chain_screen_20261030_v1":
        raise ValueError("Expected previously verified, fixed five-link screen")
    links=screen["links"]
    assert len(links)==5
    targets={(names_to_keys[l["butterfly"]],e["island"])
             for l in links for e in l["endpoints"]}
    key_years=defaultdict(set)
    for row in occurrence_year_rows:
        k=(str(row["speciesKey"]),str(row["landmass_id"]))
        if k not in targets:
            continue
        y=int(row["year"])
        if y<1900 or y>2026:
            raise ValueError(f"Unreasonable indexed occurrence year {k}, {y}")
        key_years[k].add(y)

    results=[]
    for l in links:
        species=l["butterfly"]
        if species not in names_to_keys:
            raise ValueError(f"Missing source species key {species}")
        butterfly_key=names_to_keys[species]
        ends=[]
        for prior in l["endpoints"]:
            island=prior["island"]
            years=sorted(key_years.get((butterfly_key,island),[]))
            initial=years[0] if years else None
            if (initial is not None and prior["first_year_in_2011_2017_followup_table"] is not None
                    and initial != prior["first_year_in_2011_2017_followup_table"]):
                raise ValueError(f"Source-year disagreement for {species} on {island}")
            ends.append({
                "island":island,
                "first_indexed_year_from_original_year_table":initial,
                "all_indexed_years_from_original_year_table":years,
                "year_state":"KNOWN_INDEXED_YEAR" if years else "NO_ROW_IN_PROVIDED_YEAR_TABLE",
                "missing_row_is_not_biological_absence":True,
            })
        have=all(e["first_indexed_year_from_original_year_table"] is not None for e in ends)
        ordered=have and ends[0]["first_indexed_year_from_original_year_table"] != ends[1]["first_indexed_year_from_original_year_table"]
        upstream=None
        prior_re_record=None
        years_apart=None
        if ordered:
            earlier,later=sorted(ends,key=lambda x:x["first_indexed_year_from_original_year_table"])
            upstream=earlier["island"]
            t0=earlier["first_indexed_year_from_original_year_table"]
            t1=later["first_indexed_year_from_original_year_table"]
            years_apart=t1-t0
            prior_re_record=any(t0<y<t1 and y<=t0+5
                                for y in earlier["all_indexed_years_from_original_year_table"])
        results.append({
            "butterfly":species,"known_host_before_window":l["catalogue_compatible_prior_host"],
            "island_port_representative_distance_km":l["port_representative_great_circle_km"],
            "endpoints":ends,
            "both_original_first_years_available":have,
            "first_record_order_known_and_distinct":ordered,
            "earlier_indexed_island":upstream,
            "other_island_recorded_after_years":years_apart,
            "upstream_re_recorded_before_later_island_first_record":prior_re_record,
            "not_a_measured_vessel_or_butterfly_route":True,
        })
    return {
        "schema":"chocho_full_original_year_prior_host_link_resolution_v1",
        "input_five_link_screen_schema":screen["schema"],
        "only_five_preselected_hosts_not_a_colonization_test":True,
        "rows":results,
        "fully_resolved_pair_count":sum(x["both_original_first_years_available"] for x in results),
        "distinct_order_and_prior_re_record_count":
            sum(x["upstream_re_recorded_before_later_island_first_record"] is True for x in results),
        "not_biological_reproduction_or_actual_island_relay":True,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-year-parquet",type=Path,required=True)
    parser.add_argument("--source-root",type=Path,default=ROOT)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    f=args.original_year_parquet
    if not f.is_file():
        parser.error(f"Original local species×island×year parquet unavailable: {f}; no date imputation")
    try:
        import pandas as pd
    except ImportError:
        parser.error("pandas with parquet engine is required for original local source")
    data=pd.read_parquet(f,columns=["speciesKey","landmass_id","year"])
    screen=json.loads((args.source_root/SCREEN).read_text(encoding="utf-8"))
    mapping=species_keys_from_followup(args.source_root)
    result=full_year_resolution(screen,mapping,data.to_dict("records"))
    h=hashlib.sha256()
    with f.open("rb") as reader:
        while True:
            b=reader.read(1<<20)
            if not b:
                break
            h.update(b)
    result["local_original_year_parquet_sha256"]=h.hexdigest()
    target=args.output
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "output":str(target),
        "resolved":result["fully_resolved_pair_count"],
        "ordered_with_upstream_re_record":result["distinct_order_and_prior_re_record_count"],
    }))


if __name__=="__main__":
    main()
