#!/usr/bin/env python3
"""Strict importer for the EXACT privately supplied original year Parquet.

The original bytes are NOT in this Git repository, and are NOT uploaded by this
script. The CLI requires an authorized local path and verifies the known
SHA-256, data schema, 83,573 original source rows and all original 225x51
index margins before producing a DERIVED index-only temporal risk table.

Source non-indexing is NEVER interpreted as ecological absence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_SOURCE_SHA256="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
EXPECTED_COLUMNS={
    "speciesKey","landmass_id","year","records",
    "precise1km_records","landmass_area_km2",
}
EXPORT_FIELDS=[
    "species","taxon_key","landmass_id","original_first_indexed_year",
    "original_recorded_year_count_2001_2024",
    "pre_2011_first_index_state",
    "new_first_index_2011_2017_state",
    "new_first_index_2018_2022_state",
    "post_2022_first_index_state",
    "not_indexed_2001_2024_state",
    "SOURCE_OUTCOME_NOT_BIOLOGICAL_PRESENCE",
]


def sha256_file(path:Path)->str:
    sha=hashlib.sha256()
    with path.open("rb") as handle:
        for buf in iter(lambda:handle.read(1024*1024),b""):
            sha.update(buf)
    return sha.hexdigest()


def source_spec_check(year_records:list[dict], *, full_original:bool)->dict:
    if not year_records:
        raise ValueError("Empty original-year source")
    seen=set()
    species=set()
    islands=set()
    years=set()
    record_sum=0
    for r in year_records:
        if set(r)!=EXPECTED_COLUMNS:
            raise ValueError("Unexpected original year source columns")
        try:
            k=int(r["speciesKey"]); lid=int(r["landmass_id"]); y=int(r["year"])
            count=int(r["records"]); precise=int(r["precise1km_records"])
            area=float(r["landmass_area_km2"])
        except (ValueError,TypeError,OverflowError) as e:
            raise ValueError("Non-numeric original source ID/record/area field") from e
        if (str(k)!=str(r["speciesKey"]) or
            str(lid)!=str(r["landmass_id"]) or
            str(y)!=str(r["year"])):
            # Supports pandas numeric values, while rejecting silent IDs
            # containing suffixes, leading-zero changes or group symbols.
            if not (isinstance(r["speciesKey"],int)
                    and isinstance(r["landmass_id"],int)
                    and isinstance(r["year"],int)):
                raise ValueError("Original source ID/year key is not canonical")
        if not (k>0 and lid>0 and 2001<=y<=2024 and count>0 and
                precise==count and math.isfinite(area) and area>0):
            raise ValueError("Original source row violates year, count, precision or area")
        key=(k,lid,y)
        if key in seen:
            raise ValueError("Duplicate species×island×year original record")
        seen.add(key);species.add(k);islands.add(lid);years.add(y)
        record_sum+=count
    stats={
        "source_rows":len(year_records),"species":len(species),
        "island_groups":len(islands),"year_min":min(years),
        "year_max":max(years),"occurrence_records_sum":record_sum,
        "all_unique":True,
    }
    if full_original and stats!={"source_rows":83573,
                                  "species":224,
                                  "island_groups":3384,
                                  "year_min":2001,
                                  "year_max":2024,
                                  "occurrence_records_sum":14425536,
                                  "all_unique":True}:
        raise ValueError("Parquet matches its declared SHA but original row count/identity statistics changed")
    return stats


def first_year_index(year_records:list[dict])->tuple[dict[tuple[str,str],int],dict[tuple[str,str],int]]:
    first={}
    count=Counter()
    for row in year_records:
        key=(str(int(row["speciesKey"])),str(int(row["landmass_id"])))
        y=int(row["year"])
        count[key]+=1
        if key not in first or y<first[key]:
            first[key]=y
    return first,dict(count)


def build_original_visibility(year_records:list[dict], species:list[dict], islands:list[dict],
                               positive_2011_2017:dict[tuple[str,str],int],
                               selected_late_4:dict[tuple[str,str],int],
                               *, require_original_population:bool=True)->tuple[list[dict],dict]:
    stats=source_spec_check(year_records,full_original=require_original_population)
    first,year_counts=first_year_index(year_records)
    species_by_key={}
    for sp in species:
        name=sp["species"]
        key=str(sp["taxon_key"])
        if key in species_by_key:
            raise ValueError("Duplicate accepted butterfly key in source roster")
        species_by_key[key]=name
    island_ids={str(z["landmass_id"]) for z in islands}
    if len(island_ids)!=len(islands):
        raise ValueError("Duplicate original landmass ID")
    rows=[]
    pre=Counter();new1=Counter();new2=Counter()
    first_after_2022=0
    by_sp_pre=Counter();by_sp_new=Counter();by_is_pre=Counter();by_is_new=Counter()
    for k,name in sorted(species_by_key.items(),key=lambda z:z[1]):
        for island in sorted(island_ids,key=int):
            combo=(k,island)
            year=first.get(combo)
            old=int(year is not None and year<=2010)
            y11=int(year is not None and 2011<=year<=2017)
            y18=int(year is not None and 2018<=year<=2022)
            ypost=int(year is not None and year>=2023)
            if old:by_sp_pre[name]+=1;by_is_pre[island]+=1
            if y11 or y18:by_sp_new[name]+=1;by_is_new[island]+=1
            pre[old]+=1;new1[y11]+=1;new2[y18]+=1;first_after_2022+=ypost
            rows.append({
                "species":name,"taxon_key":k,"landmass_id":island,
                "original_first_indexed_year":year if year is not None else "",
                "original_recorded_year_count_2001_2024":year_counts.get(combo,0),
                "pre_2011_first_index_state":"INDEXED_IN_SOURCE" if old else "NOT_PRE2011_INDEXED_NOT_ECOLOGICAL_ABSENCE",
                "new_first_index_2011_2017_state":"FIRST_INDEXED" if y11 else "NOT_FIRST_INDEXED_IN_THIS_WINDOW",
                "new_first_index_2018_2022_state":"FIRST_INDEXED" if y18 else "NOT_FIRST_INDEXED_IN_THIS_WINDOW",
                "post_2022_first_index_state":"FIRST_INDEXED" if ypost else "NOT_FIRST_INDEXED_IN_THIS_WINDOW",
                "not_indexed_2001_2024_state":"NO_SOURCE_INDEX_NOT_ECOLOGICAL_ABSENCE" if year is None else "INDEXED_IN_SOURCE",
                "SOURCE_OUTCOME_NOT_BIOLOGICAL_PRESENCE":"TRUE",
            })
            if (name,island) in positive_2011_2017:
                if year!=positive_2011_2017[(name,island)]:
                    raise ValueError("Original Parquet contradicts committed positive 2011–2017 source year")
            if (name,island) in selected_late_4:
                if year!=selected_late_4[(name,island)]:
                    raise ValueError("Original Parquet contradicts independently verified four 2018–2022 source years")
    row_total=len(species_by_key)*len(island_ids)
    if require_original_population:
        if (len(species)!=225 or len(islands)!=51 or row_total!=11475
                or pre[1]!=524 or new1[1]!=302 or new2[1]!=128):
            raise ValueError("Original butterfly×island source margins do not match verified records")
        for sp in species:
            if (by_sp_pre[sp["species"]]!=int(sp["early_indexed_islands_2001_2010"]) or
                by_sp_new[sp["species"]]!=int(sp["new_indexed_islands_2011_2022"])):
                raise ValueError("Original butterfly marginal source contradicts Parquet")
        for isle in islands:
            id=str(isle["landmass_id"])
            if (by_is_pre[id]!=int(isle["early_butterfly_species_recorded_2001_2010"]) or
                by_is_new[id]!=int(isle["new_butterfly_species_indexed_2011_2022"])):
                raise ValueError("Original island marginal source contradicts Parquet")
    support={
        "schema":"chocho_NC_exact_original_parquet_strict_recovery_v1",
        "status":"DERIVED_ORIGINAL_GBIF_INDEX_TABLE_NOT_ECOLOGICAL_OCCUPANCY",
        "source_statistics":stats,
        "original_species":len(species_by_key),
        "original_islands":len(island_ids),
        "source_grain":"species×island×first indexed year among 2001–2024 occurrence-index records",
        "species_island_combinations":row_total,
        "source_2001_2010_indexed":pre[1],
        "source_2011_2017_first_indexed":new1[1],
        "source_2018_2022_first_indexed":new2[1],
        "source_first_index_post_2022":first_after_2022,
        "source_never_indexed_2001_2024_within_225x51":sum(z["original_first_indexed_year"]=="" for z in rows),
        "original_2011_2017_known_source_years_verified":len(positive_2011_2017),
        "original_selected_2018_2022_known_source_years_verified":len(selected_late_4),
        "source_nonindex_is_not_confirmed_butterfly_absence":True,
        "source_future_host_dataset_NOT_integrated":True,
        "requires_detection_effort_for_record_event_model":True,
    }
    return rows,support


def read_visibility_inputs():
    spec=importlib.util.spec_from_file_location("ncomms_risk_source",ROOT/"scripts/build_ncomms_island_risk_visibility.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Source visibility module not present")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    original,orig_pins=module.load_original()
    v,ii,ss,audit=module.build_visibility(original,orig_pins)
    early={(z["species"],z["landmass_id"]):int(z["first_index_year_IF_2011_2017_case"])
           for z in v if z["butterfly_first_index_2011_2017_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX"}
    chain=ROOT/"provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/full_original_island_chain_years_observed_20261030.json"
    obj=json.loads(chain.read_text())
    selected={}
    for link in obj["rows"]:
        for endpoint in link["endpoints"]:
            y=int(endpoint["first_indexed_year_from_original_year_table"])
            if y>=2018 and y<=2022:
                selected[(link["butterfly"],str(endpoint["island"]))]=y
    if len(early)!=302 or len(selected)!=4:
        raise ValueError("Original indexed source identities have changed")
    return ii,ss,early,selected


def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--original-year-parquet",type=Path,required=True,
                     help="Existing AUTHORISED local original byte path; never uploads raw data.")
    arg.add_argument("--output-dir",type=Path,required=True)
    args=arg.parse_args()
    source=args.original_year_parquet
    if not source.is_file():
        raise FileNotFoundError("Authorized original file is NOT currently mounted: "+str(source))
    h=sha256_file(source)
    if h!=ORIGINAL_SOURCE_SHA256:
        raise ValueError("Original source bytes have an unexpected SHA256; no conversion attempted")
    try:
        import pandas as pd
    except ImportError as e:
        raise RuntimeError("Install pandas and pyarrow to read the Parquet; do not change data source") from e
    data=pd.read_parquet(source)
    if set(data.columns)!=EXPECTED_COLUMNS:
        raise ValueError("Original Parquet columns have changed")
    records=data.to_dict(orient="records")
    islands,species,early,late=read_visibility_inputs()
    out,report=build_original_visibility(records,species,islands,early,late)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    p=args.output_dir/"original225x51_exact_index_status_REQUIRES_EFFORT_FOR_MODEL.csv"
    with p.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,EXPORT_FIELDS,lineterminator="\n");w.writeheader();w.writerows(out)
    report["verified_raw_parquet_sha256"]=h
    report["export_csv_sha256"]=hashlib.sha256(p.read_bytes()).hexdigest()
    (args.output_dir/"strict_original_parquet_reconstruction_receipt.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":report["status"],"source_sha256_verified":True,
                      "cells":report["species_island_combinations"],
                      "export":str(p)},indent=2))


if __name__=="__main__":
    main()
