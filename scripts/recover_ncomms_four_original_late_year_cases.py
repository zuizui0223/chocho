#!/usr/bin/env python3
"""Recover four original-year 2018–2022 butterfly index cases from verified five-link source.

The originally supplied 83,573-row Parquet is NOT in tracked Git. Rather than
invent its contents, recover exactly those endpoint-year observations already
independently verified and archived by the original scientific analysis.

This is a partial source-positive overlay only; no regression or ecological
absence claims, and no host effect from selection on five prior-host links.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CHAIN=Path("provenance/reviewer_defenses/results/maritime_host_steppingstone_20261030/full_original_island_chain_years_observed_20261030.json")
CHAIN_GIT_SHA1="68de667c1f93c4e7d99f25ffba54b341d7f1f697"
PARENT=ROOT/"scripts/build_ncomms_island_risk_visibility.py"
SOURCE_OLD_PRIVATE_SHA256="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"

EXTRA_COLS=[
  "original_year_first_index_from_SELECTED_five_links",
  "source_named_prior_host_from_SELECTED_five_links",
  "original_year_evidence_source",
]
LATE_COLS=[
  "butterfly", "landmass_id", "first_index_year", "selected_prior_host",
  "previous_2018_2022_event_identity_state", "verified_source", 
]
SOURCE_URI=CHAIN.as_posix()


def load_builder():
    spec=importlib.util.spec_from_file_location("ncomms_visibility_for_chain",PARENT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Original source visibility builder is unavailable")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_chain(root:Path=ROOT)->tuple[dict,dict]:
    raw=(root/CHAIN).read_bytes()
    sha1=hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+b"\0"+raw).hexdigest()
    if sha1!=CHAIN_GIT_SHA1:
        raise ValueError("Original source-chain Git blob SHA changed")
    obj=json.loads(raw)
    if obj.get("schema")!="chocho_full_original_year_prior_host_link_resolution_v1":
        raise ValueError("Original chronology schema changed")
    if obj.get("local_original_year_parquet_sha256")!=SOURCE_OLD_PRIVATE_SHA256:
        raise ValueError("Original raw Parquet identity changed")
    if not obj.get("only_five_preselected_hosts_not_a_colonization_test"):
        raise ValueError("Archive's five-link ascertainment disclaimer is missing")
    if not obj.get("not_biological_reproduction_or_actual_island_relay"):
        raise ValueError("Original source inference guard has been altered")
    if len(obj.get("rows",[]))!=5 or obj.get("fully_resolved_pair_count")!=5:
        raise ValueError("Original explicitly preselected five links changed")
    return obj,{"repo_path":SOURCE_URI, "git_blob_sha1":sha1, 
                "raw_sha256":hashlib.sha256(raw).hexdigest(),
                "upstream_local_original_parquet_sha256_NOT_CURRENTLY_AVAILABLE":SOURCE_OLD_PRIVATE_SHA256}


def recover_visibility(original_rows:list[dict],species_rows:list[dict],
                       island_rows:list[dict],original_audit:dict,chain:dict,
                       source_receipt:dict)->tuple[list[dict],list[dict],dict]:
    pairs={(x["species"],x["landmass_id"]):x for x in original_rows}
    if len(pairs)!=225*51 or len(pairs)!=len(original_rows):
        raise ValueError("Source 225×51 grid is incomplete or duplicates a pair")
    spp={x["species"]:x for x in species_rows}
    isl={x["landmass_id"]:x for x in island_rows}
    original_positive_2011_2017=sum(
      z["butterfly_first_index_2011_2017_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX"
      for z in original_rows
    )
    if original_positive_2011_2017!=302:
        raise ValueError("Original exact early positive record population changed")
    for r in original_rows:
        r.update({
            "original_year_first_index_from_SELECTED_five_links":"",
            "source_named_prior_host_from_SELECTED_five_links":"",
            "original_year_evidence_source":"",
        })
    selected={}
    late=[]
    first_early=Counter()
    first_late=Counter()
    for link in chain["rows"]:
        name=link["butterfly"]
        prior=link["known_host_before_window"]
        if not prior or not isinstance(prior,str):
            raise ValueError("Selected historical host identity was lost")
        if len(link["endpoints"])!=2 or not link["both_original_first_years_available"]:
            raise ValueError("Unverified source event in five-link chronology")
        if not link["first_record_order_known_and_distinct"]:
            raise ValueError("Source link chronology changed")
        for end in link["endpoints"]:
            id=str(end["island"])
            k=(name,id)
            if k in selected or k not in pairs:
                raise ValueError("Repeated or foreign original species×island endpoint: "+str(k))
            if name not in spp or id not in isl:
                raise ValueError("Selected endpoint not in fixed source margins")
            year=int(end["first_indexed_year_from_original_year_table"])
            yrs=end["all_indexed_years_from_original_year_table"]
            if not yrs or sorted(set(yrs))!=yrs or yrs[0]!=year or year < 2011 or year > 2022:
                raise ValueError("Invalid saved source first year and indexed followup chronology")
            if end["year_state"]!="KNOWN_INDEXED_YEAR" or not end["missing_row_is_not_biological_absence"]:
                raise ValueError("Original year evidence status is not independently validated")
            row=pairs[k]
            row["original_year_first_index_from_SELECTED_five_links"]=str(year)
            row["source_named_prior_host_from_SELECTED_five_links"]=prior
            row["original_year_evidence_source"]=SOURCE_URI
            if year<=2017:
                if row["first_index_year_IF_2011_2017_case"]!=str(year):
                    raise ValueError("Selected original 2011–2017 source year does not match earlier CSV")
                if row["butterfly_first_index_2011_2017_at_THIS_species_island"]!="CONFIRMED_FIRST_INDEX":
                    raise ValueError("Earlier positive butterfly case was overwritten")
                first_early[k]+=1
            else:
                if row["first_index_year_IF_2011_2017_case"]:
                    raise ValueError("Late original first index was already an early original case")
                if row["butterfly_first_index_2011_2017_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX":
                    raise ValueError("Late original first index conflicts with archived case")
                original_state=row["full_index_2011_2022_at_THIS_species_island"]
                if not original_state.startswith("NOT_RECOVERED_"):
                    raise ValueError("Late original first case was already known; no new recovery")
                row["butterfly_pre2011_index_at_THIS_species_island"]="KNOWN_NOT_PRE2011"
                row["butterfly_first_index_2011_2017_at_THIS_species_island"]="KNOWN_NOT_FIRST_2011_2017"
                row["first_index_2018_2022_at_THIS_species_island"]="CONFIRMED_FIRST_INDEX_ORIGINAL_YEAR_SOURCE"
                row["full_index_2011_2022_at_THIS_species_island"]="CONFIRMED_FIRST_INDEX_2018_2022_ORIGINAL_YEAR_SOURCE"
                late.append({
                    "butterfly":name, "landmass_id":id, "first_index_year":year,
                    "selected_prior_host":prior,
                    "previous_2018_2022_event_identity_state":original_state,
                    "verified_source":SOURCE_URI,
                })
                first_late[k]+=1
            selected[k]=year
    if len(selected)!=10 or sum(first_early.values())!=6 or sum(first_late.values())!=4:
        raise ValueError("Original five-link selected source must yield six known and four NEW late years")
    if len({(z["butterfly"],z["landmass_id"]) for z in late})!=4:
        raise ValueError("Selected late positives contain duplicate ids")
    if {(z["butterfly"],z["landmass_id"],z["first_index_year"]) for z in late}!={
        ("Colias croceus","169",2019),
        ("Colias hyale","169",2019),
        ("Pieris rapae","951",2019),
        ("Nymphalis antiopa","715",2021)
    }:
        raise ValueError("Recovered case identities do not match original validated endpoint screen")
    for k in selected:
        if k not in pairs:
            raise AssertionError("Unreachable selected endpoint")
    late_species=Counter(z["butterfly"] for z in late)
    late_island=Counter(z["landmass_id"] for z in late)
    new_species_unassigned={}
    for z in species_rows:
        rest=int(z["unassigned_new_islands_2018_2022"])-late_species[z["species"]]
        if rest<0:
            raise ValueError("Newly recovered late case exceeds original species margin")
        new_species_unassigned[z["species"]]=rest
    new_island_unassigned={}
    for z in island_rows:
        rest=int(z["unassigned_new_species_2018_2022"])-late_island[z["landmass_id"]]
        if rest<0:
            raise ValueError("Newly recovered late case exceeds original island margin")
        new_island_unassigned[z["landmass_id"]]=rest
    if sum(new_species_unassigned.values())!=124 or sum(new_island_unassigned.values())!=124:
        raise ValueError("Recovered four source cases do not reconcile original 128 unknown late identities")
    known_late=sum(z["first_index_2018_2022_at_THIS_species_island"]=="CONFIRMED_FIRST_INDEX_ORIGINAL_YEAR_SOURCE" for z in original_rows)
    if known_late!=4:
        raise ValueError("Late source overlay did not preserve four and only four unique cases")
    report={
       "schema":"chocho_NC_verified_four_original_year_2018_2022_cases_v1",
       "status":"FOUR_REAL_ORIGINAL_YEAR_POSITIVE_CASES_RECOVERED_FROM_ARCHIVED_SELECTED_SOURCE",
       "source_recording_period":"2001–2024 original record index; late positives first seen 2018–2022",
       "source_original_private_parquet_NOT_available_for_full_rerun":True,
       "original_five_selected_chains_git_receipt":source_receipt,
       "input_visibility_grid_original_positive_2011_2017":original_positive_2011_2017,
       "prior_known_2018_2022_case_IDs":0,
       "new_original_source_year_2018_2022_case_IDs":len(late),
       "all_year_resolved_original_source_endpoints":len(selected),
       "previous_2018_2022_positive_identity_deficit":128,
       "remaining_2018_2022_positive_identity_deficit":124,
       "source_225_candidate_species_and_51_island_margins_preserved":True,
       "recovered_2018_2022_events":sorted(late,key=lambda z:(z["first_index_year"],z["butterfly"])),
       "source_rows_not_identified_rowwise_AS_negative_outcomes":True,
       "selection_bias":"Only preselected five prior-host links were available; recovered cases cannot be used to estimate resource effects",
       "source_named_prior_host_is_NOT_a_valid_complete_population_exposure":True,
       "no_ship_movement_or_butterfly_colonisation_established":True,
       "remaining_unassigned_by_species_positive_sum":sum(new_species_unassigned.values()),
       "remaining_unassigned_by_island_positive_sum":sum(new_island_unassigned.values()),
    }
    return original_rows,late,report


def save_outputs(out:Path,grid:list[dict],new:list[dict],receipt:dict,builder):
    out.mkdir(parents=True,exist_ok=True)
    fgrid=out/"species_island_225x51_with_4_original_late_positive_years.csv"
    flate=out/"original_four_2018_2022_first_index_cases.csv"
    with fgrid.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=builder.MATRIX_FIELDS+EXTRA_COLS,lineterminator="\n")
        w.writeheader();w.writerows(grid)
    with flate.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=LATE_COLS,lineterminator="\n")
        w.writeheader();w.writerows(new)
    receipt["export_sha256"]={
        fgrid.name:hashlib.sha256(fgrid.read_bytes()).hexdigest(),
        flate.name:hashlib.sha256(flate.read_bytes()).hexdigest(),
    }
    (out/"original_four_late_source_audit.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8"
    )
    return {"grid":str(fgrid),"late":str(flate),"audit":str(out/"original_four_late_source_audit.json")}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    args=parser.parse_args()
    builder=load_builder()
    source,pins=builder.load_original()
    grid,island_rows,species_rows,old_audit=builder.build_visibility(source,pins)
    chain,receipt=load_chain()
    updated,new,checks=recover_visibility(grid,species_rows,island_rows,old_audit,chain,receipt)
    paths=save_outputs(args.out_dir,updated,new,checks,builder)
    print(json.dumps({
        "status":checks["status"],
        "recovered_additional_year_identified_cases":len(new),
        "remaining_late_positive_cases_with_unknown_identity":checks["remaining_2018_2022_positive_identity_deficit"],
        "files":paths
    },indent=2))


if __name__=="__main__":
    main()
