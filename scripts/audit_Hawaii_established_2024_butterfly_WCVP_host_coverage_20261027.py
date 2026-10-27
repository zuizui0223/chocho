#!/usr/bin/env python3
"""Hawaii positive alien butterfly establishments vs original WCVP host status.

2024 Zenodo includes 10 independently recorded alien-ESTABLISHED original
chocho butterflies in Hawaii. Is any accepted globally reported HOST plant
actually WCVP native, introduced wild or unlisted wild IN HAWAII (HAW)?
WCVP wild location is NOT local cultivated garden host, consumer feeding
observation, or first arrival. Do not treat 229 other source butterflies as
failed establishment/absent from Hawaii.
"""
from __future__ import annotations

import argparse
import json
import hashlib
from collections import Counter,defaultdict
from pathlib import Path

from analyze_hosts_location_wcvp_origin_20261027 import (
    accepted_pairs, read_rows, sha256, wcvp_range
)

SCHEMA="chocho_zenodo_2024_10_aliens_Hawaii_original_WCVP_hosts_20261027_v01"
SOURCE_ZIP_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
SOURCE_ZENODO_OCC_SHA="05483a2a9a7c926fb5105c44058ee52a4598039dc82c5bb10d1db284b3b66fe7"
HAWAII_WCVP="HAW"
STATUS=("NATIVE_HAW_WILD","INTRODUCED_HAW_WILD","NOT_RECORDED_HAW_WILD")

def classify_plants(plants,native,current,region=HAWAII_WCVP):
    labels=Counter()
    for pid in plants:
        if region in native.get(pid,()):
            status="NATIVE_HAW_WILD"
        elif region in current.get(pid,()):
            status="INTRODUCED_HAW_WILD"
        else:
            status="NOT_RECORDED_HAW_WILD"
        labels[status]+=1
    if sum(labels.values())!=len(plants):
        raise ValueError("Original accepted source host assignments do not conserve")
    return labels


def classify_source_positive_species(hawaii,hostsets,native,current):
    out=[]
    for item in sorted(hawaii,key=lambda x:x["species"]):
        sp=item["species"]
        plants=hostsets.get(sp,set())
        counts=classify_plants(plants,native,current)
        nat=counts["NATIVE_HAW_WILD"]
        intro=counts["INTRODUCED_HAW_WILD"]
        if nat and intro:
            response="BOTH_NATIVE_AND_INTRODUCED_WILD_HOSTS"
        elif intro:
            response="ONLY_WCVP_INTRODUCED_WILD_HOSTS"
        elif nat:
            response="ONLY_WCVP_NATIVE_WILD_HOSTS"
        else:
            response="NO_HOST_IN_WCVP_WILD_HAWAII_MAP"
        reference=int(item["original_239_host_species_breadth"])
        if len(plants)!=reference:
            raise ValueError("Original accepted source per species host link count changed: "+sp)
        years=[
            rec["year_source_first_report"]
            for rec in item["original_source_region_and_year_records"]
            if rec["recipient_zone"]=="Hawaii"
        ]
        flags=[
            {
                "intentional_release":r["source_intentional_release_flag"],
                "eradicated":r["source_eradication_flag"],
                "indoor_outdoor":r["source_established_indoor_outdoor_flag"],
            }
            for r in item["original_source_region_and_year_records"]
            if r["recipient_zone"]=="Hawaii"
        ]
        suspicious=any(y and y.isdigit() and int(y)<=1800 for y in years)
        out.append({
            "species":sp,
            "known_original_accepted_larval_host_species":len(plants),
            "WCVP_Hawaii_wild_host_plants_native":nat,
            "WCVP_Hawaii_wild_host_plants_introduced":intro,
            "known_host_species_without_WCVP_wild_presence_Hawaii":
                counts["NOT_RECORDED_HAW_WILD"],
            "accepted_host_support_class":response,
            "original_global_added_resource_units":item["original_239_source_added_WCVP_resource_units"],
            "2024_source_establishment_year_for_Hawaii":years,
            "source_original_status_flags":flags,
            "chronology_flag":"YEAR_<=1800_CANNOT_UNCRITICALLY_BE_INTERPRETED_AS_FIRST_INVASION" if suspicious else
                  "HISTORICAL_FIRST_YEAR_UNVALIDATED_EXTERNAL_ARCHIVAL_SOURCE",
            "does_NOT_identify_actual_Hawaii_larval_host_use":True
        })
    return out


def run(args):
    initial=json.loads(args.zenodo_positive_json.read_text())
    if initial.get("schema")!="chocho_239_butterfly_original_Zenodo_positive_alien_establishment_20261027_v01":
        raise ValueError("Not original source-locked 2024 established butterflies result")
    if initial["source_zenodo"]["record_id"]!=11519262 or \
       initial["source_zenodo"]["sources"]["occurrence_table_20230317.csv"]["sha256"]!=SOURCE_ZENODO_OCC_SHA:
        raise ValueError("Zenodo original positive source SHA/version changed")
    if initial["source_original239"]["original_239_butterflies"]!=239 or \
       initial["source_positive_matches"][
           "original_239_butterflies_matched_as_established_non_native_in_at_least_one_source_region"]!=18:
        raise ValueError("Published original 239 vs Zenodo 18 source matches changed")
    rows=initial["source_positive_matches"]["source_positive_original239_species"]
    hawaii=[v for v in rows if "Hawaii" in v["source_documented_non_native_recipient_areas"]]
    if len(hawaii)!=10 or len({v["species"] for v in hawaii})!=10:
        raise ValueError("Original confirmed Hawaii alien established butterfly cohort changed")
    root=args.author_source
    manifest=json.loads((root/"SOURCE_MANIFEST.json").read_text())
    if manifest.get("inputs_count")!=11:
        raise ValueError("Not 11 member original author HOSTS/WCVP source")
    for name,digest in manifest["sha256"].items():
        if sha256(root/name)!=digest:
            raise ValueError("Original author-only source constituent hash mismatch")
    regions=json.loads((root/"inputs/level3.geojson").read_text())["features"]
    hawaii_geometry=[f for f in regions if
        (f.get("properties") or {}).get("LEVEL3_COD")==HAWAII_WCVP]
    if len(hawaii_geometry)!=1:
        raise ValueError("Hawaii original WCVP botanical region code not one-to-one validated")
    native,current=wcvp_range(
        root/"inputs/original/native_extant_nondoubtful_wgsrpd3.csv",
        root/"inputs/original/extant_nondoubtful_wgsrpd3.csv")
    original_species={r["species"] for r in
                      read_rows(root/"inputs/anthropogenic_species_metrics.csv")}
    if len(original_species)!=239:
        raise ValueError("Original butterfly source cohort changed")
    accepted,qc,_=accepted_pairs(
        root/"inputs/original/insect_host_accepted.csv",original_species)
    if qc["source_ambiguous_accepted_id_names"]:
        raise ValueError("Source accepted plant taxonomy ambiguous")
    hosts=defaultdict(set)
    for (sp,name),pid in accepted.items():
        hosts[sp].add(pid)
    if sum(len(v) for v in hosts.values())!=2596:
        raise ValueError("Source 2596 original accepted butterfly-host links changed")
    details=classify_source_positive_species(hawaii,hosts,native,current)
    classes=Counter(r["accepted_host_support_class"] for r in details)
    n=Counter()
    for r in details:
        for key in ("WCVP_Hawaii_wild_host_plants_native",
                    "WCVP_Hawaii_wild_host_plants_introduced",
                    "known_host_species_without_WCVP_wild_presence_Hawaii"):
            n[key]+=r[key]
    if sum(classes.values())!=10 or sum(n.values())!=sum(
        r["known_original_accepted_larval_host_species"] for r in details):
        raise ValueError("Source Hawaii positive butterfly host conservation failed")
    report={
        "schema":SCHEMA,
        "status":"SOURCE_HAWAII_ESTABLISHED_POSITIVE_HOST_GEOGRAPHY_ONLY_NOT_LOCAL_HOST_USE",
        "source":{
            "positive_Zenodo_2024_record_id":11519262,
            "positive_2024_occurrence_source_sha256":SOURCE_ZENODO_OCC_SHA,
            "positive_2024_18_match_original_json_sha256":
                sha256(args.zenodo_positive_json),
            "original_author_only_WCVP_source_zip_sha256":SOURCE_ZIP_SHA,
            "original_accepted_WCVP_pair_count":2596,
            "original_WGSRPD3_Hawaii_area_code":"HAW",
            "original_Hawaii_region_geometry_verified_once":True,
        },
        "documented_2024_Hawaii_aliens_from_original_239":10,
        "positive_source_Hawaii_butterfly_host_coverage_classes":dict(classes),
        "total_source_globally_accepted_host_plant_species_by_Hawaii_WCVP_wild_status":dict(n),
        "individual_source_positive_Hawaii_butterflies":details,
        "logical_nonnegative":"NO WCVP wild plant record does NOT imply no cultivated hosts or species cannot establish",
        "scientific_limits":[
            "The 2024 Zenodo source is a POSITIVE alien established fauna list, not original arrival trials or negative records.",
            "The plant WCVP HAW source range is WILD extant non-doubtful, not all nurseries/gardens/greenhouses/urban agriculture.",
            "A globally documented accepted butterfly-host pair plus WCVP HAW source plant does NOT demonstrate wild larval feeding on Hawaii.",
            "Source species can be documented as established but eradicated/introduced intentionally/historical or indoor; flags retained.",
            "This is descriptive host geography for independently established source butterflies, no ecological arrival coupling effect or orphan proportion.",
            "The 229 other original chocho species are UNLISTED at Hawaii, not proven nonestablished or without introduction attempts.",
            "No raw copyrighted HOSTS accepted botanical link plant IDs exported."
        ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":report["status"],
        "source_10_Hawaii_positive_butterflies":10,
        "source_host_coverage_classes":dict(classes),
        "plant_status":dict(n),
        "species_source_status":details,
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--zenodo-positive-json",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
