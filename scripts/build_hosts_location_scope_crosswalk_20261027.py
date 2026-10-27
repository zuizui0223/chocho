#!/usr/bin/env python3
"""Lossless original HOSTS Location-label geography scope crosswalk.

Only a conservative pre-WCVP name/type screen. No native/introduced host
status, no inferred biological location for biogeographic realm labels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

SCHEMA = "chocho_original_HOSTS_location_scope_dictionary_v01"

COUNTRY = {
    "USA":"US","India":"IN","Canada":"CA","Japan":"JP","Finland":"FI",
    "Brazil":"BR","Thailand":"TH","Taiwan":"TW","Sri Lanka":"LK",
    "Australia":"AU","Puerto Rico":"PR","Indonesia":"ID",
    "Cuba":"CU","Spain":"ES","South Africa":"ZA","Russia":"RU",
    "Mexico":"MX","China":"CN","Pakistan":"PK","Hungary":"HU",
    "Philippines":"PH","New Zealand":"NZ","Costa Rica":"CR","Nigeria":"NG",
    "Kenya":"KE","Argentina":"AR","Uganda":"UG","Ghana":"GH",
    "Singapore":"SG","Burma":"MM","Ecuador":"EC","Guam":"GU",
    "Papua New Guinea":"PG","Chile":"CL","El Salvador":"SV",
    "Germany":"DE","France":"FR","Colombia":"CO","Turkey":"TR",
    "Cameroun":"CM","Rumania":"RO","Mauritius":"MU","Israel":"IL",
    "Italy":"IT","Malawi":"MW","Guyana":"GY","Madagascar":"MG",
    "Iraq":"IQ","Jamaica":"JM","Barbados":"BB","Nepal":"NP","Brunei":"BN",
    "Somalia":"SO","Peru":"PE","Lebanon":"LB","Bangladesh":"BD",
    "Cyprus":"CY","Zimbabwe":"ZW","Tanzania":"TZ","Egypt":"EG",
    "Zambia":"ZM","Austria":"AT","Morocco":"MA","Swaziland":"SZ",
    "Laos":"LA","Turkmenistan":"TM","Venezuela":"VE",
    "Lithuania":"LT","Greece":"GR","Switzerland":"CH","Mozambique":"MZ",
    "Netherlands":"NL","Bhutan":"BT","Afghanistan":"AF","Uruguay":"UY",
    "Seychelles":"SC","Senegal":"SN","Poland":"PL","Kazakhstan":"KZ",
    "Ethiopia":"ET","Algeria":"DZ","Bulgaria":"BG","Bermuda":"BM",
    "Bolivia":"BO","Burkina Faso":"BF","Ivory Coast":"CI","Libya":"LY",
    "Fiji":"FJ","Solomon Islands":"SB","Dominican Republic":"DO",
    "Malaysia":"MY","Hong Kong":"HK",
}

SUBNATIONAL = {
    "Hawaii":"US", "West Malaysia":"MY", "East Malaysia":"MY",
    "Sabah":"MY","Sarawak":"MY","Ryukyu Islands":"JP",
    "Andaman Islands":"IN","Andaman & Nicobar Is":"IN",
    "Canary Islands":"ES","Sardinia":"IT","Java":"ID",
}
REGIONAL = {
    "Nearctic","Palaearctic","East Palaearctic","West Palaearctic",
    "Oriental","Oriental/Palaearctic","Neotropical","Afrotropical",
    "Holarctic","Cosmopolitan","Pantropical","Indo-Australian",
    "Old World","New World","Old World tropics","East Africa",
    "West Africa","North Africa","Southern Africa","Central Asia",
    "South-East Asia","Europe","Australasia","Pacific",
    "Indian subregion","West Indies","Micronesia",
    "Mediterranean","Middle East","Borneo","British Isles",
}
HISTORICAL_AMBIGUOUS = {
    "Czechoslovakia","former USSR","Yugoslavia","Rhodesia","Sudan",
    "Zaire","Korea",
}


def classify(label):
    if "[interception]" in label.lower():
        return "INTERCEPTION_NOT_FIELD_HOST_USE", None
    if "(prov)" in label.lower():
        return "PROVISIONAL_LOCATION_NOT_PRIMARY", None
    if label in HISTORICAL_AMBIGUOUS:
        return "HISTORICAL_OR_MULTINATIONAL_NOT_PRIMARY", None
    if label in COUNTRY:
        return "SINGLE_COUNTRY_OR_TERRITORY_NAME", COUNTRY[label]
    if label in SUBNATIONAL:
        return "SUBNATIONAL_NAME_NEEDS_GEOMETRY", SUBNATIONAL[label]
    if label in REGIONAL:
        return "BIOGEOGRAPHIC_OR_MULTICOUNTRY_NOT_PRIMARY", None
    return "UNREVIEWED_LABEL_NOT_PRIMARY", None


def from_counts(counts, sha):
    rows=[]
    categories=Counter()
    for label,n in sorted(counts.items(),key=lambda kv:(-kv[1],kv[0])):
        status, country=classify(label)
        rows.append({
            "original_location_label":label,
            "original_source_species_level_host_rows":n,
            "mapping_scope":status,
            "iso2_country_or_parent_hint":country,
            "wgsrpd3_source_region_code":None,
            "requires_original_geo_overlay":True,
            "plant_native_or_introduced_status":"NOT_CLASSIFIED",
        })
        categories[status]+=n
    return {
        "schema": SCHEMA,
        "status":"ALL_LABELS_SCOPED_SOURCE_ONLY_NEITHER_HOST_USE_NOR_WCVP_NATIVITY_INFERRED",
        "original_raw_hosts_resource_sha256":sha,
        "source_species_level_annotated_rows":sum(counts.values()),
        "distinct_location_labels":len(counts),
        "row_counts_by_scope":dict(categories),
        "source_verified_label_entries":rows,
        "limitations":[
            "Exact country names may span multiple WGSRPD3 areas with different accepted plant native/introduced statuses.",
            "Subnational names need exact country×botanical region geometry, never parent country aggregation.",
            "Interceptions, provisional locality, historical countries and biogeographic realms are excluded from primary geographical host-use inference.",
            "This is a label dictionary only: even single-country Location does not prove larval feeding at one exact coordinate.",
        ],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-hosts-resource",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    bylabel=Counter()
    with args.original_hosts_resource.open(encoding="utf-8-sig",newline="") as f:
        records=csv.DictReader(f)
        fields={"Location","Insect Genus","Insect Species",
                "Hostplant Genus","Hostplant Species"}
        if not fields<=set(records.fieldnames or []):
            raise ValueError("Original HOSTS Location fields missing")
        for r in records:
            if not all((r.get(name) or "").strip() for name in fields):
                continue
            bylabel[r["Location"].strip()]+=1
    if len(bylabel)!=257 or sum(bylabel.values())!=86826:
        raise ValueError("Original 257/86826 HOSTS Location frequencies drifted")
    sha=hashlib.sha256(args.original_hosts_resource.read_bytes()).hexdigest()
    result=from_counts(bylabel,sha)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({
        "source_exact_labels":result["distinct_location_labels"],
        "source_species_rows":result["source_species_level_annotated_rows"],
        "scope_counts":result["row_counts_by_scope"],
    },ensure_ascii=False))


if __name__=="__main__":
    main()
