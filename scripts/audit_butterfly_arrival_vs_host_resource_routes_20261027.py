#!/usr/bin/env python3
"""Host-resource geography and plant propagation are different channels.

Original 239 butterfly x accepted WCVP host resource contribution reconstruction;
ten POSTHOC top-source-plant entries independently linked to propagation and
nursery-occurrence SOURCES. No transborder shipment histories or transported
insect life stages are observed. Original 115 outside-native butterfly source
units are ALL observed, including the 49 not resource-recovered units.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

SCHEMA="chocho_arrival_resource_decoupling_top10_posthoc_feasibility_20261027_v01"
SOURCE_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
TOP10=("Medicago sativa","Poa pratensis","Oxalis corniculata",
       "Cynodon dactylon","Dactylis glomerata","Zea mays",
       "Avena sativa","Urtica urens","Oryza sativa","Digitaria sanguinalis")
MODES=(
    "documented_seed_planting_or_seed_spread",
    "documented_live_sod_sprigs_or_rooted_material",
    "documented_container_nursery_weed_association",
    "documented_local_seedling_transplantation",
)
TRADE=("crossborder_commodity_transit_by_this_form",
       "butterfly_egg_larva_transport_attributed_to_this_commodity")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def verify_source(root):
    path=root/"SOURCE_MANIFEST.json"
    meta=json.loads(path.read_text())
    if meta.get("inputs_count")!=11:
        raise ValueError("Not the original author's 11-member source release")
    for name,digest in meta["sha256"].items():
        if sha(root/name)!=digest:
            raise ValueError("Author-only source file SHA altered: "+name)
    species=read_rows(root/"inputs/anthropogenic_species_metrics.csv")
    if len(species)!=239 or len({s["species"] for s in species})!=239:
        raise ValueError("Wrong 239 source butterflies")
    return species


def plant_footprints(rows):
    dest=defaultdict(set)
    for row in rows:
        pid=str(row["accepted_plant_name_id"]).strip()
        region=str(row["area_code_l3"]).strip()
        if pid and region:
            dest[pid].add(region)
    return dest


def fractional_plant_credits(species, links, native, contemporary):
    """Allocate one added butterfly-region resource unit among its hosts."""
    names=defaultdict(set)
    families=defaultdict(set)
    source_hosts=defaultdict(set)
    focal={s["species"] for s in species}
    for row in links:
        if row["insect_species"] not in focal:
            continue
        pid=row["accepted_plant_name_id"]
        source_hosts[row["insect_species"]].add(pid)
        names[pid].add(row["accepted_name"])
        families[pid].add(row.get("family",""))
    if sum(map(len,source_hosts.values()))!=2596:
        raise ValueError("Original 239 butterflies × accepted plant WCVP links altered")
    credits=Counter()
    resource=Counter()
    expanded=0
    for sp in sorted(focal):
        h=source_hosts.get(sp,set())
        natives=set()
        contributors=defaultdict(list)
        for pid in h:
            natives.update(native.get(pid,set()))
            for area in contemporary.get(pid,set()):
                contributors[area].append(pid)
        current=set(contributors)
        added=current-natives
        if len(added)!=int(next(s for s in species if s["species"]==sp)["introduced_added_units"]):
            raise ValueError("Old source range-enlargement unit for butterfly differs "+sp)
        resource[sp]=len(added)
        expanded+=int(bool(added))
        for area in added:
            pids=contributors[area]
            if not pids:
                raise ValueError("An added source area lacks any documented host")
            weight=1/len(pids)
            for pid in pids:
                credits[pid]+=weight
    if len(credits)!=670 or sum(resource.values())!=14553 or expanded!=206:
        raise ValueError("Original 239/206/14553/670 HOSTS/WCVP resource reconstruction changed")
    if not math.isclose(sum(credits.values()),14553,abs_tol=1e-6):
        raise ValueError("Fractional host source unit conservation failed")
    ranked=sorted(credits.items(),key=lambda x:(-x[1],x[0]))
    top10=sum(value for _,value in ranked[:10])/14553.
    top38=sum(value for _,value in ranked[:38])/14553.
    if not math.isclose(top10,.25101051361138516,abs_tol=1e-8) or \
       not math.isclose(top38,.5056624499965091,abs_tol=1e-8):
        raise ValueError("Original fractional host plant credit concentration changed")
    top=[]
    for i,(pid,score) in enumerate(ranked[:10],1):
        if len(names[pid])!=1 or len(families[pid])!=1:
            raise ValueError("Original accepted source plant identity ambiguous")
        top.append({
            "original_rank":i,
            "accepted_botanical_name":next(iter(names[pid])),
            "family":next(iter(families[pid])),
            "source_fractional_added_resource_credit":float(score),
            "fraction_of_all_14553_original_added_units":float(score/14553.),
        })
    if tuple(x["accepted_botanical_name"] for x in top)!=TOP10:
        raise ValueError("Cited top-ten propagation mode ledger must be aligned to original output")
    return top,{
        "source_species":len(focal),
        "source_original_accepted_butterfly_host_pairs":sum(map(len,source_hosts.values())),
        "original_expanded_butterfly_species":expanded,
        "original_added_butterfly_region_resource_units":sum(resource.values()),
        "source_distinct_contributing_plant_species":len(credits),
        "top_10_host_fractional_credit_fraction":float(top10),
        "top_38_host_fractional_credit_fraction":float(top38),
    }


def load_route_pilot(path):
    rows=read_rows(path)
    if len(rows)!=10:
        raise ValueError("Ten source-top-ranked original host species required")
    if tuple(r["accepted_botanical_name"] for r in rows)!=TOP10:
        raise ValueError("Route labels were not fixed on the source top ten species")
    for row in rows:
        if not row["source_URL"].startswith("https://"):
            raise ValueError("Propagation source citation needed for every taxon")
        for field in MODES:
            if row[field] not in ("0","1"):
                raise ValueError("Original botanical practice requires binary documented flag")
        for field in TRADE:
            if row[field]!="NOT_VERIFIED":
                raise ValueError("This version contains NO species×route×country crossing observation")
        if not row["supported_practice_scope_not_transborder_trade"]:
            raise ValueError("Propagation source does not license trade vector inference")
    return rows


def audit115(receipt_path,species_table):
    read=json.loads(Path(receipt_path).read_text())
    if (read["outside_native_species_x_units"],
        read["outside_native_species_x_units_recovered_by_introduced_hosts"])!=(115,66):
        raise ValueError("Original contemporaneous butterfly observed-only units source altered")
    rows=read_rows(species_table)
    total=sum(int(r["observed_outside_native"]) for r in rows)
    recovered=sum(int(r["outside_native_units_explained_by_introduced_host_ranges"]) for r in rows)
    if total!=115 or recovered!=66 or len(rows)!=32:
        raise ValueError("Original 32-species field validation source margins changed")
    return {
        "original_panel_species":32,
        "outside_native_butterfly_observed_units":115,
        "observed_units_with_introduced_resource_overlap":66,
        "observed_units_not_recovered_by_any_reconstructed_current_host_range":49,
        "butterfly_occurrence_outcome_by_original_115_units":"ALL_115_OBSERVED",
        "no_unobserved_units_supplied_for_arrival_presence_absence_comparison":True,
        "critical_correction":"The other 49 are butterfly PRESENT but not resource-explained. They are not butterfly absences.",
    }


def join_top10(top,ledger,source):
    score_by_name={p["accepted_botanical_name"]:p for p in top}
    detail=[]
    total_mode=Counter()
    for row in ledger:
        taxon=row["accepted_botanical_name"]
        host=score_by_name[taxon]
        entry={**host,**{s:(row[s]=="1") for s in MODES}}
        entry.update({
            "plant_propagation_or_nursery_relevance_evidence_url":row["source_URL"],
            "route_annotation":row["supported_practice_scope_not_transborder_trade"],
            "actual_species_level_international_plant_shipment_observed":None,
            "butterfly_egg_or_larva_imported_with_this_plant_documented":None,
        })
        detail.append(entry)
        for mode in MODES:
            if entry[mode]:
                total_mode[mode]+=host["fraction_of_all_14553_original_added_units"]
    return {
        "original_top10":detail,
        "mode_coverage_of_14553_original_added_resource_units_not_mutually_exclusive":dict(total_mode),
        "top10_fraction_of_original_resource_gain":source["top_10_host_fractional_credit_fraction"],
        "not_a_multinomial_mode_partition":True,
        "unclassified_remaining_660_host_species_resource_fraction":
            1-source["top_10_host_fractional_credit_fraction"],
        "all_top10_international_shipment_exposures":None,
        "all_top10_butterfly_stage_interceptions_on_named_host_commodity":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--route-ledger",type=Path,required=True)
    p.add_argument("--original-observed115-receipt",type=Path,required=True)
    p.add_argument("--original-observed32-species-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    root=args.author_source
    species=verify_source(root)
    base=root/"inputs/original"
    hosts=read_rows(base/"insect_host_accepted.csv")
    native=plant_footprints(read_rows(base/"native_extant_nondoubtful_wgsrpd3.csv"))
    contemporary=plant_footprints(read_rows(base/"extant_nondoubtful_wgsrpd3.csv"))
    if any(not v<=contemporary.get(k,set()) for k,v in native.items()):
        raise ValueError("Source native WCVP accepted plant range does not nest within contemporary")
    top,credits=fractional_plant_credits(species,hosts,native,contemporary)
    mode=join_top10(top,load_route_pilot(args.route_ledger),credits)
    old=audit115(args.original_observed115_receipt,args.original_observed32_species_csv)
    report={
        "schema":SCHEMA,
        "status":"SOURCE_VERIFIED_RESOURCE_CONTRIBUTIONS_WITHOUT_OBSERVED_ARRIVAL_VECTOR",
        "original_author_source_zip_sha256":SOURCE_SHA,
        "source":credits,
        "corrected_original_occurrence_panel":old,
        "source_route_pilot":mode,
        "interpretive_stop":[
            "Original 115 observed butterfly region records cannot be used for butterfly-present versus butterfly-absent arrival regression.",
            "Plant-seeding, sodding, transplants and nursery weed association are practices/biological occurrences, not international commodity shipments.",
            "Source propagation routes can overlap: sums of mode percentages MUST NOT be used as a global share of butterfly arrival vectors.",
            "Source top-ten plant modes were annotated posthoc after rank discovery and exclude 660 other contributing host plants.",
            "This original source includes 239 focal butterflies; Cycas-feeding Chilades pandava and Pelargonium-feeding Cacyreus marshalli are absent from the original 239 cohort.",
            "The 97 European alien Lepidoptera review includes moths; its 52.5 percent park/garden habitat association is NOT a measured fraction transported on living nursery plants.",
            "No direct import shipment, eggs/larvae interception, first arrival, probability of establishment or source wild breeding data are available in this pilot.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
        "source":credits,
        "corrected_original_115":old,
        "route_mode_contribution":mode["mode_coverage_of_14553_original_added_resource_units_not_mutually_exclusive"],
        "actual_arrival_data":None,
    },sort_keys=True))


if __name__=="__main__":
    main()
