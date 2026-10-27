#!/usr/bin/env python3
"""Check whether 5 original-239 invasive butterfly taxa retain usable HOSTS location reports.

Reports only aggregated counts: no individual copyrighted HOSTS association
rows or accepted host identities are redistributed.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

SPECIES = (
    "Pieris rapae", "Thymelicus lineola", "Pieris brassicae",
    "Danaus plexippus", "Papilio demoleus",
)
RECIPIENT_LABELS = {
    "Pieris rapae": ("Canada", "USA", "Australia", "New Zealand", "Hawaii"),
    "Thymelicus lineola": ("Canada", "USA"),
    "Pieris brassicae": ("Chile",),
    "Danaus plexippus": ("Hawaii", "Australia", "New Zealand"),
    "Papilio demoleus": ("Dominican Republic", "Puerto Rico", "Cuba",
                        "Jamaica", "West Indies"),
}
SOURCE = {
    "hosts_ref": "globalbioticinteractions/HOSTS@808e0b869f9ec1adf8efff87cf6a395adda103e0",
    "accepted_source_zip_sha256": "408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323",
}


def count(records_path: Path, pairs_path: Path, source_species_path: Path):
    with source_species_path.open(encoding="utf-8", newline="") as f:
        panel={r["species"] for r in csv.DictReader(f)}
    if len(panel)!=239 or not set(SPECIES)<=panel:
        raise ValueError("The original 239 source butterfly cohort changed")
    matched_pairs = defaultdict(set)
    pair_rows = Counter()
    with pairs_path.open(encoding="utf-8",newline="") as f:
        rows=csv.DictReader(f)
        expected={"input_host_name","insect_species","accepted_plant_name_id"}
        if not expected<=set(rows.fieldnames or []):
            raise ValueError("Original accepted WCVP host pair schema missing")
        for r in rows:
            sp=r["insect_species"]
            if sp in SPECIES:
                key=(sp,r["input_host_name"])
                matched_pairs[key].add(r["accepted_plant_name_id"])
                pair_rows[sp]+=1
    if not matched_pairs:
        raise ValueError("No fixed original accepted HOSTS source names")
    invalid={(sp,pair):len(ids) for (sp,pair),ids in matched_pairs.items()
             if len(ids)!=1}
    if invalid:
        raise ValueError("Ambiguous source exact HOSTS to WCVP accepted plant joins")
    with records_path.open(encoding="utf-8-sig",newline="") as f:
        original=csv.DictReader(f)
        expected={"Location","Insect Genus","Insect Species",
                  "Hostplant Genus","Hostplant Species"}
        if not expected<=set(original.fieldnames or []):
            raise ValueError("Source raw HOSTS Location not present")
        before=Counter()
        after=Counter()
        locations=defaultdict(Counter)
        accepted_hosts=defaultdict(lambda:defaultdict(set))
        raw_pair=defaultdict(set)
        for row in original:
            insect=" ".join((row["Insect Genus"] or "").strip().split())+" "+\
                   " ".join((row["Insect Species"] or "").strip().split())
            if insect not in SPECIES:
                continue
            before[insect]+=1
            plant=" ".join((row["Hostplant Genus"] or "").strip().split())+" "+\
                  " ".join((row["Hostplant Species"] or "").strip().split())
            key=(insect,plant)
            if key not in matched_pairs:
                continue
            location=(row["Location"] or "").strip()
            after[insect]+=1
            locations[insect][location or "<EMPTY>"]+=1
            plant_id=next(iter(matched_pairs[key]))
            accepted_hosts[insect][location or "<EMPTY>"].add(plant_id)
            raw_pair[insect].add(plant_id)
    taxa=[]
    for sp in SPECIES:
        stat=locations[sp]
        taxa.append({
            "species":sp,
            "raw_source_rows_with_species":before[sp],
            "raw_rows_with_exact_original_WCVP_accepted_pair":after[sp],
            "distinct_WCVP_accepted_plant_species_in_mapped_HOSTS_rows":len(raw_pair[sp]),
            "distinct_location_strings_among_accepted_host_reports":len(stat),
            "native_status_not_yet_classified":True,
            "top_30_location_labels":[
                {"label":loc, "raw_rows":n,
                 "distinct_accepted_host_species":len(accepted_hosts[sp][loc])}
                for loc,n in stat.most_common(30)],
            "prespecified_invasion_location_source_report_counts":[
                {"label":loc, "raw_rows":stat[loc],
                 "distinct_accepted_host_species":len(accepted_hosts[sp][loc])}
                for loc in RECIPIENT_LABELS[sp]],
        })
    return {
        "schema":"chocho_5_candidate_invasive_butterflies_location_evidence_20261027_v01",
        "status":"SOURCE_EXACT_LOCATION_COVERAGE_AUDIT_NO_HOST_SHIFT_CLAIM",
        "source":{**SOURCE,
                  "original_hosts_resource_sha256":hashlib.sha256(
                      records_path.read_bytes()).hexdigest(),
                  "original_accepted_pairs_sha256":hashlib.sha256(
                      pairs_path.read_bytes()).hexdigest()},
        "eligible_species":5,
        "rows_found_for_all_species":all(t["raw_source_rows_with_species"]>0
                                      for t in taxa),
        "taxa":taxa,
        "scientific_limits":[
            "A missing geographic HOSTS record is not evidence that no host interaction exists in nature.",
            "Source list is larval host association catalogues; reported Location often broad or nonstandard.",
            "Original accepted HOSTS/WCVP pair is fixed globally, not direct verification of specific local insect feeding.",
            "This result does not label WCVP plant native/introduced region or identify species-level host adoption.",
            "Original source host pair IDs and raw association records were not exported.",
        ],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-resource",type=Path,required=True)
    p.add_argument("--accepted-original-host-pairs",type=Path,required=True)
    p.add_argument("--source-original-species",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    r=count(args.hosts_resource,args.accepted_original_host_pairs,
            args.source_original_species)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({
        "status":r["status"],
        "taxa":[{"species":x["species"],
                 "original_species_raw_records":x["raw_source_rows_with_species"],
                 "original_matched_source_rows":x["raw_rows_with_exact_original_WCVP_accepted_pair"],
                 "documented_location_labels":x["distinct_location_strings_among_accepted_host_reports"],
                 "invasion_region_reports":x["prespecified_invasion_location_source_report_counts"]}
                for x in r["taxa"]]
    },ensure_ascii=False))


if __name__=="__main__":
    main()
