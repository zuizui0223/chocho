#!/usr/bin/env python3
"""Source-only taxonomy feasibility: identify butterflies among observed caterpillars.

Only exact accepted butterfly binomials and species-rank IDs. No ecological model
or source individual/observer record included in result.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path

SOURCE_SHA="a5d5ff967ddfc79a6e11abb7bab9fae4354e11ca"
FOCAL_SHA="9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"


def git_digest(file):
    raw=file.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()


def run(source,base):
    expert=source/"data/2025-06-17_ExpertIdentification.csv"
    focal=base/"data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    if git_digest(expert)!=SOURCE_SHA or git_digest(focal)!=FOCAL_SHA:
        raise ValueError("Photo expert or source accepted butterfly taxa changed")
    with focal.open(newline="",encoding="utf-8") as f:
        butterfly={r["species"].strip() for r in csv.DictReader(f)}
    if len(butterfly)!=239:
        raise ValueError("Original accepted species cohort changed")
    groups=Counter()
    ranks=Counter()
    taxon_exact=defaultdict(set)
    conflicts=0
    total=0
    original_sightings=set()
    larval_sightings=set()
    species_identified_sightings=set()
    determined_butterfly={}
    observed_species=Counter()
    exact_other_species=Counter()
    with expert.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f)
        if not {"ArthropodSightingFK","TaxonName","Rank","Group"}<=set(reader.fieldnames or []):
            raise ValueError("Original expert survey taxonomic fields missing: "+repr(reader.fieldnames))
        for row in reader:
            total+=1
            animal=str(row.get("Group") or "").strip()
            rank=str(row.get("Rank") or "").strip().lower()
            name=str(row.get("TaxonName") or "").strip()
            sid=str(row.get("ArthropodSightingFK") or "").strip()
            groups[animal or "[empty]"]+=1
            ranks[rank or "[empty]"]+=1
            if sid:
                original_sightings.add(sid)
            if animal!="caterpillar" or not sid:
                continue
            larval_sightings.add(sid)
            if rank!="species" or not name:
                continue
            species_identified_sightings.add(sid)
            taxon_exact[sid].add(name)
            if name in butterfly:
                observed_species[name]+=1
                if sid in determined_butterfly and determined_butterfly[sid]!=name:
                    conflicts+=1
                determined_butterfly[sid]=name
            else:
                exact_other_species[name]+=1
    exact_unique=len(determined_butterfly)
    types=len(set(determined_butterfly.values()))
    return {
        "schema":"chocho_caterpillar_photo_expert_butterfly_exact_taxon_gate_20261026",
        "status":(
            "SPECIES_RESOLVED_BUTTERFLY_LARVAL_LINK_FEASIBLE" if
            exact_unique>=20 and types>=10 and conflicts==0 else
            "NOT_ENOUGH_SPECIES_RESOLVED_BUTTERFLY_CATERPILLAR_RECORDS"
        ),
        "original_research_source_commit":"2c5690ac7c6dba7eca506da886c49781c37617bd",
        "expert_file_git_blob":SOURCE_SHA,
        "original_chocho_butterfly_list_git_blob":FOCAL_SHA,
        "original_photo_identification_rows":total,
        "original_photo_distinct_sightings":len(original_sightings),
        "original_group_row_counts":dict(groups.most_common(15)),
        "original_taxon_rank_row_counts":dict(ranks.most_common(15)),
        "caterpillar_photo_distinct_sightings":len(larval_sightings),
        "species_rank_caterpillar_distinct_sightings":len(species_identified_sightings),
        "identified_species_rank_caterpillar_taxa":len({a for x in taxon_exact.values() for a in x}),
        "expert_sighting_id_with_multiple_species_rank_names":sum(len(x)>1 for x in taxon_exact.values()),
        "exact_239_butterfly_taxon_species_count":types,
        "exact_239_butterfly_unique_caterpillar_photo_sightings":exact_unique,
        "taxon_name_assignment_conflict_count":conflicts,
        "individual_sighting_identifications_owner_names_original_location_not_exported":True,
        "source_info":"Original Caterpillars Count expert photo identifications; broader unidentified moth/butterfly larvae cannot be recoded as butterflies",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=run(args.original_repo,Path(__file__).resolve().parents[1])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
