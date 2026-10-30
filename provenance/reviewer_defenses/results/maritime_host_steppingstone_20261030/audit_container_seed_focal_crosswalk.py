#!/usr/bin/env python3
"""Cross-check four documented container-borne plant taxa and known butterfly hosts.

A genus match does NOT validate the same transported species or route.
Study: Lucardi et al. 2020, Sci Rep 10:15017, doi:10.1038/s41598-020-71954-3.
Port inventory: retrospective 2001–2024 species counts, different ports.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import argparse

ROOT = Path(__file__).resolve().parents[4]
BASE = Path("provenance/reviewer_defenses/results")
SPECIES = BASE / "port_alien_broad_inventory_20261030/port_alien_host_broad_species.csv"
EARLY = BASE / "island_dated_host_match_20261030/island_dated_host_match_links.csv"
PINS = {
    str(SPECIES): "9712b8e22a7e74f7e8809e7212ed6a8cfa2158fb",
    str(EARLY): "88db773c74dea9cd2d58bde0b67a394d93226869",
}
# Only four *named focal taxa* from the paper's principal germination work.
# "spp." names cannot be interpreted as species-level detections.
PAPER_FOCAL = (
    ("Saccharum spontaneum", "Saccharum", "species"),
    ("Phragmites spp.", "Phragmites", "genus"),
    ("Typha spp.", "Typha", "genus"),
    ("Andropogon spp.", "Andropogon", "genus"),
)


def checked_csv(root: Path, path: Path):
    body=(root/path).read_bytes()
    obj=hashlib.sha1(b"blob "+str(len(body)).encode()+b"\0"+body).hexdigest()
    if obj != PINS[str(path)]:
        raise ValueError(f"Unexpected source version: {path}")
    return list(csv.DictReader(io.StringIO(body.decode("utf-8"))))


def crosswalk(root: Path = ROOT):
    species=checked_csv(root, SPECIES)
    early=checked_csv(root, EARLY)
    if len(species) != 687 or len(early) != 5:
        raise ValueError("Old source panels changed")
    mapping={}
    for r in species:
        name=r["host_name"]
        if name in mapping:
            raise ValueError("Duplicate host in original port inventory")
        mapping[name]=int(r["ports"])
    result=[]
    for named, genus, resolution in PAPER_FOCAL:
        matches=sorted([{"candidate_host_species":name,"ports_with_current_alien_candidate":
                         n} for name,n in mapping.items() if name.split(" ",1)[0] == genus],
                       key=lambda x:x["candidate_host_species"])
        exact=named if resolution=="species" and named in mapping else None
        result.append({
            "container_study_taxon": named,
            "paper_reported_resolution":resolution,
            "exact_species_level_overlap":exact,
            "port_species_with_same_genus":matches,
            "same_genus_is_not_verified_same_transport_species":True,
        })
    historical=sorted({r["host_name"] for r in early})
    return {
        "schema":"chocho_shipping_container_focal_taxon_crosswalk_20261030_v1",
        "role":"published_crosswalk_not_focal_species_transport_evidence",
        "container_study_doi":"10.1038/s41598-020-71954-3",
        "container_study_scope":"US Port of Savannah refrigerated container intake grilles, sampled in 2015-2017",
        "container_study_taxa":"four focal taxa only, not all 22 sorted taxonomic groups or at least 30 possible species",
        "source_port_inventory_time_window":"2001-2024 pooled",
        "source_port_inventory_origin_class":"contemporary regional alien-host candidate, not proven introduction at time of transport",
        "source_sha_pins":PINS,
        "four_focal_taxon_crosswalk":result,
        "verified_species_level_identity_for_container_plant_vs_port_host":0,
        "genus_level_overlapping_focal_groups":
            sum(bool(r["port_species_with_same_genus"]) for r in result),
        "early_host_species_port_inventory_counts":[
            {"source_historically_early_host":name,
             "ports_in_2001_2024_retrospective_alien_candidate_inventory":mapping.get(name,0),
             "paper_named_focal_transport_study_confirmation":False}
            for name in historical],
        "no_actual_vessel_transport_or_host_plant_origin_assigned":True,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--output",type=Path)
    opts=ap.parse_args()
    res=json.dumps(crosswalk(opts.root),indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if opts.output:
        target=opts.output if opts.output.is_absolute() else opts.root/opts.output
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(res,encoding="utf-8")
        print(target)
    else:
        print(res)


if __name__=="__main__":
    main()
