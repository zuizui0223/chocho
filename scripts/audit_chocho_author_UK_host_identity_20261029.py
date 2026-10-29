#!/usr/bin/env python3
"""Inspect original authorized GEB source archive host-species links for 27 UKBMS butterflies.

This is an author-only source audit, not replication or redistribution of raw
HOSTS/WCVP resources. Export only original filenames, schema, counts and names
of already public focal butterfly taxa.
"""
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path

ROSTER=27
EXPECTED_ROOT_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"

def audit(root):
    root=Path(root)
    manifest=root/"SOURCE_MANIFEST.json"
    if not manifest.is_file():
        raise ValueError("Original SHA-verified author-source root lacks manifest")
    source_files=[p for p in (root/"inputs").rglob("*") if p.is_file()]
    names=[str(p.relative_to(root)) for p in source_files]
    expected=[p for p in source_files if p.name=="insect_host_accepted.csv"]
    if not expected:
        return {"schema":"chocho_author_only_UK_host_identity_coverage_20261029_v01",
                "status":"NO_ACCEPTED_HOST_IDENTITY_TABLE_IN_ORIGINAL_AUTHORIZED_ARCHIVE",
                "source_files":names[:100],"host_file_count":0,
                "raw_original_HOSTS_rows_exported":False}
    if len(expected)!=1:raise ValueError("Conflicting original host identity sources")
    p=expected[0]
    counts=defaultdict(set)
    duplicates=0;n=0;accepted_names=0
    with p.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        cols=reader.fieldnames or []
        required={"insect_species","accepted_name","accepted_plant_name_id"}
        if not required<=set(cols):raise ValueError("Original WCVP host identity columns changed")
        from audit_uk_ltmn_butterfly_detection_panel_20261028 import CROSSWALK
        for row in reader:
            n+=1
            insect=row["insect_species"].strip()
            name=row["accepted_name"].strip()
            plant=row["accepted_plant_name_id"].strip()
            if insect in CROSSWALK and name and plant:
                counts[insect].add((plant,name))
    return {"schema":"chocho_author_only_UK_host_identity_coverage_20261029_v01",
            "status":"AUTHOR_ONLY_WCVP_ACCEPTED_UK_HOST_EDGES_RESOLVED",
            "source_file":str(p.relative_to(root)),
            "source_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
            "source_total_rows":n,
            "source_columns":cols,
            "original_uk_source_butterfly_taxa":len(CROSSWALK),
            "uk_taxa_with_at_least_one_resolved_host":len(counts),
            "per_source_uk_butterfly_accepted_host_species_count":{k:len(v) for k,v in sorted(counts.items())},
            "uk_original_unique_host_binomial_names":len({name for v in counts.values() for plant,name in v}),
            "no_source_plant_or_insect_individual_rows_exported":True,
            "no_original_ukbms_station_ids_exported":True}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    out=audit(args.source_root)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(out,sort_keys=True))
    if out["status"]!="AUTHOR_ONLY_WCVP_ACCEPTED_UK_HOST_EDGES_RESOLVED":
        raise SystemExit(2)
