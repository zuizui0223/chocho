#!/usr/bin/env python3
"""Source-only census of ALL 26 butterfly-positive added-host country cells.

No GBIF plant outcomes, no juvenile-stage selection, no local host-use claim.
Reconstruct exact archived HOSTS/WCVP/BCE source in runner; export only derived
accepted plant IDs/binomials and one request per unique host/country.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

SCHEMA="chocho_all26_positive_country_source_hosts_v1"
OLD="chocho_host_plant_candidates_5_juvenile_added_country_cases_20261021_v1"
EXPECTED_ZIP_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"
FOCAL=frozenset({("Anteos clorinde","CU"),("Colobura dirce","VE"),
 ("Papilio demoleus","CU"),("Papilio demoleus","DO"),("Papilio demoleus","HT")})


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def name_gate(names):
    if len(names)!=1:
        return "AMBIGUOUS_ACCEPTED_ID_NAME"
    name=names[0].strip()
    if not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][a-z-]+",name):
        return "NOT_EXACT_SPECIES_BINOMIAL"
    return "EXACT_SOURCE_BINOMIAL_NEEDS_GBIF_TAXON_RESOLUTION"


def census(model,regions,links,positive,names,old_five,*,locked=True):
    if positive.get("status")!="COMPLETE_48":
        raise ValueError("Original 48-positive GBIF source incomplete")
    if old_five.get("schema")!=OLD or old_five.get("focal_candidate_plant_links")!=20:
        raise ValueError("Five original source cases missing")
    index={sp:i for i,sp in enumerate(model["species"])}
    if len(index)!=len(model["species"]) or (locked and len(index)!=239):
        raise ValueError("Original butterfly identities unexpectedly duplicate/drift")
    original={(r["species"],r["country"]):r for r in positive["observations"]}
    if len(original)!=len(positive["observations"]):
        raise ValueError("Repeated positive butterfly-country observation")
    if locked and (len(original),len({x for x,_ in original}),
                   len({c for _,c in original}),
                   sum(bool(x["added_resource_country"]) for x in original.values()))!=(48,20,23,26):
        raise ValueError("Original positive 48 / species 20 / country 23 / added 26 drifted")
    rows=[]
    for (sp,ct),cell in sorted(original.items()):
        if sp not in index:
            raise ValueError("Original butterfly absent from pinned HOSTS")
        available=links[ct]
        material=available["material"]
        guard=available["guard"]
        if not material or not guard:
            raise ValueError(f"Missing country material/guard area for {sp}/{ct}")
        j=index[sp]
        native=model["native_species"][j]
        if any(native>>i&1 for i in guard):
            raise ValueError(f"Original native-host country exclusion failed {sp}/{ct}")
        proposed=[]
        for host in model["species_hosts"][j]:
            intro=model["introduced_sets"].get(host,set())
            overlapping=sorted(regions[i] for i in material & intro)
            if overlapping:
                n=names.get(host,[])
                proposed.append({
                    "wcvp_accepted_id":str(host),
                    "accepted_names":n,
                    "source_name_gate":name_gate(n),
                    "source_introduced_WGSRPD3":overlapping,
                })
        proposed.sort(key=lambda r:(r["wcvp_accepted_id"],r["accepted_names"]))
        if int(bool(proposed))!=cell["added_resource_country"]:
            raise ValueError(f"Source botanical host opportunity disagrees with frozen original cell {sp}/{ct}")
        rows.append({"butterfly":sp,"country":ct,
                     "added_host_resource_country":cell["added_resource_country"],
                     "source_indexed_butterfly_records":cell["record_count"],
                     "source_juvenile_stage_labels_present":bool(
                         cell["facets"].get("life_stage",{}).get("categories",{})),
                     "candidate_source_plant_count":len(proposed),
                     "candidate_accepted_plants":proposed})
    five_source={(r["species"],r["country"]):r for r in old_five["focal_five_cases"]}
    if set(five_source)!=FOCAL:
        raise ValueError("Five historic source case identities drifted")
    for r in rows:
        key=(r["butterfly"],r["country"])
        if key not in FOCAL:continue
        old=five_source[key]
        old_pairs={(str(p["wcvp_accepted_id"]),tuple(p["botanical_names"]),
                    tuple(p["country_material_introduced_wgsrpd3"]))
                   for p in old["candidate_plants"]}
        new_pairs={(p["wcvp_accepted_id"],tuple(p["accepted_names"]),
                    tuple(p["source_introduced_WGSRPD3"]))
                   for p in r["candidate_accepted_plants"]}
        if old_pairs!=new_pairs:
            raise ValueError(f"Historical five botanical source identities changed {key}")
    added=[x for x in rows if x["added_host_resource_country"]]
    if locked and (len(rows),len(added),sum(bool(x["candidate_source_plant_count"]) for x in rows))!=(48,26,26):
        raise ValueError("Full added-resource positive 26 selection changed")
    links_count=sum(x["candidate_source_plant_count"] for x in added)
    query_groups=defaultdict(lambda:{"linked_butterflies":set(),"area_codes":set(),"name_statuses":set(),"names":set()})
    seen_names=defaultdict(set)
    for row in added:
        for plant in row["candidate_accepted_plants"]:
            gid=(plant["wcvp_accepted_id"],row["country"])
            group=query_groups[gid]
            group["linked_butterflies"].add(row["butterfly"])
            group["area_codes"].update(plant["source_introduced_WGSRPD3"])
            group["name_statuses"].add(plant["source_name_gate"])
            group["names"].update(plant["accepted_names"])
            seen_names[plant["wcvp_accepted_id"]].update(plant["accepted_names"])
    unique=[
        {"wcvp_accepted_id":gid[0],"country":gid[1],
         "accepted_names":sorted(v["names"]),
         "source_name_gate":name_gate(sorted(v["names"])),
         "butterflies_linked":sorted(v["linked_butterflies"]),
         "source_WGSRPD3_intro_regions":sorted(v["area_codes"])}
        for gid,v in sorted(query_groups.items())
    ]
    histogram=dict(sorted(Counter(str(r["candidate_source_plant_count"]) for r in added).items(),
                          key=lambda p:int(p[0])))
    return {
        "schema":SCHEMA,"status":"SOURCE_ONLY_ALL26_NO_GBIF_PLANT_RECORDS",
        "frozen_48_positive_species_country_cells":len(rows),
        "frozen_added_host_cells":len(added),
        "non_added_source_cells":len(rows)-len(added),
        "candidate_butterfly_country_plant_links":links_count,
        "distinct_source_plant_ids":len(seen_names),
        "distinct_plant_country_source_queries":len(unique),
        "source_name_gates":dict(sorted(Counter(v["source_name_gate"] for v in unique).items())),
        "plant_candidate_count_per_added_cell":histogram,
        "source_only_true_absence_not_inferred":True,
        "all_48_positive_country_cells":rows,
        "unique_plant_country_query_basis":unique
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("--source-root","--crosswalk","--positive48","--original-five","--output"):
        p.add_argument(flag,required=True,type=Path)
    a=p.parse_args()
    # Import sibling source code only during real analysis. Offline tests
    # can import the pure core without loading numpy or author source.
    from analyze_bce_augmented_plant_network_climate_null_20261019 import source_augmented_model
    from derive_gbif_positive48_introduced_hosts_20261021 import load_country_links,plant_names
    model,regions,climate,prov=source_augmented_model(a.source_root)
    links=load_country_links(a.crosswalk,regions)
    result=census(model,regions,links,
                  json.loads(a.positive48.read_text(encoding="utf8")),
                  plant_names(a.source_root),
                  json.loads(a.original_five.read_text(encoding="utf8")))
    result["executed_utc"]=datetime.now(timezone.utc).isoformat()
    result["provenance_sha256"]={
        "original_archive_inner_zip":EXPECTED_ZIP_SHA,
        "crosswalk":checksum(a.crosswalk),
        "positive48":checksum(a.positive48),
        "previous_five_botanical_subset":checksum(a.original_five),
        "source_manifest":checksum(a.source_root/"SOURCE_MANIFEST.json"),
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("ALL26_SOURCE_INVENTORY_SUMMARY="+json.dumps(
        {k:v for k,v in result.items() if k not in ("all_48_positive_country_cells","unique_plant_country_query_basis")},
        separators=(",",":")),flush=True)
    print("ALL26_SOURCE_INVENTORY_FULL="+json.dumps(result,separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
