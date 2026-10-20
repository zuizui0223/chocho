#!/usr/bin/env python3
"""Replay archived GADM1/life-stage facets; NEVER infer verified host use."""
from __future__ import annotations
import argparse,hashlib,json,statistics
from collections import Counter
from pathlib import Path

SOURCE_HASHES={
  "source_strata_sha256":"6eaa91001b282a8bb8bb3f78847386ef5c4fbefcf273f4b09d173d6a691100ec",
  "source_taxon_manifest_sha256":"18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2",
  "source_country_csv_sha256":"d9f5207cd611e737279b7b579fe8f2347132d68bab36fc78a56e727f498c0359"}

def category_totals(cells,field):
    result=Counter()
    for cell in cells:
        facet=cell["facets"].get(field)
        if facet is not None:
            n=sum(facet["categories"].values())
            if n!=facet["assigned_records"] or n+facet["unassigned_records"]!=cell["record_count"]:
                raise ValueError(f"Invalid {field} facet totals")
            result.update(facet["categories"])
    return dict(sorted(result.items()))

def summarise(source):
    if source.get("status")!="COMPLETE_48" or source.get("pilot") is not False:
        raise ValueError("The pinned all-48 complete source must be used")
    if any(source["source"].get(k)!=v for k,v in SOURCE_HASHES.items()):
        raise ValueError("Source stratum/manifest/country hashes drifted")
    rows=source["observations"]
    if len(rows)!=48 or source.get("source_count_exact_cells")!=48 or source.get("source_count_drifted_cells")!=0:
        raise ValueError("Frozen 48 source cells incomplete")
    seen=set();loc=[]
    for r in rows:
        key=(r["species"],r["country"])
        if key in seen or r["status"]!="COMPLETE_GADM1" or r["source_pooled_delta"]!=0:
            raise ValueError("Duplicate or drifted source cell")
        seen.add(key)
        gadm=r["facets"]["gadm1"]
        n=r["record_count"]
        mapped=sum(gadm["categories"].values())
        if n<=0 or mapped+gadm["unassigned_records"]!=n or mapped!=r["mapped_records"]:
            raise ValueError("Invalid or missing GADM mapping")
        k=len(gadm["categories"])
        if k!=r["mapped_gadm_units"] or k<1:
            raise ValueError("Missing or inconsistent GADM1 subdivision bins")
        loc.append({"species":r["species"],"country":r["country"],
                    "added_resource_country":r["added_resource_country"],
                    "record_count":n,"mapped_count":mapped,
                    "gadm1_units":k,
                    "dominant_fraction_all_records":round(max(gadm["categories"].values())/n,6),
                    "life_stage_metadata_count":r["facets"].get("life_stage",{}).get("assigned_records",0)})
    life=category_totals(rows,"life_stage")
    est=category_totals(rows,"establishment_means")
    if (sum(r["record_count"] for r in rows),sum(r["mapped_count"] for r in loc),
        sum(life.values()),sum(est.values())) != (756,736,500,0):
        raise ValueError("Original source totals not reproduced")
    groups={}
    for g in (0,1):
        z=[x for x in loc if x["added_resource_country"]==g]
        if len(z)!=(26 if g else 22):raise ValueError("Original resource flag counts drifted")
        groups[str(g)]={"cells":len(z),"records":sum(x["record_count"] for x in z),
                       "gadm_mapped_records":sum(x["mapped_count"] for x in z),
                       "single_gadm1_cells":sum(x["gadm1_units"]==1 for x in z)}
    ge5=[x for x in loc if x["record_count"]>=5]
    focus=[x for x in loc if (x["species"],x["country"]) in
           {("Aguna asander","DO"),("Papilio demoleus","IL")}]
    juveniles={k:v for k,v in life.items() if k.strip().lower() in
       {"larva","larvae","larval","egg","eggs","pupa","pupae","caterpillar","caterpillars","juvenile"}}
    return {"schema":"gbif_positive48_subcountry_independent_source_replay_v1",
            "status":"SOURCE_REPLAY_PASSED_NO_ESTABLISHMENT_INFERENCE",
            "source_sha256":SOURCE_HASHES,
            "positive_species_country_cells":48,
            "indexed_records":756,"gadm1_assigned_records":736,
            "gadm1_unassigned_records":20,
            "life_stage_annotated_records":500,
            "establishment_means_annotated_records":0,
            "life_stage_raw_labels":life,
            "establishment_means_raw_labels":est,
            "putative_juvenile_raw_labels_unverified":juveniles,
            "cells_with_one_gadm1":sum(x["gadm1_units"]==1 for x in loc),
            "gadm1_unit_count_median":statistics.median(x["gadm1_units"] for x in loc),
            "cells_at_least_5_records":len(ge5),
            "gadm1_units_median_at_least_5":statistics.median(x["gadm1_units"] for x in ge5) if ge5 else None,
            "by_added_resource_country":groups,
            "two_source_audit_cases":focus,
            "per_species_country_safe_summary":loc,
            "interpretation":"GADM1 metadata and provider-submitted stage labels, not independent populations, verified larvae, local host use, or colonization."}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--archived-source",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    raw=a.archived_source.read_bytes()
    source=json.loads(raw)
    result=summarise(source)
    result["archived_source_file_sha256"]=hashlib.sha256(raw).hexdigest()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    compact={k:v for k,v in result.items() if k!="per_species_country_safe_summary"}
    print("GBIF_POS48_REPLAY_REPORT="+json.dumps(compact,separators=(",",":")),flush=True)
    print("GBIF_POS48_REPLAY_FULL="+json.dumps(source,separators=(",",":")),flush=True)

if __name__=="__main__":
    main()
