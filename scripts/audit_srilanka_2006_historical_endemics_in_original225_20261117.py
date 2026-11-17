#!/usr/bin/env python3
"""Historical (2006) E-marked Sri Lankan butterfly names versus original 225.

This is an exact-binomial source-selection audit, NOT modern endemicity
classification, accepted-name reconciliation or a complete audit of all endemic
taxa in the 2006 book.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PINS={
 "gbif":("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv","b6a3c98eb76864052d45d520a55b7d830816126c"),
 "historical_E":("results/butterfly_GBIF_effort_origin_gate_20261117/SRILANKA_2006_HISTORICAL_ENDEMIC_E_18_NAME_MATCH_20261117.csv","59e2c707f06e32eee04b45ffe3eda0bbcef798cf")
}

def source(k):
    path,sha=PINS[k]
    p=ROOT/path
    b=p.read_bytes()
    got=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    if got!=sha:raise ValueError("Git blob source SHA drift "+k)
    with p.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def audit(gbif,historic):
    butterflies={r["species"] for r in gbif}
    srcids={(r["speciesKey"],r["landmass_id"]) for r in gbif}
    if (len(butterflies),len(srcids),len(gbif))!=(225,11475,11475):
        raise ValueError("Original 225x51 species/island panel changed")
    if len(historic)!=18:
        raise ValueError("Historical checklist excerpt expected 18 documented E-marked taxa")
    names=set()
    pages={59,60,61,62,63,64}
    for row in historic:
        name=row["historical_2006_binomial"]
        if name in names or len(name.split())!=2 or int(row["printed_page"]) not in pages:
            raise ValueError("Unresolved source taxon duplication/name/page")
        names.add(name)
        if row["2006_checklist_annotated_endemic_E"]!="yes" or row["source_url"]!="https://portals.iucn.org/library/efiles/documents/2006-030.pdf":
            raise ValueError("Source endpoint/E marker not supported")
        if (row["current_accepted_name_reviewed"],row["current_endemicity_reviewed"])!=("no","no"):
            raise ValueError("Current classification must not be promoted from 2006 source")
        if row["original225_exact_binomial_match"]!=("no" if name not in butterflies else "yes"):
            raise ValueError("Original butterfly exact-name presence mismatches source")
    match=names&butterflies
    if match:raise ValueError("Unexpected E-marked source exact-name intersection")
    return {
      "schema":"chocho_sri_lanka_2006_IUCN_historical_E_binomials_HOSts225_exactmatch_audit_20261117_v1",
      "status":"HISTORICAL_SOURCE_SELECTION_BIAS_SCREEN_NOT_CURRENT_NATIVE_ENDEMISM",
      "source_git_blobs":{k:v[1] for k,v in PINS.items()},
      "original_butterflies":len(butterflies),
      "original_species_island_cells":len(gbif),
      "historical_2006_E_labeled_names_explicitly_checked":len(names),
      "names_matching_original_225_by_unreconciled_exact_binomial":len(match),
      "historical_2006_E_labeled_names_absent_by_exact_binomial":len(names-match),
      "names_checked":sorted(names),
      "boundary":[
        "The full published source reports 20 endemics but this table lists 18 individually verified E-marked names; this is not a complete checklist re-extraction.",
        "Historical E labels are not contemporary endemic statuses; no modern taxonomic synonymy was reviewed.",
        "No exact binomial match in 225 can reflect HOSTS selection, modern accepted name changes or nomenclatural mismatch; no claim all 18 biological taxa absent.",
        "Do not assign native, introduced or endemic status to GBIF candidate butterflies from this historical source table."
      ]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("_historical_srilanka_endemics"))
    z=p.parse_args()
    result=audit(source("gbif"),source("historical_E"))
    z.out.mkdir(parents=True,exist_ok=True)
    (z.out/"IUCN2006_SRI_LANKA_HISTORICAL_ENDEMIC_NAME_SELECTION_AUDIT.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"checked":result["historical_2006_E_labeled_names_explicitly_checked"],
         "exact_overlap":result["names_matching_original_225_by_unreconciled_exact_binomial"]}))
if __name__=="__main__":main()
