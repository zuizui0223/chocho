#!/usr/bin/env python3
"""Accepted-ID reconciliation of BCE-derived Clarke host names vs frozen HOSTS.

Executed only after full 92-species lexical source was available. The output
describes database-coverage differences, NOT new feeding observations.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path

def loadcsv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def source(path):
    v=json.loads(path.read_text(encoding="utf-8"))
    if v["schema"]!="chocho_bce_clarke_all_exact_92_source_names_v0.1":
        raise RuntimeError("Unexpected full European source")
    if v["exact_bce_chocho_intersection"] != 92:
        raise RuntimeError("Source not full 92 exact butterflies")
    if len(v["per_butterfly"]) != 92:
        raise RuntimeError("Incomplete exact intersect source")
    return v

def prepare(args):
    v=source(args.input_json)
    names=sorted({name for row in v["per_butterfly"]
                  if row["site_status"]=="SOURCE_VERIFIED"
                  for name in row["BCE_only_exact_spelling"]})
    if len(names)<10:raise RuntimeError("No botanical candidates")
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",encoding="utf-8",newline="") as f:
        wr=csv.writer(f);wr.writerow(["plant_binomial"])
        wr.writerows([[n] for n in names])
    print("BCE_ONLY_EXACT_CANDIDATE_NAMES",len(names),flush=True)

def audit(args):
    v=source(args.input_json)
    m=loadcsv(args.mapping_csv)
    mp={row["plant_binomial"]:row for row in m}
    if len(mp)!=len(m):
        raise RuntimeError("Duplicate accepted name candidates")
    canonical=defaultdict(set)
    for row in loadcsv(args.original_interaction_csv):
        b=row.get("insect_species","").strip()
        id=row.get("accepted_plant_name_id","").strip()
        if id:canonical[b].add(id)
    checked=[]
    for r in v["per_butterfly"]:
        if r["site_status"]!="SOURCE_VERIFIED":continue
        butterfly=r["species"]
        counts=Counter()
        details=[]
        for plant in r["BCE_only_exact_spelling"]:
            if plant not in mp:raise RuntimeError("Missing WCVP mapping "+plant)
            mapped=mp[plant]
            category=mapped["match_status"]
            plant_id=mapped["accepted_plant_name_id"].strip()
            if category.startswith("UNIQUE_"):
                if not plant_id:raise RuntimeError("Missing accepted ID in matched plant")
                category=("SAME_ACCEPTED_ID_ALREADY_RECORDED"
                          if plant_id in canonical[butterfly] else
                          "ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY")
            counts[category]+=1
            details.append({"plant":plant,"match_status":category,
                "accepted_plant_name_id":plant_id or None,
                "accepted_name":mapped["accepted_name"]})
        checked.append({
            "species":butterfly,
            "BCE_only_original_binomials":len(r["BCE_only_exact_spelling"]),
            "wcvp_accepted_candidate_status":dict(counts),
            "original_same_spelling_links":r["intersection_exact_spelling"],
            "details":details,
        })
    pooled=Counter()
    for r in checked:pooled.update(r["wcvp_accepted_candidate_status"])
    result={
      "schema":"chocho_bce_clarke92_wcvp_accepted_id_gap_v0.1",
      "status":"EXPLORATORY_SECONDARY_SOURCE_ID_CROSSWALK_NOT_NEW_FIELD_DATA",
      "source":"Clarke 2024 evidence-1-to-3 derived on BCE, not original Dryad raw database",
      "frozen_HOSTS":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
      "frozen_WCVP":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
      "exact_focal_butterflies":92,
      "verified_butterflies":len(checked),
      "pooled_BCE_only_taxon_resolutions":dict(pooled),
      "per_butterfly":checked,
      "caveats":[
       "The BCE publication-derived checklist is not independently collected field data.",
       "Even a unique taxonomic ID absent from frozen HOSTS may reflect unequal geographic/evidence/host-use acceptance criteria.",
       "A unique WCVP accepted botanical name does not independently certify larval feeding, development or oviposition beyond source references.",
       "No geography and no new ecological mechanism estimated. No extension from European subset to full 239 species.",
       "Do not count unresolved names as false negatives or use this output to revise the primary GEB manuscript."
      ]
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("RECONCILIATION",json.dumps({"verified_butterflies":len(checked),
           "counts":dict(pooled)}),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=["prepare","audit"])
    p.add_argument("--input-json",required=True,type=Path)
    p.add_argument("--output-csv",type=Path)
    p.add_argument("--mapping-csv",type=Path)
    p.add_argument("--original-interaction-csv",type=Path)
    p.add_argument("--output-json",type=Path)
    a=p.parse_args()
    if a.mode=="prepare":
        if not a.output_csv:p.error("prepare requires --output-csv")
        prepare(a)
    else:
        if not a.mapping_csv or not a.original_interaction_csv or not a.output_json:
            p.error("audit requires --mapping-csv, --original-interaction-csv, --output-json")
        audit(a)

if __name__=="__main__":main()
