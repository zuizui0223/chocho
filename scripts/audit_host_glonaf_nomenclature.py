#!/usr/bin/env python3
"""Diagnose name-token mismatch between fixed239 HOSTS and GloNAF plant fields.

No ecological effect, no row data, no speculative synonym matching.
Prints a small fixed set of public botanical names as a diagnostic, and
saves summary counts without rehosting original source datasets.
"""
from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path
import zipfile

from analyze_documented_host_specificity import focal_hosts,binomial
from build_lenzner_tdwg_colonial_crosswalk import load_csv_member

SOURCE=Path("_temp_historical_sources/lenzner_6640264.zip")
OUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/"
         "HOSTS_GloNAF_name_alignment_diagnosis.json")


def first_botanical_name(genus,species):
    genus=(genus or "").strip()
    species=(species or "").strip()
    if not genus or not species:
        return None
    tokens=species.split()
    if not tokens:
        return None
    if tokens[0].casefold()==genus.casefold():
        if len(tokens)<2:
            return None
        epithet=tokens[1]
    else:
        epithet=tokens[0]
    if not epithet[0].isalpha() or not all(k.isalpha() or k=="-" for k in epithet):
        return None
    return (genus+" "+epithet).casefold()


def diagnostic():
    fixed,n=focal_hosts()
    with zipfile.ZipFile(SOURCE) as z:
        allrows,damage=load_csv_member(z,"GloNAF.csv")
    old=set()
    new=set()
    source_species_toks=Counter()
    source_genus_toks=Counter()
    examples=[]
    for r in allrows:
        g=r.get("genus","")
        e=r.get("Species","")
        x=binomial(g,e)
        y=first_botanical_name(g,e)
        if x:old.add(x)
        if y:new.add(y)
        source_species_toks[min(5,len((e or "").strip().split()))]+=1
        source_genus_toks[min(5,len((g or "").strip().split()))]+=1
        if len(examples)<12 and e and g:
            examples.append({"genus":g[:70],"Species":e[:140],
                             "name":str(r.get("name",""))[:70],
                             "joined_exact":x,"joined_first_epithet":y})
    both=sorted(old&fixed)
    truncated=sorted(new&fixed)
    unmatched_known=sorted(fixed-new)[:10]
    return {
        "status":"SOURCE_NOMENCLATURE_DIAGNOSTIC_NOT_HOST_SPECIFICITY_EFFECT",
        "raw_glonaf_rows":len(allrows),
        "known_fixed_239_HOSTS":n,
        "original_glo_names":len(old),
        "first_binomial_glo_names":len(new),
        "exact_GloNAF_vs_HOSTS_overlap":len(both),
        "first_binomial_GloNAF_vs_HOSTS_overlap":len(truncated),
        "species_field_token_distribution":dict(source_species_toks),
        "genus_field_token_distribution":dict(source_genus_toks),
        "public_species_name_examples":examples,
        "first_overlap_examples":truncated[:10],
        "fixed_host_not_found_examples":unmatched_known,
        "source_text_replacements":damage,
        "do_not_infer_no_hosts_from_no_name_match":True,
    }


def main():
    r=diagnostic()
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(r,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
