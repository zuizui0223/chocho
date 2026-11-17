#!/usr/bin/env python3
"""Independently correct original GEB239 taxon selection for source-confirmed homotypic names.

The 42-species denominator and 239-species model population remain unchanged.
Never infer biological larval hosts or transfer a host association from a
widespread species to a separately recognized island endemic.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

BASE=Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
ORIGINAL=BASE/"independent_42_species_vs_frozen239.csv"
GEB=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
BIOTA=Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
TAXON_ALIAS={
    "Colias crocea":{
        "GEB239_frozen_name":"Colias croceus",
        "evidence":"NCBI Taxonomy ID 72248: Colias crocea is listed as a homotypic synonym of Colias croceus",
        "source_url":"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=72248",
        "status":"independently_verified_homotypic_synonym",
    },
}

def load(path):
    with path.open(encoding="utf-8-sig",newline="") as h:
        return list(csv.DictReader(h))

def audit(output):
    inp=load(ORIGINAL)
    model=load(GEB)
    origin=load(BIOTA/"canary_BIOTA_source_status_20.csv")+load(BIOTA/"canary_BIOTA_remaining_13.csv")
    if len(inp)!=42 or len(model)!=239 or len(origin)!=33:
        raise ValueError("Original source populations changed")
    if len({r["species"] for r in model})!=239:
        raise ValueError("Duplicate 239 model species name")
    frozen_names={r["species"] for r in model}
    if not {"Colias croceus"}<=frozen_names or "Colias crocea" in frozen_names:
        raise ValueError("Expected original GEB names absent or already corrected")
    classes={r["independent_checklist_name"]:r["canary_origin_stratum"] for r in origin}
    assert len(classes)==33
    rows=[]
    for r in inp:
        name=r["scientific_name"]
        exact=int(r["in_frozen239"])
        if exact!=int(name in frozen_names):
            raise ValueError("Original exact model selection differs from frozen audit: "+name)
        verified=TAXON_ALIAS.get(name)
        aliased=(verified["GEB239_frozen_name"] if verified else "")
        corrected=int(bool(exact or (aliased and aliased in frozen_names)))
        if name in ("Pararge xiphioides","Hipparchia bacchus","Hipparchia gomera",
                    "Pararge xiphia","Thymelicus christi"):
            if corrected:
                raise ValueError("Taxonomic concept conflicts accidentally imported into GEB population")
        rows.append({
            "independent_checklist_species":name,
            "family":r["checklist_family"],
            "archipelagos":r["archipelagos"],
            "madeira_origin":r["madeira_origin"],
            "canary_origin":classes.get(name,"not_in_canary"),
            "in_GEB239_exact_binomial":exact,
            "in_GEB239_confirmed_homotypic":corrected,
            "matched_GEB239_name":name if exact else (aliased if corrected else ""),
            "source_authority":verified["source_url"] if verified else "",
            "identification_status":"homotypic_taxon_alias_verified" if not exact and corrected else (
                "exact" if exact else "unmatched_not_absent"
            ),
        })
    exact=sum(r["in_GEB239_exact_binomial"] for r in rows)
    confirmed=sum(r["in_GEB239_confirmed_homotypic"] for r in rows)
    if (exact,confirmed)!=(12,13):
        raise ValueError(f"Unexpected source-reconciled original239 overlap: {exact},{confirmed}")
    def coverage(group):
        return {
            "n":len(group),
            "exact_original_239":sum(r["in_GEB239_exact_binomial"] for r in group),
            "source_verified_homotypic_239":sum(r["in_GEB239_confirmed_homotypic"] for r in group),
        }
    canary=[r for r in rows if r["canary_origin"]!="not_in_canary"]
    madeira_native=[r for r in rows if r["madeira_origin"]=="native_nonendemic"]
    madeira_endo=[r for r in rows if r["madeira_origin"]=="island_endemic_species"]
    canary_endo=[r for r in canary if r["canary_origin"]=="endemic_confirmed"]
    canary_native=[r for r in canary if r["canary_origin"].startswith("native_nonendemic")]
    result={
        "schema":"chocho_independent42_originalGEB239_taxon_identity_lowerbound_v2",
        "source_denom":"INBO original independent 42 butterflies, not selected from HOSTS",
        "model239_n":239,
        "source_exact_original239_n":exact,
        "source_confirmed_homotypic_original239_min":confirmed,
        "corrected_entries":[r for r in rows if r["in_GEB239_confirmed_homotypic"] and not r["in_GEB239_exact_binomial"]],
        "cohorts":{
            "all_independent42":coverage(rows),
            "Canary_all33":coverage(canary),
            "Canary_endemic13":coverage(canary_endo),
            "Canary_native_nonendemic19":coverage(canary_native),
            "Madeira_native_nonendemic7":coverage(madeira_native),
            "Madeira_endemic3":coverage(madeira_endo),
        },
        "source_verified_aliased_names":TAXON_ALIAS,
        "limitations":[
            "Only independently verified homotypic synonym links are merged for the original model-239 selection audit; 13/42 remains a lower bound rather than a fully reconciled complete name coverage.",
            "Accepted-name collisions between island endemic species and broad species are explicitly prohibited.",
            "Actual model239 sample-selection reasons (environment, plant taxonomy, trait completeness) for all other not-selected taxa remain untested.",
            "No new plant or butterfly colonization result is implied. GEB original 239-species estimated host-opportunity +54.9% unchanged."
        ]
    }
    if result["cohorts"]["Canary_endemic13"]["source_verified_homotypic_239"]!=0:
        raise ValueError("Canary endemics unexpectedly entered original GEB model")
    if result["cohorts"]["Canary_native_nonendemic19"]["source_verified_homotypic_239"]!=10:
        raise ValueError("Native source-matched model139 overlap should be 10/19")
    if result["cohorts"]["Madeira_native_nonendemic7"]["source_verified_homotypic_239"]!=7:
        raise ValueError("Madeira native 7/7 homotypic reconciliation failed")
    output.mkdir(parents=True,exist_ok=True)
    with (output/"independent42_original239_homotypic_selection.csv").open("w",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    (output/"original239_homotypic_selection_summary.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    print("ORIGINAL_239_INDEPENDENT_NAME_AUDIT_PASSED")
    print("ORIGINAL_239_INDEPENDENT_NAME_SUMMARY="+json.dumps({
        "exact":exact,"homotypic_confirmed_min":confirmed,
        "cohorts":result["cohorts"]
    },ensure_ascii=False))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    audit(a.out_dir)
