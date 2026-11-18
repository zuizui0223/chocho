#!/usr/bin/env python3
"""Replay independent living island-butterfly named larval host coverage from raw source crosswalks.

Unlike older summary-only receipts, this test rebuilds 36/41 from one row per
independently enumerated butterfly with individual source flags. Historical
extinct Pieris wollastoni is retained in table but excluded from extant rate;
its historical host reports must not be classified as absent biology.
"""
import csv,json
from pathlib import Path
from collections import defaultdict

ROOT=Path("provenance/reviewer_defenses/results")
BASE=ROOT/"independent_island_butterfly_host_coverage_20261117"
NEXT=ROOT/"island_endemic_synonym_origin_20261118"
FROZEN=BASE/"independent_42_species_vs_frozen239.csv"
OUTPUT=Path("results/source_merged_endemic_larval_host_20261118")

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def run():
    butterfly=read_csv(FROZEN)
    rescue=read_csv(NEXT/"independent_larval_host_rescue_10.csv")
    uniform=json.loads((NEXT/"uniform_source_larval_genera_24.json").read_text(encoding="utf-8"))["annotations"]
    alias=json.loads((NEXT/"VERIFIED_SOURCE_OUTCOME_20261118.json").read_text(encoding="utf-8"))["all_unique_42"]["newly_confirmed_alias_matches"]
    new=read_csv(NEXT/"E_grancanariensis_H_tamadabae_source_evidence_20261118.csv")
    canary=read_csv(NEXT/"canary_BIOTA_source_status_20.csv")+read_csv(NEXT/"canary_BIOTA_remaining_13.csv")
    origin={x["independent_checklist_name"]:x["canary_origin_stratum"] for x in canary}
    if (len(butterfly),len(rescue),len(uniform),len(alias),len(new),len(origin))!=(42,10,24,2,4,33):
        raise ValueError("Independent source denominator or literature provenance changed")
    if len({x["scientific_name"] for x in butterfly})!=42:
        raise ValueError("Independent checklist duplicate")
    unique_new={x["butterfly"] for x in new if x["plant_name_resolution"]=="species"
         and x["evidence_class"]=="field_larvae_on_named_plant_and_laboratory_pupation"}
    if unique_new!={"Euchloe grancanariensis"}:
        raise ValueError("A new larval host taxon must have direct wild larva evidence")
    literature={x["independent_checklist_species"] for x in rescue}
    if len(literature)!=10:raise ValueError("Original independent literature duplicates")
    rows=[]
    for old in sorted(butterfly,key=lambda r:r["scientific_name"]):
        name=old["scientific_name"]
        historical=int(name=="Pieris wollastoni")
        direct=int(old["hosts_named_genus"])
        synonymous=int(name in alias)
        targeted=int(name in literature)
        uniform_flag=int(name in uniform)
        new_field=int(name in unique_new)
        documented=int(bool(direct or synonymous or targeted or uniform_flag or new_field))
        row={
          "butterfly":name,"family":old["checklist_family"],
          "independent_archipelagos":old["archipelagos"],
          "Canary_independent_origin":origin.get(name,"not_Canary"),
          "Madeira_independent_origin":old["madeira_origin"],
          "historical_extinct_excluded_from_living_denominator":historical,
          "original_HOSTS_exact_host_genus":direct,
          "confirmed_Hosts_synonym_only":synonymous,
          "previous_targeted_external_named_genus":targeted,
          "uniform_Canary_source_named_genus":uniform_flag,
          "Euchloe_2010_primary_field_named_host_genus":new_field,
          "any_named_host_genus_in_audited_sources":documented,
          "living_species_named_host_genus_documented":int(documented and not historical)
        }
        if direct and synonymous:
            raise ValueError("Alias adds duplicate original strict match: "+name)
        if targeted and (direct or synonymous):
            raise ValueError("Targeted record was already source-matched: "+name)
        rows.append(row)
    if sum(r["original_HOSTS_exact_host_genus"] for r in rows)!=19:
        raise ValueError("Frozen exact HOSTS baseline changed")
    if sum(r["confirmed_Hosts_synonym_only"] for r in rows)!=2:
        raise ValueError("Frozen accepted taxon alias baseline changed")
    if sum(r["previous_targeted_external_named_genus"] for r in rows)!=10:
        raise ValueError("Original targeted additional literature mismatch")
    if sum(r["uniform_Canary_source_named_genus"] for r in rows)!=24:
        raise ValueError("Uniform 24 source-page host genus baseline changed")
    if sum(r["Euchloe_2010_primary_field_named_host_genus"] for r in rows)!=1:
        raise ValueError("New taxon rescue unexpected")
    living=[r for r in rows if not r["historical_extinct_excluded_from_living_denominator"]]
    missing=[r["butterfly"] for r in living if not r["living_species_named_host_genus_documented"]]
    expected=["Hipparchia bacchus","Hipparchia gomera","Hipparchia tamadabae","Hipparchia tilosi","Hipparchia wyssii"]
    if len(living)!=41 or sorted(missing)!=sorted(expected):
        raise ValueError("Five original Canary Hipparchia are not the only remaining extant host-genus unknowns: "+str(missing))
    living_documented=sum(r["living_species_named_host_genus_documented"] for r in rows)
    if living_documented!=36:
        raise ValueError("Named host-genus coverage of 41 living butterflies drifted")
    ca_endo=[r for r in living if r["Canary_independent_origin"]=="endemic_confirmed"]
    other=[r for r in ca_endo if not r["butterfly"].startswith("Hipparchia ")]
    hi=[r for r in ca_endo if r["butterfly"].startswith("Hipparchia ")]
    if (len(ca_endo),len(other),len(hi))!=(13,8,5):
        raise ValueError("Independently classified Canary endemics expected 8 non-Hipparchia and 5 Hipparchia")
    if sum(x["living_species_named_host_genus_documented"] for x in other)!=8:
        raise ValueError("Independent eight non-Hipparchia endemic taxa not all source-documented")
    if sum(x["living_species_named_host_genus_documented"] for x in hi)!=0:
        raise ValueError("A Hipparchia wild plant GENUS was invented")
    result={
        "schema":"chocho_independent_living_butterflies_larval_host_documentation_replay_v1",
        "original_independent_checklist_total":42,
        "extant_independent_butterflies_n":41,
        "extant_butterflies_with_named_larval_host_genus_in_audited_sources_n":living_documented,
        "extant_host_genus_documentation_fraction":living_documented/len(living),
        "Canary_independent_endemic_n":len(ca_endo),
        "Canary_independent_endemic_documented_n":8,
        "non_Hipparchia_Canary_independent_endemic_documented_n":8,
        "non_Hipparchia_Canary_independent_endemic_n":8,
        "Hipparchia_Canary_endemic_without_named_larval_host_genus_n":5,
        "Hipparchia_Canary_endemic_n":5,
        "historical_extinct_Pieris_wollastoni_retained_in_audit_but_excluded_from_current_estimate":True,
        "remaining_extant_unknowns":missing,
        "evidence_interpretation":[
            "Documented named larval-host genus after multiple nonuniform literature sources, NOT an unbiased biological true larval-host knowledge rate.",
            "2016 field Hipparchia tamadabae larvae on an unidentified Poaceae tussock grass do not yield a genus or species; in particular do not assign adjacent Hyparrhenia hirta.",
            "2014 Euchloe grancanariensis egg deposition on Descurainia preauxiana is oviposition, while 2009 wild caterpillars on Hirschfeldia incana are direct larval feeding.",
            "Original GEB239 source-selected potential host resource +54.9% globally and selected island +118.9% unchanged; new source evidence cannot extend the numerical estimand to all island endemics.",
            "All unresolved species in the current source union are a single related island radiation; do not estimate a general endemicity effect from the resulting post-hoc coverage groups."
        ]
    }
    OUTPUT.mkdir(parents=True,exist_ok=True)
    with (OUTPUT/"independent42_source_host_genus_flags_extant41.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUTPUT/"extant41_host_documentation_summary.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("SOURCE_VERIFIED_INDEPENDENT_LIVING_HOST_GENUS_DOC",json.dumps({
        "extant":len(living),"documented":living_documented,
        "canary_endemics":len(ca_endo),"non_Hipparchia_documented":8,
        "Hipparchia_genus_unknown":len(missing)
    },ensure_ascii=False))

if __name__=="__main__":
    run()
