#!/usr/bin/env python3
"""Audit source-documented contrasting butterfly host acceptance in 2010 paper.

One Fuerteventura larva accepted two Sisymbrium species and rejected offered
Hirschfeldia; 15 Gran Canaria larvae were *found in the wild* eating Hirschfeldia.
This is an observational field-vs-lab contrast, NEVER a replicated paired
host-choice experiment or population-level selection inference.
"""
import csv,json
from pathlib import Path
from collections import defaultdict

ROOT=Path("provenance/reviewer_defenses/results")
SOURCE=ROOT/"island_endemic_synonym_origin_20261118"/"Euchloe_two_island_host_use_field_vs_lab_2010.csv"
BIO=ROOT/"island_endemic_synonym_origin_20261118"/"canary_BIOTA_source_status_20.csv"
INBO=ROOT/"independent_island_butterfly_host_coverage_20261117"/"species_island_coverage.csv"
OUT=Path("results/Euchloe_island_host_comparison_20261118")

def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def main():
    source=read(SOURCE)
    statuses=read(BIO)
    inbo=read(INBO)
    allowed={(r["scientific_name"],r["island"]) for r in inbo}
    origin={r["independent_checklist_name"]:r["canary_origin_stratum"] for r in statuses}
    if len(source)!=4 or len({r["comparison_id"] for r in source})!=4:
        raise ValueError("Source evidence comparison changed")
    taxa={r["butterfly_2026_canary_species"] for r in source}
    if taxa!={"Euchloe grancanariensis","Euchloe hesperidum"}:
        raise ValueError("Cross-taxon comparison inadvertently changed")
    for e in source:
        name=e["butterfly_2026_canary_species"]
        if origin.get(name)!="endemic_confirmed":
            raise ValueError("Independent Canary endemic taxonomy not confirmed: "+name)
        if (name,"Canary Islands, "+e["source_island"]) not in allowed:
            raise ValueError("Source island not in original independent butterfly checklist")
        if "redalyc.org/pdf/455/45514996013.pdf" not in e["source_url"]:
            raise ValueError("Source not original 2010 Acosta-Fernandez")
        if e["observation_setting"] not in {"wild_larvae_found_on_named_plant","captive_larva_offered_host"}:
            raise ValueError("Unknown source observation setting")
    wild=[e for e in source if e["observation_setting"]=="wild_larvae_found_on_named_plant"]
    lab=[e for e in source if e["observation_setting"]=="captive_larva_offered_host"]
    if len(wild)!=1 or len(lab)!=3:
        raise ValueError("Wild and captive host evidence strata improperly merged")
    assert (wild[0]["butterfly_2026_canary_species"],wild[0]["offered_or_field_larval_host"],wild[0]["feeding_outcome"])==(
        "Euchloe grancanariensis","Hirschfeldia incana","feeding_confirmed_in_field")
    if int(wild[0]["n_field_larvae_collected"])!=15 or int(wild[0]["n_pupated_source_sample"])!=5:
        raise ValueError("Gran Canaria field larva/pupation source counts changed")
    if len({e["butterfly_2026_canary_species"] for e in lab})!=1 or any(
        e["butterfly_2026_canary_species"]!="Euchloe hesperidum" or
        e["sample_size_basis"]!="single_2009_larva" or
        int(e["n_field_larvae_collected"])!=1 for e in lab):
        raise ValueError("Three test plants are from ONE captive larva, not independent samples")
    actions={e["offered_or_field_larval_host"]:e["feeding_outcome"] for e in lab}
    if actions!={"Hirschfeldia incana":"rejected_offered_fruits",
                "Sisymbrium erysimoides":"accepted_feeding_fruits",
                "Sisymbrium irio":"accepted_feeding_fruits"}:
        raise ValueError("Host acceptance outcomes changed")
    res={
       "schema":"chocho_euchloe_host_use_conditional_field_vs_single_captive_larva_v1",
       "source":"Acosta-Fernandez 2010 SHILAP 38(149):121-128, source PDF DOI record https://www.redalyc.org/pdf/455/45514996013.pdf",
       "independent_butterfly_taxa":"Gran Canaria endemic E. grancanariensis and Fuerteventura endemic E. hesperidum, both recorded by INBO and independently BIOTA-endemic",
       "source_identity":"2010 article uses subspecies E. belemia grancanariensis and E. belemia hesperidum; do not project acceptance to all E. belemia",
       "grand_canaria":{
         "taxon":"Euchloe grancanariensis",
         "island":"Gran Canaria",
         "host":"Hirschfeldia incana",
         "field_larvae_found_on_host":15,
         "reported_parasitized":9,
         "reported_cause_unknown_death":1,
         "reported_pupated":5,
         "two_malformed_pupae_died":2,
         "reported_adult_emergences":None
       },
       "fuerteventura":{
         "taxon":"Euchloe hesperidum",
         "island":"Fuerteventura",
         "field_larvae_sample_n":1,
         "captive_offered_host_fruits":["Hirschfeldia incana","Sisymbrium erysimoides","Sisymbrium irio"],
         "accepted_hosts_in_captivity":["Sisymbrium erysimoides","Sisymbrium irio"],
         "rejected_host_in_captivity":["Hirschfeldia incana"],
         "pupated_after_feeding":True,
         "emerged_as_adult_confirmed":False
       },
       "inference":{
         "positive":"At least two island-endemic congeneric butterflies show different source-reported use of named local Brassicaceae species, with field feeding on H. incana in one taxon and its refusal in a single captive larva of the other.",
         "not_supported":"No replicated reciprocal host-choice experiment; cannot claim adaptive host preference divergence or general fitness comparison; a refusal by one larva is not a biological absence of H. incana host use across Fuerteventura populations.",
         "species_identity":"Host taxa must be matched as exact WCVP accepted species on Canary CNY; do not infer alien use from genus in the wider flora."
       }
    }
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"Euchloe_field_vs_lab_host_evidence.json").write_text(json.dumps(res,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("EUCHLOE_SOURCE_FIELD_VS_CULTURE_HOST_ACCEPTANCE_GUARDED",json.dumps({
        "Gran_Canaria_field_larvae":15,
        "Fuerteventura_captive_larvae":1,
        "lab_accepted_plant_species":2,
        "lab_rejected_plant_species":1,
        "replicated_experiment":False
    },ensure_ascii=False))

if __name__=="__main__":
    main()
