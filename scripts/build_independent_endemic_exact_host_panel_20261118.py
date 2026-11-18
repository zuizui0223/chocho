#!/usr/bin/env python3
"""Build a source-typed, independent island-endemic butterfly × exact plant-species panel.

Independently defined butterfly denominator comes from INBO + Canary BIOTA +
Madeira Borges 2008 + 2023 official Azores source. Do not choose butterflies
based on whether a resource expansion is positive. Panel is a *partially
ascertained subset* of 18 independently classified endemic butterflies, not a
representative 18-species sample. Species-level host names are source-typed.

The original NHM HOSTS gives only one non-local Festuca jubata record for each
Azores endemic. All other links are sourced from frozen independent literature
evidence CSVs. Host statuses and local presence are checked separately by WCVP.
"""
import argparse,csv,json,re
from pathlib import Path

BASE=Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
HERE=Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
INBO=BASE/"independent_42_species_vs_frozen239.csv"
CANARY20=HERE/"canary_BIOTA_source_status_20.csv"
CANARY13=HERE/"canary_BIOTA_remaining_13.csv"
Focal=HERE/"independent_butterfly_host_plant_species_events_v1.csv"
EGRAN=HERE/"E_grancanariensis_H_tamadabae_source_evidence_20261118.csv"
EHESP=HERE/"Euchloe_two_island_host_use_field_vs_lab_2010.csv"
AZORES={"Hipparchia azorina","Hipparchia miguelensis"}
VALID_GRADES={"high","medium","species"}
# Tier 1 uses reported actual larvae (wild or garden) and independently
# reported larval-host species (not only oviposition/captive acceptance);
# Tier 2 extends to oviposition-only/captive-only evidence.
TIER_FOR={
    "wild_larva_photo":"larval_reported",
    "garden_larva_photo":"larval_reported",
    "wild_eggs_larvae_report":"larval_reported",
    "wild_larva_report":"larval_reported",
    "field_host_report":"larval_reported",
    "literature_host_report":"larval_reported",
    "curated_host_report":"larval_reported",
    "field_larvae_on_named_plant_and_laboratory_pupation":"larval_reported",
    "direct_oviposition_and_egg_photo_named_plant":"oviposition_only",
    "accepted_feeding_fruits":"captivity_only",
    "NHM_original_nonlocal":"nonlocal_database",
}

def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
def check_binomial(name):
    if not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][a-z-]+",name or ""):
        raise ValueError("Not an exact species-level binomial: "+str(name))
def key(status):
    if status=="endemic_confirmed":return "CNY"
    return "AZO"

def run(hosts_path,out_dir):
    butterflies=read(INBO)
    canary=read(CANARY20)+read(CANARY13)
    if len(butterflies)!=42 or len(canary)!=33 or len({r["independent_checklist_name"] for r in canary})!=33:
        raise ValueError("Frozen independent butterfly origins changed")
    can_names={r["independent_checklist_name"] for r in canary if r["canary_origin_stratum"]=="endemic_confirmed"}
    md_names={r["scientific_name"] for r in butterflies if r["madeira_origin"]=="island_endemic_species"}
    endemic=can_names|md_names|AZORES
    if (len(can_names),len(md_names),len(endemic))!=(13,3,18):
        raise ValueError("Original independently classified island endemics changed")
    entries=[]
    def add(butterfly,plant,region,grade,source,scope,locator):
        if butterfly not in endemic:raise ValueError("Non-endemic butterfly in controlled panel: "+butterfly)
        if region not in ("CNY","MDR","AZO"):raise ValueError("Not an independently classified source region")
        if butterfly in AZORES and region!="AZO":raise ValueError("Azores source origin mismatch")
        if butterfly in md_names and region!="MDR":raise ValueError("Madeira origin mismatch")
        if butterfly in can_names and region!="CNY":raise ValueError("Canary origin mismatch")
        check_binomial(plant)
        if grade not in TIER_FOR:raise ValueError("Unverified host evidence stage")
        entries.append({"butterfly":butterfly,"plant_species":plant,"source_region_TDWG3":region,
            "source_evidence":grade,"evidence_tier":TIER_FOR[grade],
            "literature_url":source,"source_scope_note":scope,
            "observed_island":locator,"in_original_GEB239":0})
    for row in read(Focal):
        tax=row["butterfly_name"]
        if tax not in endemic or row["name_resolution_grade"] not in ("high","medium"):
            continue
        if row["source_evidence_type"] in ("historical_breeding","historical_host_report"):
            continue
        if not row["reported_host_binomial"].strip():continue
        region="MDR" if tax in md_names else "CNY"
        add(tax,row["reported_host_binomial"],region,row["source_evidence_type"],
            row["butterfly_host_source_url"],row["scope_note"],row["named_island"])
    for row in read(EGRAN):
        if row["butterfly"]!="Euchloe grancanariensis":continue
        if row["plant_name_resolution"]!="species":continue
        if row["evidence_class"]=="independent_secondary_hostplant_listing":continue
        add(row["butterfly"],row["host_plant"],"CNY",row["evidence_class"],
            row["source_url"],row["primary_evidence_note"],row["island"])
    for row in read(EHESP):
        if row["butterfly_2026_canary_species"]!="Euchloe hesperidum" or row["feeding_outcome"]!="accepted_feeding_fruits":
            continue
        add("Euchloe hesperidum",row["offered_or_field_larval_host"],"CNY",
            row["feeding_outcome"],row["source_url"],row["evidence_boundary"],row["source_island"])
    with hosts_path.open(encoding="utf-8-sig",newline="") as f:
        host_records=list(csv.DictReader(f))
    if len(host_records)!=140485:raise ValueError("Pinned original NHM HOSTS count changed")
    for butterfly in sorted(AZORES):
        genus,epithet=butterfly.split(" ")
        matches=[r for r in host_records
           if r.get("Insect Genus","")==genus and r.get("Insect Species","")==epithet
           and r.get("Insect Family","")=="Nymphalidae"]
        if len(matches)!=1 or (matches[0].get("Hostplant Genus"),matches[0].get("Hostplant Species"))!=("Festuca","jubata"):
            raise ValueError("Two Azores endemic host plant species differ in NHM HOSTS")
        add(butterfly,"Festuca jubata","AZO","NHM_original_nonlocal",
            "https://github.com/globalbioticinteractions/HOSTS/tree/808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "Original HOSTS butterfly host identity only; no Azores island field larvae, location or performance",
            "")
    # Distinct plant-by-butterfly evidence; an identical repeated plant
    # recorded in 2 islands counts once for global WCVP potential geography.
    groups={}
    rank={"wild_larva_photo":1,"wild_eggs_larvae_report":1,
          "field_larvae_on_named_plant_and_laboratory_pupation":1,
          "garden_larva_photo":2,"literature_host_report":3,
          "field_host_report":3,"curated_host_report":4,
          "NHM_original_nonlocal":5,"direct_oviposition_and_egg_photo_named_plant":6,
          "accepted_feeding_fruits":7}
    for row in entries:
        k=(row["butterfly"],row["plant_species"])
        if k not in groups or rank[row["source_evidence"]]<rank[groups[k]["source_evidence"]]:
            groups[k]=row.copy()
            groups[k]["all_evidence_source_islands"]={row["observed_island"]} if row["observed_island"] else set()
        else:
            if row["observed_island"]:groups[k]["all_evidence_source_islands"].add(row["observed_island"])
    panel=[]
    for _,row in sorted(groups.items()):
        row["all_evidence_source_islands"]="; ".join(sorted(row["all_evidence_source_islands"]))
        panel.append(row)
    # Eligible source-typed panel: 8 independently named butterflies (2 AZO,
    # 1 Madeira, 5 Canary), not all 18 endemic species.
    expected={
      "Hipparchia azorina":{"Festuca jubata"},
      "Hipparchia miguelensis":{"Festuca jubata"},
      "Gonepteryx maderensis":{"Rhamnus glandulosa"},
      "Pieris cheiranthi":{"Tropaeolum majus","Crambe strigosa","Descurainia millefolia"},
      "Gonepteryx cleobule":{"Rhamnus crenulata","Rhamnus glandulosa"},
      "Euchloe grancanariensis":{"Hirschfeldia incana","Descurainia preauxiana"},
      "Euchloe hesperidum":{"Sisymbrium erysimoides","Sisymbrium irio"},
      "Pararge xiphioides":{"Brachypodium sylvaticum"}
    }
    by_butterfly={}
    for row in panel:by_butterfly.setdefault(row["butterfly"],set()).add(row["plant_species"])
    if by_butterfly!=expected:
        raise ValueError("Unexpected source-documented species panel: "+repr(by_butterfly))
    if len(panel)!=13:raise ValueError("Expected 13 butterfly×host-plant species edges")
    if sum(x["evidence_tier"]=="nonlocal_database" for x in panel)!=2:
        raise ValueError("Original 2 Azores HOSTS-only source edges not isolated")
    out_dir.mkdir(parents=True,exist_ok=True)
    with (out_dir/"eight_endemic_source_host_species13.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(panel[0]))
        writer.writeheader();writer.writerows(panel)
    summary={"schema":"chocho_independent8_endemic_exact_species_host_source_panel_v1",
        "original_independent_endemic_denominator":18,
        "source_selected_exact_host_endemic_species":len(expected),
        "source_verified_host_edges":len(panel),
        "data_coverage_fraction_selected":len(expected)/18,
        "reasons_for_excluding_other_ten":"No accepted exact larval host plant species in the source-typed data examined; not proof absent hosts",
        "distinct_host_plants":len({r["plant_species"] for r in panel}),
        "source_known_butterflies":sorted(expected),
        "all_endemic_species":sorted(endemic),
        "source_region_correct":"AZO, MDR and CNY are WGSRPD3 archipelagos, NOT individual islands",
        "sensitivity_tiers":{
            "larval_reported":"Larva observed or species-level literature host reported",
            "oviposition_only":"Egg/oviposition on plant; no larval success",
            "captivity_only":"Offered-host laboratory feeding on one larva",
            "nonlocal_database":"NHM HOSTS foodplant record lacking island field context"},
        "no_performance_or_global_colonization_inferred":True}
    (out_dir/"eight_endemic_host_source_panel_summary.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    print("INDEPENDENT_EIGHT_SOURCE_TYPED_EXACT_HOST_PLANTS="+json.dumps({
        "butterflies":len(expected),"host_edges":len(panel),"host_plant_species":summary["distinct_host_plants"]
    },ensure_ascii=False))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    run(args.hosts,args.out_dir)
