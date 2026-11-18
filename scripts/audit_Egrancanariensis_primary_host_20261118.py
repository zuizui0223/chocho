#!/usr/bin/env python3
"""Independent native-host rescue for Canary endemic Euchloe grancanariensis.

Evidence rules:
- Primary source 'Euchloe belemia grancanariensis' maps to independent
  Canary BIOTA species 'Euchloe grancanariensis' only as a named regional
  taxon concept, NOT the whole widespread E. belemia complex.
- Grancanariensis eggs deposited on Descurainia preauxiana (2014) are
  oviposition, not proven successful larval feeding on that plant.
- Acosta-Fernandez 2010 wild larvae on Hirschfeldia incana (2009) are
  direct larval feeding with captive pupation counts, not adult emergence.
- Hipparchia tamadabae 2016 field larvae on UNIDENTIFIED tussock grass
  identify only Poaceae; NEVER promote Hyparrhenia hirta to host.
"""
import csv,json
from pathlib import Path
from collections import Counter

HERE=Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
INPUT=HERE/"E_grancanariensis_H_tamadabae_source_evidence_20261118.csv"
BEFORE=HERE/"VERIFIED_GENERALIZATION_GATE_20261118.json"
BIO20=HERE/"canary_BIOTA_source_status_20.csv"
OUTPUT=Path("results/independent_endemic_species_host_rescue_20261118")

def load(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def run():
    raw=load(INPUT);old=json.loads(BEFORE.read_text(encoding="utf-8"))
    biota=load(BIO20)
    if len(raw)!=4 or len({r["event_id"] for r in raw})!=4:
        raise ValueError("Independent source documents changed")
    source_status={x["independent_checklist_name"]:x for x in biota}
    for name in ("Euchloe grancanariensis","Hipparchia tamadabae"):
        if source_status[name]["canary_origin_stratum"]!="endemic_confirmed":
            raise ValueError("Cannot reclassify independent Canary endemic origin")
    unknown=set(old["butterfly_taxon_ascertainment"]["literature_undocumented"]["contemporary_Canary_endemic"])
    if len(unknown)!=6 or unknown!={"Euchloe grancanariensis",
      "Hipparchia bacchus","Hipparchia gomera","Hipparchia tamadabae",
      "Hipparchia tilosi","Hipparchia wyssii"}:
        raise ValueError("Original literature-unknown independent endemic denominator changed")
    viable_old=old["butterfly_taxon_ascertainment"]["documentation_union_confirmed_min"]
    if viable_old!=35 or old["butterfly_taxon_ascertainment"]["extant_only_independent_denominator"]!=41:
        raise ValueError("Existing source-confirmed host-genus baseline changed")
    granc=[x for x in raw if x["butterfly"]=="Euchloe grancanariensis"]
    tama=[x for x in raw if x["butterfly"]=="Hipparchia tamadabae"]
    if len(granc)!=3 or len(tama)!=1:
        raise ValueError("Curated source event identity changed")
    if {x["host_plant"] for x in granc}!={"Descurainia preauxiana","Hirschfeldia incana"}:
        raise ValueError("Euchloe documented larval host species changed")
    direct=[x for x in granc if x["evidence_class"]=="field_larvae_on_named_plant_and_laboratory_pupation"]
    if len(direct)!=1 or direct[0]["host_plant"]!="Hirschfeldia incana" or direct[0]["observation_date"]!="2009-03":
        raise ValueError("Original 2010 primary source larval feeding evidence missing")
    n=direct[0]
    if tuple(map(int,(n["n_larvae_collected"],n["n_parasitized"],n["n_pupated"],n["n_pupae_failed"])))!=(15,9,5,2):
        raise ValueError("2010 collected larvae outcomes differ")
    if n["n_emerged_adults_confirmed"]!="unknown":
        raise ValueError("Original source did not prove adult emergence")
    oviposition=[x for x in granc if x["evidence_class"]=="direct_oviposition_and_egg_photo_named_plant"]
    if len(oviposition)!=1 or oviposition[0]["host_plant"]!="Descurainia preauxiana" or oviposition[0]["observation_date"]!="2014-01-21":
        raise ValueError("Original 2014 Descurainia oviposition evidence missing")
    field_grass=tama[0]
    if field_grass["plant_name_resolution"]!="family_only" or field_grass["host_plant"]!="Poaceae sp. undetermined":
        raise ValueError("Unknown grass illegitimately assigned host genus or species")
    if "Hyparrhenia hirta" not in field_grass["primary_evidence_note"]:
        raise ValueError("Explicit nearby non-host observed grass distinction lost")
    if any(not x["source_url"].startswith("https://") for x in raw):
        raise ValueError("Source URLs must remain inspectable")
    result={
        "schema":"chocho_island_endemic_exact_foodplants_native_rescue_20261118_v1",
        "source_scope":"Independently selected Canary 33 butterflies from BIOTA/INBO; source documents are targeted, not a uniform survey",
        "baseline":{
          "unique_Macaronesia_butterflies_n":42,
          "existing_species_with_named_larval_host_genus_min":35,
          "contemporary_extant_denominator":41,
          "Canary_endemic_total":13,
          "Canary_endemic_documented_genus_before":7,
        },
        "new_native_endemic_host_source":{
          "butterfly":"Euchloe grancanariensis",
          "regional_taxonomy":"Euchloe belemia grancanariensis in 2010 primary source; separately recognized as endemic species by 2026 Canary BIOTA",
          "island":"Gran Canaria",
          "named_host_species":[
            {"species":"Hirschfeldia incana","evidence":"2009 wild larvae feeding on fruits, independent primary source Acosta-Fernandez 2010","rearing_pupation_n":5},
            {"species":"Descurainia preauxiana","evidence":"2014 direct wild egg laying observed and photographed, not directly proven post-egg larval feeding"}
          ],
          "field_larvae_collected_2009_n":15,
          "parasitized_larvae_n":9,
          "unknown_cause_larva_death_n":1,
          "pupated_n":5,
          "of_pupated_malformed_and_died_n":2,
          "pupae_reaching_observed_summer_diapause_n":3,
          "adult_emergence_confirmed":False,
          "additional_butterfly_species_with_named_larval_host_genus":1,
        },
        "still_unidentified_grass":{
          "butterfly":"Hipparchia tamadabae",
          "island":"Gran Canaria",
          "field_larvae_year":2016,
          "field_foodplant_resolution":"Poaceae family only; tussock grass unidentified",
          "do_not_assume_larval_foodplant":"Hyparrhenia hirta: author explicitly reported adjacent but not hosting larvae",
          "additional_named_host_genus_species":0,
          "evidence":"Source Pyrgus H. wyssii complex explicitly identifies Gran Canaria taxon tamadabae"
        },
        "updated_confirmed_min_documentation":{
          "documented_butterflies_with_named_larval_host_genus_n":36,
          "independent_unique_species_denominator":42,
          "minimum_fraction":36/42,
          "documented_extant_n":36,
          "extant_only_denominator":41,
          "Canary_endemic_named_host_genus_n":8,
          "Canary_endemic_independently_classified_n":13,
          "still_no_named_larval_host_genus":[
            "Hipparchia bacchus","Hipparchia gomera","Hipparchia tamadabae",
            "Hipparchia tilosi","Hipparchia wyssii"
          ],
          "historical_extinct_excluded_from_current_resource_denominator":["Pieris wollastoni"],
          "historic_taxon_note":"Pieris wollastoni has separate historical host reports but is excluded from extant host-coverage estimates; no biological absence implied"
        },
        "scientific_boundary":[
          "Documented larval host genus remains a source-sampled lower bound, NOT unbiased true biological host knowledge prevalence.",
          "Hirschfeldia incana 2009 is direct wild larval foodplant with laboratory pupation; 2014 Descurainia preauxiana is direct oviposition, not demonstrated development.",
          "Observed parasitoid infection among fifteen source larvae is NOT a regional population parasitism rate.",
          "An unidentified tussock grass cannot be assigned Hyparrhenia, other grass genera, or native/introduced origin.",
          "No endemic butterfly population response to alien plant introduction or numerical mean host-resource geography growth established.",
          "Original GEB positive +54.9% and +118.9% host-opportunity estimates remain unchanged and refer only to original source-selected species."
        ]
    }
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT/"native_endemic_host_species_field_evidence.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("INDEPENDENT_CANARY_END_EMIC_NATIVE_HOSTS_RESOLVED",json.dumps({
        "confirmed_min":36,
        "out_of":42,
        "still_missing_five":result["updated_confirmed_min_documentation"]["still_no_named_larval_host_genus"],
        "collected":15,"parasitized":9,"pupated":5,"adult_emergence_proven":False
    },ensure_ascii=False))
    return result

if __name__=="__main__":
    run()
