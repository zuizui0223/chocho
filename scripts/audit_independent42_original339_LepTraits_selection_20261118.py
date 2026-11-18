#!/usr/bin/env python3
"""Why independent Macaronesian endemic butterflies are absent from frozen GEB239.

Source denominator is 42 INBO species. Read exact original 339-source manifest
and its S1 resource descriptors, actual source LepTraits, and frozen GEB239.
A missing species from an *already frozen source panel* is not an ecological
absence of larval host; trait/source-panel selection causes external-validity
uncertainty before any host-resource geographic effect was estimated.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import defaultdict

BASE=Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
HERE=Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
INBO=BASE/"independent_42_species_vs_frozen239.csv"
CAN20=HERE/"canary_BIOTA_source_status_20.csv"
CAN13=HERE/"canary_BIOTA_remaining_13.csv"
S1=Path("data/frozen/s1_species_manifest.json")
DESC=Path("data/frozen/figure_sources/s1_resource_descriptors.csv")
GEB=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
LEP=Path("data/external/leptraits_consensus_v1.0.csv")
LEP_SHA="6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b"
DES_SHA="894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc"
ALIASES={
  "Colias crocea":{"source_name":"Colias croceus","source":"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=72248"},
  "Hipparchia azorina":{"source_name":"Hipparchia azorinus","source":"https://www.sciencedirect.com/org/science/article/pii/S131329891800040X"},
}
def csvrows(p):
    with p.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def nz(v):
    return str(v or "").strip() not in ("","NA","nan","None")
def positive_number(value):
    if not nz(value):return False
    try:return float(value)>0 and math.isfinite(float(value))
    except(ValueError,TypeError):return False
def main(output):
    original=csvrows(INBO)
    canary=csvrows(CAN20)+csvrows(CAN13)
    traits=csvrows(LEP)
    desc=csvrows(DESC)
    geb=csvrows(GEB)
    source=json.loads(S1.read_text(encoding="utf-8"))
    if sha(LEP)!=LEP_SHA or sha(DESC)!=DES_SHA:
        raise ValueError("Original LepTraits or original descriptor SHA256 diverged")
    if (len(original),len(canary),len(desc),len(geb))!=(42,33,339,239):
        raise ValueError("Original denominators changed")
    if len(set(r["scientific_name"] for r in original))!=42:
        raise ValueError("Duplicate independent butterfly name")
    if source.get("schema")!="ttf_relational_prior_S1_species_exclusion_v0.1" or len(source.get("species",[]))!=339:
        raise ValueError("Original frozen 339 species source selection manifest changed")
    dnames={r["species"] for r in desc}
    mnames={r["species"] for r in geb}
    srcnames=set(source["species"])
    if len(dnames)!=339 or dnames!=srcnames:
        raise ValueError("Original 339 descriptor list differs from actual source-manifest species")
    if len(mnames)!=239 or not mnames<=dnames:
        raise ValueError("Original GEB 239 must be a subset of the source 339 butterfly panel")
    if "Species" not in traits[0]:
        raise ValueError("Original source LepTraits lacks expected Species key")
    tnames={r["Species"]:r for r in traits if r.get("Species")}
    if not dnames<=set(tnames):
        raise ValueError("339 original source species unexpectedly missing from exact LepTraits")
    origin={r["independent_checklist_name"]:r["canary_origin_stratum"] for r in canary}
    if len(origin)!=33:raise ValueError("Original independent Canary source statuses missing")
    rows=[]
    for r in sorted(original,key=lambda x:x["scientific_name"]):
        name=r["scientific_name"]
        alt=ALIASES.get(name,{})
        candidate_names=[name]+([alt["source_name"]] if alt else [])
        choices=[x for x in candidate_names if x in tnames]
        source339=[x for x in candidate_names if x in srcnames]
        orig239=[x for x in candidate_names if x in mnames]
        if len(source339)>1 or len(orig239)>1:
            raise ValueError("Independent accepted name maps to multiple source-panel species: "+name)
        assert len(source339)==0 or len(choices)>=1
        trait=tnames[choices[0]] if choices else None
        raw_fields={
          "NumberOfHostplantFamilies":positive_number(trait.get("NumberOfHostplantFamilies")) if trait else False,
          "Voltinism":nz(trait.get("Voltinism")) if trait else False,
          "CanopyAffinity":nz(trait.get("CanopyAffinity")) if trait else False,
          "EdgeAffinity":nz(trait.get("EdgeAffinity")) if trait else False,
          "MoistureAffinity":nz(trait.get("MoistureAffinity")) if trait else False,
          "DisturbanceAffinity":nz(trait.get("DisturbanceAffinity")) if trait else False,
          "WingSize":any(positive_number(trait.get(x)) for x in (
              "WS_L","WS_U","FW_L","FW_U","WS_L_Fem","WS_U_Fem",
              "WS_L_Mal","WS_U_Mal","FW_L_Fem","FW_U_Fem","FW_L_Mal","FW_U_Mal"
          )) if trait else False,
        }
        absent_fields=sorted(k for k,v in raw_fields.items() if not v) if trait else []
        traitcomplete=int(bool(trait and not absent_fields))
        rows.append({
          "independent_species":name,
          "family":r["checklist_family"],
          "canary_origin":origin.get(name,"not_in_Canary"),
          "madeira_origin":r["madeira_origin"],
          "archipelagos":r["archipelagos"],
          "has_original_HOSTS_exact_genus":int(r["hosts_named_genus"]),
          "in_original_LepTraits_any_name":int(bool(choices)),
          "source_LepTraits_matched_name":choices[0] if choices else "",
          "source_LepTraits_major_complete":traitcomplete,
          "source_LepTraits_missing_major_fields":";".join(absent_fields) if trait else "species_name_absent_from_LepTraits",
          "source_LepTraits_raw_host_families":trait.get("NumberOfHostplantFamilies","") if trait else "",
          "in_original_frozen_339_S1_species_panel":int(bool(source339)),
          "in_original_frozen_239_final_resource_sample":int(bool(orig239)),
          "source_panel_matching_name":source339[0] if source339 else "",
          "model_final_matching_name":orig239[0] if orig239 else "",
          "nomenclatural_source_url":alt.get("source",""),
          "sampling_stage":"selected_239" if orig239 else
             "source339_not_final239" if source339 else
             "raw_leptraits_named_but_not_selected339" if choices else
             "absent_from_original_leptraits_under_approved_names"
        })
    if sum(x["in_original_frozen_339_S1_species_panel"] for x in rows)!=15:
        raise ValueError("Original 42 into 339 source panel expected 15")
    if sum(x["in_original_frozen_239_final_resource_sample"] for x in rows)!=13:
        raise ValueError("Original 42 into final239 expected 13")
    if sum(x["has_original_HOSTS_exact_genus"] for x in rows)!=19:
        raise ValueError("Original source-STRICT HOSTS 19 count changed")
    if sum(x["source_LepTraits_major_complete"] for x in rows if x["in_original_frozen_339_S1_species_panel"])!=15:
        raise ValueError("Frozen original 339 source descriptors were major-complete")
    endemic_canary=[x for x in rows if x["canary_origin"]=="endemic_confirmed"]
    endemic_madeira=[x for x in rows if x["madeira_origin"]=="island_endemic_species"]
    endemic_azores=[x for x in rows if x["independent_species"] in ("Hipparchia azorina","Hipparchia miguelensis")]
    all_endemics={x["independent_species"]:x for x in endemic_canary+endemic_madeira+endemic_azores}
    if len(endemic_canary)!=13 or len(endemic_madeira)!=3 or len(endemic_azores)!=2 or len(all_endemics)!=18:
        raise ValueError("Independent 18 island endemic butterfly species changed")
    if any(r["in_original_frozen_339_S1_species_panel"] for r in all_endemics.values()):
        raise ValueError("One island endemic was already included in the frozen 339; cannot claim zero")
    counter=defaultdict(int)
    for r in rows:
        counter[r["sampling_stage"]]+=1
    def describe(group):
        return {
            "n":len(group),
            "in_original_LepTraits":sum(r["in_original_LepTraits_any_name"] for r in group),
            "LepTraits_major_complete":sum(r["source_LepTraits_major_complete"] for r in group),
            "in_339":sum(r["in_original_frozen_339_S1_species_panel"] for r in group),
            "in_239":sum(r["in_original_frozen_239_final_resource_sample"] for r in group),
            "original_exact_HOSTS_genus":sum(r["has_original_HOSTS_exact_genus"] for r in group),
        }
    missing_major_by_group={}
    for label,group in (
        ("all42",rows),
        ("all18_endemics",list(all_endemics.values())),
        ("Canary13_endemics",endemic_canary),
        ("Madeira3_endemics",endemic_madeira),
        ("Azores2_endemics",endemic_azores),
    ):
        matched=[r for r in group if r["in_original_LepTraits_any_name"]]
        missing_major_by_group[label]={
            "source_taxa_with_LepTraits":len(matched),
            "missing_major_field_counts":{
                name:sum(name in r["source_LepTraits_missing_major_fields"].split(";") for r in matched)
                for name in ("NumberOfHostplantFamilies","Voltinism","CanopyAffinity",
                             "EdgeAffinity","MoistureAffinity","DisturbanceAffinity","WingSize")
            },
            "complete_major_within_LepTraits":sum(x["source_LepTraits_major_complete"] for x in matched),
            "taxon_level_gaps":[
                {"name":x["independent_species"],"source_taxon":x["source_LepTraits_matched_name"],
                 "missing_fields":x["source_LepTraits_missing_major_fields"]}
                for x in matched if x["source_LepTraits_missing_major_fields"]
            ]
        }
    summary={
      "schema":"chocho_independent42_originalGEB_339_stage_vs_LepTraits_audit_v2",
      "frozen_sources":{"s1_species_count":339,"s1_schema":source["schema"],
            "descriptor_SHA256":DES_SHA,"leptraits_SHA256":LEP_SHA,
            "original_final_GEB239_rows":len(mnames)},
      "independent_checklist_n":42,
      "source339_overlap_n":15,
      "final239_overlap_n":13,
      "independent_endemic_n":18,
      "independent_endemic_source339_overlap_n":0,
      "independent_endemic_final239_overlap_n":0,
      "LepTraits_total_distinct_taxa":len(tnames),
      "sampling_stages":dict(counter),
      "major_trait_missingness_by_group":missing_major_by_group,
      "cohorts":{
        "independent42":describe(rows),
        "all_18_independent_island_endemics":describe(list(all_endemics.values())),
        "Canary_13_endemics":describe(endemic_canary),
        "Madeira_3_endemics":describe(endemic_madeira),
        "Azores_2_endemics":describe(endemic_azores),
        "Canary_19_natives":describe([r for r in rows if r["canary_origin"].startswith("native_nonendemic")]),
        "Madeira_7_natives":describe([r for r in rows if r["madeira_origin"]=="native_nonendemic"]),
        "Azores_6_natives":describe([r for r in rows if r["independent_species"] in (
             "Colias crocea","Danaus plexippus","Lampides boeticus","Pieris brassicae",
             "Vanessa atalanta","Vanessa cardui")]),
      },
      "interpretation":[
        "No named island endemic butterfly in the 18 independently classified taxa was in the pre-existing frozen 339 species panel, before the final resource effect calculation. This is a source-population selection boundary, not evidence of no hosts.",
        "The 339 list is a prior response-blind fixed source panel of major-complete butterflies; its genesis is not identified by checking only host databases.",
        "If a taxon is present in raw LepTraits but absent from 339, source-panel design rather than missing LepTraits name explains that exclusion stage. Do NOT assume physiological ecology or resource shortage.",
        "All original GEB +54.9% global and +118.9% sampled island known-host potential opportunity effects remain unchanged for their original selected population.",
        "Some island statuses differ by source/dates. Archipelago×species cohorts are not independent biological evolution replicates."
      ]
    }
    output.mkdir(parents=True,exist_ok=True)
    with (output/"independent42_339_239_leptraits_selection.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (output/"independent42_339_239_selection_summary.json").write_text(
       json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("INDEPENDENT42_GEB_SELECTION_339_239_LEPTRAITS_AUDIT="+json.dumps({
       "overall":summary["cohorts"]["independent42"],
       "endemic":summary["cohorts"]["all_18_independent_island_endemics"],
       "stage_counts":summary["sampling_stages"],
       "endemic_raw_trait_missingness":summary["major_trait_missingness_by_group"]["all18_endemics"],
       "missing_original339_but_LepTraits": [
          r["independent_species"] for r in rows if r["in_original_LepTraits_any_name"]
          and not r["in_original_frozen_339_S1_species_panel"]]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    main(a.out_dir)
