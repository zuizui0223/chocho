#!/usr/bin/env python3
"""Quantify the island endemic taxa recoverable by relaxing unrelated habitat-trait completeness.

Original GEB prior S1 339 panel required full LepTraits axes; endpoint-specific
host-only questions do not logically require canopy/edge/moisture/disturbance.
This is eligibility under a *new* data design, not a corrected 339 or original
GEB statistical estimate, and not a guarantee of resolved plant host SPECIES.
"""
import csv,hashlib,json,math
from pathlib import Path
from collections import Counter
ROOT=Path("provenance/reviewer_defenses/results")
INBO=ROOT/"independent_island_butterfly_host_coverage_20261117"/"independent_42_species_vs_frozen239.csv"
CAN20=ROOT/"island_endemic_synonym_origin_20261118"/"canary_BIOTA_source_status_20.csv"
CAN13=ROOT/"island_endemic_synonym_origin_20261118"/"canary_BIOTA_remaining_13.csv"
SOURCE=Path("data/external/leptraits_consensus_v1.0.csv")
MANIFEST=Path("data/frozen/s1_species_manifest.json")
LEP_SHA="6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b"
OUT=Path("results/independent42_relaxed_host_only_traits_20261118")
ALIASES={"Colias crocea":"Colias croceus","Hipparchia azorina":"Hipparchia azorinus"}

def rows(path):
  with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def positive(v):
  try:return math.isfinite(float(v)) and float(v)>0
  except(TypeError,ValueError):return False
def present(v):
  return str(v or "").strip() not in ("","NA","None","nan")
def evaluate():
  if hashlib.sha256(SOURCE.read_bytes()).hexdigest()!=LEP_SHA:
    raise ValueError("Original LepTraits source digest drift")
  original=rows(INBO)
  biota=rows(CAN20)+rows(CAN13)
  traits=rows(SOURCE)
  source=json.loads(MANIFEST.read_text())
  if len(original)!=42 or len(biota)!=33 or len(source["species"])!=339:
    raise ValueError("Source denominators changed")
  indexed={r["Species"]:r for r in traits if r.get("Species")}
  classes={r["independent_checklist_name"]:r["canary_origin_stratum"] for r in biota}
  if len(classes)!=33:raise ValueError("Original Canary classes missing")
  species339=set(source["species"])
  host_fields=[
    "NumberOfHostplantFamilies","Voltinism","WingSize"
  ]
  extra_fields=["CanopyAffinity","EdgeAffinity","MoistureAffinity","DisturbanceAffinity"]
  wing_keys=("WS_L","WS_U","FW_L","FW_U","WS_L_Fem","WS_U_Fem",
             "WS_L_Mal","WS_U_Mal","FW_L_Fem","FW_U_Fem","FW_L_Mal","FW_U_Mal")
  results=[]
  for r in original:
    name=r["scientific_name"]
    selected=ALIASES.get(name,name)
    data=indexed.get(selected)
    if data is None:data=indexed.get(name)
    host={"NumberOfHostplantFamilies":bool(data and positive(data.get("NumberOfHostplantFamilies"))),
          "Voltinism":bool(data and present(data.get("Voltinism"))),
          "WingSize":bool(data and any(positive(data.get(k)) for k in wing_keys))}
    contexts={k:bool(data and present(data.get(k))) for k in extra_fields}
    has_raw=bool(data)
    host_eligible=has_raw and all(host.values())
    major_complete=host_eligible and all(contexts.values())
    endemic=(
      classes.get(name)=="endemic_confirmed" or
      r["madeira_origin"]=="island_endemic_species" or
      name in ("Hipparchia azorina","Hipparchia miguelensis")
    )
    results.append({
      "original_independent_species":name,
      "source_trait_species_name":data.get("Species","") if data else "",
      "original_archipelagos":r["archipelagos"],
      "Canary_origin":classes.get(name,"not_Canary"),
      "Madeira_origin":r["madeira_origin"],
      "confirmed_independent_species_level_endemic":int(endemic),
      "source_LepTraits_entry_present":int(has_raw),
      "raw_host_family_count_present":int(host["NumberOfHostplantFamilies"]),
      "voltinism_present":int(host["Voltinism"]),
      "wing_size_proxy_possible":int(host["WingSize"]),
      "habitat_canopy_present":int(contexts["CanopyAffinity"]),
      "habitat_edge_present":int(contexts["EdgeAffinity"]),
      "habitat_moisture_present":int(contexts["MoistureAffinity"]),
      "habitat_disturbance_present":int(contexts["DisturbanceAffinity"]),
      "new_host_only_three_traits_eligible":int(host_eligible),
      "old_all_major_traits_complete":int(major_complete),
      "old_frozen339_present":int(name in species339 or selected in species339),
      "named_host_genus_in_original_NHM_exact":int(r["hosts_named_genus"]),
    })
  endemic=[r for r in results if r["confirmed_independent_species_level_endemic"]]
  other=[r for r in results if not r["confirmed_independent_species_level_endemic"]]
  if len(results)!=42 or len(endemic)!=18:
    raise ValueError("Independent endemic source cohort changed")
  if sum(r["old_frozen339_present"] for r in endemic)!=0:
    raise ValueError("No original endemic should have entered previous 339")
  if sum(r["source_LepTraits_entry_present"] for r in endemic)!=10:
    raise ValueError("Original 10 independent endemics in LepTraits changed")
  if sum(r["new_host_only_three_traits_eligible"] for r in endemic)!=10:
    raise ValueError("Host-only eligibility on raw LepTraits endemic 10 changed")
  if sum(r["old_all_major_traits_complete"] for r in endemic)!=0:
    raise ValueError("Source full environment traits for endemics unexpectedly available")
  def summary(rr):
    return {
      "n":len(rr),
      "raw_LepTraits":sum(x["source_LepTraits_entry_present"] for x in rr),
      "relaxed_three_host_relevant_traits_eligible":sum(x["new_host_only_three_traits_eligible"] for x in rr),
      "old_full_traits_complete":sum(x["old_all_major_traits_complete"] for x in rr),
      "old_frozen339_selected":sum(x["old_frozen339_present"] for x in rr),
      "host_only_but_excluded_by_full_trait_criteria_names":[x["original_independent_species"] for x in rr
        if x["new_host_only_three_traits_eligible"] and not x["old_all_major_traits_complete"]],
    }
  result={
    "schema":"chocho_source_frozen_endpoint_specific_relaxed_trait_gate_v1",
    "original_LepTraits_SHA256":LEP_SHA,
    "source_selected_S1_339_prior_schema":source["schema"],
    "independent42":summary(results),
    "independently_classified_endemic18":summary(endemic),
    "independent_other24":summary(other),
    "host_only_relevant_fields":host_fields,
    "habitat_context_fields_not_required_for_host_only_data_eligibility":extra_fields,
    "limitations":[
      "A relaxed source eligibility flag is NOT a newly fitted butterfly host-resource expansion effect.",
      "Source 339 was frozen independently prior to butterfly outcome analysis, not simply the set of all major-complete species.",
      "Some endemics are absent from original LepTraits taxonomically, and some host networks remain only genus-level even when trait data are recoverable.",
      "A host-family breadth number and wing/voltinism do NOT establish accepted plant species in a particular island.",
      "The original GEB239 +54.9% global / +118.9% source-selected island host opportunity is retained unchanged."
    ]
  }
  OUT.mkdir(parents=True,exist_ok=True)
  with (OUT/"independent42_original_vs_relaxed_trait_gate.csv").open("w",newline="",encoding="utf8") as f:
    w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
  (OUT/"independent42_relaxed_host_only_trait_gate.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
  print("SOURCE_FROZEN_RELAXED_HOST_ONLY_GATE="+json.dumps({
    "all42":result["independent42"],
    "endemics18":result["independently_classified_endemic18"]
  },ensure_ascii=False))
if __name__=="__main__":evaluate()
