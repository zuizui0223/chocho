#!/usr/bin/env python3
"""Partial identification and target-population overlap audit for island-endemic butterflies.

Original independent universe: 42 INBO butterflies, of which one historical
Madeira endemic is extinct; target = 18 *extant* independently classified
endemic butterflies among 41 extant butterflies.

Falsification:
  * 0/18 endemics in original source panel339, 15/23 extant nonendemics.
    No observed support for applying original GEB239 estimand via IPW.
  * 8/18 have exact source-named host plant species, 4/8 show positive source-
    fixed global WGSRPD3 potential added area, 4/8 show zero.
  * The remaining 10 are unobserved at species-level, NOT zero effects.
  * The added-region difference is NOT monotone to discovery of new host
    species: a new native host may erase previously introduced-only regions
    from the count, while contemporary potential coverage itself cannot shrink.
Thus 4/18 is a rate of *currently evidenced* source-list positives, NOT a
mathematical lower bound on the complete biological target's true positive rate.

This audit never extrapolates field larval use, demographic establishment,
port causality, nor transfers the original +54.9%/+118.9% to all island endemics.
"""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

BASE = Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
HERE = Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
INBO = BASE / "independent_42_species_vs_frozen239.csv"
CANARY_FILES = [HERE/"canary_BIOTA_source_status_20.csv",
                HERE/"canary_BIOTA_remaining_13.csv"]
SOURCE8 = HERE/"VERIFIED_INDEPENDENT8_HOST_GEOGRAPHY_AND_FIELD_GATES_20261118.json"
GENUS = HERE/"VERIFIED_EUCHLOE_NATIVE_HOST_AND_HIPPARCHIA_GATE_20261118.json"
EXTINCT = {"Pieris wollastoni"}
AZORES_END = {"Hipparchia azorina", "Hipparchia miguelensis"}

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def comb(n,k):
    if k < 0 or k > n:return 0
    return math.comb(n,k)

def exact_conditional_source_selection_probability(table):
    """Descriptive, one-sided fixed-margin probability, NOT independent species."""
    endemic_n = table["endemic_n"]
    other_n = table["other_n"]
    selected_total = table["selected_endemic"] + table["selected_other"]
    observed = table["selected_endemic"]
    denominator=comb(endemic_n + other_n, selected_total)
    return sum(comb(endemic_n,k)*comb(other_n,selected_total-k)
       for k in range(max(0,selected_total-other_n),observed+1))/denominator

def observed_plant_coverage(hosts):
    """A source-list global plant geography metric. Uses hypothetical sets only."""
    native=set()
    contemporary=set()
    for n,i in hosts:
        native|=set(n)
        contemporary|=set(n)|set(i)
    return {"native":native,"contemporary":contemporary,
            "added":contemporary-native}

def nonmonotone_unit_test():
    a=observed_plant_coverage([({"AZO"},{"MDR"})])
    b=observed_plant_coverage([({"AZO"},{"MDR"}),({"MDR"},set())])
    c=observed_plant_coverage([({"AZO"},set())])
    d=observed_plant_coverage([({"AZO"},set()),({"AZO"},{"MDR"})])
    assert len(a["added"])==1 and len(b["added"])==0
    assert len(c["added"])==0 and len(d["added"])==1
    assert a["contemporary"]<=b["contemporary"]
    assert c["contemporary"]<=d["contemporary"]
    return {
      "hypothetical_two_region_test_only":True,
      "initial_host_native":["AZO"],"initial_host_introduced":["MDR"],
      "add_hypothetical_new_native_host_region":["MDR"],
      "initial_added_region_count":len(a["added"]),
      "after_new_native_host_added_region_count":len(b["added"]),
      "contemporary_plant_coverage_monotone_under_new_hosts":True,
      "potential_added_region_difference_nonmonotone_under_new_hosts":True,
      "reverse_zero_to_positive_possible":len(c["added"])==0 and len(d["added"])==1,
      "not_a_real_butterfly_host_species_suggestion":True,
    }

def main(selection_path, output):
    raw = read(INBO)
    selection = read(selection_path)
    canary = [r for path in CANARY_FILES for r in read(path)]
    eight = json.loads(SOURCE8.read_text(encoding="utf-8"))
    genus = json.loads(GENUS.read_text(encoding="utf-8"))
    if len(raw)!=42 or len(selection)!=42 or len(canary)!=33:
        raise ValueError("Original 42 species and BIOTA 33 origin source data drift")
    from_inbo={r["scientific_name"]:r for r in raw}
    from_selection={r["independent_species"]:r for r in selection}
    if len(from_inbo)!=42 or len(from_selection)!=42 or set(from_selection)!=set(from_inbo):
        raise ValueError("Independent origin universe and audited selection differ")
    canon={r["independent_checklist_name"]:r["canary_origin_stratum"] for r in canary}
    if len(canon)!=33:raise ValueError("Duplicate Canary source name or incomplete official origin")
    canary_endemic={name for name,status in canon.items() if status=="endemic_confirmed"}
    madeira_endemic={r["scientific_name"] for r in raw if r["madeira_origin"]=="island_endemic_species"}
    endemic=canary_endemic|madeira_endemic|AZORES_END
    if (len(canary_endemic),len(madeira_endemic),len(AZORES_END),len(endemic))!=(13,3,2,18):
        raise ValueError("Independently classified 18 endemic species denominator changed")
    if not EXTINCT<=set(from_inbo) or EXTINCT & endemic:
        raise ValueError("Historical extinct Madeira endemic incorrectly coded extant")
    living=set(from_inbo)-EXTINCT
    living_endemic=endemic & living
    living_nonendemic=living-endemic
    if (len(living),len(living_endemic),len(living_nonendemic))!=(41,18,23):
        raise ValueError("Extant 41 vs 18/23 origin strata changed")
    if any(not int(from_selection[name]["in_original_frozen_339_S1_species_panel"]) in (0,1)
           for name in living):
        raise ValueError("Source 339 stage not binary")
    selection_records=[]
    for name in sorted(living):
        sel=from_selection[name]
        member=name in living_endemic
        selection_records.append({
          "independent_species":name,
          "origin": "independent_island_endemic" if member else "extant_nonendemic",
          "original_LepTraits_any_name":int(sel["in_original_LepTraits_any_name"]),
          "original_LepTraits_major_complete":int(sel["source_LepTraits_major_complete"]),
          "in_original_GEB339":int(sel["in_original_frozen_339_S1_species_panel"]),
          "in_original_GEB239":int(sel["in_original_frozen_239_final_resource_sample"]),
          "frozen_selection_stage":sel["sampling_stage"],
        })
    endemic_sel=[r for r in selection_records if r["origin"]=="independent_island_endemic"]
    other_sel=[r for r in selection_records if r["origin"]=="extant_nonendemic"]
    source339={k:sum(x["in_original_GEB339"] for x in g)
       for k,g in (("endemic",endemic_sel),("nonendemic",other_sel))}
    source239={k:sum(x["in_original_GEB239"] for x in g)
       for k,g in (("endemic",endemic_sel),("nonendemic",other_sel))}
    if source339!={"endemic":0,"nonendemic":15} or source239!={"endemic":0,"nonendemic":13}:
        raise ValueError("Original selection overlap must remain 0 endemic / 15 or 13 nonendemic")
    if sum(x["original_LepTraits_any_name"] for x in endemic_sel)!=10 or \
       sum(x["original_LepTraits_major_complete"] for x in endemic_sel)!=0:
        raise ValueError("The exact original LepTraits missingness barrier changed")
    source_entries=eight["source_endemic_resource"]["source_butterfly_rows"]
    if len(source_entries)!=8:
        raise ValueError("Known plant-species 8-subset drift")
    recorded={r["butterfly"]:r for r in source_entries}
    if len(recorded)!=8 or not set(recorded)<=endemic:
        raise ValueError("Exact source known species are not the 18 independent endemic universe")
    if len({r["butterfly"] for r in source_entries})!=8:
        raise ValueError("Duplicate source butterfly")
    positive={name for name,r in recorded.items()
              if int(r["potential_global_added_region_units"])>0}
    zero=set(recorded)-positive
    missing=endemic-set(recorded)
    if len(positive)!=4 or len(zero)!=4 or len(missing)!=10:
        raise ValueError("Four positive, four zero, ten no exact-host case structure changed")
    if eight["source_endemic_resource"]["strict_direct_wild_named_plant_butterfly_taxa_n"]!=3:
        raise ValueError("Original strict field source count changed")
    no_genus=set(genus["extant_independent_coverage"]["five_Hipparchia_endemics_missing_named_host_genus"])
    if len(no_genus)!=5 or not no_genus<=missing:
        raise ValueError("Remaining unknown species must include exactly five Canary Hipparchia")
    known_genus_missing_species=missing-no_genus
    if len(known_genus_missing_species)!=5:
        raise ValueError("Expected five genus-known but exact host-species-unknown endemics")
    # The observed eight were selected on ability to resolve host plants.
    # Positive and zero labels are deterministic for their *frozen source lists*.
    # Under UNCHECKED additional host discoveries, both labels may change.
    toy=nonmonotone_unit_test()
    assumption_fixed_8={
      "scenario":"Assume all eight known source host lists are biologically complete and stay fixed. Allow remaining ten species to take either outcome.",
      "positive_counts_min":4,
      "positive_counts_max":14,
      "positive_fraction_min":4/18,
      "positive_fraction_max":14/18,
      "counterfactual_only":True,
      "not_identification_of_true_endemic_butterfly_effect":True
    }
    unrestricted={
      "scenario":"No additional-host completeness assumption; new native hosts can erase 'introduced-only' added units, new introduced hosts can create new ones.",
      "positive_counts_min":0,
      "positive_counts_max":18,
      "positive_fraction_min":0,
      "positive_fraction_max":1,
      "is_vacuous_mathematical_bound_not_empirical_estimate":True
    }
    # Exact conditional probability is descriptively illustrative and
    # inappropriate as an evolutionary independence test for shared clades.
    table339={"endemic_n":18,"other_n":23,"selected_endemic":0,"selected_other":15}
    table239={"endemic_n":18,"other_n":23,"selected_endemic":0,"selected_other":13}
    summary={
      "schema":"chocho_extant_endemic_butterfly_external_validity_nonidentification_v1",
      "original_independent_species":42,
      "historic_extinct_excluded":"Pieris wollastoni",
      "target_extant":41,
      "target_extant_endemic":18,
      "target_extant_nonendemic":23,
      "selection":{
        "in_original_LepTraits_endemics":10,
        "complete_major_LepTraits_endemics":0,
        "original_GEB339":source339,
        "original_GEB239":source239,
        "original339_nonendemic_selected_fraction":15/23,
        "original239_nonendemic_selected_fraction":13/23,
        "endemic_target_selection_overlap":False,
        "inverse_probability_weighting_of_endemic_effect_identifiable":False,
        "population_effect_transport_to_all_18_identifiable":False,
        "descriptive_fixed_margin_probabilities_not_phylogenetic_tests":{
          "339_selection_one_sided":exact_conditional_source_selection_probability(table339),
          "239_selection_one_sided":exact_conditional_source_selection_probability(table239)
        },
        "interpretation":"No endemic observed in source or final resource cohorts. An IPW parameter for 18 endemic target species has zero source support; cannot reweight original 239 outcomes."
      },
      "source_plant_taxonomic_coverage":{
        "original_independent_endemic_species":18,
        "exact_named_larval_host_plant_species_known":8,
        "exact_plant_species_not_known":10,
        "genus_known_but_species_not_known_5":sorted(known_genus_missing_species),
        "genus_and_species_not_known_5":sorted(no_genus),
        "documented_source_positive_global_added_region_count":4,
        "documented_source_zero_global_added_region_count":4,
        "documented_source_positives":sorted(positive),
        "documented_source_zeros":sorted(zero),
        "source_positive_fraction_among_8_source_selected":4/8,
        "count_with_positive_documented_source_envelopes_among_original18":4,
        "source_selected_8_are_nonrandom":True,
        "reference_original_replay":"https://github.com/zuizui0223/chocho/actions/runs/44221587347",
      },
      "counterfactual_effect_bounds_not_actual_prevalence":{
        "assuming_the_8_host_lists_complete":assumption_fixed_8,
        "allowing_additional_true_host_species_in_all_18":unrestricted,
        "source_of_nonmonotonicity":toy,
        "mathematical_explanation":"Let N be union of native WGSRPD3 units across all larval-host species and C be union of native+introduced units. Adding an unobserved host monotonically expands N and C; the difference C minus N is NOT monotonic. Therefore observed 4 positive source-list outcomes are NOT a guaranteed lower bound for the true full-host-list positive fraction."
      },
      "biological_evidence_limits":{
        "independently_checked_eight_butterflies_named_island_cells":22,
        "source_direct_named_plant_larval_events":4,
        "source_direct_wild_events":3,
        "source_direct_garden_events":1,
        "no_plant_species_presence_or_larva_evidence_measured_on_other_18_cells":True,
        "actual_population_establishment_inferred":False,
        "original_GEB_known_HOSTS_global_54_9_pct_unchanged":True,
        "original_GEB_selected_island_118_9_pct_unchanged":True
      },
      "scientific_decision":"Original GEB inference is valid for original source-selected cohort. Independent island endemic source-panel evidence supports heterogeneous potential host geographies and directly observed alien host use in P. cheiranthi, but no population-level mean resource effect can be generalized to 18 endemic taxa without complete plant-species host lists and comparable target-unit habitat data."
    }
    output.mkdir(parents=True,exist_ok=True)
    with (output/"extant41_independent_origin_vs_original_GEB_selection.csv").open(
        "w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(selection_records[0]))
        writer.writeheader();writer.writerows(selection_records)
    per_endemic=[]
    for name in sorted(endemic):
        if name in positive: status="source_positive_added_region"
        elif name in zero:status="source_zero_added_region_not_biological_zero"
        else:status="exact_host_species_unresolved"
        per_endemic.append({
          "independent_endemic_butterfly":name,
          "exact_source_plant_status":status,
          "has_genus_but_not_exact_species":int(name in known_genus_missing_species),
          "has_no_named_host_genus_in_reviewed_sources":int(name in no_genus),
          "global_added_WGSRPD3_source_only":recorded[name]["potential_global_added_region_units"] if name in recorded else "",
          "not_assumed_no_hosts":1
        })
    with (output/"independent18_source_plants_possibilities_not_imputed.csv").open(
        "w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(per_endemic[0]))
        writer.writeheader();writer.writerows(per_endemic)
    (output/"nonidentification_transportability_result.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ENDEMIC_TRANSPORTABILITY_NONIDENTIFICATION_AUDIT_PASSED="+json.dumps({
       "extant":41,"endemic":18,"extant_nonendemic":23,
       "original_339":source339,"original_239":source239,
       "source_positive_species":sorted(positive),
       "source_zero_species":sorted(zero),
       "source_unknown_species":sorted(missing),
       "fixed_known_8_scenario_range":[4,14],
       "unrestricted_possible_range":[0,18],
       "source_added_difference_nonmonotone":toy["potential_added_region_difference_nonmonotone_under_new_hosts"],
       "original_IPW_transportability_identifiable":False,
    },ensure_ascii=False),flush=True)
    return summary

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection-source-csv",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    main(args.selection_source_csv,args.out_dir)
