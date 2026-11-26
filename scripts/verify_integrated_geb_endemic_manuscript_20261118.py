#!/usr/bin/env python3
"""Source-denominator and inference gate for the one-paper GEB endemic integration.

No third-party raw data are reanalysed here. The source records were independently
replayed in the archived PR #177; this guard checks exactly what the paper claims.
"""
from pathlib import Path
import json
import re
import csv

ROOT=Path(__file__).resolve().parents[1]
MAN=ROOT/"manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md"
SI=ROOT/"manuscript/integrated_butterfly_resource_access_20261109/SUPPLEMENTARY_ENDEMIC_HOST_VALIDITY_20261118.md"
P=ROOT/"provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118"
FIELD=P/"VERIFIED_FOCAL_ALIEN_LARVAL_HOST_PLANT_SPECIES_20261118.json"
ENDEMIC=P/"VERIFIED_ENDEMIC18_NONIDENTIFICATION_TRANSPORT_20261118.json"
EIGHT=P/"VERIFIED_INDEPENDENT8_HOST_GEOGRAPHY_AND_FIELD_GATES_20261118.json"
FAMILY=ROOT/"results/hosts_use_vs_flora_20261116/VERIFIED_FLORISTIC_AND_FAMILY_BASELINE_20261116.json"
SP=ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv"
ISLAND=ROOT/"results/island_resource_gap_filling/by_species.csv"
CROSSED=ROOT/"results/WCVP126_crossed_native_alien_20261117/VERIFIED_SOURCE_NATIVE_ALIEN_CROSSED_20261117.json"
MARINE=ROOT/"results/island51_cross_community_marine_20261107/comparison_scope_summary.csv"

def check(condition,msg):
    if not condition: raise AssertionError(msg)

def main():
    doc=MAN.read_text(encoding="utf-8")
    si=SI.read_text(encoding="utf-8")
    field=json.loads(FIELD.read_text(encoding="utf-8"))
    endemic=json.loads(ENDEMIC.read_text(encoding="utf-8"))
    eight=json.loads(EIGHT.read_text(encoding="utf-8"))
    family=json.loads(FAMILY.read_text(encoding="utf-8"))
    with SP.open(encoding="utf-8-sig",newline="") as src:
        species=list(csv.DictReader(src))
    check(len(species)==239 and sum(int(x["introduced_added_units"])>0 for x in species)==206,
          "Frozen known-host butterfly population")
    check(sum(int(x["introduced_added_units"]) for x in species)==14553,
          "Global plant-resource additions")
    check("Coauthor review draft; unsubmitted" in doc,"Working paper status")
    crossed=json.loads(CROSSED.read_text(encoding="utf-8"))
    with ISLAND.open(encoding="utf-8-sig",newline="") as src:
        island_rows=[r for r in csv.DictReader(src) if r["group"]=="island"]
    with MARINE.open(encoding="utf-8-sig",newline="") as src:
        marine_rows=list(csv.DictReader(src))
    check(len(island_rows)==239 and sum(int(r["added"]) for r in island_rows)==5741,
          "Original island-only opportunity denominator")
    island_rows.sort(key=lambda r:(-int(r["native"]),r["species"]))
    check(int(island_rows[58]["native"])==29 and int(island_rows[59]["native"])==28,
          "Untied top-59 native island breadth stratum")
    check(sum(int(r["added"]) for r in island_rows[:59])==2271,
          "Top-59 resource geographically broad beneficiaries")
    check(crossed["original_GSHHG126_verified_model_coastal_islands"]==104
          and crossed["original_exact_common_WCVP_native_alien_and_GSHHG_pair_sample"]==5355,
          "Matched plant geometry source sample")
    assert abs(crossed["relative_MSE_positive_is_improvement"]["native"]["G_vs_B_pct"]-2.2343211230)<1e-8
    assert abs(crossed["relative_MSE_positive_is_improvement"]["alien"]["G_vs_B_pct"]-1.7750986997)<1e-8
    matches=[r for r in marine_rows if r["panel"]=="later_butterflies_retro_hosts"
             and r["pair_class"]=="all_both_plants_indexed" and r["contrast"] in ("BASE_to_AIS","BASE_to_BOTH")]
    check(len(matches)==2 and all(int(r["unique_source_pairs"])==990 for r in matches),
          "Stricter butterfly geographic transfer cohort")
    gains={r["contrast"]:float(r["MSE_gain_percent"]) for r in matches}
    check(abs(gains["BASE_to_AIS"]-.677555746745)<1e-8 and
          abs(gains["BASE_to_BOTH"]+4.69859540475)<1e-8,
          "Stricter AIS-only positive and AIS-plus-wind adverse results")
    # v0.4 has four Results subsections and the independent endemic methods at 2.6.
    # Original numerical cross-geography studies remain in its untouched S.A–S.F.
    split_key="## References (DOI identities checked; final source and style audit outstanding)"
    check(doc.count(split_key)==1, "Original reference and source supplement boundary preserved")
    main_text, appendix=doc.split(split_key,1)
    journal_si=MAN.with_name("GEB_SUPPORTING_INFORMATION_V05.md").read_text(encoding="utf8")
    captions=MAN.with_name("FIGURE_CAPTIONS_V03.md").read_text(encoding="utf8")
    precise=main_text+"\n"+journal_si+"\n"+captions
    for phrase in ("21.4", "17.7", "39.56%", "2.95", "4.239%", "−0.110%"):
        check(phrase in precise, "Exact resource/port/marine result missing from retained evidence: "+phrase)
    for phrase in ("Pieris cheiranthi", "Tropaeolum majus", "nearly flat gradient", "no direction-specific support"):
        check(phrase in main_text, "Main ecological result or qualification missing: "+phrase)
    for phrase in ("### S.F.", "McRae, B. H.", "ESR (2009)", "Vargas, P."):
        check(phrase in appendix, "Original references/sea methodological audit removed: "+phrase)
    check("no new butterfly colonization" not in doc.lower(),
          "Do not infer a population-level absence of colonization from nonidentification")
    abstract=main_text.split("## Structured abstract",1)[1].split("## 1. Introduction",1)[0]
    check("Pieris cheiranthi" in abstract and "Tropaeolum majus" in abstract,
          "Direct independent endemic larvae and introduced host must be named in the revised abstract")
    check("18 are endemic" in main_text and "none of those 18" in main_text,
          "Main text must report the 0/18 independently source-defined endemic ascertainment gap")
    results=main_text.split("## 3. Results",1)[1].split("## 4. Discussion",1)[0]
    headings=[list(re.finditer(r"^### 3\."+str(n)+r"\.",results,re.M)) for n in range(1,5)]
    check(all(len(h)==1 for h in headings) and
          all(headings[n][0].start()<headings[n+1][0].start() for n in range(3)) and
          len(re.findall(r"^### 3\.[1-9]\.",results,re.M))==4,
          "Four unique and ordered v0.4 main Results sections")
    check("### 2.6. Independently documented island-endemic" in main_text and
          "### S.E. Independent island-endemic" in appendix,
          "Independent source Methods 2.6 and preserved original Supplement S.E")
    check("SUPPLEMENTARY_ENDEMIC_HOST_VALIDITY_20261118.md" in main_text,
          "Complete companion endemic methods source linked")
    for name in ("Pieris cheiranthi","Tropaeolum majus","Tenerife","La Palma",
                 "Euchloe grancanariensis","Hirschfeldia incana","Gonepteryx cleobule"):
        check(name in doc and name in si,"Original named source observation or island missing: "+name)
    for exact in ("+118.9%","+40.9%","4.12%"):
        check(exact in precise,"Verified butterfly source positive missing from main/SI/legends: "+exact)
    for exact in ("5.08%","−3.87 points","−15.85 to +8.13"):
        check(exact in appendix,"Archived source sensitivity or adverse check missing: "+exact)
    check("newly host-present region" in doc and "home archipelago" in doc,
          "Source native-host presence and endemic island origin not distinguished")
    sc=endemic["source_scope"]
    check(sc["living_butterflies"]==41 and
          sc["independently_origin_classified_endemic_extant_butterflies"]==18 and
          sc["living_nonendemic_butterflies"]==23,
          "Independent living/endemic denominators")
    check(endemic["original_selection"]["original_final239"]["endemic_n"]==0 and
          endemic["original_selection"]["original_final239"]["nonendemic_n"]==13,
          "Original source panel has zero of 18 endemic butterflies")
    c=endemic["foodplant_source_selected_effect"]
    check(c["exact_named_larval_host_plant_species_known_butterflies_n"]==8 and
          c["source_list_observed_positive_added_WGSRPD3_n"]==4 and
          c["source_list_observed_zero_added_WGSRPD3_n"]==4,
          "Independent 8 documented species source selection")
    check(eight["origin_archipelago_not_individual_island"][
        "newly_host_present_origin_archipelago_region_units_created"]==0,
        "No newly source-archipelago host-presence cell")
    rows=eight["source_endemic_resource"]["source_butterfly_rows"]
    check(len(rows)==8 and sum(x["potential_global_added_region_units"]>0 for x in rows)==4,
          "8 species independent source-result rows")
    check(all(x["original_WCVP_native_plus_introduced_global_region_units"]-
          x["original_WCVP_native_global_region_units"]==x["potential_global_added_region_units"]
          for x in rows),"Per-species WCVP source arithmetic")
    positive=field["primary_new_positive"]
    check(positive["endemic_butterfly"]=="Pieris cheiranthi" and
          positive["introduced_named_larval_host"]=="Tropaeolum majus",
          "Accepted named field host and butterfly")
    check({x["island"] for x in positive["independently_observed_larvae"]}==
          {"Tenerife","La Palma"}, "Two named island original photos")
    check(not field["site_resident_versus_breeding"]["survival_and_recruitment_measured_in_current_audit"],
          "Do not invent host-caused adult recruitment")
    f=family["scopes"]["HOSTS_Papilionoidea__exclude_FAO_exact__min10"][
        "contrasts"]["family_flora_adjusted_excess"]
    check(f["median_difference"]<0 and f["descriptive_locality_bootstrap95"][1]>0,
          "Family-matched residual not a demonstrated positive preference")
    abstract=doc.split("## Structured abstract",1)[1].split("## 1. Introduction",1)[0]
    wc=len(abstract.replace("*","").split())
    check(wc<=300,"GEB structured abstract currently >300 words: "+str(wc))
    print(json.dumps({"status":"PASS","frozen_butterflies":239,
      "independent_living_island_butterflies":41,"independent_endemics":18,
      "original_endemic_overlap":0,"host_species_documented_endemics":8,
      "source_positive_global_resource_endemics":4,"home_archipelago_new_cells":0,
      "field_alien_host_named_islands":["Tenerife","La Palma"],
      "family_adjusted_extra_preference_supported":False,
      "abstract_word_count":wc,"island_gap_start_pp":21.4,"island_gap_end_pp":17.7,
      "top59_island_gain_fraction_pct":round(100*2271/5741,4),
      "marine_crossed_native_shortest_pct":round(crossed["relative_MSE_positive_is_improvement"]["native"]["G_vs_B_pct"],4),
      "marine_crossed_alien_shortest_pct":round(crossed["relative_MSE_positive_is_improvement"]["alien"]["G_vs_B_pct"],4),
      "stricter_AIS_gain_pct":round(gains["BASE_to_AIS"],4),
      "stricter_AIS_wind_gain_pct":round(gains["BASE_to_BOTH"],4),
      "proof_type":"source claim consistency, not population inference"},
      sort_keys=True))

if __name__=="__main__":
    main()
