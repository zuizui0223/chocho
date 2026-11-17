#!/usr/bin/env python3
"""Source-matched descriptive controls for 12 NOAA-model 2005 shore approaches.

Checks whether the 10/12 co-listed source-native → recipient-alien species
criterion is exceptional compared with geographically and botanically
comparable pairings from the SAME source island. Do not recode NOAA model
nonapproach as absence of real ocean dispersal, or contemporaneous accepted
WCVP origin flags as actual dated introduction.

All candidate matching uses only SOURCE ISLAND, botanical Level1 equality,
log geodesic distance, and log recipient alien richness — not response.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path

MAX_DISTANCE_LOG=math.log(2.0)
MAX_RECIPIENT_ALIEN_LOG=math.log(3.0)
MIN_CONTROLS=3
EXCLUDE_MODEL_POSITIVE=True

def unique(root,name):
    vals=list(root.rglob(name))
    if len(vals)!=1:raise RuntimeError(f"Original source {name} does not have exact one source artifact copy: {len(vals)}")
    return vals[0]

def rows(p):
    with p.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def backward_pair(plant,source,dest):
    if plant["source_region"]==source and plant["destination_region"]==dest:
        return int(plant["native_in_source_and_alien_in_dest"]),int(plant["alien_in_source_and_native_in_dest"])
    if plant["source_region"]==dest and plant["destination_region"]==source:
        return int(plant["alien_in_source_and_native_in_dest"]),int(plant["native_in_source_and_alien_in_dest"])
    raise RuntimeError("Invalid original two-way WCVP taxon status pair mapping")

def origin_dest_margins(p,source,dest):
    if p["source_region"]==source and p["destination_region"]==dest:
        return int(p["source_native_WCVP_accepted"]),int(p["dest_alien_WCVP_accepted"])
    if p["source_region"]==dest and p["destination_region"]==source:
        return int(p["dest_native_WCVP_accepted"]),int(p["source_alien_WCVP_accepted"])
    raise RuntimeError("Original source/dest WCVP denominators do not match")

def main():
    parser=argparse.ArgumentParser()
    for name in ("current-positive-root","full-WCVP-root","current-matrix-root","out-dir"):
        parser.add_argument("--"+name,required=True)
    a=parser.parse_args()
    from pathlib import Path
    o=Path(a.out_dir)
    positive=rows(unique(Path(a.current_positive_root),"NOAA2005_source_positive12_WCVP_flora_species_geographic_index.csv"))
    floras=rows(unique(Path(a.full_WCVP_root),"WCVP126_flora_similarity_strict_WPI2020_geographic_pairs.csv"))
    network=rows(unique(Path(a.current_matrix_root),"original126_NOAA2005_directed_source_nearshore_matrix.csv"))
    if len(positive)!=12 or len(floras)!=7875 or len(network)!=15750:
        raise RuntimeError("Original 2005 12 positive / 126 7875 plants / 15750 directed geographic sources changed")
    pairs={}
    for p in floras:
        key=tuple(sorted((p["source_region"],p["destination_region"])))
        if key in pairs:raise RuntimeError("Original WCVP unordered plant pair duplicate")
        pairs[key]=p
    currents={(r["source"],r["destination"]):r for r in network}
    if len(currents)!=15750:raise RuntimeError("Duplicate NOAA actual directed source and recipient IDs")
    cases=[]
    for z in positive:
        source,dest=z["source"],z["destination"]
        if (source,dest) not in currents:raise RuntimeError("NOAA model-positive source direction not found")
        r=currents[(source,dest)]
        if r["model_physical_source_evaluation_status"]!="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE":
            raise RuntimeError("One positive source 2005 direction not original eligible")
        if int(r["NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"])<=0:
            raise RuntimeError("Positive WCVP source not positive original NOAA modeled nearshore")
        p=pairs[tuple(sorted((source,dest)))]
        fwd,rev=backward_pair(p,source,dest)
        if fwd!=int(z["WCVP_native_in_source_alien_in_dest_accepted_taxa_species_count"]):
            raise RuntimeError("Accepted WCVP source-native to recipient-alien matching not independent")
        native_source,alien_recipient=origin_dest_margins(p,source,dest)
        pos_logd=math.log1p(float(p["centroid_geodesic_km"]))
        pos_logalien=math.log1p(alien_recipient)
        level1=p["same_plant_level1"]
        possible=[]
        for (src,other),m in currents.items():
            if src!=source or other==dest:continue
            if m["model_physical_source_evaluation_status"]!="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE":
                continue
            if EXCLUDE_MODEL_POSITIVE and int(m["NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"])>0:
                continue
            q=pairs.get(tuple(sorted((src,other))))
            if q is None or q["same_plant_level1"]!=level1:continue
            _,alien_t=origin_dest_margins(q,src,other)
            if (abs(math.log1p(float(q["centroid_geodesic_km"]))-pos_logd)>MAX_DISTANCE_LOG or
                abs(math.log1p(alien_t)-pos_logalien)>MAX_RECIPIENT_ALIEN_LOG):
                continue
            c_fw,c_rev=backward_pair(q,src,other)
            possible.append({"other_region":other,"source_native_to_recipient_alien":c_fw,
              "reverse_native_to_alien":c_rev,
              "source_alien_recipient":alien_t,
              "geodesic_km":float(q["centroid_geodesic_km"]),
              "botanical_level1_same":q["same_plant_level1"]})
        control_freq=(sum(c["source_native_to_recipient_alien"]>0 for c in possible)/len(possible)
                      if possible else None)
        control_normalized=(statistics.median(
          [c["source_native_to_recipient_alien"]/(1+c["source_alien_recipient"]) for c in possible])
          if possible else None)
        cases.append({
            "source":source,"destination":dest,
            "original_noaa_model_nearshore_count":int(r["NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"]),
            "geodesic_km":float(p["centroid_geodesic_km"]),
            "same_TDWG_Level1":level1=="1",
            "source_native_in_recipient_alien_accepted_taxa":fwd,
            "reverse_recipient_native_in_source_alien_accepted_taxa":rev,
            "source_native_taxa":native_source,"recipient_alien_taxa":alien_recipient,
            "same_source_matched_NOAA_model_nonapproach_control_n":len(possible),
            "source_matched_control_same_native_to_alien_positive_fraction":control_freq,
            "source_matched_control_median_native_to_alien_per1p_recipient_alien":
              control_normalized,
            "minimum_three_matched_controls_met":len(possible)>=MIN_CONTROLS,
            "model_nonapproach_is_no_natural_seed_dispersal":False,
            "positive_current_model_is_actual_seed_immigration":False,
        })
    if len({(c["source"],c["destination"]) for c in cases})!=12:
        raise RuntimeError("Original 12 source directions not unique")
    sufficient=[c for c in cases if c["minimum_three_matched_controls_met"]]
    total_pos=sum(c["source_native_in_recipient_alien_accepted_taxa"]>0 for c in cases)
    if total_pos!=10:raise RuntimeError("Original positive 10/12 taxon co-listing changed")
    observed_sufficient=sum(c["source_native_in_recipient_alien_accepted_taxa"]>0 for c in sufficient)
    freq=[c["source_matched_control_same_native_to_alien_positive_fraction"] for c in sufficient]
    result={
      "schema":"chocho_WCVP126_2005_OSCAR_12_current_positive_directional_source_matched_botanical_controls_v1",
      "status":"SOURCE_MATCHED_NOAA_MODEL_NONAPPROACH_CONTROLS_DESCRIPTIVE_NOT_CAUSAL_PLANT_TRANSPORT",
      "original_2005_NOAA_particle_CI":44012093616,
      "original_directional_WCVP_taxon_crosswalk_CI":44012684281,
      "original_126island_WCVP_all_species_pair_CI":44010307551,
      "original_2005_source_positive_directed_pairs":len(cases),
      "original_2005_source_native_to_recipient_alien_taxon_positive_pairs":total_pos,
      "same_source_distance_factor_max":2.0,
      "same_source_recipient_alien_richness_factor_max":3.0,
      "same_TDWG_Level1_boolean_required":True,
      "source_nonapproach_controls_do_not_represent_no_actual_dispersal":True,
      "minimum_nonapproach_controls_per_case":MIN_CONTROLS,
      "n_current_model_positive_directions_with_3plus_matched_controls":len(sufficient),
      "observed_positive_co_listing_among_matched_supported_cases":observed_sufficient,
      "expected_co_listing_positive_count_from_matched_control_frequencies_sum":
         float(sum(freq)),
      "source_controls_co_listing_fraction_median":statistics.median(freq) if freq else None,
      "n_source_current_model_positive_directed_pairs_with_reverse_co_listing":
          sum(c["reverse_recipient_native_in_source_alien_accepted_taxa"]>0 for c in cases),
      "case_details":cases,
      "no_formal_effect_model_or_hypothesis_test_fitted":True,
      "limitations":[
        "The case-control design was chosen posthoc on 12 current-positive modeled 2005 source→recipient directions, no independent replicated sea drift or prespecified transport-null hypothesis.",
        "Real source currents were decimated to 1deg from 1/3deg NOAA old OSCAR, with artificial launch up to 200km from shoreline and no actual seed observations.",
        "The model-nonapproach eligible source directions are NOT independently proven natural non-dispersal controls.",
        "Source native in one area and recipient introduced in WCVP is a present-day plant distribution co-report, not a dated local seed movement.",
        "Source-matched constraints select based on distance, recipient alien richness and broad region only; habitat, floristic detection, historical trade and plate tectonic histories remain uncontrolled.",
        "If fewer than 3 same-source matched control directions exist, no case-specific contrast is interpreted.",
        "This independent botanical taxon co-listing is not a butterfly HOSTS larval feeding or a new insular population.",
      ]
    }
    o.mkdir(parents=True,exist_ok=True)
    (o/"NOAA2005_native_to_alien_12_model_nearshore_source_matched_controls.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("NOAA2005_DIRECTIONAL_WCVP_NATIVE_TO_ALIEN_SOURCE_MATCHED_CONTROLS",json.dumps({
      "status":result["status"],
      "12_model_positive_pairs":len(cases),
      "10_source_native_to_recipient_alien_cases":total_pos,
      "number_with_three_or_more_source_distance_richness_controls":len(sufficient),
      "actual_positive_co_listed_cases_with_controls":observed_sufficient,
      "control_expected_positive_co_listed_cases":result["expected_co_listing_positive_count_from_matched_control_frequencies_sum"],
      "reverse_native_to_alien_co_listing_positive":result["n_source_current_model_positive_directed_pairs_with_reverse_co_listing"],
      "matched_control_candidate_n_per_case":[{"source":x["source"],"dest":x["destination"],"n":x["same_source_matched_NOAA_model_nonapproach_control_n"],
        "control_prevalence":x["source_matched_control_same_native_to_alien_positive_fraction"]}
       for x in cases]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
