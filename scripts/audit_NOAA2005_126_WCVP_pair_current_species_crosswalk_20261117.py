#!/usr/bin/env python3
"""Source-verified 126-island WCVP flora × direction-aware NOAA 2005 source gate.

This joins already ACTUALLY RUN old NOAA OSCAR 90-day modeled coastal
PROXIMITY to the original accepted WCVP native/alien taxon similarity at
identical WGSRPD3 region IDs. It deliberately does not convert 12 simulated
directed nearshore pairs or four raster contacts into real seed import,
ecological selection, butterfly establishment or real 2020 shipping OD.

The source is sparse and deliberately FAILS a modelling feasibility gate if
current-positive directed region pairs <30 or exact other-TDWG contacts <15.
These descriptive minimums are a post-source diagnostic, not registered
inferential power requirements and not a test of no real dispersal.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path

SHIP_CI=44010307551
OSCAR_CI=44012093616
SOURCE_2005_LIMIT="2005 three start dates, 90d, old NOAA 1degree subsampled current"

def read(path):
    with path.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def finite_or_none(x):
    if x in ("",None):return None
    z=float(x)
    return z if math.isfinite(z) else None

def median_or_none(x):
    return float(statistics.median(x)) if x else None

def main():
    p=argparse.ArgumentParser()
    for a in ("current-audit","current-matrix","current-trajectories",
              "wcvp-pair-geography","sea-geometry","out-dir"):
        p.add_argument("--"+a,type=Path,required=True)
    a=p.parse_args()
    status=json.loads(a.current_audit.read_text())
    if (status["status"]!="REAL_2005_GLOBAL_NOAA_LAGRANGIAN_1DEG_FEASIBILITY_SOURCE_ONLY_NOT_HISTORIC_PLANT_SEED_DELIVERY" or
        status["original_directed_126_minus_self_matrix_rows"]!=15750 or
        status["source_with_at_least_one_original_NOAA_2005_fixed_release"]!=94):
        raise RuntimeError("Original 2005 NOAA modelled source/region sample drift")
    mat=read(a.current_matrix)
    tracks=json.loads(a.current_trajectories.read_text())
    plants=read(a.wcvp_pair_geography)
    geom=read(a.sea_geometry)
    if len(mat)!=15750 or len(plants)!=7875 or len(geom)!=7875 or len(tracks)!=280:
        raise RuntimeError("Original NOAA/WCVP/geometry source cohort not same")
    orig={(r["source_region"],r["destination_region"]):r for r in plants}
    src={(r["source_region"],r["destination_region"]):r for r in geom}
    if len(orig)!=7875 or len(src)!=7875:
        raise RuntimeError("WCVP botanical pair source duplicated")
    if {tuple(sorted([a,b])) for a,b in orig}!={tuple(sorted([a,b])) for a,b in src}:
        raise RuntimeError("WCVP and GSHHG source islands don't pair exactly")
    seen=set();positive=[]
    actual_four=[]
    for r in mat:
        i,j=r["source"],r["destination"]
        if (i,j) in seen or i==j:raise RuntimeError("Original direction matrix duplicated or diagonal")
        seen.add((i,j))
        if r["model_physical_source_evaluation_status"]!="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE":
            if r["NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"] not in ("",None):
                raise RuntimeError("NOAA original source missingness recoded as a positive or zero")
            continue
        n=int(r["NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"])
        ncontact=int(r["NOAA_2005_fixed90d_releases_contacting_another_exact_TDWG3_polygon"])
        if n>0:positive.append((i,j,r,n,ncontact))
        if ncontact>0:actual_four.append((i,j,ncontact))
    if len(seen)!=15750 or len(positive)!=12 or len(actual_four)!=4:
        raise RuntimeError(f"Original directional 12 NOAA nearshore and 4 exact raster contact source pairs drifted: {len(seen)}, {len(positive)}, {len(actual_four)}")
    route_samples=defaultdict(list)
    for track in tracks:
        for dest,details in track["nearshore_regions"].items():
            route_samples[(track["source"],dest)].append({
             "date":track["source_actual_date_utc"],
             "first_day":int(details["earliest_day"]),
             "first_displacement_from_fixed_source_km":float(details["first_displacement_from_source_km"]),
             "minimum_reported_nearshore_distance_km":float(details["closest_km"]),
             "original_track_stopped":track["status"],
             "actual_90d_trajectory_finished":track["status"]=="COMPLETED_90_DAYS",
             })
    # Convert a source→destination link to the ORIGINAL ordered WCVP
    # matrix. Source-native/destination-alien is only taxonomic
    # origin status in two contemporary regions, not dated seed arrival.
    detail=[]
    for i,j,r,n,ncontact in positive:
        k=tuple(sorted((i,j)))
        # WCVP input is ascending source/dest originally, validate source.
        w=orig.get(k)
        if w is None:raise RuntimeError("No accepted WCVP botanical pair for actual NOAA-positive direction")
        source_is_a=(w["source_region"]==i and w["destination_region"]==j)
        if not source_is_a and not (w["source_region"]==j and w["destination_region"]==i):
            raise RuntimeError("Source native/alien taxonomic direction order not recoverable")
        x=src[k]
        has_Ocean_geometry=x["source_sea_geometry_coverage_status"]=="BOTH_ANCHORED"
        species_native_to_alien=int(w["native_in_source_and_alien_in_dest"]
                        if source_is_a else w["alien_in_source_and_native_in_dest"])
        direct_ship=int(w["AIS2020_direct_WPI_source_trips_A_to_B"] if source_is_a
                         else w["AIS2020_direct_WPI_source_trips_B_to_A"])
        same_level1=(w["same_plant_level1"]=="1")
        trajectories=route_samples.get((i,j),[])
        if len(trajectories)!=n:
            raise RuntimeError(f"NOAA trajectories and direction matrix positive links disagree {i}->{j}: {n} vs {len(trajectories)}")
        detail.append({
          "source":i,"destination":j,
          "NOAA_2005_fixed_release_nearshore_count":n,
          "NOAA_2005_source_release_count":int(r["original_source_releases_available"]),
          "NOAA_2005_first_nearshore_day":min(x["first_day"] for x in trajectories),
          "NOAA_2005_min_first_nearshore_displacement_from_source_km":
              round(min(x["first_displacement_from_fixed_source_km"] for x in trajectories),4),
          "NOAA_2005_has_any_nearshore_first_event_after_day7_and_displacement50km":
              any(x["first_day"]>=7 and x["first_displacement_from_fixed_source_km"]>=50 for x in trajectories),
          "NOAA_2005_has_model_GSHHG_mask_contact_inside_another_exact_TDWG":bool(ncontact),
          "geodesic_centroid_km":float(w["centroid_geodesic_km"]),
          "native_Jaccard":finite_or_none(w["native_Jaccard"]),
          "alien_Jaccard":finite_or_none(w["alien_Jaccard"]),
          "native_shared_accepted_WCVP_species":int(w["native_shared_accepted_taxa"]),
          "alien_shared_accepted_WCVP_species":int(w["alien_shared_accepted_taxa"]),
          "WCVP_native_in_source_alien_in_dest_accepted_taxa_species_count":species_native_to_alien,
          "both_islands_same_TDWG_Level1":same_level1,
          "GSHHG_geography_distance_eligible":has_Ocean_geometry,
          "modern_WPI2023_port_proxy_log1p_product":float(w["modern_port_pair_proxy_log1p_product"]),
          "strict_actual_AIS2020_port_OD_source_episode_count":direct_ship,
          "observed_AIS2020_vessel_seed_cargo":False,
          "is_confirmed_real_biological_seed_dispersal":False
        })
    if len(detail)!=12:raise RuntimeError("Original positive directional route report incomplete")
    noaa_late=sum(r["NOAA_2005_has_any_nearshore_first_event_after_day7_and_displacement50km"]
                  for r in detail)
    num_WCVP_ni_ali=sum(r["WCVP_native_in_source_alien_in_dest_accepted_taxa_species_count"]>0 for r in detail)
    overlap2020=sum(r["strict_actual_AIS2020_port_OD_source_episode_count"]>0 for r in detail)
    evidence={
      "schema":"chocho_true2005_global_NOAA_COARSE_directional_proximity_WCVP126_origin_matched_botanical_source_gate_v1",
      "status":"DIRECTED_2005_CURRENT_AND_WCVP_WHOLE_FLORA_CROSSWALK_COMPLETED_EFFECT_MODEL_UNIDENTIFIABLE",
      "NOAA_original_real_flow_source_CI":OSCAR_CI,
      "WCVP_126_native_alien_species_source_CI":44010011288,
      "WPI2020_source_coverage_CI":SHIP_CI,
      "GSHHG126_sea_geometry_baseline_CI":44010673580,
      "original_WCVP_islands":126,
      "original_undirected_botanical_region_pairs":7875,
      "original_source_2005_Lagrangian_coarse_1deg_trajectories":len(tracks),
      "actual_original_NOAA_2005_source_regions":94,
      "original_directional_pair_status_count":dict(Counter(r["model_physical_source_evaluation_status"] for r in mat)),
      "original_modelled_evaluable_source_destination_direction_pairs":7313,
      "original_directed_nearshore_source_positive_pairs":len(detail),
      "original_exact_GSHHG_other_TDWG_raster_contact_positive_pairs":len(actual_four),
      "NOAA_source_positive_dir_pairs_with_first_nearshore_day7_or_later_and_50km_source_displacement":noaa_late,
      "source_positive_dir_pairs_with_WCVP_source_native_recipient_alien_named_taxa":num_WCVP_ni_ali,
      "source_positive_dir_pairs_also_with_2020_actual_WPI_vessel_OD":overlap2020,
      "source_positive_native_Jaccard_median":median_or_none([r["native_Jaccard"] for r in detail if r["native_Jaccard"] is not None]),
      "source_positive_alien_Jaccard_median":median_or_none([r["alien_Jaccard"] for r in detail if r["alien_Jaccard"] is not None]),
      "original_source_precondition_before_adding_global_directional_current_to_flora_prediction":{
        "minimum_positive_directional_island_region_pairs_for_exploratory_signal_model":30,
        "minimum_exact_distinct_target_GSHHG_raster_mask_region_hits":15,
        "this_cohort_meets_gate":False,
        "explanation":"Only 12 source-positive modeled directions and four exact-target land-mask contacts. Both fail exploratory support gates; no current×origin/species-specific or shipping comparison regression was fitted."
      },
      "no_native_vs_alien_current_ship_regression_fitted":True,
      "source_only_2005_historiographic_2005_to_WCVP_vs_AIS2020_not_causal":True,
      "important_limits":[
        "A NOAA 1deg-grid old-OSCAR 90-day model near a shore does not observe any seed and does not measure physical colonisation probability.",
        "Original 2005 280 fixed source simulations cover 94/126 botanical island regions; missing 32 origins are not biologically unreachable.",
        "Only 12 directed source-eligible region pairs had modelled nearshore pass among 7313 evaluated; this does not indicate that all other ocean pair pathways lack natural dispersion at other years/horizons/resolutions.",
        "Four model raster contacts inside a different exact TDWG shape are not four viable seed landings.",
        "Contemporary WCVP accepted source-native to dest-introduced taxa are a geographic status combination, not observed plant origin/date or transport.",
        "Verified 2020 AIS port links cover only seven unordered original 126 botanical pairs, and could not directly cause 2005 modelled currents.",
        "Even a matched current-positive and ship-positive link does not demonstrate that a species was transported by either route.",
        "No confirmatory ecological cross-hypothesis model can be fit from the sparse independently observed transport source; keep earlier 22-island EuDiS +2.95pp as biased trait composition descriptive evidence."
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"NOAA2005_source_positive12_WCVP_flora_species_geographic_index.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(detail[0]));w.writeheader();w.writerows(detail)
    (a.out_dir/"NOAA2005_126_WCVP_species_source_crosswalk_feasibility.json").write_text(
      json.dumps(evidence,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("NOAA2005_WCVP126_ORIGIN_LINK_SOURCE_CROSSWALK",json.dumps({
      "status":evidence["status"],"current_source_positive_directions":len(detail),
      "exact_terrestrial_raster_contacts":len(actual_four),
      "late_nontrivial_first_nearshore_direction_pairs":noaa_late,
      "positive_directions_with_source_native_recipient_alien_WCVP_taxa":num_WCVP_ni_ali,
      "positive_directions_with_matched_actual_2020_AIS_source_OD":overlap2020,
      "native_Jaccard_median":evidence["source_positive_native_Jaccard_median"],
      "alien_Jaccard_median":evidence["source_positive_alien_Jaccard_median"],
      "effect_model_fitted":False
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
