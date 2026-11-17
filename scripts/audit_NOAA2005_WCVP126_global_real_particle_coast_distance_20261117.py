#!/usr/bin/env python3
"""Fixed-source QA of actual 2005 NOAA global 1deg-decimated particle starts.

Original paper has WCVP 126 accepted-flora islands and 7875 same plant pairs.
Original NOAA real global pilot 94 source islands and 280 90-day maximum paths:
some numerical GSHHG land contacts and nearshore approaches (NOT arrivals).
This posthoc source-distance filter removes artificial releases far > shoreline
and reports which potential directed paths remain before any floral association.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path

THRESH=(25,50,80,100,150,200)
RELEASES=(907,919,931)

def unique(root,name):
    r=list(root.rglob(name))
    if len(r)!=1:raise RuntimeError(f"expected original unique artifact {name}, found {len(r)}")
    return r[0]

def load(root,name):
    return json.loads(unique(root,name).read_text(encoding="utf8"))

def csvrows(p):
    with p.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser()
    for arg in ("particle-root","source-preflight-root","coast-root","flora-root","out-dir"):
        p.add_argument("--"+arg,required=True,type=Path)
    a=p.parse_args()
    tracks=load(a.particle_root,"original126_NOAA2005_first90day_fixed_particle_tracks.json")
    receipt=load(a.particle_root,"NOAA2005_WCVP126_directed_source_coverage_and_interpretation.json")
    source=load(a.source_preflight_root,"NOAA2005_global126_coastal_current_source_eligibility.json")
    coast=load(a.coast_root,"GSHHG126_sea_geometry_source_coverage_audit.json")
    pair_file=unique(a.flora_root,"WCVP126_flora_similarity_strict_WPI2020_geographic_pairs.csv")
    pairs=csvrows(pair_file)
    if (len(pairs)!=7875 or len(tracks)!=280 or
        receipt["original_directed_126_minus_self_matrix_rows"]!=15750 or
        receipt["source_with_at_least_one_original_NOAA_2005_fixed_release"]!=94 or
        receipt["source_pairs_evaluable_at_least_one_30day_current_track"]!=7313):
        raise RuntimeError("Original source 126 island / 280 particle / 7313 directional-eligible cohort changed")
    if (coast["source_126_regions_with_at_least_one_verified_coastal_water_anchor"]!=104 or
        source["source_islands_valid_noaa_all_three_release_dates"]!=93):
        raise RuntimeError("Original NOAA2005 coastal source eligibility or land mask changed")
    loc={x["region"]:x for x in source["126_source_region_release_eligibility"]}
    recipients={k for k,v in coast["coastal_anchors_by_TDWG_region"].items() if v}
    if len(recipients)!=104 or len(loc)!=126:raise RuntimeError("Original 104 recipients/126 sources changed")
    floran={}
    for row in pairs:
        k=tuple(sorted((row["source_region"],row["destination_region"])))
        if k in floran:raise RuntimeError("Duplicated source-verified pair")
        floran[k]=row
    valid=[]
    for tr in tracks:
        code=tr["source"]
        time=tr["original_NOAA_source_time_index"]
        if time not in RELEASES:raise RuntimeError("Original source season indices changed")
        idx=RELEASES.index(time)
        d=loc[code]["nearest_valid_offshore_source_km_per_season"][idx]
        ind=loc[code]["original_NOAA_1deg_eligible_source_lon_lat_index_per_season"][idx]
        if d is None or ind is None:raise RuntimeError("Actual particle lacks matching original coastal source")
        if [int(v) for v in ind]!=tr["source_noaa_grid_row_col"]:
            raise RuntimeError("Particle original NOAA launch indices differ from independently source verified coastal starts")
        z=dict(tr)
        z["original_shore_to_NOAA_start_offset_km"]=float(d)
        valid.append(z)
    expected=Counter({
      "ORIGINAL_GSHHG_RASTER_LAND_CONTACT":45,
      "COMPLETED_90_DAYS":129,
      "NOAA_VELOCITY_MISSING_OR_OUTSIDE_FIELD":27,
      "NOAA_MIDPOINT_VELOCITY_MISSING_OR_OUTSIDE_FIELD":79})
    if Counter(t["status"] for t in valid)!=expected:raise RuntimeError("Original particle outcome categories changed")

    report_by_threshold={}
    rows=[]
    for km in THRESH:
        chosen=[x for x in valid if x["original_shore_to_NOAA_start_offset_km"]<=km]
        by_src=defaultdict(list)
        for x in chosen:by_src[x["source"]].append(x)
        qualifying_sources={s for s,xs in by_src.items() if any(x["complete_elapsed_days"]>=30 for x in xs)}
        all3={s for s,xs in by_src.items() if len(xs)==3}
        pairs_evaluable={(src,dst) for src in qualifying_sources for dst in recipients if src!=dst}
        near=set()
        contacts=set()
        for x in chosen:
            src=x["source"]
            for dst in x["nearshore_regions"]:
                if dst in recipients and dst!=src:
                    near.add((src,dst))
            p=x["original_GSHHG_raster_contact"]
            if p and p["is_other_TDWG3_raster_contact"]:
                dst=p["TDWG3_destination_exact"]
                if dst in recipients:contacts.add((src,dst))
        near30={pair for pair in near
                if any(x["complete_elapsed_days"]>=30 and x["source"]==pair[0]
                       and pair[1] in x["nearshore_regions"] for x in chosen)}
        near_eval=sorted(near & pairs_evaluable)
        contact_eval=sorted(contacts & pairs_evaluable)
        actual_evidence=[]
        for pair in near_eval:
            flora=floran.get(tuple(sorted(pair)))
            if flora is None:raise RuntimeError("Original flora source pair missing")
            actual_evidence.append({
                "source":pair[0],"recipient":pair[1],
                "native_Jaccard_WCVP":flora["native_Jaccard"],
                "alien_Jaccard_WCVP":flora["alien_Jaccard"],
                "interstate_geodesic_km":flora["centroid_geodesic_km"],
                "verified_source_2020_ships_between_pair":
                    int(flora["AIS2020_symmetric_direct_source_rows"])>0
            })
        assert len(actual_evidence)==len(near_eval)
        count={
          "nearshore_model_source_max_offshore_offset_km":km,
          "source_accepted_paths":len(chosen),
          "sources_with_at_least_one_path":len(by_src),
          "sources_with_three_release_dates":len(all3),
          "source_regions_with_at_least_one_30day_path":len(qualifying_sources),
          "original_recipients_with_GSHHG_coast_source":len(recipients),
          "evaluable_source_destination_directed_pairs":len(pairs_evaluable),
          "source_paths_90days_completed":sum(x["status"]=="COMPLETED_90_DAYS" for x in chosen),
          "source_paths_landmask_contacts":sum(x["status"]=="ORIGINAL_GSHHG_RASTER_LAND_CONTACT" for x in chosen),
          "source_paths_noaa_missing_or_field_end":sum("MISSING" in x["status"] for x in chosen),
          "all_paths_within50km_of_another_region_raw_directional_pairs":len(near),
          "evaluable_directed_pairs_with_other_recipient_within50km":len(near_eval),
          "evaluable_directed_pairs_with_other_recipient_within50km_and_one_30day_path":len(near30 & pairs_evaluable),
          "evaluable_other_region_exact_contact_source_pairs":len(contact_eval),
          "source_real_current_nearshore_positive_pairs_evaluable":actual_evidence,
          "no_empirical_biological_arrival_estimated":True
        }
        report_by_threshold[str(km)]=count
        rows.append({k:v for k,v in count.items() if k!="source_real_current_nearshore_positive_pairs_evaluable"})
    if report_by_threshold["200"]["source_accepted_paths"]!=280:
        raise RuntimeError("Original source 280 trajectories missing after 200km threshold")
    if len(report_by_threshold["200"]["source_real_current_nearshore_positive_pairs_evaluable"])!=12:
        raise RuntimeError("Original 12 pair nearshore pilot source validation changed")
    for x,y in zip(THRESH,THRESH[1:]):
        if report_by_threshold[str(x)]["source_accepted_paths"]>report_by_threshold[str(y)]["source_accepted_paths"]:
            raise RuntimeError("Offshore source eligibility must be monotonic")
    r={
      "schema":"chocho_original_126_WCVP_island_noaa_2005_1deg_source_50km_shore_distance_restriction_v1",
      "status":"SOURCE_VERIFIED_NOAA_GLOBAL_CURRENT_LAGRANGIAN_OFFSHORE_DISTANCE_CENSITIVITY_NOT_REAL_SEED_ARRIVAL",
      "source_original_global_53step_particle_CI":44012093616,
      "source_original_NOAA_three_release_date_coastal_preflight_CI":44011612861,
      "source_original_GSHHG126_104_coast_anchors_CI":44010673580,
      "source_original_126island_7875_flora_pairs_CI":44010307551,
      "source_126_island_full_directed_matrix":15750,
      "source_200km_280_particle_source_pilot_verified":True,
      "source_200km_12_nearshore_pair_source_pilot_verified":True,
      "offshore_distance_thresholds_km":list(THRESH),
      "results":report_by_threshold,
      "limits":[
        "Coastal source distance 50km restriction is a posthoc sensitivity; source starts remain artificially chosen from 1-degree sampled 15m NOAA ocean pixels, not real floating seed release locations.",
        "The 90day particle paths have 2005 three date historical sources and cannot reconstruct contemporary WCVP plant introduction timing.",
        "A 50km approach to TDWG source geography is a modeled offshore proximity diagnostic, never seed landing nor butterfly range expansion.",
        "The original WCVP undirected native/introduced Jaccard cannot identify direction or age of any botanical migration; matched flora scores on few model-positive paths are descriptive.",
        "Only 7 exact WPI source-linked AIS2020 actual observed ship OD pairs among all 126 botanical islands, no matched historical cargo.",
        "Missing original source and recipient shoreline, censored flows or missing field are UNKNOWN, never disallowing natural transport."
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"NOAA126_global_current_2005_source50km_restriction_audit.json").write_text(
      json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (a.out_dir/"NOAA126_global_current_source_distance_quality_counts.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print("NOAA2005_GLOBAL126_50KM_DIRECTIONAL_COAST_PAIR_IDENTITIES",json.dumps({
        "max_source_offshore_km":50,
        "positive_directed_pairs_WCVP_flora_overlap":
            report_by_threshold["50"]["source_real_current_nearshore_positive_pairs_evaluable"]
    },ensure_ascii=False),flush=True)
    print("NOAA2005_GLOBAL126_TRUE_CURRENTS_OFFSHORE_RESTRICTION",json.dumps({
       "status":r["status"],
       "by_distance":{k:{
        "n_model_paths":v["source_accepted_paths"],
        "n_current_evaluable_directions":v["evaluable_source_destination_directed_pairs"],
        "n_valid_nearshore_directional_pairs":v["evaluable_directed_pairs_with_other_recipient_within50km"],
        "n_valid_other_region_mask_contact_pairs":v["evaluable_other_region_exact_contact_source_pairs"],
        "n_original_sea_30day_source_regions":v["source_regions_with_at_least_one_30day_path"]
       } for k,v in report_by_threshold.items()}
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
