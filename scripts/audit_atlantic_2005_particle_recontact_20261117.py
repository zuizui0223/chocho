#!/usr/bin/env python3
"""Independent source-artifact diagnostic for the fixed 2005 NOAA Lagrangian pilot.

Distinguishes near-origin land-mask recontact, later coastline collisions,
missing velocities, and valid duration. Source WPI proxy NEARNESS != landing.
All scenarios are deterministic source starts, NOT independent seed samples.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from collections import Counter,defaultdict
from pathlib import Path

from run_atlantic_noaa_2005_midpoint_particles_20261117 import latlon_dist

EXPECTED={"OBSERVED_NOAA_FLOW":{"GSHHG_FINE_MASK_LAND_CONTACT":19,
            "COMPLETED_90_DAYS":7,"SOURCE_VELOCITY_MISSING":1},
          "MATHEMATICAL_REVERSED_FLOW":{"GSHHG_FINE_MASK_LAND_CONTACT":10,
            "COMPLETED_90_DAYS":17},
          "STATIONARY_NO_CURRENT":{"COMPLETED_90_DAYS":27}}

def median_or_none(v):
    return round(float(statistics.median(v)),3) if v else None

def unique_result(root,name):
    xs=list(root.rglob(name))
    if len(xs)!=1:raise RuntimeError(f"Exactly one source artifact file {name} expected; found {len(xs)}")
    return xs[0]

def analyze(root,out):
    tracks=json.loads(unique_result(root,"particles_diagnostic_tracks.json").read_text())
    qa=json.loads(unique_result(root,"lagrangian_source_quality.json").read_text())
    cube=json.loads(unique_result(root,"source_uv_cube_receipt.json").read_text())
    if qa["status"]!="ACTUAL_2005_LEGACY_OSCAR_MIDPOINT_LAGRANGIAN_SOURCE_PILOT_NOT_SEED_ARRIVAL" or len(tracks)!=81:
        raise RuntimeError("Original executed 81 source trajectories missing or status changed")
    if len(qa["fixed_origins"])!=3 or sum(cube["grid_shape"])!=380:
        raise RuntimeError("Original NOAA 2005 source/fixed WPI cohort changed")
    scenarios=tuple(EXPECTED)
    by_scenario={s:[r for r in tracks if r.get("scenario")==s] for s in scenarios}
    if any(Counter(r["status"] for r in by_scenario[s])!=EXPECTED[s] for s in scenarios):
        raise RuntimeError("Original source categories changed (do not reinterpret missing as zero)")
    diagnostics={}
    for s,rs in by_scenario.items():
        by_status=Counter(r["status"] for r in rs)
        contacts=[r for r in rs if r["status"]=="GSHHG_FINE_MASK_LAND_CONTACT"]
        near_origin=[]
        for r in contacts:
            origin=qa["fixed_origins"][r["source"]]
            contact=r["trajectory"][-1]
            d=latlon_dist((float(origin["latitude"]),float(origin["longitude"])),
                          (float(contact["latitude"]),float(contact["longitude"])))
            near_origin.append({ "source":r["source"],"start":r["start_time_utc"],
                "failure_day":int(r["failure_day"]),
                "from_source_WPI_km":round(d,3),
                "land_contact_is_named_destination":False,
            })
        comp=[r for r in rs if r["status"]=="COMPLETED_90_DAYS"]
        diagnostics[s]={
           "total_tracks":len(rs),"statuses":dict(by_status),
           "survived_to_30days":sum(int(r["days_computed"])>=30 for r in rs),
           "survived_to_60days":sum(int(r["days_computed"])>=60 for r in rs),
           "survived_to_90days":sum(int(r["days_computed"])>=90 for r in rs),
           "median_launch_offset_from_source_WPI_km":median_or_none([r["source_offshore_offset_km"] for r in rs]),
           "median_displacement_completed_90days_km":median_or_none([r["displacement_km"] for r in comp]),
           "median_path_length_completed_90days_km":median_or_none([r["path_length_km"] for r in comp]),
           "proxy_25km_other_WPI_count":sum(r["within_25km_of_other_WPI_port_proxy"] for r in rs),
           "proxy_50km_other_WPI_count":sum(r["within_50km_of_other_WPI_port_proxy"] for r in rs),
           "coast_contacts":len(contacts),
           "coast_contacts_in_first_7days":sum(c["failure_day"]<=7 for c in near_origin),
           "coast_contacts_in_first_30days":sum(c["failure_day"]<=30 for c in near_origin),
           "coast_contacts_within_50km_of_source_WPI":sum(c["from_source_WPI_km"]<=50 for c in near_origin),
           "coast_contacts_within_100km_of_source_WPI":sum(c["from_source_WPI_km"]<=100 for c in near_origin),
           "median_coast_contact_day":median_or_none([c["failure_day"] for c in near_origin]),
           "median_land_contact_distance_from_source_WPI_km":median_or_none([c["from_source_WPI_km"] for c in near_origin]),
           "coast_contact_detail":near_origin
        }
    # source time and category counts must be paired across 3 scenarios, not separate biological experiments.
    by_src={}
    for src in ("AZO","CNY","CVI"):
        by_src[src]={s:dict(Counter(r["status"] for r in by_scenario[s] if r["source"]==src))
                     for s in scenarios}
    result={
      "schema":"chocho_atlantic_2005_true_oscar_particles_original_source_recontact_diagnostic_v1",
      "status":"ACTUAL_OSCAR_PARTICLE_SOURCE_RECONTACT_CLASSIFICATION_NOT_ISLAND_ARRIVAL",
      "parent_source_CI_run":43902616033,
      "source_replay_checks":{"tracks":len(tracks),"scenario_counts_match":True,
        "source_cube_sha256":qa["source_NetCDF_file_sha256"],
        "GSHHG_ocean_mask_SHA256":qa["source_GSHHG_fine_ocean_mask_SHA256"],
        "source_dates": [cube["time_start_utc"],cube["time_end_utc"]],
        "source_cube_missing_uv_cell_times":cube["missing_uv_cell_time"]
      },
      "scenario_diagnostics":diagnostics,"source_by_scenario_status":by_src,
      "decision_rules":[
        "Early GSHHG land contact is geographic mask contact; source recontacts cannot be promoted to new destination arrivals.",
        "The 0/27 WPI-destination proxy at 50km in this 90-day pilot is NOT a statement that ocean dispersal is absent.",
        "A stationary control and mathematically reversed velocity are implementation baselines, not independent randomization nulls.",
        "Observed source trajectories have 3 close NOAA offshore grid releases per port×3 release dates, and are highly nonindependent.",
        "GSHHG 0.0625deg coast mask and NOAA 1/3deg currents have unresolved coastal topology.",
        "WPI ports are geographic reference points, not all-island coastline or historical vessel connectivity.",
        "No seed viability, breeding butterflies, or plant establishment was measured."
      ]
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("ATLANTIC_PARTICLE_RECONTACT_DIAGNOSTIC",json.dumps({
        "status":result["status"],
        "scenarios":{s:{k:v[k] for k in
          ("survived_to_30days","survived_to_60days","survived_to_90days",
            "coast_contacts_in_first_7days","coast_contacts_in_first_30days",
            "coast_contacts_within_50km_of_source_WPI",
            "median_coast_contact_day","median_land_contact_distance_from_source_WPI_km",
            "proxy_50km_other_WPI_count")} for s,v in diagnostics.items()}
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    analyze(a.source_root,a.out)
