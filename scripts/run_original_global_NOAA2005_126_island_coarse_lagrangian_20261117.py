#!/usr/bin/env python3
"""Directional coarse historical NOAA old-OSCAR particle reachability pilot.

Start from original verified TDWG/GSHHG coastal anchors (NOT WPI ports or WCVP
flora observations), using the original NOAA global 1° DECIMATION of true
1/3° old OSCAR gridded 2005 15-m actual velocities, three source dates.
RK2 daily Eulerian displacement for at most 90 days. Stop at original GSHHG
raster land contact, NOAA field/velocity missing, or time limit. Also compute
50km proximity to independently pinned TDWG3 recipient shore polygons.

The matrix is directional MODELED PROXIMITY over fixed pilot releases, NOT
ocean-current physical arrival probability, natural seed buoyancy/float time,
historical ship OD, biological plant exchange, butterfly establishment, or
a global 126-island measured colonisation mechanism.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
from pyproj import Geod
from shapely.geometry import Point,shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree

from run_atlantic_noaa_2005_midpoint_particles_20261117 import (
 uv_interpolate,integrate_vector,crosses_land,latlon_dist)

G=Geod(ellps="WGS84")
MAX_DAYS=90
DAY=86400
RELEASE_INDICES=(907,919,931)
START_TIME_INDEX=898
GRAPH_SHA="79a25cb1491c30bc5542333bd4020535bc88772c98f8aa8fa1fe9c12f7b8e906"
SOURCE_PREP_CI=44011612861
SOURCE_GSHHG126_CI=44010673580
MAX_PROXIMITY_KM=50.
MIN_SEPARATION_KM=0.
MAX_SOURCE_DATA_GAP_DAYS=30

def SHA(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def load_geo(path,region_codes):
    source=json.loads(path.read_text(encoding="utf-8"))
    reg={}
    for feat in source["features"]:
        code=str((feat.get("properties") or {}).get("LEVEL3_COD","")).strip()
        if code in region_codes:
            if code in reg:raise RuntimeError("Duplicate original TDWG3 botanical polygon")
            reg[code]=shape(feat["geometry"])
    if set(reg)!=set(region_codes):raise RuntimeError("Missing original 126 official botanical polygons")
    codes=sorted(reg)
    geoms=[reg[c] for c in codes]
    return codes,geoms,STRtree(geoms)

def candidate_nearby_regions(pt,source,region_codes,reg,tree):
    # Probe all potential shores via R-tree bounding box with width 4deg,
    # generous enough for 50km at latitude up to 80N/S. Actual geodesic
    # nearest distances (not degree units) make the final decision.
    geom=pt.buffer(4.0)
    out={}
    for ix in set(int(j) for j in tree.query(geom)):
        code=region_codes[ix]
        if code==source:continue
        pol=reg[ix]
        near=nearest_points(pt,pol)[1]
        _,_,meters=G.inv(pt.x,pt.y,near.x,near.y)
        dk=float(abs(meters))/1000
        if dk<=MAX_PROXIMITY_KM:
            out[code]=min(out.get(code,float("inf")),dk)
    return out

def exact_land_region(pt,region_codes,reg,tree):
    contained=[region_codes[int(ix)] for ix in tree.query(pt)
               if reg[int(ix)].covers(pt)]
    if len(contained)==1:return contained[0],"EXACT_TDWG3_POLYGON_CONTAINS_CONTACT"
    if len(contained)>1:return None,"MULTIPLE_TDWG3_REGIONS_AMBIGUOUS"
    return None,"GSHHG_MASK_CONTACT_OUTSIDE_TDWG3_REGIONS"

def run_path(code,idx,original_source,source,mask,region_codes,geoms,tree):
    time_idx=idx-START_TIME_INDEX
    latAxis=source["latitude_deg"]
    lonAxis=source["longitude_deg_0_420"]
    rr=int(original_source[0]);cc=int(original_source[1])
    if not (0<=rr<len(latAxis) and 0<=cc<len(lonAxis)):
        raise RuntimeError("Frozen NOAA launch indices outside source grid")
    lat=float(latAxis[rr])
    lon=float(((lonAxis[cc]+180)%360)-180)
    t=float(source["time_seconds_utc"][time_idx])
    original=(lat,lon)
    current=original
    nearshore={} # code -> earliest day, shortest distance and max track displacement
    trajectory=[{"day":0,"lat":lat,"lon":lon}]
    status="COMPLETED_90_DAYS"
    stopped_day=None
    contact=None
    traveled=0.
    for day in range(1,MAX_DAYS+1):
        sec=t+(day-1)*DAY
        vel=uv_interpolate(source["time_seconds_utc"],latAxis,lonAxis,
               source["u_m_s"],source["v_m_s"],sec,current[0],current[1])
        if vel is None:
            status="NOAA_VELOCITY_MISSING_OR_OUTSIDE_FIELD"
            stopped_day=day
            break
        halfway=integrate_vector(*current,*vel,DAY/2)
        mid=uv_interpolate(source["time_seconds_utc"],latAxis,lonAxis,
               source["u_m_s"],source["v_m_s"],sec+DAY/2,*halfway)
        if mid is None:
            status="NOAA_MIDPOINT_VELOCITY_MISSING_OR_OUTSIDE_FIELD"
            stopped_day=day
            break
        nxt=integrate_vector(*current,*mid,DAY)
        hit,loc=crosses_land(mask,current,nxt)
        if hit:
            status="ORIGINAL_GSHHG_RASTER_LAND_CONTACT"
            stopped_day=day
            p=Point(float(loc[1]),float(loc[0]))
            dest,land_state=exact_land_region(p,region_codes,geoms,tree)
            contact={"contact_day":day,"lon":p.x,"lat":p.y,
                     "TDWG3_destination_exact":dest,
                     "classification":land_state,
                     "is_other_TDWG3_raster_contact":dest is not None and dest!=code}
            break
        traveled+=latlon_dist(current,nxt)
        current=nxt
        # Evaluate destination shoreline proximity ONLY after valid source
        # advection, never treat a source's initial offshore proximity as
        # a new landfall or a shipborne introduction event.
        q=Point(current[1],current[0])
        near=candidate_nearby_regions(q,code,region_codes,geoms,tree)
        for dest,km in near.items():
            v=nearshore.setdefault(dest,{"earliest_day":day,
               "closest_km":float(km),"first_displacement_from_source_km":
                latlon_dist(original,current),"minimum_day_after_source_start":day})
            v["closest_km"]=min(v["closest_km"],float(km))
            v["earliest_day"]=min(v["earliest_day"],day)
        if day%5==0 or day in (30,60,90):
            trajectory.append({"day":day,"lat":float(current[0]),"lon":float(current[1])})
    return {
      "source":code,"original_NOAA_source_time_index":idx,
      "source_actual_date_utc":datetime.fromtimestamp(t,tz=timezone.utc).isoformat(),
      "start_lat":original[0],"start_lon":original[1],
      "source_noaa_grid_row_col":[rr,cc],
      "status":status,"complete_elapsed_days":MAX_DAYS if stopped_day is None else stopped_day-1,
      "termination_attempt_day":stopped_day,
      "original_GSHHG_raster_contact":contact,
      "net_displacement_km":round(latlon_dist(original,current),4),
      "valid_water_path_length_km":round(traveled,4),
      "number_of_other_botanical_TDWG_regions_neared_within50km":len(nearshore),
      "nearshore_regions":nearshore,
      "trajectory_every_5days":trajectory,
      "is_confirmed_propagule_seed_or_butterfly_arrival":False
    }

def main():
    ap=argparse.ArgumentParser()
    for a in ("noaa-cube","noaa-source-receipt","noaa-preflight","original-GSHHG-graph",
              "126-coast-audit","tdwg3","original-WCVP-islands","out-dir"):
        ap.add_argument("--"+a,type=Path,required=True)
    a=ap.parse_args()
    if SHA(a.original_GSHHG_graph)!=GRAPH_SHA:
        raise RuntimeError("Original global GSHHG source SHA changed")
    previous=json.loads(a.noaa_preflight.read_text())
    source_receipt=json.loads(a.noaa_source_receipt.read_text())
    if (previous["source_islands_valid_noaa_at_least_one_release_date"]!=94 or
        source_receipt["source_previous_coastal_dates_preflight_CI"]!=SOURCE_PREP_CI or
        SHA(a.noaa_cube)!=source_receipt["source_NPZ_sha256"]):
        raise RuntimeError("Original NOAA global source-preflight or 53-grid source SHA changed")
    anchors=json.loads(a.coast_audit.read_text())
    if anchors["source_126_regions_with_at_least_one_verified_coastal_water_anchor"]!=104:
        raise RuntimeError("Original source coastline geometry source changed")
    geo_rows=list(csv.DictReader(a.original_WCVP_islands.open()))
    region_codes=sorted(r["code"] for r in geo_rows)
    if len(region_codes)!=126:raise RuntimeError("Original accepted 126 WCVP region IDs changed")
    region_codes,geoms,tree=load_geo(a.tdwg3,region_codes)
    original_anchor_available=set(c for c,x in anchors["coastal_anchors_by_TDWG_region"].items() if x)
    with np.load(a.original_GSHHG_graph) as g:
        mask=g["fine_cell_to_node"]
        if mask.shape!=(2720,5760):raise RuntimeError("Original GSHHG land-mask shape changed")
    with np.load(a.noaa_cube) as q:
        src={key:q[key] for key in (
             "time_seconds_utc","latitude_deg","longitude_deg_0_420","u_m_s","v_m_s")}
    if src["u_m_s"].shape!=(53,161,401) or src["v_m_s"].shape!=(53,161,401):
        raise RuntimeError("Wrong original 2005 old OSCAR global source cube")
    prep={x["region"]:x for x in previous["126_source_region_release_eligibility"]}
    if set(prep)!=set(region_codes):raise RuntimeError("Original source 126 launch count changed")
    paths=[]
    for code in region_codes:
        p=prep[code]
        locs=p.get("original_NOAA_1deg_eligible_source_lon_lat_index_per_season")
        if not locs:continue
        for idx,loc in zip(RELEASE_INDICES,locs):
            if loc is None:continue
            z=run_path(code,idx,loc,src,mask,region_codes,geoms,tree)
            paths.append(z)
    expected=sum(x["fixed_season_validities"] for x in prep.values())
    if len(paths)!=expected:
        raise RuntimeError("Original independent source-time 126 NOAA trajectory population drift")
    totals=Counter(p["status"] for p in paths)
    sources={x["source"] for x in paths}
    if len(sources)!=94:raise RuntimeError("Original global NOAA preflight source islands eligible changed")
    dest_record=defaultdict(lambda:{"near_release_ids":set(),"exact_contact_ids":set(),
                                    "first_nearshore_day":None})
    for r in paths:
        source=r["source"];reln=str(r["original_NOAA_source_time_index"])
        for dest,x in r["nearshore_regions"].items():
            k=(source,dest)
            dest_record[k]["near_release_ids"].add(reln)
            y=dest_record[k]["first_nearshore_day"]
            dest_record[k]["first_nearshore_day"]=min(y,x["earliest_day"]) if y is not None else x["earliest_day"]
        contact=r["original_GSHHG_raster_contact"]
        if contact and contact["is_other_TDWG3_raster_contact"]:
            dest_record[(source,contact["TDWG3_destination_exact"])]["exact_contact_ids"].add(reln)
    rows=[]
    for srcname in region_codes:
        tries=prep[srcname]["fixed_season_validities"]
        pv=[r for r in paths if r["source"]==srcname]
        adequate=sum(r["complete_elapsed_days"]>=MAX_SOURCE_DATA_GAP_DAYS for r in pv)
        for dest in region_codes:
            if srcname==dest:continue
            audited_target=dest in original_anchor_available
            if not tries:stat="NO_VERIFIED_NOAA_SOURCE_START"
            elif not audited_target:stat="UNRESOLVED_GSHHG_DESTINATION_SHORE"
            elif adequate==0:stat="ALL_SOURCE_TRACKS_CENSORED_BEFORE_30D"
            else:stat="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE"
            rec=dest_record.get((srcname,dest))
            nnear=len(rec["near_release_ids"]) if rec else 0
            ncontact=len(rec["exact_contact_ids"]) if rec else 0
            rows.append({
              "source":srcname,"destination":dest,
              "original_source_releases_available":tries,
              "original_source_releases_with_valid_current_at_least30d":adequate,
              "original_target_GSHHG_raster_shore_valid":audited_target,
              "model_physical_source_evaluation_status":stat,
              "NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km":nnear if stat=="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE" else "",
              "NOAA_2005_fixed90d_releases_contacting_another_exact_TDWG3_polygon":ncontact if stat=="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE" else "",
              "model_50km_proximity_count_fraction_valid_start":round(nnear/tries,6) if stat=="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE" else "",
              "physical_current_is_old_OSCAR_1deg_source_decimation":True,
              "no_observed_biological_arrival_or_shipping":True
            })
    if len(rows)!=126*125:raise RuntimeError("Original directed 126 island source/dest matrix shape differs")
    evaluable=[r for r in rows if r["model_physical_source_evaluation_status"]=="COARSE_NOAA_2005_DIRECTIONAL_PROXIMITY_EVALUABLE"]
    status=Counter(r["model_physical_source_evaluation_status"] for r in rows)
    if not evaluable:raise RuntimeError("No original NOAA2025 historic physical source-pair outcome recoverable")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"original126_NOAA2005_first90day_fixed_particle_tracks.json").write_text(
       json.dumps(paths,ensure_ascii=False,indent=2)+"\n")
    with (a.out_dir/"original126_NOAA2005_directed_source_nearshore_matrix.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report={
      "schema":"chocho_WCVP126_original_NOAA_old_OSCAR_2005_decimated_real_currents_directional_source_proximity_v1",
      "status":"REAL_2005_GLOBAL_NOAA_LAGRANGIAN_1DEG_FEASIBILITY_SOURCE_ONLY_NOT_HISTORIC_PLANT_SEED_DELIVERY",
      "source_GSHHG_original_graph_SHA256":GRAPH_SHA,
      "source_previous_preflight_CI":SOURCE_PREP_CI,
      "source_original_GSHHG_126_polygon_geometry_CI":SOURCE_GSHHG126_CI,
      "source_NetCDF_53_timestep_sha256":source_receipt["source_NPZ_sha256"],
      "NOAA_grid_shape":[53,161,401],
      "NOAA_not_NASA_daily_OSCAR_v2":True,
      "source_decimation_of_original_NOAA_1over3_degree":3,
      "source_original_126_WCVP_regions":126,
      "source_original_GSHHG_coastal_polygon_anchorable_regions":104,
      "source_with_at_least_one_original_NOAA_2005_fixed_release":len(sources),
      "fixed_release_true_NOAA_timestamps":source_receipt["fixed_release_dates"],
      "actual_source_particle_tracks":len(paths),
      "original_source_velocity_censored_and_90day_stops":dict(totals),
      "original_directed_126_minus_self_matrix_rows":len(rows),
      "directed_pair_source_coverage_status":dict(status),
      "source_pairs_evaluable_at_least_one_30day_current_track":len(evaluable),
      "evaluable_source_dest_pairs_with_any_50km_offshore_nearshore_pass":
          sum(int(r["NOAA_2005_fixed90d_releases_near_TDWG_shore_within50km"])>0 for r in evaluable),
      "evaluable_source_dest_pairs_with_any_GSHHG_land_contact_inside_other_exact_TDWG3":
          sum(int(r["NOAA_2005_fixed90d_releases_contacting_another_exact_TDWG3_polygon"])>0 for r in evaluable),
      "source_coastal_geometry_coverage_cannot_assume_missing_equal_no_access":True,
      "scientific_limits":[
        "Real historical source NOAA OSCAR 1/3degree sampled every 3rd original cell to 1degree. This coarsening changes coastal currents and is NOT a full 1/3degree run.",
        "Source launch sites may be up to 200km offshore from the original geographically anchored coast; these are artificial model starts, not real floating seed locations.",
        "Only three fixed 2005 dates, maximum90d deterministic RK2 tracks, not season/year dispersion probabilities or known viable-seed float durations.",
        "True observed u/v current at 15m is not the actual plant buoyant surface drift speed and omits tides, waves, Stokes drift, windage, mortality and local estuarine flow.",
        "Proximity to an official TDWG region within50km is NOT a verified onshore landing, plant arrival or butterfly colonisation.",
        "Original GSHHG 0.0625degree mask collision at another TDWG polygon is a model geometry event, not viable propagule landing.",
        "No NOAA source records for 32/126 islands or for destination shorelines absent in the GSHHG model should NEVER be imputed to source-directed connectivity zero.",
        "The computed source/destination fractions are 1–3 deterministic simulation starts, NOT empirical population-scale probabilities.",
        "The contemporary WCVP native/alien shared flora records pre- or post-date the 2005 current sample unpredictably; current vs modern plants is spatial association only.",
        "Real 2020 AIS strict port origin-destination exists for only7 botanical island-pairs and postdates 2005; no causal ship-current comparative effect can be established."
      ]
    }
    (a.out_dir/"NOAA2005_WCVP126_directed_source_coverage_and_interpretation.json").write_text(
       json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("NOAA2005_126_ISLAND_REAL_CURRENT_LAGRANGIAN_PILOT",json.dumps({
       "status":report["status"],"source_regions":len(sources),
       "fixed_source_release_trajectories":len(paths),
       "trajectory_outcomes":dict(totals),
       "directed_source_pairs_evaluable":len(evaluable),
       "source_pair_coverage_states":dict(status),
       "source_pair_current_nearshore_positive":report["evaluable_source_dest_pairs_with_any_50km_offshore_nearshore_pass"],
       "source_pair_exact_another_TDWG_mask_contact_positive":report["evaluable_source_dest_pairs_with_any_GSHHG_land_contact_inside_other_exact_TDWG3"],
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
