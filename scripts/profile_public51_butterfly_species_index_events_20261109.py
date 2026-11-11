#!/usr/bin/env python3
"""Which butterflies are indexed later and how far from prior conspecific indices?

This descriptive butterfly-first evidence panel reuses the exact verified
225 x 51 public 2026-11-09 GBIF record index, frozen source geography, and
LepTraits S1. Candidate prior-indexed island is NOT the actual origin and
first GBIF indexed record is NOT first butterfly arrival, dispersal,
establishment, breeding, or larval host use.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
import math
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
ANALYSIS=ROOT/"scripts/run_public51_butterfly_trait_access_20261109.py"
PANEL=ROOT/"data/exploratory/public51_butterfly_index_20261109"
PERIODS=("2011_2017","2018_2022")

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise ValueError("Missing original source model")
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def rank_higher_is_closer(x,distribution):
    values=np.asarray(distribution,dtype=float)
    # Half ties: ranks valid even when all at-risk islands are equidistant.
    return float((np.count_nonzero(values>x)+.5*np.count_nonzero(values==x))/len(values))

def percentile_between(x,distribution):
    values=np.asarray(distribution,dtype=float)
    return float((np.count_nonzero(values<x)+.5*np.count_nonzero(values==x))/len(values))

def write_csv(path,rows):
    if not rows:raise ValueError("Empty source evidence")
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def run(out):
    m=load_module(ANALYSIS,"source51")
    panel,receipt=m.load_source(PANEL)
    src,trait,data,source_pins,bridge,bridge_summary=m.data_bridge()
    ids={r["landmass_id"] for r in panel}
    clusters=m.geographic_blocks(data,ids)
    pair={}
    for r in data["ocean"]:
        key=tuple(sorted((r["island_a"],r["island_b"]),key=int))
        pair[key]=(float(r["distance_km_verified"]),
                   float(r["uniform_distributed_port_resistance"]))
    if len(pair)!=1275:
        raise ValueError("Not the original 51x51 Circuitscape matrix")
    source_index={p:{ } for p in PERIODS}
    for period in PERIODS:
        by_species=defaultdict(list)
        for z in panel:
            old=int(z["pre2001_2010_indexed_records"])
            mid=int(z["future2011_2017_indexed_records"])
            if old>0 or (period=="2018_2022" and mid>0):
                by_species[z["speciesKey"]].append(z["landmass_id"])
        source_index[period]={sp:set(islands) for sp,islands in by_species.items()}
    all_species=[]
    all_cases=[]
    meta={}
    for period in PERIODS:
        rows,cohort=m.make_rows(panel,period,data,bridge,clusters)
        assert len(rows)=={"2011_2017":4374,"2018_2022":4933}[period]
        assert sum(z["record_index_case"] for z in rows)=={"2011_2017":187,"2018_2022":83}[period]
        traits_complete={s:z for s,z in bridge.items() if z.get("trait_available")}
        at_risk=defaultdict(list)
        distance_all=np.array([math.expm1(r["features"]["min_log_geodesic_km"]) for r in rows])
        distances_by_species=defaultdict(list)
        for r in rows:
            d=float(math.expm1(r["features"]["min_log_geodesic_km"]))
            distances_by_species[r["species"]].append(d)
        for r in rows:
            sp=r["species"];key=r["speciesKey"];target=r["landmass_id"]
            candidate_sources=source_index[period].get(key,set())-{target}
            if not candidate_sources:
                raise ValueError("Selected prior-source risk without prior source")
            near=min(candidate_sources,key=lambda other:pair[tuple(sorted((target,other),key=int))][0])
            dist,resistance=pair[tuple(sorted((target,near),key=int))]
            if not np.isclose(dist,math.expm1(r["features"]["min_log_geodesic_km"]),rtol=1e-10,atol=1e-8):
                raise ValueError("Nearest indexed source and original model distance disagree")
            z=bridge[sp]
            item={
                "period":period,"species":sp,"speciesKey":key,
                "target_landmass_id":target,
                "indexed_before_window_elsewhere":len(candidate_sources),
                "candidate_nearest_prior_index_landmass_id":near,
                "candidate_nearest_distance_km":round(float(dist),5),
                "min_prior_source_uniform_sea_resistance":round(
                    float(math.expm1(r["features"]["min_log_uniform_effective_resistance"])),7),
                "candidate_is_actual_biological_origin":"FALSE",
                "GBIF_first_index_not_colonization":"TRUE",
                "wing_size_proxy":z["wing_size_proxy"],
                "wing_median_group":z["wing_median_group"],
                "disturbance_class":int(z["disturbance_category"]),
                "at_risk_geodesic_proximity_percentile":round(
                    rank_higher_is_closer(dist,distances_by_species[sp]),6),
                "initial_index_event":int(r["record_index_case"]),
                "source_index_effort_log":round(float(
                    r["features"]["log_target_before_window_record_effort"]),7),
            }
            at_risk[sp].append(item)
            if r["record_index_case"]:
                all_cases.append(item)
        meta[period]={"risk_rows":len(rows),"indexed_events":sum(r["record_index_case"] for r in rows),
                      "event_species":len({z["species"] for z in all_cases if z["period"]==period}),
                      "median_closest_source_distance_km_of_events":
                          float(np.median([z["candidate_nearest_distance_km"]
                                          for z in all_cases if z["period"]==period]))}
        for sp, rr in at_risk.items():
            pos=[x for x in rr if x["initial_index_event"]]
            if not pos:continue
            z=traits_complete[sp]
            all_species.append({
                "period":period,"species":sp,"trait_wing_size_proxy":z["wing_size_proxy"],
                "wing_median_group":z["wing_median_group"],
                "disturbance_class":int(z["disturbance_category"]),
                "candidate_risk_islands":len(rr),"new_indexed_islands":len(pos),
                "median_nearest_prior_index_distance_km_in_risk":round(
                    float(np.median([x["candidate_nearest_distance_km"] for x in rr])),4),
                "median_nearest_prior_index_distance_km_at_events":round(
                    float(np.median([x["candidate_nearest_distance_km"] for x in pos])),4),
                "mean_event_within_species_proximity_percentile":round(
                    float(np.mean([x["at_risk_geodesic_proximity_percentile"] for x in pos])),5),
                "species_event_sample_is_observational_not_dispersal":"TRUE",
            })
    if len(all_cases)!=270:raise ValueError("Wrong original species-by-island case total")
    late=[z for z in all_species if z["period"]=="2018_2022"]
    late_sorted=sorted(late,key=lambda z:(-z["new_indexed_islands"],
                                           -z["mean_event_within_species_proximity_percentile"],z["species"]))
    # Tail distances are relative to other FIRST-INDEXED species-island records,
    # NOT a threshold for true dispersal events or directed island hops.
    late_distances=[z["candidate_nearest_distance_km"] for z in all_cases
                    if z["period"]=="2018_2022"]
    threshold=float(np.quantile(late_distances,.9))
    tail=[z for z in all_cases if z["period"]=="2018_2022"
          and z["candidate_nearest_distance_km"]>=threshold]
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/"butterfly_first_index_candidate_source_events_270.csv",all_cases)
    write_csv(out/"butterfly_species_indexed_distance_profiles.csv",all_species)
    write_csv(out/"butterfly_late_distance_tail_source_candidates.csv",tail)
    audit={
        "schema":"chocho_butterfly_first_index_nearest_prior_record_source_profiles_20261109",
        "data_source":"independently reconstructed 2026-11-09 GBIF occurrence-index facets, NOT a biological dispersal census",
        "original_51_source_sha256":receipt["merged_table_sha256"],
        "trait_source_git_blob":trait.TRAITS[1],
        "original_source_matrix_git_blob":src.SOURCES["ocean"][1],
        "record_index_cohort":meta,
        "per_species_event_profiles":len(all_species),
        "late_new_indexed_distances_90th_percentile_km":threshold,
        "late_90th_percentile_event_count":len(tail),
        "late_most_recorded_species":late_sorted[:15],
        "late_long_distance_candidate_events":tail,
        "limits":[
            "Nearest earlier GBIF-indexed other island is merely a candidate source, not a historical point of origin",
            "Modern port coordinates approximate pair distances; no actual ship or butterfly route",
            "The minimum source distance is among 51 source-selected GSHHG landmasses, not all geographical landmasses",
            "GBIF non-indexed islands are not demonstrated biological absences",
            "Wing size is a LepTraits composite proxy, not independently measured travel capacity",
            "Species with multiple target islands are not independent biological dispersal events",
            "This descriptive rank does not control target GBIF survey effort or island biogeography",
            "Any apparent geographically long-distance source pair can be caused by missing indexing closer to target",
            "No historical larval host record at the butterfly first-index site is checked",
        ]}
    (out/"butterfly_candidate_source_event_audit.json").write_text(
        json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":"OBSERVED_GBIF_INDEX_NOT_DISPERSAL",
        "species_profile_rows":len(all_species),"event_rows":len(all_cases),
        "meta":meta,"late_top5":late_sorted[:5],
        "late_event_long_distance_threshold_km":round(threshold,2),
        "late_tail_count":len(tail)},ensure_ascii=False,indent=2))

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out-dir",type=Path,required=True)
    v=a.parse_args()
    run(v.out_dir)
