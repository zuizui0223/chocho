#!/usr/bin/env python3
"""Offline, source-artifact-only bias/footprint audit of WPI79 matched botany.

No fresh GBIF calls and no re-selection after observing flora. The only
inputs are the frozen coastal geometry and its paired botanical observations.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


def haversine_km(lat1,lon1,lat2,lon2):
    p1,p2=math.radians(lat1),math.radians(lat2)
    dp=p2-p1
    dl=math.radians(lon2-lon1)
    t=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 6371.0088*2*math.asin(math.sqrt(min(1,max(0,t))))


def group_island_mean(pairs,key):
    values=defaultdict(list)
    for x in pairs:
        if x[key] is not None:
            values[x["island_id"]].append(x[key])
    return {k:statistics.mean(v) for k,v in values.items()}


def bootstrap_islands(means,*,draws=9999,seed=81273):
    v=list(means.values())
    if not v: return None
    rng=random.Random(seed)
    estimates=sorted(statistics.mean([rng.choice(v) for _ in v]) for _ in range(draws))
    return [estimates[round((len(estimates)-1)*.025)],
            estimates[round((len(estimates)-1)*.975)]]


def summarize_period(geometry,observations,period,*,radius_km=5):
    selected=[x for x in observations if x["period"]==period]
    sites=defaultdict(dict)
    for x in selected:
        key=(str(x["island_id"]),int(x["port_id"]))
        if x["site_type"] in sites[key]:
            raise ValueError("Duplicated site-period")
        sites[key][x["site_type"]]=x
    pairs=[]
    for key,item in sorted(sites.items(),key=lambda a:(int(a[0][0]),a[0][1])):
        if set(item)!={"port","nonport"}: raise ValueError("Missing paired observation "+str(key))
        if key not in geometry: raise ValueError("Unknown coastline site ID "+str(key))
        g=geometry[key]
        distance=haversine_km(float(g["port_sample_lat"]),float(g["port_sample_lon"]),
                              float(g["nonport_lat"]),float(g["nonport_lon"]))
        no_overlap=distance>=2*radius_km
        a,b=item["port"],item["nonport"]
        general_quality=all(x["all_plant_records"]>=30 for x in (a,b))
        trait_quality=general_quality and all(x["candidate_hosts_with_N_indexed_species"]>=2 for x in (a,b))
        ratio=(min(a["all_plant_records"],b["all_plant_records"])/
               max(a["all_plant_records"],b["all_plant_records"])) if max(
                    a["all_plant_records"],b["all_plant_records"])>0 else 0
        pairs.append({
            "island_id":key[0],"port_id":key[1],
            "site_centres_geodesic_proxy_separation_km":distance,
            "five_km_radius_circular_buffers_disjoint":no_overlap,
            "general_plant_effort_both_30plus":general_quality,
            "trait_usable_both_2plus":trait_quality,
            "general_record_count_ratio_min_over_max":ratio,
            "port_general_species_key_count":a["all_plant_species_keys"],
            "nonport_general_species_key_count":b["all_plant_species_keys"],
            "port_N_host_species_count":a["candidate_hosts_with_N_indexed_species"],
            "nonport_N_host_species_count":b["candidate_hosts_with_N_indexed_species"],
            "port_minus_nonport_mean_N":(
                a["candidate_hosts_mean_N_unweighted_by_occurrences"]-
                b["candidate_hosts_mean_N_unweighted_by_occurrences"]
                if trait_quality else None),
            "log_general_plant_species_ratio":(
                math.log1p(a["all_plant_species_keys"])-math.log1p(b["all_plant_species_keys"])
                if general_quality else None),
            "log_N_host_species_ratio":(
                math.log1p(a["candidate_hosts_with_N_indexed_species"])-
                math.log1p(b["candidate_hosts_with_N_indexed_species"])
                if general_quality else None),
        })
    eligible=[x for x in pairs if x["trait_usable_both_2plus"]]
    nonoverlap=[x for x in eligible if x["five_km_radius_circular_buffers_disjoint"]]
    effort=[x for x in nonoverlap if x["general_record_count_ratio_min_over_max"]>=.2]
    def summary(sample):
        m=group_island_mean(sample,"port_minus_nonport_mean_N")
        vals=list(m.values())
        return {
            "eligible_pairs":len(sample),
            "eligible_independent_islands":len(m),
            "island_weighted_mean_port_minus_nonport_N":statistics.mean(vals) if vals else None,
            "median_island_mean_port_minus_nonport_N":statistics.median(vals) if vals else None,
            "island_bootstrap_descriptive_95":bootstrap_islands(m),
            "positive_islands":sum(v>0 for v in vals),
            "negative_islands":sum(v<0 for v in vals),
            "no_hypothesis_test_or_causal_inference":True,
        }
    general=[x for x in pairs if x["general_plant_effort_both_30plus"] and x["five_km_radius_circular_buffers_disjoint"]]
    general_species=group_island_mean(general,"log_general_plant_species_ratio")
    host_species=group_island_mean(general,"log_N_host_species_ratio")
    return {
        "all_pairs":len(pairs),"all_islands":len({x["island_id"] for x in pairs}),
        "pairs_with_overlapping_5km_buffers":sum(not x["five_km_radius_circular_buffers_disjoint"] for x in pairs),
        "pairs_with_both_sites_at_least_30_general_plant_records":sum(x["general_plant_effort_both_30plus"] for x in pairs),
        "original_reported_threshold_only":summary(eligible),
        "after_spatial_nonoverlap":summary(nonoverlap),
        "after_nonoverlap_and_5fold_effort_comparability":summary(effort),
        "background_species_effort_screen": {
            "eligible_pairs":len(general),
            "general_plant_species_ratio_island_mean_log":statistics.mean(general_species.values()) if general_species else None,
            "candidate_N_trait_host_species_ratio_island_mean_log":statistics.mean(host_species.values()) if host_species else None,
            "not_adjusted_for_observer_or_repeated_species":True,
        },
        "paired_site_audit":pairs,
    }


def run(geometry_path,botany_path):
    geometry=json.loads(geometry_path.read_text(encoding="utf8"))
    botany=json.loads(botany_path.read_text(encoding="utf8"))
    if geometry["summary"]["paired_sites"]!=79 or geometry["summary"]["paired_islands"]!=36:
        raise ValueError("Wrong source GSHHG/WPI cohort")
    if botany["sample_design"]["selected_pairs"]!=79 or botany["sample_design"]["botanical_search_radius_km"]!=5:
        raise ValueError("Wrong botanical 79-site / 5km census")
    if len(botany["observations"])!=316 or botany["match_audit"]["gbif_strict_exact_accepted_species_resolved"]<170:
        raise ValueError("Expected source complete record and strict taxonomy gate")
    sites={(str(x["island_id"]),int(x["port_id"])):x for x in geometry["rows"]}
    if len(sites)!=79:raise ValueError("Duplicate or missing WPI IDs")
    result={
        "schema":"chocho_WPI79_5km_botanical_offline_geometry_effort_sensitivity_v1",
        "status":"OBSERVATIONAL_SENSITIVITY_ONLY_NO_CAUSAL_NUTRIENT_GATE",
        "source_GBIF_checked":botany["checked_utc"],
        "source_geometry_gshhg_sha256":geometry["source_gshhg_sha256"],
        "source_site_geometry_sha256_claim":botany["site_geometry_sha256"],
        "source_sample_79_WPI_ports_36_islands":True,
        "periods":{p:summarize_period(sites,botany["observations"],p) for p in
                   ("early_2001_2010","later_2011_2024")},
        "limitations":[
            "The difference is measured in detected trait-bearing candidate hosts, not host abundance, leaf N or local sewage loading.",
            "Nonport sites were newly reselected, not the original 2021 WorldCover sites.",
            "1km pilot was under-recorded; this is a separate 5km spatial scale.",
            "GBIF recording effort remains different despite minimum count and ratio screens.",
            "Only Germany-derived Ellenberg N-scored species are represented; no plant-trait match imputation.",
            "2019 WPI ports are a retrospective proxy, not verified historic shipping origins.",
            "The same island contains many ports in some cases; inference is clustered at island level.",
            "Local plant alien status, period-specific colonization and butterfly utilization have not been established.",
            "Overlapping buffers, if any, are excluded in strict spatial sensitivity.",
        ]
    }
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geometry",required=True,type=Path)
    p.add_argument("--botany",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    result=run(args.geometry,args.botany)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("COAST_BOTANY_OFFLINE_BIAS_AUDIT",json.dumps({
        k:{name:value for name,value in v.items() if name!="paired_site_audit"}
        for k,v in result["periods"].items()},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
