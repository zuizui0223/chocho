#!/usr/bin/env python3
"""Explore whether nutrient affinity among *butterfly hosts* is port-specific
or reflects the whole recorded plant community at matched island coasts.

Source-locked: independent disjoint WPI78 site selection, GBIF occurrence
speciesKey facets for 2001–2010 and 2011–2024, historic German Ellenberg N.
Non-target plants can themselves be butterfly hosts; they are an ecological
background, NOT verified non-host plants. No N fluxes or host abundance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import statistics
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from exploratory_port_nonport_botany import query, GBIF_MATCH, facet_species
from exploratory_port_nutrient_ellenberg import parse_value

PERIODS = ("early_2001_2010", "later_2011_2024")
MAX_WORKERS = 5


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_trait_names(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    result = {}
    for r in rows:
        name = (r.get("Name") or "").strip()
        score = parse_value(r.get("N") or "")
        if score is None or not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][a-z-]+", name):
            continue
        key = name.casefold()
        if key in result and result[key]["N"] != score:
            raise ValueError(f"Conflicting Ellenberg N: {name}")
        result[key] = {"name": name, "N": score}
    return [result[k] for k in sorted(result)]


def resolve_one(row: dict) -> dict:
    name = row["name"]
    d = query(GBIF_MATCH, {"name": name, "rank": "SPECIES", "strict": "true"})
    canonical = str(d.get("canonicalName") or "").strip().casefold()
    eligible = (d.get("rank") == "SPECIES" and d.get("matchType") == "EXACT"
                and d.get("status") == "ACCEPTED" and canonical == name.casefold()
                and d.get("usageKey"))
    return {"name": name, "N": row["N"],
            "gbif_speciesKey": int(d["usageKey"]) if eligible else None,
            "gbif_status": str(d.get("status") or ""),
            "gbif_matchType": str(d.get("matchType") or ""),
            "gbif_canonical": d.get("canonicalName")}


def create_mapping(ellenberg: Path, focal: Path, output: Path, workers: int) -> dict:
    traitrows = read_trait_names(ellenberg)
    prev = read_json(focal)
    focalrows = prev["species_record_inputs"]
    if len(focalrows) != 182:
        raise ValueError("Frozen independent global host source changed")
    byname = {x["name"].casefold(): x for x in focalrows}
    # The 182 previously GBIF-validated candidate host identities are pinned.
    reuse, need = [], []
    for r in traitrows:
        known = byname.get(r["name"].casefold())
        if known and known["N"] == r["N"]:
            reuse.append({"name": r["name"], "N": r["N"],
                          "gbif_speciesKey": int(known["gbif_speciesKey"]),
                          "source": "GLOBAL_GBIF_182_LOCKED"})
        else:
            need.append(r)
    rows = list(reuse)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        fs = {pool.submit(resolve_one, r): r["name"] for r in need}
        for fut in as_completed(fs):
            val = fut.result()
            val["source"] = "LIVE_GBIF_ACCEPTED_EXACT"
            rows.append(val)
    rows.sort(key=lambda x:x["name"].casefold())
    # Multiple synonyms with identical speciesKey are not safe for unweighted N means.
    bykey = defaultdict(list)
    for r in rows:
        if r["gbif_speciesKey"] is not None:
            bykey[r["gbif_speciesKey"]].append(r)
    collisions = {str(k): v for k,v in bykey.items() if len(v)>1}
    if collisions:
        keys = set(map(int, collisions))
        for row in rows:
            if row["gbif_speciesKey"] in keys:
                row["ambiguous_key_collision"] = True
    usable = {str(k):v[0] for k,v in bykey.items() if len(v)==1}
    hostmap = {str(x["gbif_speciesKey"]):x for x in focalrows}
    missing_focal = [k for k in hostmap if k not in usable]
    if missing_focal:
        raise ValueError(f"Missing locked candidate N host key in all plant trait mapping: {missing_focal[:12]}")
    counts={"source_ellenberg_valid_binomial_N":len(traitrows),
            "reused_locked_host_species":len(reuse),
            "fresh_name_matches_attempted":len(need),
            "exact_accepted_unambiguous_species":len(usable),
            "names_without_exact_accepted_key":sum(r["gbif_speciesKey"] is None for r in rows),
            "duplicate_speciesKey_collisions":len(collisions)}
    result={
        "schema":"chocho_full_ellenberg_N_exact_GBIF_key_mapping_v1",
        "status":"GBIF_EXACT_ACCEPTED_SPECIESKEY_ONLY",
        "focal_host_input_run":"https://github.com/zuizui0223/chocho/actions/runs/41403516863",
        "ellenberg_sha256":hashlib.sha256(ellenberg.read_bytes()).hexdigest(),
        "counts":counts,"accepted_unambiguous_by_speciesKey":usable,
        "focal_host_by_speciesKey":hostmap,
        "unresolved_names":[r["name"] for r in rows if r["gbif_speciesKey"] is None],
        "duplicate_key_collision_detail":collisions,
        "no_trait_imputation":True
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("ALL_PLANTS_N_MAPPING_GATE",json.dumps(counts,sort_keys=True),flush=True)
    return result


def selected_sites(source: dict) -> list[dict]:
    if (len(source["rows"]),source["summary"]["focal_ports"],
        source["summary"]["paired_sites"],source["summary"]["paired_islands"]) != (79,79,78,35):
        raise ValueError("Expected immutable corrected WPI79 / eligible 78 paired sample")
    rows = [r for r in source["rows"] if r["status"]=="PAIRED"]
    if any(r["five_km_radius_buffers_disjoint"] is not True for r in rows):
        raise ValueError("Five-km port/non-port sample footprints overlap")
    return sorted(rows,key=lambda r:(int(r["island_id"]),r["port_id"]))


def summarize_facet(facet: dict, mapping: dict) -> dict:
    plant = {int(k):int(v) for k,v in facet["faceted_specieskey_counts"].items()}
    scored = {int(k):v for k,v in mapping["accepted_unambiguous_by_speciesKey"].items()}
    focal = {int(k):v for k,v in mapping["focal_host_by_speciesKey"].items()}
    matched = {k:row for k,row in scored.items() if k in plant}
    target = [row["N"] for k,row in matched.items() if k in focal]
    background = [row["N"] for k,row in matched.items() if k not in focal]
    scores = [row["N"] for row in matched.values()]
    score_mean=lambda x:statistics.mean(x) if x else None
    counts=facet["faceted_specieskey_counts"]
    return {
       "plant_records":facet["indexed_plant_records"],
       "all_plant_species_keys":len(counts),
       "all_plant_species_faceted_record_coverage":(
           facet["records_with_speciesKey"]/facet["indexed_plant_records"]
           if facet["indexed_plant_records"] else None),
       "N_scored_plant_species":len(scores),
       "N_scored_fraction_of_all_indexed_plant_species":len(scores)/len(counts) if counts else None,
       "N_scored_candidate_host_species":len(target),
       "N_scored_other_plant_species":len(background),
       "mean_N_all_scored_plants":score_mean(scores),
       "mean_N_candidate_hosts":score_mean(target),
       "mean_N_other_scored_plants":score_mean(background),
       "mean_N_hosts_minus_other":(
           statistics.mean(target)-statistics.mean(background)
           if target and background else None),
       "candidate_hosts_high_N_7_9":sum(x>=7 for x in target),
       "other_plants_high_N_7_9":sum(x>=7 for x in background),
       "no_butterfly_larval_usage_or_plant_abundance_measured":True,
    }


def records_from_cache(cache: Path) -> list[dict]:
    if not cache.exists():return []
    data=[]
    for line in cache.read_text(encoding="utf-8").splitlines():
        if line.strip():
            data.append(json.loads(line))
    keys=[(r["island_id"],r["port_id"],r["side"],r["period"]) for r in data]
    if len(keys)!=len(set(keys)):
        raise ValueError("Duplicate site facets in raw cache")
    return data


def collect_observations(geometry: Path, traits_map: Path, source: Path, output: Path,
                         cache: Path, workers: int, radius_km: int) -> dict:
    geo=read_json(geometry)
    mapping=read_json(traits_map)
    prior=read_json(source)
    if prior["sample_design"]["selected_pairs"]!=78 or prior["sample_design"]["botanical_search_radius_km"]!=5:
        raise ValueError("Previous 78-site observation source ID invalid")
    sites=selected_sites(geo)
    if hashlib.sha256(geometry.read_bytes()).hexdigest()!=prior["site_geometry_sha256"]:
        raise ValueError("Not the same exact geometry bytes as original host-only GBIF comparison")
    ids=[(r["island_id"],r["port_id"],s,p) for r in sites
         for s in ("port","nonport") for p in PERIODS]
    prior_map={(r["island_id"],r["port_id"],r["site_type"],r["period"]):r for r in prior["observations"]}
    if len(ids)!=312 or len(prior_map)!=312 or set(ids)!=set(prior_map):
        raise ValueError("Site IDs differ from strict previous 312-observation source")
    cache.parent.mkdir(parents=True,exist_ok=True)
    existing=records_from_cache(cache)
    done={(r["island_id"],r["port_id"],r["side"],r["period"]):r for r in existing}
    if not set(done).issubset(ids):raise ValueError("Cache has foreign sites")
    site_lookup={(r["island_id"],r["port_id"]):r for r in sites}
    remaining=[x for x in ids if x not in done]

    def acquire(item):
        island,port,side,period=item
        p=site_lookup[(island,port)]
        lat,lon=(p["port_sample_lat"],p["port_sample_lon"]) if side=="port" else (
            p["nonport_lat"],p["nonport_lon"])
        facet=facet_species(float(lat),float(lon),period,radius_km=radius_km)
        # Store original all-plant species facets in the private workflow artifact.
        return {"island_id":island,"port_id":port,"side":side,"period":period,
                "lat":lat,"lon":lon,"facet":facet}

    with cache.open("a",encoding="utf-8") as f, ThreadPoolExecutor(max_workers=workers) as pool:
        fs={pool.submit(acquire,i):i for i in remaining}
        for fut in as_completed(fs):
            val=fut.result()
            f.write(json.dumps(val,sort_keys=True,separators=(",",":"))+"\n")
            f.flush()
            key=(val["island_id"],val["port_id"],val["side"],val["period"])
            done[key]=val
            print("ALL_PLANT_SPECIESKEY_SITE",*key,
                  "records",val["facet"]["indexed_plant_records"],
                  "species",val["facet"]["indexed_species_keys"],flush=True)
    if len(done)!=len(ids):
        raise ValueError("Incomplete full 312 site×period facet retrieval")
    obs=[]
    for key in ids:
        row=done[key]
        summary=summarize_facet(row["facet"],mapping)
        historic=prior_map[key]
        obs.append({**{k:row[k] for k in ("island_id","port_id","side","period")},
                    **summary,
                    "historical_host_only_source_records":historic["all_plant_records"],
                    "historical_host_only_source_species":historic["all_plant_species_keys"],
                    "GBIF_index_record_drift":summary["plant_records"]-historic["all_plant_records"],
                    "GBIF_index_species_drift":summary["all_plant_species_keys"]-historic["all_plant_species_keys"],
                    "historical_host_only_N_species":historic["candidate_hosts_with_N_indexed_species"],
                    "candidate_host_N_drift":summary["N_scored_candidate_host_species"]-historic["candidate_hosts_with_N_indexed_species"]})
    result={
       "schema":"chocho_WPI78_5km_all_flora_N_background_v1",
       "status":"OBSERVATIONAL_GBIF_TRAIT_COMPOSITION_NOT_NITROGEN_OR_ESTABLISHMENT",
       "checked_at_utc":datetime.now(timezone.utc).isoformat(),
       "original_geometry_sha256":prior["site_geometry_sha256"],
       "original_host_only_index_checked_utc":prior["checked_utc"],
       "all_plant_mapping_sha256":hashlib.sha256(traits_map.read_bytes()).hexdigest(),
       "raw_facets_jsonl_sha256":hashlib.sha256(cache.read_bytes()).hexdigest(),
       "raw_facets_not_in_Git":True,
       "predeclared_source_pairs":78,
       "independent_islands":35,
       "site_period_facets":len(obs),
       "mean_N_focal_hosts_is_not_local_abundance":True,
       "candidate_host_set_is_globally_not_native_in_catalogue_not_island_introduction_status":True,
       "other_scored_plants_are_NOT_verified_nonhosts":True,
       "drift_counts":{"records_changed_sites":sum(x["GBIF_index_record_drift"]!=0 for x in obs),
                       "species_changed_sites":sum(x["GBIF_index_species_drift"]!=0 for x in obs),
                       "host_N_changed_sites":sum(x["candidate_host_N_drift"]!=0 for x in obs)},
       "observations":obs
    }
    result["paired_analysis"]=analyze(obs)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("ALL_FLORA_N_PAIRED_GATES",json.dumps({
        "drift":result["drift_counts"],"output":result["paired_analysis"]
    },sort_keys=True),flush=True)
    return result


def island_bootstrap(vals: list[float], draws=4999, seed=18711) -> list[float]|None:
    if not vals:return None
    rng=random.Random(seed)
    estimates=sorted(statistics.mean(rng.choices(vals,k=len(vals))) for _ in range(draws))
    return [estimates[round((len(estimates)-1)*.025)],
            estimates[round((len(estimates)-1)*.975)]]


def island_mean(records, field):
    groups=defaultdict(list)
    for r in records:
        if r[field] is not None:
            groups[str(r["island_id"])].append(r[field])
    return {k:statistics.mean(v) for k,v in groups.items()}


def one_metric(pairs,key):
    groups=island_mean(pairs,key)
    v=list(groups.values())
    return {"pairs":sum(r[key] is not None for r in pairs),
            "islands":len(groups),
            "island_equal_weight_mean":statistics.mean(v) if v else None,
            "island_equal_weight_median":statistics.median(v) if v else None,
            "island_bootstrap_95_descriptive":island_bootstrap(v),
            "islands_positive":sum(x>0 for x in v),
            "islands_negative":sum(x<0 for x in v)}


def pair_rows(obs, min_records=30, min_hosts=2, min_other=5):
    lookup=defaultdict(dict)
    for r in obs:
        lookup[(r["island_id"],r["port_id"],r["period"])][r["side"]]=r
    pairs=[]
    for (island,port,period),samples in sorted(lookup.items(),key=lambda kv:(kv[0][2],int(kv[0][0]),kv[0][1])):
        if set(samples)!={"port","nonport"}:raise ValueError("Unpaired site-period")
        p,n=samples["port"],samples["nonport"]
        eligible=all(x["plant_records"]>=min_records and
                     x["N_scored_candidate_host_species"]>=min_hosts and
                     x["N_scored_other_plant_species"]>=min_other for x in (p,n))
        pair={"island_id":island,"port_id":port,"period":period,
              "both_sites_general_records_30plus":min(p["plant_records"],n["plant_records"])>=30,
              "host_N_and_other_N_supported":eligible,
              "observed_general_plant_effort_ratio_minmax":(
                  min(p["plant_records"],n["plant_records"])/max(p["plant_records"],n["plant_records"])
                  if max(p["plant_records"],n["plant_records"]) else 0),
              "both_sites_scored_plant_species_coverage_min":min(
                  p["N_scored_fraction_of_all_indexed_plant_species"] or 0,
                  n["N_scored_fraction_of_all_indexed_plant_species"] or 0),
              "ports_minus_nonport_background_N":None,
              "ports_minus_nonport_target_host_N":None,
              "difference_in_differences_hosts_minus_other_N":None}
        if eligible:
            pair["ports_minus_nonport_background_N"]=p["mean_N_other_scored_plants"]-n["mean_N_other_scored_plants"]
            pair["ports_minus_nonport_target_host_N"]=p["mean_N_candidate_hosts"]-n["mean_N_candidate_hosts"]
            pair["difference_in_differences_hosts_minus_other_N"]=(
                pair["ports_minus_nonport_target_host_N"]-pair["ports_minus_nonport_background_N"])
        pairs.append(pair)
    return pairs


def analyze(obs):
    pairs=pair_rows(obs)
    results={"selection":{"site_periods":len(obs),"paired_port_periods":len(pairs),
                          "predeclared_min_record_per_site":30,
                          "predeclared_min_target_N_species_per_site":2,
                          "predeclared_min_other_N_species_per_site":5},
             "periods":{}}
    for period in PERIODS:
        raw=[r for r in pairs if r["period"]==period]
        eligible=[r for r in raw if r["host_N_and_other_N_supported"]]
        effort=[r for r in eligible if r["observed_general_plant_effort_ratio_minmax"]>=0.2]
        def summary(sel):
            return {"eligible_pairs":len(sel),
                    "eligible_islands":len({x["island_id"] for x in sel}),
                    "background_N":one_metric(sel,"ports_minus_nonport_background_N"),
                    "candidate_host_N":one_metric(sel,"ports_minus_nonport_target_host_N"),
                    "host_minus_background_N":one_metric(sel,"difference_in_differences_hosts_minus_other_N")}
        results["periods"][period]={
            "original_wpi_pairs":len(raw),"both_sites_30plus_records":sum(x["both_sites_general_records_30plus"] for x in raw),
            "base_quality":summary(eligible),"max_5fold_recording_effort_difference":summary(effort)}
    a={(p["island_id"],p["port_id"]):p for p in pairs if p["period"]==PERIODS[0] and p["host_N_and_other_N_supported"]}
    b={(p["island_id"],p["port_id"]):p for p in pairs if p["period"]==PERIODS[1] and p["host_N_and_other_N_supported"]}
    shared=sorted(set(a)&set(b),key=lambda x:(int(x[0]),x[1]))
    temporal=[]
    for island,port in shared:
        row={"island_id":island,"port_id":port}
        for k in ("ports_minus_nonport_background_N","ports_minus_nonport_target_host_N",
                  "difference_in_differences_hosts_minus_other_N"):
            row[k]=b[(island,port)][k]-a[(island,port)][k]
        temporal.append(row)
    results["identical_pair_period_difference"]={
        "shared_pairs":len(shared),
        "shared_islands":len({k[0] for k in shared}),
        "background_N":one_metric(temporal,"ports_minus_nonport_background_N"),
        "candidate_host_N":one_metric(temporal,"ports_minus_nonport_target_host_N"),
        "host_minus_background_N":one_metric(temporal,"difference_in_differences_hosts_minus_other_N")}
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage",choices=("map","observe"),required=True)
    p.add_argument("--ellenberg",type=Path)
    p.add_argument("--focal",type=Path)
    p.add_argument("--mapping",type=Path)
    p.add_argument("--geometry",type=Path)
    p.add_argument("--original-botany",type=Path)
    p.add_argument("--raw-cache",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=5)
    args=p.parse_args()
    if args.stage=="map":
        assert args.ellenberg and args.focal
        create_mapping(args.ellenberg,args.focal,args.out,args.workers)
    else:
        assert args.geometry and args.mapping and args.original_botany and args.raw_cache
        collect_observations(args.geometry,args.mapping,args.original_botany,args.out,args.raw_cache,args.workers,5)


if __name__=="__main__":
    main()
