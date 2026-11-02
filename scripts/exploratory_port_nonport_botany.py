#!/usr/bin/env python3
"""Source-locked GBIF plant-facet pilot at independently reselected island coasts.

Descriptive observations only. Non-records are not biological absences. No
N deposition, sewage or butterfly feeding or establishment is measured.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from exploratory_port_nutrient_ellenberg import load_data

ROOT = SCRIPT_DIR.parent
HOSTS = ROOT / "provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_species.csv"
GBIF_MATCH = "https://api.gbif.org/v1/species/match"
GBIF_SEARCH = "https://api.gbif.org/v1/occurrence/search"
HEADERS = {"User-Agent": "chocho-20261102-coastal-botanical-source-locked-pilot/1.0",
           "Accept": "application/json"}
PERIODS = {"early_2001_2010": "2001,2010", "later_2011_2024": "2011,2024"}


def query(url: str, params: dict, *, attempts: int = 4) -> dict:
    link = url + "?" + urllib.parse.urlencode(params)
    for i in range(attempts):
        request = urllib.request.Request(link, headers=HEADERS)
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                obj = json.load(response)
            if not isinstance(obj, dict) or "error" in obj:
                raise ValueError("Unexpected upstream response: " + str(obj)[:350])
            return obj
        except urllib.error.HTTPError as ex:
            if ex.code not in (408, 425, 429, 500, 502, 503, 504) or i == attempts - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if i == attempts - 1:
                raise
        time.sleep(min(10, 1.5 * 2**i))
    raise RuntimeError("Unreachable")


def map_gbif_one(row: dict) -> dict:
    name = row["name"]
    data = query(GBIF_MATCH, {"name": name, "rank": "SPECIES", "strict": "true"})
    canonical = str(data.get("canonicalName") or "").casefold().strip()
    exact = data.get("matchType") == "EXACT"
    accepted = data.get("status") == "ACCEPTED"
    species = data.get("rank") == "SPECIES"
    ok = exact and accepted and species and canonical == name.casefold().strip()
    key = data.get("usageKey") if ok else None
    return {"name": name, "N": row["N"], "gbif_speciesKey": int(key) if key else None,
            "matchType": data.get("matchType"), "status": data.get("status"),
            "rank": data.get("rank"), "canonicalName": data.get("canonicalName"),
            "status_label": "RESOLVED_EXACT_ACCEPTED" if key else "UNRESOLVED"}


def map_gbif(hosts: list[dict], *, workers: int = 6) -> tuple[dict[int,dict], list[dict]]:
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(map_gbif_one,row):row["name"] for row in hosts}
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda x:x["name"])
    by_key = defaultdict(list)
    for row in rows:
        if row["gbif_speciesKey"]:
            by_key[row["gbif_speciesKey"]].append(row)
    collision = [k for k,v in by_key.items() if len(v)>1]
    if collision:
        raise ValueError(f"Multiple matched species names collapsed onto GBIF key: {collision[:10]}")
    return {key:vals[0] for key,vals in by_key.items()}, rows


def choose_pilot(sites: list[dict], n: int) -> list[dict]:
    paired = [x for x in sites if x["status"] == "PAIRED"]
    if n < 0:
        return sorted(paired, key=lambda x:(int(x["island_id"]),int(x["port_id"])))
    by_island = defaultdict(list)
    for x in paired:
        by_island[x["island_id"]].append(x)
    candidate = []
    for island, records in by_island.items():
        candidate.append(max(records, key=lambda x:(x["prior_2001_2010_5km_indexed_botanical_records"], -x["port_id"])))
    candidate.sort(key=lambda x:(x["prior_2001_2010_5km_indexed_botanical_records"],int(x["island_id"])))
    if n<=0:
        return candidate
    n=min(n,len(candidate))
    # Equal sampling across 4 predetermined quantile strata of old recording effort.
    groups = [candidate[i*len(candidate)//4:(i+1)*len(candidate)//4] for i in range(4)]
    quota = [n//4 + (1 if i<n%4 else 0) for i in range(4)]
    selected = []
    for group, q in zip(groups,quota):
        ranked = sorted(group,key=lambda r:((int(r["island_id"])*2654435761)%2**32,r["port_id"]))
        selected.extend(ranked[:q])
    return sorted(selected,key=lambda x:(int(x["island_id"]),x["port_id"]))


def facet_species(lat: float, lon: float, period: str, radius_km: int = 1) -> dict:
    params = {
        "taxonKey":6,
        "year":PERIODS[period],
        "geoDistance":f"{lat:.7f},{lon:.7f},{radius_km}km",
        "hasCoordinate":"true",
        "hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT",
        "facet":"speciesKey",
        "facetLimit":1000,
        "limit":0
    }
    offset = 0
    keys = {}
    total = None
    while True:
        params["facetOffset"] = offset
        data = query(GBIF_SEARCH, params)
        n = int(data["count"])
        if total is not None and n != total:
            raise ValueError("GBIF index changed during facet pagination")
        total = n
        facets = next((x.get("counts",[]) for x in data.get("facets",[])
                       if x.get("field","").upper() == "SPECIES_KEY"),[])
        for rec in facets:
            key = int(rec["name"])
            if key in keys:
                raise ValueError("Species facet overlap across pages")
            keys[key] = int(rec["count"])
        if len(facets)<1000:
            break
        offset += 1000
        if offset>40000:
            raise ValueError("Unreasonably high botanical facet count")
    if sum(keys.values())>total:
        raise ValueError("Species facet exceeds overall indexed plants")
    return {"indexed_plant_records":total,"indexed_species_keys":len(keys),
            "records_with_speciesKey":sum(keys.values()),"faceted_specieskey_counts":keys,
            "facet_pages":offset//1000+1}


def plant_summary(facet: dict, key_to_trait: dict[int,dict]) -> dict:
    values = [(key_to_trait[k]["N"],v) for k,v in facet["faceted_specieskey_counts"].items() if k in key_to_trait]
    n = len(values)
    return {
        "all_plant_records":facet["indexed_plant_records"],
        "all_plant_species_keys":facet["indexed_species_keys"],
        "all_plant_species_key_facet_coverage":(
            facet["records_with_speciesKey"]/facet["indexed_plant_records"] if facet["indexed_plant_records"] else None),
        "candidate_hosts_with_N_indexed_species":n,
        "candidate_hosts_N_high_7_9":sum(val>=7 for val,_ in values),
        "candidate_hosts_N_low_1_3":sum(val<=3 for val,_ in values),
        "candidate_hosts_mean_N_unweighted_by_occurrences":(
            sum(val for val,_ in values)/n if n else None),
        "candidate_hosts_N_high_fraction":(
            sum(val>=7 for val,_ in values)/n if n else None),
        "candidate_hosts_species_key_matches":sorted(k for k in facet["faceted_specieskey_counts"] if k in key_to_trait),
        "not_observed_does_not_mean_absent": True,
    }


def summarise(receipt: dict, *, min_plant_records: int, min_trait_species: int) -> dict:
    groups = defaultdict(dict)
    for row in receipt["observations"]:
        groups[(row["island_id"],row["port_id"],row["period"])][row["site_type"]] = row
    paired = []
    for (island, port, period), items in groups.items():
        if set(items) != {"port","nonport"}:
            continue
        a,b=items["port"],items["nonport"]
        sufficient = all(
            x["all_plant_records"] >= min_plant_records and
            x["candidate_hosts_with_N_indexed_species"] >= min_trait_species
            for x in (a,b))
        diff = (a["candidate_hosts_mean_N_unweighted_by_occurrences"]-
                b["candidate_hosts_mean_N_unweighted_by_occurrences"]) if sufficient else None
        paired.append({
            "island_id":island, "port_id":port, "period":period,
            "port_all_plant_records":a["all_plant_records"],
            "nonport_all_plant_records":b["all_plant_records"],
            "port_N_trait_candidates":a["candidate_hosts_with_N_indexed_species"],
            "nonport_N_trait_candidates":b["candidate_hosts_with_N_indexed_species"],
            "data_adequate_for_trait_composition_compare":sufficient,
            "port_minus_nonport_candidate_host_mean_N":diff
        })
    stats={}
    for period in PERIODS:
        q=[x for x in paired if x["period"]==period]
        valid=[x["port_minus_nonport_candidate_host_mean_N"] for x in q if x["data_adequate_for_trait_composition_compare"]]
        by_island=defaultdict(list)
        for x in q:
            if x["data_adequate_for_trait_composition_compare"]:
                by_island[x["island_id"]].append(x["port_minus_nonport_candidate_host_mean_N"])
        island_values={island:sum(values)/len(values) for island,values in by_island.items()}
        site_coverage=lambda key: sum(x[key]>0 for x in q)
        stats[period]={
            "pairs_observed":len(q), "pairs_trait_comparable":len(valid),
            "independent_islands_trait_comparable":len(island_values),
            "paired_mean_N_differences":valid,
            "mean_of_paired_differences_unweighted_by_port":sum(valid)/len(valid) if valid else None,
            "island_weighted_mean_N_difference":(
                sum(island_values.values())/len(island_values) if island_values else None),
            "positive_paired_differences":sum(x>0 for x in valid),
            "islands_positive_mean_N_differences":sum(x>0 for x in island_values.values()),
            "pairs_both_sides_any_general_botanical_record":sum(
                x["port_all_plant_records"]>0 and x["nonport_all_plant_records"]>0 for x in q),
            "pairs_both_sides_30plus_general_plant_records":sum(
                x["port_all_plant_records"]>=30 and x["nonport_all_plant_records"]>=30 for x in q),
            "not_a_significance_test":True,
            "data_gate_passed":len(island_values)>=6,
        }
    return {"periods":stats,
            "pair_rows":sorted(paired,key=lambda x:(x["period"],int(x["island_id"])))}


def run(sites_path: Path, traits_path: Path, *, n_pilot: int, radius_km: int,
        min_plant_records: int, min_trait_species: int) -> dict:
    full = json.loads(sites_path.read_text(encoding="utf8"))
    if (len(full["rows"]),full["summary"]["focal_islands"]) != (79,36):
        raise ValueError("Site sample must be from WPI79 locked geometry gate")
    if radius_km not in (1,5):
        raise ValueError("Only predeclared 1km primary pilot or 5km equal-footprint comparison allowed")
    selected=choose_pilot(full["rows"],n_pilot)
    hosts, counts=load_data(HOSTS,traits_path)
    resolved, mapped=map_gbif(hosts)
    def acquire_item(item):
        pair,site_type,period=item
        lat=pair["port_sample_lat"] if site_type=="port" else pair["nonport_lat"]
        lon=pair["port_sample_lon"] if site_type=="port" else pair["nonport_lon"]
        if lat is None or lon is None:
            raise ValueError("Unpaired geometry in selected sites")
        facet=facet_species(float(lat),float(lon),period,radius_km=radius_km)
        summary=plant_summary(facet,resolved)
        return {"island_id":pair["island_id"],"port_id":pair["port_id"],
                "site_type":site_type,"period":period,"latitude":lat,
                "longitude":lon,**summary}
    items=[(pair,site_type,period) for pair in selected
           for site_type in ("port","nonport") for period in PERIODS]
    observations=[]
    # Bounded I/O concurrency; deterministic sorted data after all requests.
    with ThreadPoolExecutor(max_workers=4) as pool:
        future_map={pool.submit(acquire_item,item):item for item in items}
        for future in as_completed(future_map):
            row=future.result()
            observations.append(row)
            print("BOTANICAL_SITE",row["island_id"],row["port_id"],
                  row["site_type"],row["period"],
                  "plants",row["all_plant_records"],
                  "species",row["all_plant_species_keys"],
                  "N_hosts",row["candidate_hosts_with_N_indexed_species"],flush=True)
    observations.sort(key=lambda x:(int(x["island_id"]),int(x["port_id"]),x["site_type"],x["period"]))
    result={
        "schema":"chocho_port_nonport_botanical_pilot_source_20261102_v1",
        "status":"EXPLORATORY_GBIF_INDEXED_PLANT_COMPOSITION_NOT_NITROGEN_EFFECT",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "site_geometry_sha256":hashlib.sha256(sites_path.read_bytes()).hexdigest(),
        "site_selection_method":"blind to focal botanical and butterfly responses; one port per island across four prior-effort quartiles",
        "sample_design":{"selected_islands":len(selected),"selected_pairs":len(selected),
                         "intended_observations":len(selected)*2*len(PERIODS),
                         "period_windows":PERIODS,"botanical_search_radius_km":radius_km,
                         "minimum_general_plant_records":min_plant_records,
                         "minimum_N_scored_host_species":min_trait_species},
        "match_audit":{"source_host_counts":counts,
                       "gbif_strict_exact_accepted_species_resolved":len(resolved),
                       "gbif_unresolved_count":len(mapped)-len(resolved),
                       "gbif_unresolved_names":[x["name"] for x in mapped if x["gbif_speciesKey"] is None]},
        "observations":observations,
        "caveats":[
            "GBIF records are detection/recording indices, not independent occurrence or plant abundance.",
            "Non-records are not biological absence; sparse sites are marked incomparable.",
            "Matched N trait is from Germany-centric Ellenberg ecology; ~3/4 of global candidate hosts lack exact scores.",
            "Contemporary catalogued alien status is regional and does not identify historical introductions on individual islands.",
            "WPI 2019 port positions are retrospective proxies; no vessel histories or sewage loads measured.",
            "Coordinates are independently reconstructed; original nonport shoreline controls were not recovered.",
            "Not yet stratified by urbanized versus natural nonport or climate/coastline geometry.",
        ]
    }
    result["paired_summary"]=summarise(result,min_plant_records=min_plant_records,
                                      min_trait_species=min_trait_species)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites",type=Path,required=True)
    p.add_argument("--ellenberg",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--pilot-pairs",type=int,default=12,
                   help="-1 means all independently paired WPI79 sites; 0 means one per island")
    p.add_argument("--radius-km",type=int,choices=(1,5),default=1)
    p.add_argument("--min-general-plant-records",type=int,default=30)
    p.add_argument("--min-N-trait-species",type=int,default=2)
    args=p.parse_args()
    out=run(args.sites,args.ellenberg,n_pilot=args.pilot_pairs,
            radius_km=args.radius_km,
            min_plant_records=args.min_general_plant_records,
            min_trait_species=args.min_N_trait_species)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("BOTANICAL_PAIRED_GATE",json.dumps({
        "design":out["sample_design"],
        "GBIF_taxa_resolved":out["match_audit"]["gbif_strict_exact_accepted_species_resolved"],
        "periods":out["paired_summary"]["periods"]},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
