#!/usr/bin/env python3
"""Exploratory global-prevalence test of candidate non-native butterfly hosts at ports.

Tests whether site nutrient-affinity (German Ellenberg N) relates to the number
of WPI port buffer catalogues *conditional on worldwide GBIF reporting*.
Global GBIF record count/country coverage are imperfect geography and
observation-intensity proxies, not objective range sizes. No sewage effect.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from exploratory_port_nutrient_ellenberg import DEFAULT_HOSTS, load_data, ranks, spearman
from exploratory_port_nonport_botany import GBIF_SEARCH, map_gbif, query

PERIOD="2001,2024"

def get_worldwide_one(gbif_species_key: int) -> dict:
    params={"taxonKey":str(gbif_species_key),"year":PERIOD,
            "hasCoordinate":"true","hasGeospatialIssue":"false",
            "occurrenceStatus":"PRESENT","facet":"country",
            "facetLimit":300,"limit":0}
    try:
        data=query(GBIF_SEARCH,params)
    except Exception as ex:
        raise RuntimeError(
            f"Global GBIF facet failed: taxonKey={gbif_species_key}; "
            f"query_fields={params}; original_error={type(ex).__name__}: {ex}"
        ) from ex
    countries=next((x.get("counts",[]) for x in data.get("facets",[])
                    if x.get("field","").upper()=="COUNTRY"),[])
    if len(countries)>=300:
        raise ValueError("Country facets possibly truncated (300 entries)")
    n_records=int(data["count"])
    counted=sum(int(x["count"]) for x in countries)
    if counted>n_records:
        raise ValueError("Country facets exceed global indexed records")
    return {"gbif_speciesKey":gbif_species_key,"global_indexed_records_2001_2024":n_records,
            "number_countries_with_indexed_observations":len(countries),
            "count_with_valid_country_in_facets":counted}


def get_worldwide(keys:list[int],workers:int=6)->dict[int,dict]:
    results={}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(get_worldwide_one,k):k for k in keys}
        for future in as_completed(futures):
            k=futures[future]
            row=future.result()
            if k in results:raise ValueError("Duplicate taxon result")
            results[k]=row
    return results


def centered(values):
    avg=statistics.mean(values)
    return [v-avg for v in values]


def pearson_centered(x,y):
    if len(x)!=len(y) or len(x)<4:return math.nan
    a,b=centered(x),centered(y)
    dot=lambda z,w:sum(u*v for u,v in zip(z,w))
    sx,sy=dot(a,a),dot(b,b)
    return dot(a,b)/math.sqrt(sx*sy) if sx>0 and sy>0 else math.nan


def residualize(y,controls):
    y=centered(y)
    if not controls:return y
    if len(controls)>2:raise ValueError("At most two predeclared global report controls")
    x=[centered(c) for c in controls]
    dot=lambda a,b:sum(u*v for u,v in zip(a,b))
    if len(x)==1:
        d=dot(x[0],x[0])
        beta=(dot(x[0],y)/d) if d>1e-12 else 0
        return [v-beta*w for v,w in zip(y,x[0])]
    a=dot(x[0],x[0]);b=dot(x[0],x[1]);d=dot(x[1],x[1])
    e=dot(x[0],y);f=dot(x[1],y)
    det=a*d-b*b
    if det<=max(1e-12,1e-10*a*d):
        return residualize(y,[x[0]])  # predeclared conservative rank fallback
    c0=(e*d-f*b)/det
    c1=(f*a-e*b)/det
    return [u-c0*v-c1*w for u,v,w in zip(y,x[0],x[1])]


def partial_spearman(rows,controls):
    if len(rows)<8:return math.nan
    y=ranks([r["ports"] for r in rows])
    x=ranks([r["N"] for r in rows])
    c=[ranks([r[name] for r in rows]) for name in controls]
    ry=residualize(y,c)
    rx=residualize(x,c)
    return pearson_centered(rx,ry)


def bootstrap_by_genus(rows, *, n_draws=2999, seed=12121):
    by_genus=defaultdict(list)
    for row in rows:by_genus[row["genus"]].append(row)
    genera=sorted(by_genus)
    rng=random.Random(seed)
    keys=["raw","partial_countries","partial_counts","partial_both"]
    draws={k:[] for k in keys}
    for _ in range(n_draws):
        chosen=[item for _ in genera for item in by_genus[rng.choice(genera)]]
        vals=associations(chosen)
        for k,v in vals.items():
            if math.isfinite(v):draws[k].append(v)
    out={}
    for k in keys:
        arr=sorted(draws[k])
        out[k]={"draws":len(arr),"descriptive_genus_bootstrap_95":[
            arr[round((len(arr)-1)*.025)],arr[round((len(arr)-1)*.975)]] if arr else None}
    return out


def associations(rows):
    return {
        "raw":spearman([r["N"] for r in rows],[r["ports"] for r in rows]),
        "partial_countries":partial_spearman(rows,["number_countries_with_indexed_observations"]),
        "partial_counts":partial_spearman(rows,["global_indexed_records_2001_2024"]),
        "partial_both":partial_spearman(rows,["number_countries_with_indexed_observations","global_indexed_records_2001_2024"]),
    }


def run(ellenberg_path:Path,n_bootstrap:int)->dict:
    hosts,meta=load_data(DEFAULT_HOSTS,ellenberg_path)
    key_traits,matching=map_gbif(hosts,workers=6)
    global_info=get_worldwide(sorted(key_traits))
    joined=[]
    for key,trait in key_traits.items():
        if key not in global_info:raise ValueError("Missing global GBIF response")
        worldwide=global_info[key]
        if worldwide["global_indexed_records_2001_2024"]<1:continue
        joined.append({**trait,**worldwide,
                       "N":trait["N"],"genus":trait["name"].split(" ")[0],
                       "ports":next(x["ports"] for x in hosts if x["name"]==trait["name"])})
    joined.sort(key=lambda x:x["name"])
    assoc=associations(joined)
    boot=bootstrap_by_genus(joined,n_draws=n_bootstrap)
    return {
        "schema":"chocho_port_nutrient_global_GBIF_country_prevalence_control_v1",
        "status":"EXPLORATORY_CORRELATION_AFTER_GBIF_WORLDWIDE_PREVALENCE_CONTROL",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "global_GBIF_year_window":PERIOD,
        "GBIF_country_counts_are_indexed_records_not_biological_distribution":True,
        "species_matches_from_pinned_German_ellenberg":{"eligible_N_host_species":meta["n_exact_species_matches_with_N"],
                "exact_GBIF_matched_species":len(key_traits),"global_observation_supported":len(joined),
                "excluded_unmatched": [x["name"] for x in matching if x["gbif_speciesKey"] is None]},
        "number_distinct_genera":len({r["genus"] for r in joined}),
        "correlations":assoc,
        "genus_cluster_bootstrap":boot,
        "global_indexed_count_median":statistics.median(x["global_indexed_records_2001_2024"] for x in joined),
        "countries_observed_median":statistics.median(x["number_countries_with_indexed_observations"] for x in joined),
        "species_record_inputs":joined,
        "limitations":[
            "Global GBIF country records and count are observation/reporting proxies, not independent native geographic range data.",
            "Records near indexed WPI ports are a subset of global GBIF count, so part-whole overlap can influence partial coefficients.",
            "Only 182/687 global port candidate introduced larval-host plant species with German nutrient indicator and accepted taxon match.",
            "Country facet does not account for unequal region area, detectability, planting, introduced status or observation intensity.",
            "Year periods pooled and retrospective 2019 port locations; does not identify historical shipping pathways.",
            "No site nitrogen loads, wastewater, butterfly larval plant utilization, or colonization observed.",
        ],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ellenberg",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=2999)
    a=p.parse_args()
    z=run(a.ellenberg,a.bootstrap)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("PORT_N_GBIF_GLOBAL_COVERAGE_ASSOCIATIONS",json.dumps({
        "n":z["species_matches_from_pinned_German_ellenberg"],
        "genera":z["number_distinct_genera"],
        "corr":z["correlations"],
        "genus_bootstrap":z["genus_cluster_bootstrap"],
        "global_indexed_count_median":z["global_indexed_count_median"],
        "countries_observed_median":z["countries_observed_median"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
