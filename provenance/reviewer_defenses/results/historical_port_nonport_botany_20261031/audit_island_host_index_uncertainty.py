#!/usr/bin/env python3
"""Offline, source-locked bootstrap on 2001–2010 matched port/offport plants.

The historical 1-km botanical site comparison uses the ten top GEB larval
host-opportunity plant names, NOT a validated island-specific introduced
plant inventory. This is a retrospective descriptive observation, NOT
demonstrated maritime transport or butterfly settlement.

Resample independent island units, not individual ports; the original fixed
geometric design can share an island across multiple port pairs. Report no
bootstrap confidence interval where fewer than three eligible islands remain.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
CSV=HERE/"port_offport_historical_pairs.csv"
CSV_GIT_BLOB="ae7fa4986e1a5be73a3c50b89669d05bc5138ac1"
RECEIPT=HERE/"historical_port_offport_results.json"
RECEIPT_GIT_BLOB="53ef6dddbe1a8efceaf518812bf2b3c77cf285d6"
THRESHOLDS=(1,5,10,20,50)
SEED=20261031
DRAWS=5000


def verify(path,expected):
    raw=path.read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+b"\0"+raw).hexdigest()
    if actual!=expected:
        raise ValueError(f"Source Git blob changed: {path} {actual}")
    return raw


def percentile(a,p):
    arr=sorted(a)
    position=(len(arr)-1)*p
    low=int(math.floor(position)); high=int(math.ceil(position))
    return arr[low]+(arr[high]-arr[low])*(position-low)


def island_bootstrap(values,seed,draws):
    """LCG high bits, not parity-biased low-bit modulus for tiny island sets."""
    state=seed
    out=[]
    for _ in range(draws):
        total=0.0
        for _ in range(len(values)):
            state=(state*1664525+1013904223)&0xffffffff
            i=(state*len(values))//4294967296
            total+=values[i]
        out.append(total/len(values))
    return [percentile(out,.025),percentile(out,.975)]


def analyze():
    rows=list(csv.DictReader(verify(CSV,CSV_GIT_BLOB).decode("utf8").splitlines()))
    original=json.loads(verify(RECEIPT,RECEIPT_GIT_BLOB).decode("utf8"))
    if len(rows)!=38 or original["source_ports_in_design"]!=38:
        raise ValueError("38 original port-source input changed")
    if original["historical_plant_query_complete"]!=114 or original["historical_plant_query_failed"]!=0:
        raise ValueError("Complete matched historical plant query support changed")
    result={
        "schema":"chocho_pre2011_38port_17island_top10_matched_index_uncertainty_v1",
        "status":"THIN_HISTORICAL_BOTANICAL_RECORD_ASSOCIATION_NO_CONFIRMED_ALIEN_HOST_ENRICHMENT",
        "sources":{"pair_csv_git_blob":CSV_GIT_BLOB,"receipt_git_blob":RECEIPT_GIT_BLOB},
        "source_ports":38,
        "source_islands":17,
        "geometric_biological_query_sites":114,
        "all_queries_complete":True,
        "species_key_cutoffs":{},
        "bootstrap":{"method":"resample eligible island means; non-independent same-island ports averaged first; fixed observation sources, not causal or full-uncertainty CI",
                     "seed":SEED,"draws":DRAWS,
                     "index_method":"upper 32-bit LCG scale to avoid modulo-4 pseudorandom cycling"},
        "limits":[
            "Globally top-ten resource candidate plants are NOT necessarily historically introduced on each island.",
            "GBIF recorded proportion is not vegetation abundance, larval feeding, breeding or colonization.",
            "Many islands have no adequate botanical index opportunities at port and both coast controls.",
            "One island's multiple port comparisons are not independent replications.",
            "Controls are away from listed WPI harbors but unlisted old ports or urban planting remain possible.",
            "2019 port positions are retrospective proxy locators for 2001–2010 botanical observations.",
            "Bootstrap resamples a very small number of eligible islands and ignores geographic/taxon dependence.",
            "When fewer than three islands meet quality cutoffs, no interval is reported.",
            "No causal port enrichment, ship transport, or subsequent butterfly range expansion can be concluded."
        ]
    }
    for threshold in THRESHOLDS:
        key=f"port_minus_control_top10_index_fraction_{threshold}"
        flag=f"meets_min_species_{threshold}"
        selected=[r for r in rows if r[flag]=="True"]
        grouped=defaultdict(list)
        for r in selected:
            grouped[r["island_id"]].append(float(r[key]))
        per_island=[
            {"island":island,"ports":len(values),"mean":sum(values)/len(values)}
            for island,values in grouped.items()
        ]
        means=[r["mean"] for r in per_island]
        measured=sum(means)/len(means) if means else None
        bootstrap=island_bootstrap(means,SEED+threshold,DRAWS) if len(means)>=3 else None
        leave_out=[sum(v for j,v in enumerate(means) if j!=i)/(len(means)-1)
                   for i in range(len(means))] if len(means)>1 else []
        saved=original["results_by_min_background_plant_species_keys"][str(threshold)]
        if len(selected)!=saved["eligible_port_pairs"] or len(grouped)!=saved["eligible_islands"]:
            raise ValueError("Archived group supports changed")
        if measured is not None and not math.isclose(
            measured,saved["island_equal_mean_port_minus_offport_host_fraction"],
            rel_tol=0,abs_tol=1e-12
        ):
            raise ValueError("Island-equal historical host index mean mismatch")
        result["species_key_cutoffs"][str(threshold)]={
            "eligible_ports":len(selected),
            "eligible_islands":len(means),
            "island_equal_mean_port_minus_offport_top10_index_fraction":measured,
            "island_bootstrap95":bootstrap,
            "leave_one_island_out_range":[min(leave_out),max(leave_out)] if leave_out else None,
            "positive_islands":sum(v>0 for v in means),
            "negative_islands":sum(v<0 for v in means),
            "zero_islands":sum(v==0 for v in means),
            "island_support":per_island
        }
    return result


def main():
    result=analyze()
    (HERE/"historical_38port_island_uncertainty.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8"
    )
    print(json.dumps({"status":result["status"],
                      "cutoffs":{k:{p:v[p] for p in ("eligible_ports","eligible_islands",
                          "island_equal_mean_port_minus_offport_top10_index_fraction",
                          "island_bootstrap95")}
                      for k,v in result["species_key_cutoffs"].items()}},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
