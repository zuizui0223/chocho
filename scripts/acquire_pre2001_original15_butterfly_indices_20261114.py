#!/usr/bin/env python3
"""Retrieve pre-2001 GBIF indices for the 15 putative first-index island cases.

The 2026 GBIF facet panel used an arbitrary 2001 baseline. Existing butterflies
before 2001 must NOT be interpreted as 2011–2017 immigrants. Historical GBIF
nonindex remains UNKNOWN status; literature may show earlier populations.

Uses exact frozen original GSHHG 2.3.7, same speciesKeys/island IDs, and
bounded 1700–2000 GBIF georeferenced species facet source queries.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

from build_public51_butterfly_risk_20261109 import (
    GBIF, GSHHG_SHA256, get_json, load_polygons,
)

START_YEAR=1700
END_YEAR=2000


def write_csv(path, rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def audit():
    p=argparse.ArgumentParser()
    p.add_argument("--gshhg-zip",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--source-case-audit",type=Path,
                   default=Path("results/port_host_reindex_20261114/audit.json"))
    p.add_argument("--out-dir",type=Path,
                   default=Path("results/port_host_pre2001_20261114"))
    p.add_argument("--max-facets",type=int,default=50)
    p.add_argument("--wkt-max-chars",type=int,default=2300)
    args=p.parse_args()

    raw=args.source_case_audit.read_bytes()
    cases=json.loads(raw)["introduced_case_index"]
    if len(cases)!=15:
        raise ValueError(f"Expected 15 original host-indexed first-event candidates, got {len(cases)}")
    species_by_island=defaultdict(set)
    for c in cases:
        species_by_island[str(c["island"])].add(str(c["species"]))
    islands=set(species_by_island)
    geom,geom_audit,gshhg_sha=load_polygons(args.gshhg_zip,args.work,islands,args.wkt_max_chars)
    if gshhg_sha!=GSHHG_SHA256:
        raise ValueError("Unexpected polygon source")
    expected=sum(map(len,geom.values()))
    if expected>args.max_facets:
        raise ValueError(f"Planned {expected} source fragments exceed predeclared cap {args.max_facets}")
    print("PREDECLARED_1700_2000_QUERIES",expected,flush=True)

    counts=defaultdict(int)
    queries=[]
    errors=[]
    for isl in sorted(species_by_island,key=int):
        keys=sorted(species_by_island[isl],key=int)
        for i,wkt in enumerate(geom[isl],1):
            params={
                "geometry":wkt,
                "year":f"{START_YEAR},{END_YEAR}",
                "taxonKey":keys,
                "hasCoordinate":"true",
                "hasGeospatialIssue":"false",
                "occurrenceStatus":"PRESENT",
                "facet":"speciesKey","facetLimit":1000,
                "limit":0,
            }
            try:
                response=get_json(GBIF+"/occurrence/search",params)
                n=int(response["count"])
                facet={str(t["field"]).upper():t.get("counts",[]) for t in response.get("facets",[])}
                if n>0 and "SPECIES_KEY" not in facet:
                    raise ValueError("Nonempty GBIF index without speciesKey facet")
                found=defaultdict(int)
                for item in facet.get("SPECIES_KEY",[]):
                    key=str(item.get("name"))
                    if key in keys:
                        found[key]+=int(item["count"])
                if sum(found.values())>n:
                    raise ValueError("Species facet sum exceeds total")
                for key,value in found.items():
                    counts[(isl,key)]+=value
                queries.append({
                    "island":isl,"fragment":i,"taxonKeys":keys,
                    "period":f"{START_YEAR}-{END_YEAR}",
                    "total_indexed_matching_query":n,
                    "candidate_taxon_indexed":dict(found),
                    "facet_status":"COMPLETE",
                    "source_polygon_WKT_SHA256":hashlib.sha256(wkt.encode()).hexdigest(),
                    "GBIF_response_sha256":hashlib.sha256(json.dumps(response,sort_keys=True).encode()).hexdigest(),
                })
            except Exception as e:
                errors.append({"island":isl,"fragment":i,"error":repr(e)})
                queries.append({"island":isl,"fragment":i,"taxonKeys":keys,
                                "period":f"{START_YEAR}-{END_YEAR}","facet_status":"INCOMPLETE_UNKNOWN",
                                "source_polygon_WKT_SHA256":hashlib.sha256(wkt.encode()).hexdigest()})
            print("HIST_FACET",isl,i,queries[-1]["facet_status"],
                  queries[-1].get("candidate_taxon_indexed",{}),flush=True)
            time.sleep(0.1)

    records=[]
    for c in cases:
        isl,key=str(c["island"]),str(c["species"])
        ok=all(q["facet_status"]=="COMPLETE" for q in queries if q["island"]==isl)
        records.append({
            "species":c["name"],"speciesKey":key,"island":isl,
            "first_index_window":"2011-2017",
            "pre2001_indexed_records":counts[(isl,key)] if ok else "",
            "pre2001_positive":int(counts[(isl,key)]>0) if ok else "",
            "historical_GBIF_status":"INDEX_POSITIVE" if ok and counts[(isl,key)]>0
            else ("INDEX_ZERO_UNKNOWN_OCCUPANCY" if ok else "INCOMPLETE_SOURCE_UNKNOWN"),
            "documented_external_historical_source":
                "Powell 1927/Miller 1984-85 resident" if (key,isl)==("1920496","951") else "",
            "do_not_infer_true_recent_colonization":True,
        })
    args.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(args.out_dir/"original15_pre2001_species_island_indices.csv",records)
    audit={
        "status":"COMPLETED" if not errors else "PARTIAL_SOURCE_ERRORS",
        "snapshot_time_utc":datetime.now(timezone.utc).isoformat(),
        "source_cases_sha256":hashlib.sha256(raw).hexdigest(),
        "GSHHG_zip_sha256":gshhg_sha,
        "query_count_planned":expected,"query_count_executed":len(queries),
        "queries":queries,"query_errors":errors,
        "geometry_receipts":geom_audit,
        "case_records":records,
        "original15_with_positive_pre2001_index":[x for x in records if x["historical_GBIF_status"]=="INDEX_POSITIVE"],
        "pre2001_GBIF_index_zero_not_true_absence":True,
        "exclusion_counterexample":"P. rapae Catalina has independent 1927 published record even if pre2001 GBIF index is zero",
        "limitations":[
            "Only original 15 source-selected post-2011 GBIF first-index cells examined",
            "No original private 83,573-row time series analyzed",
            "Pre2001 nonindex is unknown occupancy, NOT no previous colonization",
            "Other historical sources and uncatalogued butterflies are not covered",
            "Retrospective data may drift with GBIF source indexing"
        ]
    }
    (args.out_dir/"source_audit.json").write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n")
    print("SUMMARY",json.dumps({
        "status":audit["status"],
        "cases":len(records),"queries":len(queries),"errors":len(errors),
        "positive_pre2001":len(audit["original15_with_positive_pre2001_index"]),
        "positives":[(z["species"],z["island"],z["pre2001_indexed_records"])
                     for z in audit["original15_with_positive_pre2001_index"]]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    audit()
