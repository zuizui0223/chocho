#!/usr/bin/env python3
"""Diagnose GBIF HTTP400s without repeating network calls.

Data: original Q4 500-pair historical plant source, fixed previously-failed
12-pair response artifact (retrieved from run 41412312725), official pinned
WGSRPD polygons. Distinguish WKT geometry errors from absent plant records.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter
from pathlib import Path
from shapely import wkt
from audit_q4_gbif_plant_first_year_pilot import build_geometries

def rows(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser()
    for name in ("retry-source","wgsrpd3-geojson","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    recs=rows(a.retry_source)
    if len(recs)!=12:raise RuntimeError(f"Expected exact 12 source results: {len(recs)}")
    valid_status={"HTTP_400","DATED_OCCURRENCE_FOUND","NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT"}
    if {x["status"] for x in recs}!=valid_status:
        raise RuntimeError("Unexpected GBIF retry source statuses")
    codes={r["island_code"] for r in recs}
    lookup=build_geometries(a.wgsrpd3_geojson,codes)
    details=[]
    for x in recs:
        region=x["island_code"]
        info=lookup[region]
        gwkt=info.get("wkt","")
        polygon=wkt.loads(gwkt) if gwkt else None
        pieces=len(polygon.geoms) if polygon and polygon.geom_type=="MultiPolygon" else 1 if polygon else 0
        details.append({
         "plant_id":x["plant_id"],"region":region,
         "status":x["status"],"geometry_type":polygon.geom_type if polygon else None,
         "wkt_chars":len(gwkt),"n_geographic_parts":pieces,
         "geometry_fidelity":info.get("geometry_fidelity"),
         "gbif_error_body_excerpt":x.get("http_error_body_excerpt",""),
        })
    errs=[r for r in details if r["status"]=="HTTP_400"]
    succeeded=[r for r in details if r["status"]!="HTTP_400"]
    report={
     "schema":"chocho_GBIF_historical_plant_year_HTTP400_geographic_diagnosis_v0.1",
     "source_retried_run":41412312725,
     "network_requests_made_by_this_diagnostic":0,
     "n_retried":len(recs),"n_http400":len(errs),
     "n_retried_dated":sum(r["status"]=="DATED_OCCURRENCE_FOUND" for r in details),
     "status_by_geometry_type":{
       category:dict(Counter(r["geometry_type"] for r in rr))
       for category,rr in (("HTTP_400",errs),("non_HTTP_400",succeeded))
     },
     "http400_body_excerpts":sorted({r["gbif_error_body_excerpt"] for r in errs if r["gbif_error_body_excerpt"]}),
     "per_region":details,
     "diagnosis_is_causal":False,
     "notes":[
       "The geometries are the original query geometries, not simplified bounding boxes or location guesses.",
       "An HTTP400 response is a rejected GBIF API request. Never mark plant as absent.",
       "If MULTIPOLYGON geometry is a consistent failure mode, split-query logic requires overlap-safe aggregation and quality checks before a full rerun.",
       "Avoid treating a diagnostic comparison of a few WKT geometries as proof that all 25 original 400 queries share a single cause.",
       "The successful dated 4 cases must be preserved, not folded into the original 219 response-blind record count without a new version."
     ]
    }
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in report.items() if k!="per_region"},ensure_ascii=False,indent=2))
if __name__=="__main__":main()
