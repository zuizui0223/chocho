#!/usr/bin/env python3
"""GBIF pre-2001 historical index for *every* frozen 225 species x 51 islands.

Source-locked temporal eligibility, not true absence/arrival. Unlike the
positive-event-only 317 screen, this uses identical historical GBIF queries on
ALL 11,475 frozen butterfly x original-island cells, so a well-defined
longitudinal index risk set can be reconstructed without response-conditioned
historical filtering.

Requires archived GSHHG 2.3.7 polygons and original frozen 225 GBIF keys.
One bounded speciesKey facet query per disjoint island polygon fragment. Request
facets must contain exclusively the frozen speciesKey panel or fail closed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

from build_public51_butterfly_risk_20261109 import GBIF,GSHHG_SHA256,get_json,load_polygons

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"

def rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def emit(path,values):
    if not values:raise ValueError("Missing expected source cells")
    with path.open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(values[0]))
        writer.writeheader(); writer.writerows(values)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gshhg-zip",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--batch-id",required=True)
    p.add_argument("--islands",required=True)
    p.add_argument("--out-dir",required=True,type=Path)
    p.add_argument("--max-queries",type=int,default=120)
    p.add_argument("--wkt-max-chars",type=int,default=2300)
    a=p.parse_args()

    source=rows(SOURCE)
    if len(source)!=11475 or len({(x["speciesKey"],x["landmass_id"]) for x in source})!=11475:
        raise ValueError("Frozen original 225 x 51 source changed")
    valid={x["landmass_id"] for x in source}
    ids=set(a.islands.split(","))
    if not ids or not ids<=valid:raise ValueError("Batch contains invalid original GSHHG IDs")
    species={x["speciesKey"] for x in source}
    if len(species)!=225:raise ValueError("Original taxonKeys changed")
    selection=[x for x in source if x["landmass_id"] in ids]
    if len(selection)!=len(ids)*225:raise ValueError("Incomplete original island species cells")

    polygons,shapereceipts,digest=load_polygons(
        a.gshhg_zip,a.work,ids,a.wkt_max_chars)
    if digest!=GSHHG_SHA256:raise ValueError("Wrong GSHHG archive")
    n=sum(len(x) for x in polygons.values())
    if n>a.max_queries:raise ValueError(f"{n} original fragments exceed cap {a.max_queries}")
    print("FULL_HISTORY_BATCH_START",json.dumps({
        "batch":a.batch_id, "original_islands":len(ids),"panel_species":225,
        "total_cells":len(selection), "original_fragments":n}),flush=True)
    source_key_list=sorted(species,key=int)
    counts=defaultdict(int)
    receipts=[]
    failures=[]
    for island in sorted(ids,key=int):
        for part,polygon in enumerate(polygons[island],1):
            receipt={
                "island":island,"part":part,"taxonKey_count":225,
                "source_fragment_WKT_sha256":hashlib.sha256(polygon.encode()).hexdigest(),
                "year":"1700-2000",
                "source_status":"UNKNOWN",
            }
            try:
                result=get_json(GBIF+"/occurrence/search",{
                    "geometry":polygon,
                    "year":"1700,2000",
                    "taxonKey":source_key_list,
                    "hasCoordinate":"true",
                    "hasGeospatialIssue":"false",
                    "occurrenceStatus":"PRESENT",
                    "facet":"speciesKey","facetLimit":1000,
                    "limit":0,
                })
                total=int(result["count"])
                facet={str(x.get("field","")).upper():x.get("counts",[])
                    for x in result.get("facets",[])}
                if total>0 and "SPECIES_KEY" not in facet:
                    raise ValueError("Nonzero GBIF count without SPECIES_KEY facet")
                observed=defaultdict(int)
                for z in facet.get("SPECIES_KEY",[]):
                    key=str(z["name"])
                    if key not in species:raise ValueError("GBIF returned speciesKey not frozen in panel "+key)
                    observed[key]+=int(z["count"])
                if sum(observed.values())>total:
                    raise ValueError("GBIF facet sum exceeds total")
                for key,count in observed.items():
                    counts[island,key]+=count
                receipt.update({
                    "source_status":"COMPLETE",
                    "indexed_total":total,
                    "candidate_indexed":sum(observed.values()),
                    "source_response_sha256":hashlib.sha256(
                        json.dumps(result,sort_keys=True,ensure_ascii=False).encode()
                    ).hexdigest(),
                    "source_faceted_species":len(observed),
                })
            except Exception as exc:
                receipt.update(source_status="INCOMPLETE_UNKNOWN",error=repr(exc))
                failures.append(receipt)
            receipts.append(receipt)
            print("HISTORICAL_ALL225",a.batch_id,island,part,receipt["source_status"],
                  receipt.get("candidate_indexed"),flush=True)
            time.sleep(0.08)

    complete={isl:all(r["source_status"]=="COMPLETE"
            for r in receipts if r["island"]==isl) for isl in ids}
    cases=[]
    for x in selection:
        isl,key=x["landmass_id"],x["speciesKey"]
        before=counts[isl,key] if complete[isl] else None
        old=int(x["pre2001_2010_indexed_records"])
        early=int(x["future2011_2017_indexed_records"])
        later=int(x["future2018_2022_indexed_records"])
        source_pre=before is not None and before>0
        cases.append({
            "species":x["species"],"speciesKey":key,
            "landmass_id":isl,
            "pre2001_indexed_records":before if before is not None else "",
            "pre2001_source_complete":int(complete[isl]),
            "pre2001_any_index":int(source_pre) if before is not None else "",
            "gbif_index_2001_2010":old,
            "gbif_index_2011_2017":early,
            "gbif_index_2018_2022":later,
            "original_2011_2017_first_index":int(x["first_indexed_in_2011_2017"]),
            "original_2018_2022_first_index":int(x["first_indexed_in_2018_2022"]),
            "eligible_2011_2017_longer_history":int(before==0 and old==0)
                if before is not None else "",
            "indexed_first_2011_2017_longer_history":int(before==0 and old==0 and early>0)
                if before is not None else "",
            "eligible_2018_2022_longer_history":int(before==0 and old==0 and early==0)
                if before is not None else "",
            "indexed_first_2018_2022_longer_history":int(before==0 and old==0 and early==0 and later>0)
                if before is not None else "",
            "index_nonobservation_not_biological_absence":True,
            "indexed_first_not_verified_arrival":True,
        })
    a.out_dir.mkdir(parents=True,exist_ok=True)
    emit(a.out_dir/f"full_history_batch_{a.batch_id}.csv",cases)
    audit={
        "batch_id":a.batch_id,
        "source_status":"COMPLETE" if not failures else "PARTIAL",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "source_11_475_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "GSHHG_sha256":digest,
        "original_islands":sorted(ids,key=int),
        "original_source_case_count":len(cases),
        "source_fragment_queries":len(receipts),
        "source_fragment_failures":failures,
        "source_polygon_geometry":shapereceipts,
        "source_query_receipts":receipts,
        "pre2001_positive_all_cells":sum(r["pre2001_any_index"]==1 for r in cases),
        "original_early_first_count":sum(r["original_2011_2017_first_index"]==1 for r in cases),
        "original_late_first_count":sum(r["original_2018_2022_first_index"]==1 for r in cases),
        "new_early_first_count":sum(r["indexed_first_2011_2017_longer_history"]==1 for r in cases),
        "new_late_first_count":sum(r["indexed_first_2018_2022_longer_history"]==1 for r in cases),
        "not_inferred":[
            "true prior butterfly absence","actual immigration","breeding",
            "food plant use", "shipborne propagation","wind-driven routes",
        ],
    }
    (a.out_dir/f"full_history_receipt_batch_{a.batch_id}.json").write_text(
        json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("FULL_HISTORY_BATCH_SUMMARY",json.dumps({
        k:v for k,v in audit.items()
        if k in ("batch_id","source_status","original_source_case_count",
                 "source_fragment_queries","pre2001_positive_all_cells",
                 "original_early_first_count","original_late_first_count",
                 "new_early_first_count","new_late_first_count")
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
