#!/usr/bin/env python3
"""Low-rate diagnostic retry of fixed Q4 GBIF plant year-facet HTTP 429 failures.

Response-blind deterministic first 12 unresolved host-island pairs only, one
request at a time, ≥4s between requests. After two successive 429 failures,
STOP rather than increase load on public GBIF API. This is a SOURCE RETRY
diagnostic, not an added species sample and not plant arrival dating.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time,sys
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_q4_gbif_plant_first_year_pilot import build_geometries

API="https://api.gbif.org/v1/occurrence/search"
USER_AGENT="chocho-fixed-Q4-HTTP429-conservative-retry (https://github.com/zuizui0223/chocho)"

def readcsv(path):
    with path.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))
def csvwrite(path, data):
    if not data:return
    fields=list(data[0])
    for row in data[1:]:
        for name in row:
            if name not in fields: fields.append(name)
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(data)
def query_once(taxon,geo):
    params={"taxonKey":str(taxon),"geometry":geo,"hasCoordinate":"true",
            "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT","year":"1700,2026",
            "facet":"year","facetLimit":"350","limit":"0"}
    url=API+"?"+urlencode(params)
    if len(url)>15000:raise RuntimeError("URL_TOO_LONG")
    try:
        with urlopen(Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json",
                                          "Accept-Encoding":"identity"}),timeout=35) as f:
            if f.status!=200:return {"status":f"HTTP_{f.status}"}
            raw=f.read(10_000_000)
            obj=json.loads(raw)
        count=int(obj.get("count") or 0)
        facets=obj.get("facets") or []
        yr=next((x.get("counts") for x in facets if str(x.get("field","")).upper()=="YEAR"),None)
        if yr is None:return {"status":"FACET_MISSING","n_records":count}
        year_counts={}
        for x in yr:
            try:
                year=int(x["name"]);n=int(x["count"])
                if 1700<=year<=2026:year_counts[year]=n
            except (ValueError,TypeError,KeyError):continue
        if count>0 and not year_counts:return {"status":"FACET_EMPTY_WITH_RECORDS","n_records":count}
        if sum(year_counts.values())>count:return {"status":"FACET_INCONSISTENT","n_records":count}
        return {"status":"DATED_OCCURRENCE_FOUND" if year_counts else "NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT",
                "n_records":count,"earliest_indexed_year":min(year_counts) if year_counts else "",
                "latest_indexed_year":max(year_counts) if year_counts else "",
                "distinct_indexed_years":len(year_counts),"pre1900_records":sum(v for y,v in year_counts.items() if y<1900),
                "year_counts_sum":sum(year_counts.values()),"all_year_facets_complete":sum(year_counts.values())==count}
    except HTTPError as e:
        try: body=e.read(800).decode("utf-8",errors="replace")[:500]
        except Exception: body=""
        return {"status":f"HTTP_{e.code}","message":str(e)[:160],
                "http_error_body_excerpt":body}
    except (URLError,TimeoutError,json.JSONDecodeError) as e:
        return {"status":"NETWORK_ERROR","message":str(e)[:160]}

def main():
    ap=argparse.ArgumentParser()
    for k in ("source-pairs","original-audit","wgsrpd-geojson","out-dir"):
        ap.add_argument("--"+k,type=Path,required=True)
    ap.add_argument("--max-pairs",type=int,default=12)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    rows=readcsv(a.source_pairs)
    audit=json.loads(a.original_audit.read_text(encoding="utf-8"))
    if len(rows)!=500 or audit["pilot_dated_pairs"]!=219:
        raise RuntimeError("Source fixed 500-pair GBIF pilot drift")
    unqueried=[r for r in rows if r["q4_pilot_status"]=="GBIF_QUERY_FAILED" and
               "HTTP Error 429" in r.get("error","") and r["gbif_species_key"]]
    if len(unqueried)!=69:
        raise RuntimeError(f"Expected 69 true 429 failures; got {len(unqueried)}")
    if not 1<=a.max_pairs<=12:raise RuntimeError("Public API diagnostic must be <=12 queries")
    unqueried.sort(key=lambda r:hashlib.sha256(
        (r["plant_id"]+"|"+r["island_code"]).encode("utf-8")).hexdigest())
    selected=unqueried[:a.max_pairs]
    geometries=build_geometries(a.wgsrpd_geojson,{r["island_code"] for r in selected})
    if any(geometries[r["island_code"]]["status"]!="GEOMETRY_READY" for r in selected):
        raise RuntimeError("Original 429 retry has changed official region WKT resolution")
    results=[]
    consecutive_429=0
    halted=False
    for r in selected:
        if halted:break
        g=geometries[r["island_code"]]["wkt"]
        data=query_once(r["gbif_species_key"],g)
        if data["status"]=="HTTP_429":
            consecutive_429+=1
            # Stop at second successive 429; do not hammer service.
            if consecutive_429>=2:halted=True
        else:
            consecutive_429=0
        results.append({
            "plant_id":r["plant_id"],"plant_name":r["wcvp_name"],
            "island_code":r["island_code"],"gbif_species_key":r["gbif_species_key"],
            "original_status":r["q4_pilot_status"],
            **data
        })
        print("Q4_CONSERVATIVE_RETRY",len(results),a.max_pairs,r["island_code"],data["status"],flush=True)
        if not halted:time.sleep(4)
    if not results:raise RuntimeError("No GBIF requests issued")
    csvwrite(a.out_dir/"http429_original_pilot_retry_subset.csv",results)
    statuses=Counter(x["status"] for x in results)
    summary={
      "schema":"chocho_q4_fixed500_public_GBIF_HTTP429_conservative_retry_v0.1",
      "status":"LOW_RATE_TECHNICAL_RECOVERY_DIAGNOSTIC_NOT_NEW_COLONIZATION_DATA",
      "original_fixed_500_pairs":500,"original_dated_pairs":219,
      "original_http429_failed_pairs":69,"selected_for_retry":len(selected),
      "attempted":len(results),"not_attempted_in_selected":len(selected)-len(results),
      "stopped_after_two_consecutive_429":halted,
      "status_counts":dict(statuses),
      "source_selection":"SHA256(plant_id|WGSRPD3_code), fixed before API retry, no selection on butterfly presence",
      "new_dated_pairs_in_original_failure_cohort":statuses["DATED_OCCURRENCE_FOUND"],
      "zero_indexed_plant_records_on_retry":statuses["NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT"],
      "no_automatic_imputation_of_failed_queries":True,
      "limits":[
       "Technical retries can improve record-coverage estimates but not detect original plant introduction or naturalization year.",
       "This is only up to 12 of 69 original HTTP 429 failures; success fraction must not be extrapolated to all missing pairs.",
       "Fixed original 500 response-blind source and successful outcomes remain unchanged.",
       "The speciesKey/WKT request reuses the same original GBIF year facet query, 1700–2026.",
       "HTTP 400 source errors are not retried or recoded as absence.",
       "Two consecutive GBIF 429 throttles halt queries to protect API service."
      ]
    }
    (a.out_dir/"http429_retry_audit.json").write_text(
      json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()
