#!/usr/bin/env python3
"""Validate GBIF repeated POLYGON geometry against a known control, then rescue
two fixed original HTTP400 host×island pairs by querying disjoint polygon
chunks. No missing portions are ever imputed absent and no first indexed
year is interpreted as plant introduction year.
"""
from __future__ import annotations
import argparse,csv,json,time
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError
from shapely import wkt
from audit_q4_gbif_plant_first_year_pilot import build_geometries

MAX_URL=3500
MAX_CALLS=20
PAUSE=4.5

def records(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def savecsv(path, data):
    if not data:return
    fields=[]
    for row in data:
        for name in row:
            if name not in fields: fields.append(name)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(data)

def geo_parts(region,lookup):
    full=wkt.loads(lookup[region]["wkt"])
    if full.geom_type=="Polygon":polys=[full]
    elif full.geom_type=="MultiPolygon":polys=list(full.geoms)
    else:raise RuntimeError("Not a polygon-based WGSRPD region "+region)
    if not full.is_valid or any(not p.is_valid for p in polys):
        raise RuntimeError("Invalid polygon geometry "+region)
    areasum=sum(p.area for p in polys)
    if areasum<=0 or abs(areasum-full.area)>max(1e-9,areasum*1e-9):
        raise RuntimeError("Overlapping polygon areas; cannot aggregate without deduplication")
    return [wkt.dumps(p,rounding_precision=4,trim=True) for p in polys]

def makeurl(key,polys):
    items=[("taxonKey",str(key)),("hasCoordinate","true"),("hasGeospatialIssue","false"),
           ("occurrenceStatus","PRESENT"),("year","1700,2026"),
           ("facet","year"),("facetLimit","350"),("limit","0")]
    items.extend(("geometry",p) for p in polys)
    return "https://api.gbif.org/v1/occurrence/search?"+urlencode(items)

def split(key,polys):
    groups=[]
    chunk=[]
    for poly in polys:
        if len(makeurl(key,[poly]))>MAX_URL:
            raise RuntimeError("Single component too long; skip rather than drop part")
        if chunk and len(makeurl(key,chunk+[poly]))>MAX_URL:
            groups.append(chunk)
            chunk=[]
        chunk.append(poly)
    if chunk:groups.append(chunk)
    return groups

def search(key,polys):
    url=makeurl(key,polys)
    if len(url)>MAX_URL:raise RuntimeError("URL exceeds fixed cap")
    try:
        with urlopen(Request(url,headers={"User-Agent":"chocho-GbifYearPolygonRescue/0.1",
                                          "Accept":"application/json"}),timeout=45) as r:
            obj=json.loads(r.read(6_000_000))
        n=int(obj.get("count") or 0)
        facet=next((x.get("counts") for x in (obj.get("facets") or [])
                    if str(x.get("field","")).upper()=="YEAR"),None)
        if facet is None:return {"status":"YEAR_FACET_MISSING","records":n}
        annual={}
        for item in facet:
            try:
                year=int(item["name"]); v=int(item["count"])
                if 1700<=year<=2026 and v>=0:annual[year]=v
            except (KeyError,TypeError,ValueError):
                continue
        if n>0 and not annual:return {"status":"YEAR_FACET_EMPTY_WITH_RECORDS","records":n}
        if sum(annual.values())>n:return {"status":"YEAR_FACET_INCONSISTENT","records":n}
        return {"status":"OK","records":n,"annual":annual}
    except HTTPError as e:
        try:body=e.read(260).decode("utf-8",errors="replace")
        except Exception:body=""
        return {"status":"HTTP_"+str(e.code),"body":body}
    except (URLError,TimeoutError,json.JSONDecodeError) as e:
        return {"status":"NETWORK_ERROR","error":str(e)[:200]}

def run_case(row,lookup,shared,kind):
    reg=row["island_code"]
    polygons=geo_parts(reg,lookup)
    chunks=split(row["gbif_species_key"],polygons)
    if shared["calls"]+len(chunks)>MAX_CALLS:
        return {"kind":kind,"region":reg,"status":"TOO_MANY_REQUIRED_REQUESTS","chunks":len(chunks)},[]
    byyear=defaultdict(int)
    counts=0
    details=[]
    for i,poly in enumerate(chunks,1):
        if shared["stopped"]:break
        if shared["calls"]:time.sleep(PAUSE)
        response=search(row["gbif_species_key"],poly)
        shared["calls"]+=1
        code=response["status"]
        if code=="HTTP_429":shared["consecutive429"]+=1
        else:shared["consecutive429"]=0
        if shared["consecutive429"]>=2:shared["stopped"]=True
        details.append({"region":reg,"plant_id":row["plant_id"],"kind":kind,
                        "chunk":i,"chunks_total":len(chunks),"polygons":len(poly),
                        "url_chars":len(makeurl(row["gbif_species_key"],poly)),
                        "status":code,"records":response.get("records",""),
                        "message":response.get("body",response.get("error",""))[:200]})
        print("YEAR_POLYGON_QUERY",reg,i,len(chunks),len(poly),code,flush=True)
        if code!="OK":break
        counts+=response["records"]
        for year,n in response["annual"].items():byyear[year]+=n
    complete=(len(details)==len(chunks) and all(d["status"]=="OK" for d in details))
    result={"kind":kind,"region":reg,"plant_id":row["plant_id"],"plant_name":row["plant_name"],
            "source_status":row["status"],"source_component_polygons":len(polygons),
            "total_chunks":len(chunks),"successful_chunks":sum(x["status"]=="OK" for x in details),
            "complete_full_WGSRPD_region":complete,
            "status":("DATED_OCCURRENCE_FOUND" if byyear else
                     "NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT") if complete else "PARTIAL_REGION_UNKNOWN",
            "total_records_if_complete":counts if complete else "",
            "earliest_indexed_year_if_complete":min(byyear) if complete and byyear else "",
            "latest_indexed_year_if_complete":max(byyear) if complete and byyear else "",
            "dated_yearfacet_records_if_complete":sum(byyear.values()) if complete else ""}
    return result,details

def main():
    ap=argparse.ArgumentParser()
    for name in ("retry-source","wgsrpd3-geojson","out-dir"):
        ap.add_argument("--"+name,required=True,type=Path)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    source=records(a.retry_source)
    if len(source)!=12:raise RuntimeError("Fixed retry source should have exactly 12 records")
    failed=[r for r in source if r["status"]=="HTTP_400"]
    control=[r for r in source if r["plant_id"]=="223906" and r["island_code"]=="MRQ"]
    if len(failed)!=7 or len(control)!=1:
        raise RuntimeError("Source HTTP400 and reference species were changed")
    lookup=build_geometries(a.wgsrpd3_geojson,{r["island_code"] for r in failed+control})
    # 2 unique island regions, sorted by increasing source geometry complexity.
    failed.sort(key=lambda r:(len(geo_parts(r["island_code"],lookup)),r["plant_id"]))
    unique=[]
    used=set()
    for r in failed:
        if r["island_code"] not in used:
            used.add(r["island_code"]);unique.append(r)
        if len(unique)==2:break
    rows_out=[]
    detail=[]
    shared={"calls":0,"consecutive429":0,"stopped":False}
    ref,log=run_case(control[0],lookup,shared,"VALIDATION")
    rows_out.append(ref);detail.extend(log)
    validated=bool(ref["complete_full_WGSRPD_region"] and
                   ref["total_records_if_complete"]==6 and ref["earliest_indexed_year_if_complete"]==2004)
    if validated:
        for r in unique:
            if shared["stopped"]:break
            row,log=run_case(r,lookup,shared,"RESCUE")
            rows_out.append(row);detail.extend(log)
    else:
        print("VALIDATION_DID_NOT_REPLICATE_SOURCE_6_RECORDS_FIRST2004_STOP",flush=True)
    savecsv(a.out_dir/"split_wgsrpd_years.csv",rows_out)
    savecsv(a.out_dir/"split_polygon_chunk_query_receipt.csv",detail)
    recovered=[r for r in rows_out if r["kind"]=="RESCUE" and r["complete_full_WGSRPD_region"]]
    audit={"schema":"chocho_GBIF_Q4_repeatable_polygon_split_year_rescue_v0.1",
           "status":"VALIDATED_SPLIT_POLYGON_GBIF_QUERY" if validated else "STOP_REFERENCE_NOT_EQUIVALENT",
           "reference_MRQ_Canna_indica_original_records":6,
           "reference_MRQ_Canna_indica_original_first_year":2004,
           "positive_control_match":validated,
           "selected_http400_islands":list(x["island_code"] for x in unique),
           "original_http400_failures":7,
           "total_GBIF_queries":shared["calls"],"max_queries":MAX_CALLS,
           "stopped_after_two_consecutive_429":shared["stopped"],
           "fully_recovered_original_400_pairs":len(recovered),
           "new_dated_host_island_pairs":sum(x["status"]=="DATED_OCCURRENCE_FOUND" for x in recovered),
           "fully_covered_no_plant_index_pairs":sum(x["status"]=="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT" for x in recovered),
           "per_case":rows_out,
           "limits":[
             "Geographically disjoint source polygons, but boundary points can be counted in adjacent components.",
             "If any chunk fails, first year for the full region is unknown. Partial search is not absence.",
             "GBIF indexed plant years are not introduction, establishment or residence time.",
             "Two original 400 failures are selected by source polygon complexity, not by ecological outcome.",
             "Public API rate limits and per-query URL length cap bound the diagnostic.",
             "No butterfly occurrences or historical transport routes are inferred."
           ]}
    (a.out_dir/"split_polygon_year_audit.json").write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(audit,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
