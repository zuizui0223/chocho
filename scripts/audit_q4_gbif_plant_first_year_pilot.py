#!/usr/bin/env python3
"""GBIF historical plant-index coverage for a response-blind fixed Q4 pilot.

Queries first GBIF indexed PLANT occurrence year by exact accepted species key
within documented WGSRPD3 geometry. This is NOT the introduction year, and a
blank or zero count is NEVER treated as botanical absence. Per-island annual
all-plant indexed record effort is queried with identical geometry.

Reproducibility: pinned WCVP IDs/names + GBIF query metadata, source checks,
precommitted island-host pilot and failure/censoring statuses.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,math,statistics,time
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError,URLError
from shapely.geometry import shape
from shapely.geometry.polygon import orient
from shapely.ops import unary_union,transform
from shapely import wkt,make_valid
from pyproj import Transformer

API="https://api.gbif.org/v1"
AGENT="chocho-Q4-host-detection-lag-audit/0.1 (https://github.com/zuizui0223/chocho)"
YEARS="1700,2026"
MAX_WKT=7500
MAX_RESPONSE_BYTES=10_000_000
STATUS_NO_EVIDENCE="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT"

def readcsv(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))
def outcsv(path,rr):
    if not rr:
        return
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]))
        w.writeheader();w.writerows(rr)
def desc(x):
    arr=sorted(x)
    return {"min":arr[0],"median":statistics.median(arr),"max":arr[-1]} if arr else {}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def request(path,params,retries=3):
    url=API+path+"?"+urlencode(params)
    if len(url)>15000:raise ValueError(f"Query URL too long ({len(url)})")
    last=None
    for k in range(retries):
        try:
            with urlopen(Request(url,headers={
                "User-Agent":AGENT,"Accept":"application/json",
                "Accept-Encoding":"identity"}),timeout=24) as r:
                if r.status!=200:raise RuntimeError(f"GBIF HTTP {r.status}")
                raw=r.read(MAX_RESPONSE_BYTES+1)
                if len(raw)>MAX_RESPONSE_BYTES:raise RuntimeError("Response >10MB")
                data=json.loads(raw)
                if not isinstance(data,dict) or data.get("error"):
                    raise RuntimeError(f"GBIF error: {str(data.get('error'))[:300]}")
                return data
        except (HTTPError,TimeoutError,URLError,json.JSONDecodeError,RuntimeError) as e:
            last=e
            if isinstance(e,HTTPError) and e.code not in (429,500,502,503,504):
                break
            time.sleep(2**k)
    raise RuntimeError(f"{type(last).__name__}: {str(last)[:180]}")

def key_from_gbif(name):
    try:
        obj=request("/species/match",{"name":name,"rank":"SPECIES","strict":"true"})
    except Exception as e:
        return {"status":"API_ERROR","key":"","canonical":"","confidence":"",
                "match_type":"","detail":str(e)}
    canonical=str(obj.get("canonicalName") or "").strip()
    rank=str(obj.get("rank") or "").upper()
    matched=str(obj.get("matchType") or "").upper()
    confidence=int(obj.get("confidence") or 0)
    accepted=int(obj.get("speciesKey") or obj.get("acceptedUsageKey") or obj.get("usageKey") or 0)
    if not accepted or rank!="SPECIES" or matched!="EXACT":
        return {"status":"NO_EXACT_ACCEPTED_SPECIES","key":"","canonical":canonical,
                "confidence":confidence,"match_type":matched,"detail":str(obj.get("status") or "")}
    # exact name equality with authorship removed / accepted concept synonym;
    # do not accept higher-rank or accepted binomial mismatch without evidence.
    canon_norm=" ".join(canonical.lower().split())
    input_norm=" ".join(str(name).lower().split())
    if canon_norm!=input_norm:
        return {"status":"CANONICAL_NAME_MISMATCH_REVIEW","key":"",
                "canonical":canonical,"confidence":confidence,
                "match_type":matched,"detail":str(obj.get("status") or "")}
    if confidence<90:
        return {"status":"LOW_CONFIDENCE_EXACT_REVIEW","key":"",
                "canonical":canonical,"confidence":confidence,
                "match_type":matched,"detail":str(obj.get("status") or "")}
    return {"status":"EXACT_SPECIES_MATCH","key":accepted,"canonical":canonical,
            "confidence":confidence,"match_type":matched,"detail":str(obj.get("status") or "")}

def normalize_wkt(g):
    if g.geom_type=="Polygon":return orient(g,sign=1.0)
    if g.geom_type=="MultiPolygon":
        return type(g)([orient(x,sign=1.0) for x in g.geoms])
    return g

def build_geometries(geojson, requested):
    source=json.loads(geojson.read_text(encoding="utf-8"))
    coord_by_code=defaultdict(list)
    for f in source["features"]:
        code=str((f.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        if code in requested and f.get("geometry"):
            geom=make_valid(shape(f["geometry"]))
            if geom.geom_type in ("Polygon","MultiPolygon"):coord_by_code[code].append(geom)
    trans=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform
    out={}
    for code in sorted(requested):
        if code not in coord_by_code:
            out[code]={"status":"MISSING_OFFICIAL_GEOMETRY","wkt":"","geometry_fidelity":None}
            continue
        orig=make_valid(unary_union(coord_by_code[code]))
        lonmin,latmin,lonmax,latmax=orig.bounds
        if lonmax-lonmin>180:
            out[code]={"status":"ANTIMERIDIAN_GEOMETRY_REVIEW","wkt":"",
                       "geometry_fidelity":None,"original_bounds":[lonmin,latmin,lonmax,latmax]}
            continue
        selected=None
        for tol in (0,0.001,0.003,0.0075,0.015,0.03,0.05,0.10,0.20):
            candidate=normalize_wkt(make_valid(orig.simplify(tol,preserve_topology=True)))
            if candidate.geom_type not in ("Polygon","MultiPolygon"):continue
            value=wkt.dumps(candidate,rounding_precision=4,trim=True)
            if len(value)>MAX_WKT:continue
            try:
                oo=make_valid(transform(trans,orig))
                cc=make_valid(transform(trans,candidate))
                denom=oo.union(cc).area
                quality=1.0-(oo.symmetric_difference(cc).area/denom) if denom>0 else 0
            except Exception:quality=0
            if quality>=0.97:
                selected={"status":"GEOMETRY_READY","wkt":value,
                          "geometry_fidelity":quality,"wkt_length":len(value),
                          "simplify_tolerance_degrees":tol}
                break
        out[code]=selected or {"status":"WKT_COMPLEX_OR_INACCURATE","wkt":"",
                                "geometry_fidelity":None,
                                "original_bounds":[lonmin,latmin,lonmax,latmax]}
    return out

def facet_year(key,geometry):
    params={"taxonKey":str(key),"geometry":geometry,
            "hasCoordinate":"true","hasGeospatialIssue":"false",
            "occurrenceStatus":"PRESENT","year":YEARS,
            "facet":"year","facetLimit":"350","limit":"0"}
    obj=request("/occurrence/search",params)
    count=int(obj.get("count") or 0)
    facets=obj.get("facets") or []
    facets_year=next((f.get("counts") or [] for f in facets
        if str(f.get("field") or "").upper()=="YEAR"),None)
    if facets_year is None:raise RuntimeError("GBIF did not return year facet")
    annual={}
    for item in facets_year:
        try:
            year=int(item["name"])
            n=int(item["count"])
        except (ValueError,KeyError,TypeError):continue
        if 1700<=year<=2026 and n>=0:annual[year]=n
    if count>0 and not annual:
        raise RuntimeError("GBIF positive count but no year facets")
    if sum(annual.values())>count:
        raise RuntimeError("GBIF year facets exceed count")
    return {"status":"DATED_OCCURRENCE_FOUND" if annual else STATUS_NO_EVIDENCE,
            "gbif_record_count":count,
            "record_sum_year_facets":sum(annual.values()),
            "earliest_indexed_year":min(annual) if annual else "",
            "latest_indexed_year":max(annual) if annual else "",
            "distinct_recorded_years":len(annual),
            "pre1900_records":sum(v for y,v in annual.items() if y<1900),
            "pre1950_records":sum(v for y,v in annual.items() if y<1950),
            "pre2001_records":sum(v for y,v in annual.items() if y<2001),
            "pre2011_records":sum(v for y,v in annual.items() if y<2011),
            "recorded_last_decade":sum(v for y,v in annual.items() if y>=2016),
            "year_facet_complete_within_1700_2026":sum(annual.values())==count}

def runner(pairs,names,geos,workers,outdir):
    distinct_plants=sorted({r["plant_id"] for r in pairs})
    matched={}
    def taxon_task(pid):
        v=key_from_gbif(names[pid])
        return pid,v
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future={pool.submit(taxon_task,pid):pid for pid in distinct_plants}
        for i,f in enumerate(as_completed(future),1):
            pid,result=f.result()
            matched[pid]=result
            if i%40==0:print("GBIF_TAXON_MATCH",i,len(distinct_plants),flush=True)
    taxonrows=[{"plant_id":pid,"wcvp_accepted_name":names[pid],
                **matched[pid]} for pid in distinct_plants]
    outcsv(outdir/"gbif_taxon_resolution.csv",taxonrows)
    counts=Counter(row["status"] for row in taxonrows)
    print("GBIF_TAXON_STATUSES",json.dumps(counts),flush=True)
    # One background botanical year-facet query per island where the fixed
    # geometry is a fidelity-checked WKT, regardless of host detection.
    q=[]
    for isle,geom in geos.items():
        if geom["status"]=="GEOMETRY_READY":
            q.append(("EFFORT",isle,"6",geom["wkt"]))
    for row in pairs:
        isle,pid=row["island_code"],row["plant_id"]
        if geos[isle]["status"]=="GEOMETRY_READY" and matched[pid]["status"]=="EXACT_SPECIES_MATCH":
            q.append(("PLANT",pid+"|"+isle,str(matched[pid]["key"]),geos[isle]["wkt"]))
    answers={}
    def facet_task(item):
        kind,key,taxon,geom=item
        try:r=facet_year(taxon,geom)
        except Exception as e:r={"status":"GBIF_QUERY_FAILED","gbif_record_count":"",
                     "earliest_indexed_year":"","error":str(e)}
        return kind,key,r
    with ThreadPoolExecutor(max_workers=workers) as pool:
        fut={pool.submit(facet_task,item):item for item in q}
        for i,f in enumerate(as_completed(fut),1):
            kind,key,data=f.result()
            answers[(kind,key)]=data
            if i%40==0:print("GBIF_YEAR_FACETS",i,len(q),
                json.dumps(Counter(v.get("status") for v in answers.values())),flush=True)
    effort=[]
    for isle,geo in sorted(geos.items()):
        d=answers.get(("EFFORT",isle),{"status":geo["status"]})
        effort.append({"island_code":isle,
            "gbif_all_plants_index_status":d["status"],
            "gbif_all_plants_first_record_year":d.get("earliest_indexed_year",""),
            "gbif_all_plants_pre1900_count":d.get("pre1900_records",""),
            "gbif_all_plants_pre1950_count":d.get("pre1950_records",""),
            "gbif_all_plants_pre2001_count":d.get("pre2001_records",""),
            "gbif_all_plants_total_records":d.get("gbif_record_count",""),
            "geometry_fidelity":geo.get("geometry_fidelity"),
            "error":d.get("error","")})
    outcsv(outdir/"island_plant_effort_year_facets.csv",effort)
    rows=[]
    for row in pairs:
        isle,pid=row["island_code"],row["plant_id"]
        taxon=matched[pid]
        geo=geos[isle]
        if taxon["status"]!="EXACT_SPECIES_MATCH":response={"status":taxon["status"]}
        elif geo["status"]!="GEOMETRY_READY":response={"status":geo["status"]}
        else:response=answers.get(("PLANT",pid+"|"+isle),{"status":"MISSING_RESPONSE"})
        rows.append({"plant_id":pid,"wcvp_name":names[pid],
                     "gbif_species_key":taxon.get("key",""),
                     "gbif_taxon_match_status":taxon["status"],
                     "island_code":isle,"level1":row["level1"],
                     "q4_pilot_status":response["status"],
                     "earliest_gbif_indexed_year":response.get("earliest_indexed_year",""),
                     "latest_gbif_indexed_year":response.get("latest_indexed_year",""),
                     "gbif_plant_record_count":response.get("gbif_record_count",""),
                     "gbif_plant_pre1900_count":response.get("pre1900_records",""),
                     "gbif_plant_pre1950_count":response.get("pre1950_records",""),
                     "gbif_plant_pre2001_count":response.get("pre2001_records",""),
                     "distinct_indexed_years":response.get("distinct_recorded_years",""),
                     "year_facet_complete":response.get("year_facet_complete_within_1700_2026",""),
                     "island_gbif_all_plants_total":answers.get(("EFFORT",isle),{}).get("gbif_record_count",""),
                     "island_gbif_all_plants_pre1900":answers.get(("EFFORT",isle),{}).get("pre1900_records",""),
                     "geometry_fidelity":geo.get("geometry_fidelity"),
                     "error":response.get("error","")})
    outcsv(outdir/"gbif_first_year_pilot_pairs.csv",rows)
    found=[r for r in rows if r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND"]
    valid=[int(r["earliest_gbif_indexed_year"]) for r in found]
    counts=Counter(r["q4_pilot_status"] for r in rows)
    date_bins={
     "1700_1849":sum(1700<=y<1850 for y in valid),
     "1850_1899":sum(1850<=y<1900 for y in valid),
     "1900_1949":sum(1900<=y<1950 for y in valid),
     "1950_1999":sum(1950<=y<2000 for y in valid),
     "2000_2010":sum(2000<=y<2011 for y in valid),
     "2011_2026":sum(2011<=y<=2026 for y in valid)
    }
    by_level={}
    for level in sorted({r["level1"] for r in rows}):
        vv=[r for r in rows if r["level1"]==level]
        by_level[level]={"pilot_pairs":len(vv),
            "dated":sum(r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND" for r in vv),
            "exact_taxon_and_queryable_geometry":sum(
                r["gbif_taxon_match_status"]=="EXACT_SPECIES_MATCH" and r["geometry_fidelity"] is not None for r in vv)}
    summary={
      "schema":"chocho_q4_gbif_earliest_plant_observation_pilot_v0.1",
      "status":"GBIF_INDEXED_FIRST_YEAR_COVERAGE_NOT_BIOLOGICAL_INTRODUCTION",
      "source_inputs":{"fixed_host_island_pairs":len(pairs),
                       "distinct_wcvp_plants_in_pilot":len(distinct_plants),
                       "source_fixed_pair_csv_sha256":sha(a.pairs),
                       "source_wcvp_name_map_csv_sha256":sha(a.wcvp_names),
                       "wgsrpd_geojson_source_pinned":"52da7828aba9d461dd133c27b3bd7a4407161f54"},
      "gbif_index":{"api":API,"searched_years":YEARS,
                    "retrieved_on_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    "query":"GBIF occurrence/search taxonKey+WKT, valid coordinates, no geospatial issues, present, year facets",
                    "max_polygon_wkt_bytes":MAX_WKT},
      "taxon_match_status":dict(counts if False else Counter(r["status"] for r in taxonrows)),
      "geometry_status":dict(Counter(v["status"] for v in geos.values())),
      "pilot_status_counts":dict(counts),
      "pilot_dated_pairs":len(found),
      "pilot_fraction_all_pairs_with_observed_plant_year":len(found)/len(rows),
      "pilot_fraction_queried_pairs_with_observed_plant_year":len(found)/sum(
         r["q4_pilot_status"] in ("DATED_OCCURRENCE_FOUND",STATUS_NO_EVIDENCE) for r in rows)
         if any(r["q4_pilot_status"] in ("DATED_OCCURRENCE_FOUND",STATUS_NO_EVIDENCE) for r in rows) else None,
      "earliest_indexed_year_summary":desc(valid),
      "earliest_indexed_year_count_bands":date_bins,
      "islands_with_any_botanical_effort_year":sum(r["gbif_all_plants_index_status"]=="DATED_OCCURRENCE_FOUND" for r in effort),
      "islands_with_any_pre1900_botanical_record":sum(int(r["gbif_all_plants_pre1900_count"] or 0)>0 for r in effort),
      "by_level1":by_level,
      "caution":[
        "GBIF first recorded plant year is a left-censored INDEX DATE, not the historical introduction or establishment date.",
        "Plants can arrive without GBIF presence, and old specimen dates may be unrelated to the introduction date.",
        "The same island's all-plant record count has sampling composition, digitization, provider and date-varying bias.",
        "Some island geometry is not safely representable in short WKT; these pairs are marked unqueried, not absent.",
        "The pilot is a response-blind 500-pair sample, not the entire 5,545 pair population.",
        "Taxon mapping statuses identify matches and ambiguity, and unresolvable names are not taken as zero.",
        "No use of island butterflies, colonial history or port counts in pilot sampling.",
        "No biological estimate of consumer lag or half-arrival time until observation effort, native butterfly origins, residence and temporal risk sets are evaluated."
      ]
    }
    (outdir/"gbif_first_year_pilot_audit.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("Q4_PILOT_RESULT",json.dumps(summary,ensure_ascii=False),flush=True)

def main():
    global a
    ap=argparse.ArgumentParser()
    for name in ("pairs","wcvp-names","wgsrpd-geojson","out-dir"):
        ap.add_argument("--"+name,required=True,type=Path)
    ap.add_argument("--max-workers",type=int,default=6)
    a=ap.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    pairs=readcsv(a.pairs)
    name_rows=readcsv(a.wcvp_names)
    names={r["plant_id"]:r["wcvp_accepted_name"] for r in name_rows}
    if len(pairs)!=500 or len(name_rows)!=560 or len(names)!=560:
        raise RuntimeError("Not the prespecified 500-pair Q4 pilot or 560 WCVP accepted names")
    if len({(r["plant_id"],r["island_code"]) for r in pairs})!=500:
        raise RuntimeError("Duplicate source host-island pair")
    if any(r["plant_id"] not in names for r in pairs):
        raise RuntimeError("Some WCVP species missing names")
    if not 1<=a.max_workers<=8:raise RuntimeError("Too much GBIF request concurrency")
    geos=build_geometries(a.wgsrpd_geojson,{r["island_code"] for r in pairs})
    print("GBIF_WKT_GEO",json.dumps(Counter(v["status"] for v in geos.values())),flush=True)
    runner(pairs,names,geos,a.max_workers,a.out_dir)

if __name__=="__main__":
    main()
