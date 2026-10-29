#!/usr/bin/env python3
"""Bounded, geography-selected true-island GBIF recording panel.

No colonisation/absence claims and no host/port effect fitting. Exact dissolved
GSHHG source geometries and original API bytes remain in the local cache.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ENDPOINT = "https://api.gbif.org/v1/occurrence/search/predicate"
CHECKLIST = "7ddf754f-d193-4cc9-b351-99906754a03b"
FIELDS = ("SPECIES_KEY", "BASIS_OF_RECORD", "DATASET_KEY", "YEAR")
WINDOWS = ((2010, 2017), (2018, 2025))
ROOT = Path(__file__).resolve().parents[1]
NETWORK_LOCK = Lock()
LAST_NETWORK_START = 0.0


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def haversine(a, b):
    x1,y1,x2,y2 = map(math.radians, (float(a["representative_lon"]),float(a["representative_lat"]),float(b["representative_lon"]),float(b["representative_lat"])))
    h = math.sin((y2-y1)/2)**2+math.cos(y1)*math.cos(y2)*math.sin((x2-x1)/2)**2
    return 2*6371.0088*math.asin(math.sqrt(min(1,max(0,h))))


def stratum(row):
    area = float(row["area_km2"])
    if int(row["antarctica_could_be_nearer"]) or not -55 <= float(row["representative_lat"]) <= 55:
        return None
    if not 20 <= area < 20000:
        return None
    area_bin = "20_to_100" if area < 100 else "100_to_1000" if area < 1000 else "1000_to_20000"
    lo,hi = float(row["mainland_coast_gap_lower_km"]),float(row["mainland_coast_gap_upper_km"])
    gap = "under100" if hi < 100 else "atleast500" if lo >= 500 else None
    plo,phi = float(row["nearest_wpi_coast_distance_lower_km"]),float(row["nearest_wpi_coast_distance_upper_km"])
    port = "within50" if phi <= 50 else "beyond100" if plo > 100 else None
    if not gap or not port:
        return None
    return ":".join((area_bin,gap,port))


def select_panel(rows):
    groups = {}
    seen = set()
    for row in rows:
        if row["island_key"] in seen:
            raise ValueError("Duplicate source island identity")
        seen.add(row["island_key"])
        key = stratum(row)
        if key:
            groups.setdefault(key,[]).append(row)
    selected = []
    for key, candidates in sorted(groups.items()):
        order = sorted(candidates,key=lambda x:(-float(x["area_km2"]),x["island_key"]))
        first = order[0]
        chosen = [first]
        if len(order)>1:
            chosen.append(sorted(order[1:],key=lambda x:(-haversine(first,x),x["island_key"]))[0])
        for rank,row in enumerate(chosen,1):
            selected.append({**row,"selection_stratum":key,"selection_rank":rank,"stratum_candidates":len(order)})
    return selected


def collect_facets(fetch_page, fields=FIELDS, page_size=100, max_pages=100):
    """Fetch every positive facet bin, checking totals and the first page twice."""
    bins = {field:{} for field in fields}
    done = set()
    first = None
    total = None
    for number in range(max_pages):
        page = fetch_page(number*page_size)
        if type(page.get("count")) is not int or page["count"]<0 or page.get("results") != []:
            raise ValueError("Invalid count-only response")
        if first is None:
            first = page
            total = page["count"]
        if page["count"] != total:
            raise ValueError("Live GBIF total changed across facet pages")
        facets = page.get("facets")
        if not isinstance(facets,list):
            raise ValueError("Missing facets array")
        values = {}
        for facet in facets:
            field = facet.get("field")
            if field in values:
                raise ValueError("Duplicate facet field")
            values[field] = facet.get("counts")
        for field in fields:
            if field in done:
                continue
            entries = values.get(field,[])
            if not isinstance(entries,list) or len(entries)>page_size:
                raise ValueError("Malformed facet page")
            for entry in entries:
                name,count = entry.get("name"),entry.get("count")
                if not isinstance(name,str) or not name or type(count) is not int or count<=0:
                    raise ValueError("Invalid positive facet bin")
                if name in bins[field]:
                    raise ValueError("Duplicate facet identifier across pages")
                bins[field][name] = count
            if len(entries)<page_size:
                done.add(field)
        if len(done)==len(fields):
            break
    else:
        raise ValueError("Facet pagination exceeded bounded page budget")
    again = fetch_page(0,recheck=True)
    def signature(p):
        return {"count":p.get("count"),"results":p.get("results"),"facets":sorted(p.get("facets",[]),key=lambda x:x["field"])}
    if signature(again) != signature(first):
        raise ValueError("Live GBIF first page changed during acquisition")
    residual = {field:total-sum(counts.values()) for field,counts in bins.items()}
    if any(n<0 for n in residual.values()):
        raise ValueError("Facet margin exceeds total")
    if any(residual.get(field,0) for field in ("DATASET_KEY","YEAR")):
        raise ValueError("Incomplete dataset/year facet margin")
    return {"count":total,"facets":bins,"unassigned_facet_records":residual,
            "facet_pages":number+1,"first_page_recheck_passed":True}


def query_safely(work):
    try:
        return {"status":"COMPLETE",**work()}
    except Exception as exc:
        return {"status":"MISSING_QUERY_FAILED","count":None,"error_type":type(exc).__name__,"error":str(exc)[:350]}


def payload(wkt, window, offset=0, page_size=100):
    start,end = window
    filters = [
        {"type":"equals","key":"TAXON_KEY","value":"5G9","checklistKey":CHECKLIST},
        {"type":"greaterThanOrEquals","key":"YEAR","value":str(start)},
        {"type":"lessThanOrEquals","key":"YEAR","value":str(end)},
        {"type":"equals","key":"HAS_COORDINATE","value":"true"},
        {"type":"equals","key":"HAS_GEOSPATIAL_ISSUE","value":"false"},
        {"type":"equals","key":"OCCURRENCE_STATUS","value":"PRESENT"},
        {"type":"within","geometry":wkt},
    ]
    return {"predicate":{"type":"and","predicates":filters},
            "parameters":{"CHECKLIST_KEY":[CHECKLIST]},"limit":0,"offset":0,
            "facets":list(FIELDS),"facetLimit":page_size,"facetOffset":offset,
            "facetMinCount":1,"facetMultiSelect":False}


def acquire_one(row, window, wkt, cache):
    win = f"{window[0]}_{window[1]}"
    folder = cache/row["island_key"]/win
    folder.mkdir(parents=True,exist_ok=True)
    receipts = []
    started = utc()
    def fetch_page(offset,recheck=False):
        request_bytes = canonical(payload(wkt,window,offset))
        prefix = f"page_{offset:05d}"+("_recheck" if recheck else "")
        (folder/(prefix+"_request.json")).write_bytes(request_bytes)
        bodypath = folder/(prefix+"_response.json")
        metapath = folder/(prefix+"_receipt.json")
        # Resume from original successful bytes, including the actual first-page recheck.
        if bodypath.exists() and metapath.exists():
            meta = json.loads(metapath.read_bytes())
            body = bodypath.read_bytes()
            if meta["request_sha256"] != sha(request_bytes) or meta["response_sha256"] != sha(body):
                raise ValueError("Cached source request/response hash mismatch")
        else:
            for attempt in range(2):
                try:
                    req = Request(ENDPOINT,data=request_bytes,headers={"Content-Type":"application/json","User-Agent":"chocho-island-source-pilot/1.0 (bounded research)"})
                    # Bound aggregate request starts even when three islands are active.
                    global LAST_NETWORK_START
                    with NETWORK_LOCK:
                        pause = max(0,1.25-(time.monotonic()-LAST_NETWORK_START))
                        if pause:
                            time.sleep(pause)
                        LAST_NETWORK_START = time.monotonic()
                    with urlopen(req,timeout=40) as response:
                        body = response.read(5_000_001)
                    if len(body)>5_000_000:
                        raise ValueError("Count-only source response exceeds bounded size")
                    break
                except HTTPError as exc:
                    (folder/(prefix+f"_error_{attempt}.txt")).write_bytes(exc.read(100000))
                    if attempt or exc.code not in (429,500,502,503,504):
                        raise RuntimeError(f"GBIF HTTP {exc.code}; exact geometry retained; no rectangle fallback") from None
                    time.sleep(2)
                except (URLError,TimeoutError,OSError):
                    if attempt:
                        raise RuntimeError("GBIF read failed twice; no count inferred") from None
                    time.sleep(1)
            meta = {"request_sha256":sha(request_bytes),"response_sha256":sha(body),"retrieved_utc":utc(),"facet_offset":offset,"is_first_page_recheck":recheck}
            bodypath.write_bytes(body)
            metapath.write_bytes(canonical(meta))
        receipts.append(meta)
        return json.loads(body)
    result = query_safely(lambda:collect_facets(fetch_page))
    return {"island_key":row["island_key"],"selection_stratum":row["selection_stratum"],"window":win,
            "started_utc":started,"finished_utc":utc(),"geometry_sha256":sha(wkt.encode()),
            "source_pages":receipts,**result}


def write_csv(path, rows, fields):
    with path.open("w",newline="",encoding="utf8") as handle:
        writer = csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def coverage_summary(results, selected):
    complete = [r for r in results if r["status"]=="COMPLETE"]
    by_island = {row["island_key"]:[r for r in complete if r["island_key"]==row["island_key"]] for row in selected}
    return {"requested_islands":len(selected),"requested_island_windows":len(selected)*2,
        "complete_island_windows":len(complete),"failed_island_windows":len(results)-len(complete),
        "complete_zero_record_island_windows":sum(r["count"]==0 for r in complete),
        "complete_positive_record_island_windows":sum(r["count"]>0 for r in complete),
        "complete_both_windows_islands":sum(len(v)==2 for v in by_island.values()),
        "islands_atleast3_species_bins_each_window":sum(len(v)==2 and all(len(r["facets"]["SPECIES_KEY"])>=3 for r in v) for v in by_island.values()),
        "islands_zero_records_both_windows":sum(len(v)==2 and all(r["count"]==0 for r in v) for v in by_island.values()),
        "islands_with_atleast_one_complete_zero_record_window":sum(any(r["count"]==0 for r in v) for v in by_island.values()),
        "complete_windows_without_species_bins":sum(not r["facets"]["SPECIES_KEY"] for r in complete),
        "positive_records_without_any_species_bin_windows":sum(r["count"]>0 and not r["facets"]["SPECIES_KEY"] for r in complete),
        "records_without_species_bin":sum(r["unassigned_facet_records"]["SPECIES_KEY"] for r in complete),
        "records_by_window":{f"{a}_{b}":sum(r["count"] for r in complete if r["window"]==f"{a}_{b}") for a,b in WINDOWS}}


def run(source_cache, cache, output, workers=3, lock_only=False):
    from shapely import orient_polygons, to_wkt, get_num_coordinates
    from measure_island_port_geography_20261029 import source_landmasses,derived_key
    cache.mkdir(parents=True,exist_ok=True)
    output.mkdir(parents=True,exist_ok=True)
    geography_path = output/"island_geography_exposures.csv"
    geography_manifest = output/"island_geography_manifest.json"
    geo = json.loads(geography_manifest.read_bytes())
    if geo["output_sha256"][geography_path.name] != sha(geography_path.read_bytes()):
        raise ValueError("Geographic exposure source differs from manifest")
    with geography_path.open() as f:
        rows = list(csv.DictReader(f))
    land = source_landmasses(source_cache)
    geometries = {derived_key(r.parent_source_ID):r.geometry for _,r in land.iterrows()}
    for row in rows:
        point = geometries[row["island_key"]].representative_point()
        row.update(representative_lon=point.x,representative_lat=point.y)
    chosen = select_panel(rows)
    if not chosen or len(chosen)>24:
        raise ValueError("Invalid bounded physical-island panel size")
    public = []
    wkts = {}
    for row in chosen:
        geom = orient_polygons(geometries[row["island_key"]],exterior_cw=False)
        if not geom.is_valid or geom.is_empty or geom.geom_type not in ("Polygon","MultiPolygon"):
            raise ValueError("Exact original island geometry is invalid; no substitution")
        wkt = to_wkt(geom,rounding_precision=-1,trim=True)
        wkts[row["island_key"]] = wkt
        public.append({**{k:v for k,v in row.items() if not k.startswith("representative_")},
                       "geometry_sha256":sha(wkt.encode()),"geometry_points":int(get_num_coordinates(geom)),
                       "orientation":"CCW exteriors, CW holes; WGS84 lon lat; no simplification"})
    support = []
    for a in ("20_to_100","100_to_1000","1000_to_20000"):
        for gap in ("under100","atleast500"):
            for port in ("within50","beyond100"):
                key = ":".join((a,gap,port))
                support.append({"stratum":key,"candidate_islands":sum(stratum(r)==key for r in rows),"selected_islands":sum(r["selection_stratum"]==key for r in chosen)})
    selection = {"schema":"chocho_geography_only_island_butterfly_panel_selection_20261029_v1",
        "geography_manifest_sha256":sha(geography_manifest.read_bytes()),
        "geography_exposures_sha256":sha(geography_path.read_bytes()),"gshhg_source":geo["gshhg"],
        "wpi_source":geo["ports"],"selection_rule":"12 area x coast-gap x WPI-proximity strata, <=2 each: greatest area then geographically farthest representative point from first; hashed key breaks ties; no outcomes queried before locking; no outcome-based replacements",
        "latitude_filter":"representative latitude within [-55,55] and Antarctica comparison flag false",
        "scope":"Recorded butterfly observations, no absences, colonisation or establishment effects",
        "support":support,"selected":public}
    lockpath = cache/"selection_lock.json"
    if lockpath.exists():
        prior = json.loads(lockpath.read_bytes())
        if prior["selection"] != selection:
            raise ValueError("Geography-only selection changed after lock; use independent new acquisition")
        lock = prior
    else:
        lock = {"locked_before_first_query_utc":utc(),"selection_sha256":sha(canonical(selection)),"selection":selection}
        lockpath.write_bytes(canonical(lock))
    (output/"island_butterfly_pilot_selection.json").write_text(json.dumps(lock,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"stage":"geography_selection_locked","selected":len(chosen),"strata":support,"selection_sha256":lock["selection_sha256"]}),flush=True)
    if lock_only:
        return
    tasks = [(row,win) for row in chosen for win in WINDOWS]
    results = []
    with ThreadPoolExecutor(max_workers=min(3,workers)) as pool:
        futures = {pool.submit(acquire_one,row,win,wkts[row["island_key"]],cache):(row,win) for row,win in tasks}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            (cache/(result["island_key"]+"_"+result["window"]+"_complete_result.json")).write_bytes(canonical(result))
            print(json.dumps({"completed":len(results),"requested":len(tasks),"island_key":result["island_key"],"window":result["window"],"status":result["status"],"count":result["count"],"error":result.get("error")}),flush=True)
    results.sort(key=lambda r:(r["island_key"],r["window"]))
    counts = []
    facet_rows = {field:[] for field in FIELDS}
    for result in results:
        n = result["count"]
        residual = result.get("unassigned_facet_records",{})
        counts.append({"island_key":result["island_key"],"stratum":result["selection_stratum"],"window":result["window"],
            "status":result["status"],"record_count":n,"recorded_species_bins":len(result["facets"]["SPECIES_KEY"]) if n is not None else None,
            "records_without_species_bin":residual.get("SPECIES_KEY"),"records_without_basis_bin":residual.get("BASIS_OF_RECORD"),
            "dataset_bins":len(result["facets"]["DATASET_KEY"]) if n is not None else None,
            "positive_year_bins":len(result["facets"]["YEAR"]) if n is not None else None,"error":result.get("error","")})
        for field in FIELDS:
            for name,count in sorted(result.get("facets",{}).get(field,{}).items()):
                facet_rows[field].append({"island_key":result["island_key"],"window":result["window"],"facet_key":name,"record_count":count})
    write_csv(output/"island_butterfly_pilot_counts.csv",counts,list(counts[0]))
    for field,items in facet_rows.items():
        write_csv(output/("island_butterfly_pilot_"+field.lower()+"_facets.csv"),items,["island_key","window","facet_key","record_count"])
    summary = coverage_summary(results,chosen)
    summary["cached_http_error_responses"] = len(list(cache.rglob("page_*_error_*.txt")))
    query_contract = payload("<exact source geometry, held in cache only>",WINDOWS[0])
    receipt = {"schema":"chocho_true_island_Papilionoidea_recording_pilot_20261029_v1",
        "finished_utc":utc(),"selection_sha256":lock["selection_sha256"],"selection_locked_utc":lock["locked_before_first_query_utc"],
        "taxonomy":{"checklistKey":CHECKLIST,"taxonKey":"5G9","name":"Papilionoidea","rank":"SUPERFAMILY","original_239_numeric_keys_not_combined":True},
        "endpoint":ENDPOINT,"request_contract_geometry_redacted":query_contract,"windows":WINDOWS,
        "interpretation":"Current retrospective GBIF recording counts; zero means no qualifying indexed records, never absence or failed arrival. Current taxonomy and post-event digitisation are not an as-of-2017 snapshot. Static WPI 2025 ports are not pre-outcome transport exposure. No host/port effect fitted.",
        "completeness":"All positive species/basis/dataset/year facet bins paged until exhausted; constant totals; duplicate checks; dataset/year sums equal total; species/basis residuals explicit; original first page reacquired and compared. Query failures remain missing.",
        "summary":summary,"results":results}
    (output/"island_butterfly_pilot_manifest.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True),flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cache",type=Path,default=Path("/tmp/chocho-island-sources"))
    parser.add_argument("--cache-dir",type=Path,default=Path("/tmp/chocho-island-gbif-cache"))
    parser.add_argument("--output-dir",type=Path,default=ROOT/"data/exploratory/island_port_20261029")
    parser.add_argument("--workers",type=int,choices=(1,2,3),default=3)
    parser.add_argument("--lock-only",action="store_true")
    args = parser.parse_args()
    run(args.source_cache,args.cache_dir,args.output_dir,args.workers,args.lock_only)
