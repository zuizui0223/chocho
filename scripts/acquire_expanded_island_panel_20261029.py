#!/usr/bin/env python3
"""Expand a fixed true-island recording panel by geography alone.

Retains all 24 original pilot islands and their original source bytes. New
islands maximize spatial separation within the original 12 geographic strata.
No host, barrier or transport effect is estimated.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import acquire_island_butterfly_panel_20261029 as pilot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/"data/exploratory/island_port_20261029"
OUTPUT = ROOT/"data/exploratory/island_expansion_20261029"


def expand_selection(rows, original_selected, per_stratum=8):
    if not 2 <= per_stratum <= 8:
        raise ValueError("Expansion must retain original two and select at most eight")
    bykey = {r["island_key"]:r for r in rows}
    if len(bykey)!=len(rows):
        raise ValueError("Duplicate candidate island identity")
    groups = {}
    for r in rows:
        key = pilot.stratum(r)
        if key:
            groups.setdefault(key,[]).append(r)
    originals = {}
    for original in original_selected:
        key = original["selection_stratum"]
        if original["island_key"] not in bykey or pilot.stratum(bykey[original["island_key"]]) != key:
            raise ValueError("Original pilot island missing or ineligible; no replacement allowed")
        originals.setdefault(key,[]).append(original)
    if set(groups)!=set(originals):
        raise ValueError("Expansion strata differ from original pilot strata")
    result = []
    for key in sorted(groups):
        ordered = sorted(originals[key],key=lambda r:r["selection_rank"])
        if len(ordered)!=2 or [r["selection_rank"] for r in ordered]!=[1,2]:
            raise ValueError("Expected the two locked pilot seeds in each stratum")
        selected = [bykey[r["island_key"]] for r in ordered]
        selected_keys = {r["island_key"] for r in selected}
        if len(selected_keys)!=2:
            raise ValueError("Duplicated original pilot seeds")
        available = [r for r in groups[key] if r["island_key"] not in selected_keys]
        while available and len(selected)<per_stratum:
            pick = min(available,key=lambda r:(-min(pilot.haversine(r,s) for s in selected),r["island_key"]))
            selected.append(pick)
            available = [r for r in available if r["island_key"]!=pick["island_key"]]
        for rank,r in enumerate(selected,1):
            result.append({**r,"selection_stratum":key,"selection_rank":rank,
                "stratum_candidates":len(groups[key]),"provenance_origin":"ORIGINAL_PILOT_REUSED" if rank<=2 else "NEW_EXPANSION"})
    return result


def read_cached(folder, prefix, request_bytes):
    body = (folder/(prefix+"_response.json")).read_bytes()
    saved_request = (folder/(prefix+"_request.json")).read_bytes()
    receipt = json.loads((folder/(prefix+"_receipt.json")).read_bytes())
    if saved_request != request_bytes or receipt["request_sha256"] != pilot.sha(request_bytes) or receipt["response_sha256"] != pilot.sha(body):
        raise ValueError("Pinned cached request/response differs")
    return json.loads(body),receipt


def lock_selection(cache, selection):
    path = cache/"expansion_selection_lock.json"
    if path.exists():
        prior = json.loads(path.read_bytes())
        if prior["selection"] != selection or prior["selection_sha256"] != pilot.sha(pilot.canonical(selection)):
            raise ValueError("Locked geography selection cannot change on resume")
        return prior
    lock = {"locked_before_new_queries_utc":pilot.utc(),"selection_sha256":pilot.sha(pilot.canonical(selection)),"selection":selection}
    path.write_bytes(pilot.canonical(lock))
    return lock


def reused_query(row, window, wkt, old_result, old_cache, original_manifest_sha):
    win = f"{window[0]}_{window[1]}"
    folder = old_cache/row["island_key"]/win
    receipts = []
    pinned = {(r["facet_offset"],r["is_first_page_recheck"]):r for r in old_result["source_pages"]}
    def fetch(offset,recheck=False):
        prefix = f"page_{offset:05d}"+("_recheck" if recheck else "")
        result,receipt = read_cached(folder,prefix,pilot.canonical(pilot.payload(wkt,window,offset)))
        if receipt != pinned.get((offset,recheck)):
            raise ValueError("Original pilot receipt no longer matches locked original manifest")
        receipts.append(receipt)
        return result
    verified = pilot.query_safely(lambda:pilot.collect_facets(fetch))
    if verified["status"]=="COMPLETE" and any(verified[k]!=old_result[k] for k in ("count","facets","unassigned_facet_records")):
        raise ValueError("Reparsed original pilot differs from published source")
    return {"island_key":row["island_key"],"selection_stratum":row["selection_stratum"],"window":win,
        "provenance_origin":"ORIGINAL_PILOT_REUSED","original_manifest_sha256":original_manifest_sha,
        "geometry_sha256":pilot.sha(wkt.encode()),"source_pages":receipts,
        "started_utc":old_result["started_utc"],"finished_utc":old_result["finished_utc"],**verified}


class SourceClient:
    """Sequential public reads; at least 1.5 seconds between request starts."""
    def __init__(self):
        self.last_start = 0.0
        self.new_requests = 0

    def request(self, data, folder, prefix):
        for attempt in range(2):
            time.sleep(max(0,1.5-(time.monotonic()-self.last_start)))
            self.last_start = time.monotonic()
            self.new_requests += 1
            try:
                req = Request(pilot.ENDPOINT,data=data,headers={"Content-Type":"application/json","User-Agent":"chocho-geography-only-expansion/1.0 bounded research"})
                with urlopen(req,timeout=45) as response:
                    body = response.read(5_000_001)
                if len(body)>5_000_000:
                    raise ValueError("Source facet response exceeds bounded size")
                return body
            except HTTPError as exc:
                (folder/(prefix+f"_http{exc.code}_attempt{attempt}.txt")).write_bytes(exc.read(100000))
                if attempt or exc.code not in (429,500,502,503,504):
                    raise RuntimeError(f"GBIF HTTP {exc.code}; count missing, exact geometry not replaced") from None
                retry = exc.headers.get("Retry-After","")
                delay = max(60 if exc.code==429 else 5,float(retry) if retry.isdigit() else 0)
                print(json.dumps({"stage":"source_retry_cooldown","http_status":exc.code,"seconds":delay}),flush=True)
                deadline = time.monotonic()+delay
                while time.monotonic()<deadline:
                    time.sleep(min(30,max(0,deadline-time.monotonic())))
            except (URLError,TimeoutError,OSError):
                if attempt:
                    raise RuntimeError("GBIF public read failed twice; query remains missing") from None
                time.sleep(5)


def new_query(row, window, wkt, cache, client):
    win = f"{window[0]}_{window[1]}"
    folder = cache/row["island_key"]/win
    folder.mkdir(parents=True,exist_ok=True)
    receipts = []
    started = pilot.utc()
    def fetch(offset,recheck=False):
        prefix = f"page_{offset:05d}"+("_recheck" if recheck else "")
        request = pilot.canonical(pilot.payload(wkt,window,offset))
        if (folder/(prefix+"_response.json")).exists():
            result,receipt = read_cached(folder,prefix,request)
        else:
            (folder/(prefix+"_request.json")).write_bytes(request)
            body = client.request(request,folder,prefix)
            receipt = {"request_sha256":pilot.sha(request),"response_sha256":pilot.sha(body),
                "retrieved_utc":pilot.utc(),"facet_offset":offset,"is_first_page_recheck":recheck}
            (folder/(prefix+"_response.json")).write_bytes(body)
            (folder/(prefix+"_receipt.json")).write_bytes(pilot.canonical(receipt))
            result = json.loads(body)
        receipts.append(receipt)
        return result
    result = pilot.query_safely(lambda:pilot.collect_facets(fetch))
    return {"island_key":row["island_key"],"selection_stratum":row["selection_stratum"],"window":win,
        "provenance_origin":"NEW_EXPANSION","started_utc":started,"finished_utc":pilot.utc(),
        "geometry_sha256":pilot.sha(wkt.encode()),"source_pages":receipts,**result}


def publish(output, lock, results, client, final=False):
    chosen = lock["selection"]["selected"]
    results = sorted(results,key=lambda r:(r["island_key"],r["window"]))
    summary = pilot.coverage_summary(results,chosen)
    summary.update(pending_island_windows=len(chosen)*2-len(results),new_requests_this_invocation=client.new_requests,
                   original_pilot_reused_windows=sum(r["provenance_origin"]=="ORIGINAL_PILOT_REUSED" for r in results))
    receipt = {"schema":"chocho_geography_only_true_island_expansion_20261029_v1","status":"FINISHED" if final else "RUNNING",
        "updated_utc":pilot.utc(),"selection_sha256":lock["selection_sha256"],
        "selection_locked_before_new_queries_utc":lock["locked_before_new_queries_utc"],
        "taxonomy":{"checklistKey":pilot.CHECKLIST,"taxonKey":"5G9","name":"Papilionoidea","original_239_keys_not_combined":True},
        "windows":pilot.WINDOWS,"endpoint":pilot.ENDPOINT,"original_pilot":lock["selection"]["original_pilot"],
        "query_contract_geometry_redacted":pilot.payload("<exact original oriented WKT, cached only>",pilot.WINDOWS[0]),
        "interpretation":"Current retrospective recording data. Zero is no qualifying indexed record, never absence, first arrival or establishment failure. Original24 reused from prior known pilot; expansion72 geographically locked before NEW queries. No geography, port or host effect fitted. Static WPI2025 does not establish pre2018 port exposure.",
        "summary":summary,"results":results}
    (output/"expanded_island_panel_manifest.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    counts=[]
    facets={f:[] for f in pilot.FIELDS}
    for r in results:
        ok = r["status"]=="COMPLETE"
        counts.append({"island_key":r["island_key"],"stratum":r["selection_stratum"],"window":r["window"],
            "provenance_origin":r["provenance_origin"],"status":r["status"],"record_count":r["count"],
            "recorded_species_bins":len(r["facets"]["SPECIES_KEY"]) if ok else None,
            "records_without_species_bin":r.get("unassigned_facet_records",{}).get("SPECIES_KEY"),
            "records_without_basis_bin":r.get("unassigned_facet_records",{}).get("BASIS_OF_RECORD"),
            "dataset_bins":len(r["facets"]["DATASET_KEY"]) if ok else None,
            "positive_year_bins":len(r["facets"]["YEAR"]) if ok else None,"error":r.get("error","")})
        for f in pilot.FIELDS:
            for key,n in sorted(r.get("facets",{}).get(f,{}).items()):
                facets[f].append({"island_key":r["island_key"],"window":r["window"],"facet_key":key,"record_count":n})
    if counts:
        pilot.write_csv(output/"expanded_island_panel_counts.csv",counts,list(counts[0]))
    for f,rows in facets.items():
        pilot.write_csv(output/("expanded_island_panel_"+f.lower()+"_facets.csv"),rows,["island_key","window","facet_key","record_count"])
    return receipt


def run(source_cache, cache, original_cache, output, lock_only=False):
    from shapely import orient_polygons,to_wkt,get_num_coordinates
    from measure_island_port_geography_20261029 import source_landmasses,derived_key
    output.mkdir(parents=True,exist_ok=True)
    cache.mkdir(parents=True,exist_ok=True)
    old_lock_path=BASE/"island_butterfly_pilot_selection.json"
    old_manifest_path=BASE/"island_butterfly_pilot_manifest.json"
    original=json.loads(old_lock_path.read_bytes())
    old_manifest=json.loads(old_manifest_path.read_bytes())
    if original["selection_sha256"]!=pilot.sha(pilot.canonical(original["selection"])) or old_manifest["selection_sha256"]!=original["selection_sha256"]:
        raise ValueError("Original pilot selection provenance differs")
    exposure=BASE/"island_geography_exposures.csv"
    geo_manifest=BASE/"island_geography_manifest.json"
    if pilot.sha(exposure.read_bytes())!=original["selection"]["geography_exposures_sha256"] or pilot.sha(geo_manifest.read_bytes())!=original["selection"]["geography_manifest_sha256"]:
        raise ValueError("Geographic source changed since original pilot")
    with exposure.open() as f:
        rows=list(csv.DictReader(f))
    land=source_landmasses(source_cache)
    geometries={derived_key(r.parent_source_ID):r.geometry for _,r in land.iterrows()}
    for r in rows:
        p=geometries[r["island_key"]].representative_point()
        r.update(representative_lon=p.x,representative_lat=p.y)
    chosen=expand_selection(rows,original["selection"]["selected"])
    wkts={}
    public=[]
    old_geometry={r["island_key"]:r["geometry_sha256"] for r in original["selection"]["selected"]}
    for r in chosen:
        geom=orient_polygons(geometries[r["island_key"]],exterior_cw=False)
        if not geom.is_valid or geom.is_empty or geom.geom_type not in ("Polygon","MultiPolygon"):
            raise ValueError("Invalid exact source geometry; no substitution")
        wkt=to_wkt(geom,rounding_precision=-1,trim=True)
        digest=pilot.sha(wkt.encode())
        if r["island_key"] in old_geometry and old_geometry[r["island_key"]]!=digest:
            raise ValueError("Original pilot exact geometry bytes differ")
        wkts[r["island_key"]]=wkt
        public.append({**r,"representative_point_role":"Derived selection/clustering covariate only, never coastline-distance exposure",
            "geometry_sha256":digest,"geometry_points":int(get_num_coordinates(geom)),"orientation":"WGS84 lon lat, CCW exteriors/CW holes, exact unsimplified source"})
    support=[{"stratum":key,"candidate_islands":sum(pilot.stratum(r)==key for r in rows),"selected_islands":sum(r["selection_stratum"]==key for r in chosen),"original_retained":2} for key in sorted({r["selection_stratum"] for r in chosen})]
    selection={"schema":"chocho_expanded_island_geography_selection_20261029_v1",
        "original_pilot":{"selection_sha256":original["selection_sha256"],"selection_file_sha256":pilot.sha(old_lock_path.read_bytes()),
            "manifest_sha256":pilot.sha(old_manifest_path.read_bytes()),"original_locked_utc":original["locked_before_first_query_utc"],
            "original_pilot_islands":24,"original_pilot_outcomes_were_already_known":True},
        "geography_exposures_sha256":pilot.sha(exposure.read_bytes()),"geography_manifest_sha256":pilot.sha(geo_manifest.read_bytes()),
        "gshhg_source":original["selection"]["gshhg_source"],"wpi_source":original["selection"]["wpi_source"],
        "rule":"Retain original2 geography-selected pilot islands in each of12 originalstrata; add maximin representative-point haversine distance to chosen withinstratum; stable island_key ties; atmost8perstratum; geographyonly; no record-dependent replacements.",
        "support":support,"selected":public}
    lock=lock_selection(cache,selection)
    (output/"expanded_island_panel_selection.json").write_text(json.dumps(lock,indent=2,sort_keys=True)+"\n")
    pilot.write_csv(output/"expanded_island_panel_selection.csv",public,list(public[0]))
    pilot.write_csv(output/"expanded_island_panel_candidate_support.csv",support,list(support[0]))
    new_islands=sum(r["provenance_origin"]=="NEW_EXPANSION" for r in chosen)
    print(json.dumps({"stage":"EXPANSION_SELECTION_LOCKED","selected":len(chosen),"original_retained":len(chosen)-new_islands,
        "new_islands":new_islands,"new_island_windows":new_islands*2,"minimum_new_count_requests":new_islands*4,
        "minimum_spacing_runtime_seconds":max(0,new_islands*4-1)*1.5,"selection_sha256":lock["selection_sha256"],"support":support}),flush=True)
    if lock_only:
        return lock
    client=SourceClient()
    results=[]
    old_results={(r["island_key"],r["window"]):r for r in old_manifest["results"]}
    for row in chosen:
        for window in pilot.WINDOWS:
            win=f"{window[0]}_{window[1]}"
            if row["provenance_origin"]=="ORIGINAL_PILOT_REUSED":
                result=reused_query(row,window,wkts[row["island_key"]],old_results[(row["island_key"],win)],original_cache,selection["original_pilot"]["manifest_sha256"])
            else:
                result=new_query(row,window,wkts[row["island_key"]],cache,client)
            results.append(result)
            publish(output,lock,results,client)
            print(json.dumps({"completed":len(results),"total":len(chosen)*2,"island_key":row["island_key"],"window":win,
                "origin":result["provenance_origin"],"status":result["status"],"count":result["count"],"new_requests":client.new_requests,"error":result.get("error")}),flush=True)
    receipt=publish(output,lock,results,client,final=True)
    print(json.dumps(receipt["summary"],sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-cache",type=Path,default=Path("/tmp/chocho-island-sources"))
    p.add_argument("--cache-dir",type=Path,default=Path("/tmp/chocho-island-expansion-gbif-cache"))
    p.add_argument("--original-cache",type=Path,default=Path("/tmp/chocho-island-gbif-cache"))
    p.add_argument("--output-dir",type=Path,default=OUTPUT)
    p.add_argument("--lock-only",action="store_true")
    a=p.parse_args()
    run(a.source_cache,a.cache_dir,a.original_cache,a.output_dir,a.lock_only)
