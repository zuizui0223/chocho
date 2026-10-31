#!/usr/bin/env python3
"""Reconstruct *outcome-unselected* butterfly × island × period GBIF index risk.

This independent public source replaces neither the private author-original
species×island×year Parquet nor the GEB manuscript. First GBIF index != true
arrival, source country and larval-host availability. GBIF nonindex != absence.

225 butterfly candidate species and 17 islands fixed before subsequent outcomes,
GSHHG 2.3.7 polygons pin original landmass membership, 2001–2010 baseline,
2011–2017 first posterior index and 2018–2022 next index. A Papilionoidea
superfamily filter is used to avoid truncation of the all-taxa species facet.

The current GBIF index is not a historic data download: source queries are
persisted, and counts can drift when GBIF is reindexed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
SOURCES={
 "butterflies":("provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv",
                "39fe84981cab2b0e47063399f0a5c0c9731e800d"),
 "islands":("provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json",
            "54c4515947d2a3ee7545e4007985ef82a2af6868"),
}
GSHHG_SHA256="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
GSHHG_SOURCE="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GBIF="https://api.gbif.org/v1"
UA="chocho-butterfly-225x17-temporal-risk-index/20261031"
PERIODS={"pre2001_2010":(2001,2010),
         "first2011_2017":(2011,2017),
         "holdout2018_2022":(2018,2022)}

def read_source(name):
    p,expected=SOURCES[name]
    data=(ROOT/p).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if actual!=expected:raise ValueError("Git source changed: "+p+" "+actual)
    if p.endswith(".csv"):
        return list(csv.DictReader(io.StringIO(data.decode("utf8"))))
    return json.loads(data)

def get_json(endpoint,params=None):
    address=endpoint+("?"+urlencode(params,doseq=True) if params else "")
    last=None
    for attempt in range(4):
        try:
            req=Request(address,headers={"Accept":"application/json","User-Agent":UA})
            with urlopen(req,timeout=65) as f:
                if f.status!=200:raise ValueError("GBIF HTTP status "+str(f.status))
                result=json.load(f)
            if not isinstance(result,dict) or result.get("error"):
                raise ValueError("Unusable GBIF response: "+str(result.get("error"))[:200])
            return result
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as ex:
            last=ex
            if isinstance(ex,HTTPError) and ex.code not in (408,425,429,500,502,503,504):
                break
            if attempt<3:time.sleep(min(8,0.7*2**attempt))
    raise RuntimeError(f"GBIF request failed endpoint={endpoint}, error={type(last).__name__}: {str(last)[:160]}")

def resolve_superfamily():
    match=get_json(GBIF+"/species/match",{
        "name":"Papilionoidea","rank":"SUPERFAMILY","kingdom":"Animalia","strict":"true"})
    name=str(match.get("canonicalName") or match.get("scientificName") or "")
    rank=str(match.get("rank","")).upper()
    matchtype=str(match.get("matchType","")).upper()
    key=match.get("acceptedUsageKey") or match.get("usageKey")
    if name.lower()!="papilionoidea" or rank!="SUPERFAMILY" or matchtype!="EXACT" or not key:
        raise ValueError("GBIF Papilionoidea resolution ambiguous: "+json.dumps(match)[:500])
    check=get_json(GBIF+"/species/"+str(key))
    if str(check.get("rank","")).upper()!="SUPERFAMILY":
        raise ValueError("Resolved GBIF Papilionoidea taxon rank changed: "+repr(check))
    return {"superfamily_name":"Papilionoidea","GBIF_taxonKey":int(key),
            "source_match":{k:match.get(k) for k in ("matchType","rank","canonicalName","status","confidence")},
            "resolved_taxon_name":check.get("canonicalName")}

def load_polygons(zipfile,working,ids,max_chars):
    import geopandas as gpd
    import shapely
    from shapely.geometry import Polygon,MultiPolygon
    from shapely.ops import transform
    from shapely import wkt
    from pyproj import Transformer
    raw=zipfile.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=GSHHG_SHA256:raise ValueError("Original GSHHG SHA256 mismatch "+digest)
    working.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=[p for p in z.namelist() if p.endswith((".shp",".shx",".dbf",".prj")) and
                 "/h/GSHHS_h_L1." in "/"+p]
        if len(members)<3:raise ValueError("GSHHG original shapefile missing")
        for member in members:z.extract(member,working)
    shp=next(working.rglob("GSHHS_h_L1.shp"))
    land=gpd.read_file(shp)
    land["landmass_id"]=land.sibling_id.astype(int).astype(str)
    if not ids <= set(land.landmass_id):raise ValueError("missing one or more GSHHG original sibling IDs")
    projection=Transformer.from_crs(4326,6933,always_xy=True)
    polygons={}
    audit=[]
    for island in sorted(ids,key=int):
        x=land.loc[land.landmass_id.eq(island),"geometry"]
        original=shapely.union_all(x.to_numpy())
        if original.is_empty or not original.is_valid:
            raise ValueError("Invalid source GSHHG island polygon "+island)
        original_area=transform(projection.transform,original).area
        if isinstance(original,MultiPolygon):
            parts=list(original.geoms)
            poly=max(parts,key=lambda y:transform(projection.transform,y).area)
        elif isinstance(original,Polygon):poly=original
        else:raise ValueError("Unexpected source island geometry "+type(original).__name__)
        largest_area=transform(projection.transform,poly).area
        if largest_area/max(original_area,1e-10)<.985:
            raise ValueError("MULTIPOLYGON island loses >1.5% with largest component "+island)
        chosen=None; selected_tolerance=None
        for tol in (0.,.00025,.0005,.001,.002,.003,.005,.008,.012,.02,.035,.06,.10):
            simplified=poly.simplify(tol,preserve_topology=True) if tol else poly
            if simplified.is_empty or not simplified.is_valid:continue
            text=wkt.dumps(simplified,rounding_precision=5,trim=True)
            revised_area=transform(projection.transform,simplified).area
            retained=revised_area/max(original_area,1e-10)
            if len(text)<=max_chars and .96<=retained<=1.04:
                chosen=text;selected_tolerance=tol;break
        if not chosen:
            audit.append({"island_id":island,"status":"NO_VALID_AREA_PRESERVING_SHORT_WKT",
                          "original_area_km2":original_area/1e6,"polygon_wkt_length":None})
            raise ValueError(f"Cannot build <= {max_chars}-character area-preserving GBIF WKT for island {island}")
        polygons[island]=chosen
        audit.append({"island_id":island,"status":"VALID_POLYGON",
                      "original_area_km2":original_area/1e6,
                      "largest_component_frac":largest_area/original_area,
                      "simplify_degrees":selected_tolerance,
                      "simplified_area_retained":transform(projection.transform,
                         wkt.loads(chosen)).area/max(original_area,1e-10),
                      "polygon_wkt_length":len(chosen),
                      "polygon_sha256":hashlib.sha256(chosen.encode()).hexdigest()})
    return polygons,audit,digest

def indexed_butterflies_in_period(island,period,geometry,taxon):
    years=PERIODS[period]
    basic={
        "taxonKey":str(taxon),"geometry":geometry,"year":f"{years[0]},{years[1]}",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT","facet":"speciesKey","facetLimit":1000,
        "limit":0,
    }
    species={}
    total=None;offset=0;pages=0
    for _ in range(40):
        result=get_json(GBIF+"/occurrence/search",{**basic,"facetOffset":offset})
        n=int(result["count"])
        if total is not None and n!=total:
            raise ValueError("GBIF indexing count changed within facet pagination")
        total=n
        facets={str(x.get("field","")).upper():x.get("counts",[])
            for x in result.get("facets",[])}
        if total>0 and "SPECIES_KEY" not in facets:
            raise ValueError("Nonzero butterfly index but GBIF SPECIES_KEY facet missing "+str(sorted(facets)))
        page=facets.get("SPECIES_KEY",[])
        for q in page:
            key=str(q.get("name",""))
            count=int(q.get("count",0))
            if not key.isdecimal() or count<0:
                raise ValueError("Invalid GBIF butterfly speciesKey facet "+repr(q))
            if key in species:raise ValueError("Overlapping GBIF species facets in pages")
            species[key]=count
        pages+=1
        if len(page)<1000:break
        offset+=len(page)
    else:raise ValueError("More than 40 butterfly taxon facet pages: "+island)
    if sum(species.values())>total:
        raise ValueError("Butterfly species indexed counts exceed all butterfly records")
    return species,{"island_id":island,"period":period,
      "all_superfamily_papilionoidea_indexed_occurrences":total,
      "gbif_specieskey_count":len(species),"species_facet_pages":pages,
      "taxonKey":str(taxon),"GBIF_indexed_not_biological_absence":True}

def build_risk(islands,table,facts):
    # This scope is fixed before any outcomes. Explicit zeros mean only
    # no indexed GBIF butterfly records under the same query, NOT occupancy.
    expected=225*len(islands)
    out=[]
    for i in sorted(islands,key=int):
        for sp in table:
            key=sp["taxon_key"];name=sp["species"]
            pre=int(facts.get((i,"pre2001_2010"),{}).get(key,0))
            f1=int(facts.get((i,"first2011_2017"),{}).get(key,0))
            f2=int(facts.get((i,"holdout2018_2022"),{}).get(key,0))
            out.append({
              "species":name,"speciesKey":key,"landmass_id":i,
              "pre2001_2010_indexed_records":pre,
              "future2011_2017_indexed_records":f1,
              "future2018_2022_indexed_records":f2,
              "eligible_2011_2017":int(pre==0),
              "first_indexed_in_2011_2017":int(pre==0 and f1>0),
              "eligible_2018_2022":int(pre==0 and f1==0),
              "first_indexed_in_2018_2022":int(pre==0 and f1==0 and f2>0),
              "source_not_future_selected":True,
              "index_not_true_absence_or_colonization":True
            })
    assert len(out)==expected
    assert len({(r["speciesKey"],r["landmass_id"]) for r in out})==expected
    assert sum(x["eligible_2018_2022"] for x in out)<=sum(x["eligible_2011_2017"] for x in out)
    return out

def save_csv(path,rows):
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gshhg-zip",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--islands",default="all",
                   help="all or comma-separated original GSHHG sibling IDs; subset for preliminary source pilot")
    p.add_argument("--wkt-max-chars",type=int,default=4300)
    args=p.parse_args()
    candidates=read_source("butterflies")
    island_source=read_source("islands")
    all_islands={str(x["island_id"]) for x in island_source["islands"]}
    assert len(candidates)==225 and len({x["taxon_key"] for x in candidates})==225
    assert len(all_islands)==17
    selected=all_islands if args.islands=="all" else set(args.islands.split(","))
    if not selected or not selected <= all_islands:
        raise ValueError("Selected islands not all in original 17 source")
    taxon=resolve_superfamily()
    polygons,geometry_receipt,gshhg=load_polygons(args.gshhg_zip,args.work,selected,args.wkt_max_chars)
    HERE.mkdir(parents=True,exist_ok=True)
    facts={}
    queries=[]
    for island in sorted(selected,key=int):
        for period in PERIODS:
            scores,receipt=indexed_butterflies_in_period(island,period,polygons[island],
                                                        taxon["GBIF_taxonKey"])
            facts[(island,period)]=scores
            queries.append(receipt)
            print("BUTTERFLY_PERIOD",island,period,receipt["all_superfamily_papilionoidea_indexed_occurrences"],
                "all_species_facets",len(scores),"candidate225",sum(k in scores for k in
                       (x["taxon_key"] for x in candidates)),flush=True)
            time.sleep(.12)
    outcome=build_risk(selected,candidates,facts)
    fname="full17" if len(selected)==17 else "pilot_"+str(len(selected))
    save_csv(HERE/(fname+"_butterfly_island_225_atrisk.csv"),outcome)
    summary={
       "schema":"chocho_public_GBIF_unselected_225_butterfly_x_original_island_window_v1",
       "status":"PUBLIC_REINDEXED_OCCURRENCE_HISTORY_NOT_TRUE_COLONIZATION_OR_HOST_MEDIATION",
       "snapshot_timestamp_utc":datetime.now(timezone.utc).isoformat(),
       "data_source":GBIF+"/occurrence/search",
       "source_git_blob_pins":{k:{"path":v[0],"sha":v[1]} for k,v in SOURCES.items()},
       "gshhg_shape_sha256":gshhg,"gshhg_source":GSHHG_SOURCE,
       "butterfly_taxon_filter":taxon,"islands_tested":sorted(selected,key=int),
       "species_candidates_frozen_before_future_225":len(candidates),
       "all_17_original_islands":len(all_islands),
       "geometry_audit":geometry_receipt,
       "GBIF_index_query_receipts":queries,
       "risk_cells":len(outcome),
       "eligible_2011_2017":sum(x["eligible_2011_2017"] for x in outcome),
       "first_indexed_2011_2017":sum(x["first_indexed_in_2011_2017"] for x in outcome),
       "eligible_2018_2022":sum(x["eligible_2018_2022"] for x in outcome),
       "first_indexed_2018_2022":sum(x["first_indexed_in_2018_2022"] for x in outcome),
       "limits":[
           "This GBIF 2026 index snapshot differs from the private original 2001–2024 precise coordinate source.",
           "225 candidate butterfly species are fixed before looking at later island first indexing; no 60-species outcome preselection.",
           "Simplified original GSHHG island polygons preserve >=96% equal-area size but can cause nearshore spatial errors.",
           "All 225 species × island cells are included; many are geographically infeasible, so all-species crude risk rates are not transport effects.",
           "Query negatives are GBIF nonindexes, not biological absences, occupancy, or proof of no colonization.",
           "The GBIF Papilionoidea superfamily extent, taxonomy and indexing may drift in later downloads.",
           "Three window-specific speciesKey facets lack exact first year/individual observation providers and are not true sampling histories.",
           "Current source has no plant × butterfly host-edge supply, transport vector or independent regional establishment observations.",
           "Do not reuse 2001–2024 pooled port plant inventory as a pre-outcome foodplant exposure."
       ]
    }
    (HERE/(fname+"_query_receipt.json")).write_text(
         json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:summary[k] for k in ("status","risk_cells","eligible_2011_2017","first_indexed_2011_2017",
            "eligible_2018_2022","first_indexed_2018_2022")},ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
