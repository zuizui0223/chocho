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

def resolve_butterfly_taxon_filter():
    # GBIF's currently configured name matching returned NONE on strict
    # Papilionoidea SUPERFAMILY queries in the actual 2026-10-31 CI source
    # test. Use six established butterfly families instead, as OR filters;
    # accept a family only when the name itself is not an ambiguous homonym.
    match=get_json(GBIF+"/species/match",{
        "name":"Papilionoidea","rank":"SUPERFAMILY","kingdom":"Animalia","strict":"true"})
    name=str(match.get("canonicalName") or match.get("scientificName") or "")
    rank=str(match.get("rank","")).upper()
    typ=str(match.get("matchType","")).upper()
    key=match.get("acceptedUsageKey") or match.get("usageKey")
    if name.lower()=="papilionoidea" and rank=="SUPERFAMILY" and typ=="EXACT" and key:
        check=get_json(GBIF+"/species/"+str(key))
        if str(check.get("rank","")).upper()=="SUPERFAMILY":
            return {"method":"GBIF_SUPERFAMILY_TAXONKEY",
                "superfamily_name":"Papilionoidea","GBIF_taxonKey":str(key),
                "source_match":{k:match.get(k) for k in ("matchType","rank","canonicalName","status")}}
    families=["Nymphalidae","Pieridae","Lycaenidae","Papilionidae","Hesperiidae","Riodinidae"]
    by_family={}
    for family in families:
        d=get_json(GBIF+"/species/match",{"name":family,"rank":"FAMILY",
                    "kingdom":"Animalia","strict":"true"})
        rank=str(d.get("rank","")).upper()
        name=str(d.get("canonicalName") or d.get("scientificName") or "")
        key=d.get("acceptedUsageKey") or d.get("usageKey")
        if rank!="FAMILY" or name.lower().split(" ")[0]!=family.lower() or not key:
            # Some GBIF checklist deployments return NONE to a strict name
            # match. Searching candidates still requires exact family rank and
            # accepted Animalia family, never accepting arbitrary fuzzy hits.
            data=get_json(GBIF+"/species/search",{"q":family,"rank":"FAMILY","limit":40})
            candidates=[r for r in data.get("results",[])
              if str(r.get("rank","")).upper()=="FAMILY"
              and str(r.get("canonicalName") or r.get("scientificName") or "").lower().split(" ")[0]==family.lower()
              and str(r.get("status","")).upper()=="ACCEPTED"
              and str(r.get("kingdom","Animalia")).lower()=="animalia"]
            if len(candidates)!=1:
                raise ValueError("Cannot uniquely resolve butterfly family: "+family+
                                 ", matches="+str(len(candidates))+", strict="+json.dumps(d)[:220])
            key=candidates[0].get("key")
            name=str(candidates[0].get("canonicalName") or family)
        verify=get_json(GBIF+"/species/"+str(key))
        if (str(verify.get("rank","")).upper()!="FAMILY"
            or str(verify.get("canonicalName") or verify.get("scientificName") or "").lower().split(" ")[0]!=family.lower()):
            raise ValueError("GBIF family resolved to inconsistent taxon "+family+" "+str(key))
        by_family[family]=str(key)
    assert len(set(by_family.values()))==len(families)
    return {"method":"REPEATED_GBIF_FAMILY_TAXONKEY_OR_FILTER",
        "family_taxonkeys":by_family,
        "papilionoidea_strict_match_failed":True,
        "papilionoidea_source_match":match,
        "unverified_if_any_candidate_225_outside_six_families":True}

def load_polygons(zip_path,working,ids,max_chars):
    """Split intricate original GSHHG islands into area-preserving disjoint WKT."""
    import geopandas as gpd
    import shapely
    from shapely.geometry import Polygon,MultiPolygon,GeometryCollection,box
    from shapely.geometry.polygon import orient
    from shapely.ops import transform
    from shapely import wkt
    from pyproj import Transformer
    raw=zip_path.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=GSHHG_SHA256:raise ValueError("Original GSHHG SHA256 mismatch "+digest)
    working.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=[p for p in z.namelist() if p.endswith((".shp",".shx",".dbf",".prj")) and
                 "/h/GSHHS_h_L1." in "/"+p]
        if len(members)<3:raise ValueError("GSHHG original shapefile missing")
        for member in members:z.extract(member,working)
    land=gpd.read_file(next(working.rglob("GSHHS_h_L1.shp")))
    land["landmass_id"]=land.sibling_id.astype(int).astype(str)
    if not ids<=set(land.landmass_id):raise ValueError("Source GSHHG island IDs missing")
    projection=Transformer.from_crs(4326,6933,always_xy=True)

    def polys(geom):
        if geom.is_empty:return []
        if isinstance(geom,Polygon):return [geom]
        if isinstance(geom,(MultiPolygon,GeometryCollection)):
            return [p for g in geom.geoms for p in polys(g)]
        return []
    def area(geom):
        return transform(projection.transform,geom).area
    def shorten(poly,threshold):
        original=area(poly)
        if original<=0 or not poly.is_valid:return None
        for tol in (0.,.00015,.00025,.0005,.001,.002,.003,.005,.008,.012,
                    .02,.035,.06,.10):
            q=poly.simplify(tol,preserve_topology=True) if tol else poly
            if not isinstance(q,Polygon) or q.is_empty or not q.is_valid:continue
            q=orient(q,sign=1.0)
            text=wkt.dumps(q,rounding_precision=5,trim=True)
            if len(text)<=max_chars:
                ratio=area(wkt.loads(text))/original
                if threshold<=ratio<=2-threshold:
                    return text,tol,original,ratio
        return None
    def fragments(poly,depth=0):
        short=shorten(poly,.96 if depth==0 else .985)
        if short is not None:return [short]
        if depth>=10:raise ValueError("Cannot preserve coastline in 10 geographic splits")
        minx,miny,maxx,maxy=poly.bounds
        kmx=(maxx-minx)*111*max(.05,abs(math.cos(math.radians((miny+maxy)/2))))
        kmy=(maxy-miny)*111
        if kmx<=0 or kmy<=0:raise ValueError("Degenerate geographic island extent")
        if kmx>=kmy:
            mid=(minx+maxx)/2
            sides=(box(minx-1,miny-1,mid,maxy+1),box(mid,miny-1,maxx+1,maxy+1))
        else:
            mid=(miny+maxy)/2
            sides=(box(minx-1,miny-1,maxx+1,mid),box(minx-1,mid,maxx+1,maxy+1))
        part=[g for side in sides for g in polys(poly.intersection(side)) if area(g)>1.0]
        if len(part)<2:raise ValueError("Island component could not be divided into disjoint geographic tiles")
        out=[]
        for g in part:
            out.extend(fragments(g,depth+1))
            if len(out)>100:raise ValueError("More than 100 island fragments needed")
        return out

    result={};audit=[]
    for island in sorted(ids,key=int):
        original=shapely.union_all(land.loc[land.landmass_id.eq(island),"geometry"].to_numpy())
        if original.is_empty or not original.is_valid:
            raise ValueError("Invalid original GSHHG island "+island)
        original_area=area(original)
        components=polys(original)
        pieces=[x for part in components for x in fragments(part)]
        if not pieces or len(pieces)>100:
            raise ValueError("Island has zero or too many geometry pieces "+island)
        values=[x[0] for x in pieces]
        ratio=sum(area(wkt.loads(w)) for w in values)/max(original_area,1e-10)
        if not .97<=ratio<=1.03:
            raise ValueError("Island fragment area error exceeds 3%: "+island+" "+str(ratio))
        # Each intersection is made with exactly adjacent (not overlapping)
        # boxes, but simplification can create small overlap. Audit explicitly.
        source_parts=[wkt.loads(w) for w in values]
        for i in range(len(values)):
            for j in range(i+1,len(values)):
                overlap=source_parts[i].intersection(source_parts[j])
                if area(overlap)>max(2.,original_area*1e-8):
                    raise ValueError("Overlapping GBIF island fragments "+island)
        result[island]=values
        audit.append({
            "island_id":island,"status":"AREA_PRESERVED_COUNTERCLOCKWISE_DISJOINT_FRAGMENT_GEOMETRY",
            "original_area_km2":original_area/1e6,
            "original_polygon_components":len(components),
            "query_fragments":len(values),"simplified_area_retained":ratio,
            "max_fragment_WKT_chars":max(map(len,values)),
            "WKT_sha256_fragments":[hashlib.sha256(v.encode()).hexdigest() for v in values],
            "source_GSHHG_zip_SHA256":digest,
            "fragmentation_invoked":len(values)>len(components)})
        print("ISLAND_GEOMETRY",island,"fragments",len(values),
              "area_retained",round(ratio,6),
              "max_WKT_chars",max(map(len,values)),flush=True)
    return result,audit,digest


def indexed_butterflies_in_period(island,period,geometry,taxon_filter):
    years=PERIODS[period]
    basic={"geometry":geometry,"year":f"{years[0]},{years[1]}",
        "hasCoordinate":"true","hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT","facet":"speciesKey","facetLimit":1000,
        "limit":0,
    }
    if taxon_filter["method"]=="GBIF_SUPERFAMILY_TAXONKEY":
        basic["taxonKey"]=taxon_filter["GBIF_taxonKey"]
    else:
        # 'family=Nymphalidae' is silently ignored by the actual GBIF
        # occurrence/search endpoint (country + impossible-family returned
        # identical counts in source diagnosis). ONLY accepted family taxonKeys.
        basic["taxonKey"]=list(taxon_filter["family_taxonkeys"].values())
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
      "butterfly_filter_type":taxon_filter["method"],
      "GBIF_indexed_not_biological_absence":True}

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
    taxon=resolve_butterfly_taxon_filter()
    polygons,geometry_receipt,gshhg=load_polygons(args.gshhg_zip,args.work,selected,args.wkt_max_chars)
    HERE.mkdir(parents=True,exist_ok=True)
    facts={}
    queries=[]
    for island in sorted(selected,key=int):
        for period in PERIODS:
            counts=defaultdict(int)
            parts=[]
            for ix,polygon in enumerate(polygons[island],1):
                found,receipt=indexed_butterflies_in_period(island,period,polygon,taxon)
                for species_key,n in found.items():
                    counts[species_key]+=n
                parts.append({
                   "geometry_fragment":ix,
                   "WKT_sha256":hashlib.sha256(polygon.encode()).hexdigest(),
                   "indexed_occurrences":receipt["all_superfamily_papilionoidea_indexed_occurrences"],
                   "species_keys_indexed":len(found),
                   "facet_pages":receipt["species_facet_pages"]})
                time.sleep(.10)
            facts[(island,period)]=dict(counts)
            receipt={"island_id":island,"period":period,
                "all_butterfly_indexed_occurrences_disjoint_polygon_sum":
                    sum(p["indexed_occurrences"] for p in parts),
                "gbif_specieskey_count":len(counts),
                "geometry_fragments":len(parts),"fragment_receipts":parts,
                "index_is_not_biological_occupancy":True}
            queries.append(receipt)
            print("BUTTERFLY_PERIOD",island,period,receipt["all_butterfly_indexed_occurrences_disjoint_polygon_sum"],
                "all_species_facets",len(counts),"candidate225",sum(k in counts for k in
                       (x["taxon_key"] for x in candidates)),flush=True)
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
           "GSHHG disjoint fragment polygons retain >=97% island area; simplified coastlines can change nearshore assignment.",
           "Zero-overlap fragment interiors can share borders; GBIF points exactly on a split border may be double-indexed.",
           "All 225 species × island cells are included; many are geographically infeasible, so all-species crude risk rates are not transport effects.",
           "Query negatives are GBIF nonindexes, not biological absences, occupancy, or proof of no colonization.",
           "Accepted family taxonKey filters and the 225 source speciesKey concepts may drift under GBIF checklist upgrades; audit exact overlap.",
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
