#!/usr/bin/env python3
"""Repair only HTTP400 GBIF year-facet queries caused by complex Level3 WKT.

GBIF's supported occurrence-search geometry is POLYGON (not MULTIPOLYGON).
Split official WGSRPD3 union into geographic POLYGON parts, simplify with high
area fidelity, and use repeatable geometry params in bounded batches. If no
date is returned from a PARTIAL source geometry, retain unknown (not botanical
non-occurrence). The rescued year is an indexed plant documentation year,
never a naturalization or transport date.
"""
from __future__ import annotations
import argparse,csv,json,time,math
from collections import Counter,defaultdict
from pathlib import Path
from shapely.geometry import shape,Polygon,MultiPolygon
from shapely.ops import unary_union,transform
from shapely import make_valid,wkt
from shapely.geometry.polygon import orient
from pyproj import Transformer
import audit_q4_gbif_plant_first_year_pilot as gbif

MAX_WKT_BATCH=5300
MAX_COMPONENTS_PER_REGION=20
MAX_NEW_REQUESTS=100
STATUS_GOOD="DATED_OCCURRENCE_FOUND"
STATUS_NO="NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT"

def readcsv(p):
    with p.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def savecsv(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader();writer.writerows(data)

def polygons(geojson,codes):
    source=json.loads(geojson.read_text(encoding="utf-8"))
    target=defaultdict(list)
    for feat in source["features"]:
        code=str((feat.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        if code in codes and feat.get("geometry"):
            geom=make_valid(shape(feat["geometry"]))
            if not geom.is_empty:target[code].append(geom)
    t=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform
    out={}
    for code in codes:
        if code not in target:raise RuntimeError("official polygons missing "+code)
        whole=make_valid(unary_union(target[code]))
        def components(g):
            if g.geom_type=="Polygon":return [g]
            if g.geom_type=="MultiPolygon":return list(g.geoms)
            if g.geom_type=="GeometryCollection":
                return [p for o in g.geoms for p in components(o)]
            return []
        raw=components(whole)
        geo_area=lambda g:abs(transform(t,g).area)
        raw=sorted(raw,key=geo_area,reverse=True)
        total=sum(geo_area(x) for x in raw)
        selected=[]
        for poly in raw[:MAX_COMPONENTS_PER_REGION]:
            original_area=geo_area(poly)
            best=None
            for tol in (0,.0003,.0008,.0015,.003,.006,.012,.025,.05):
                shp=make_valid(poly.simplify(tol,preserve_topology=True))
                if shp.geom_type!="Polygon" or shp.is_empty:continue
                text=wkt.dumps(orient(shp,sign=1.0),rounding_precision=4,trim=True)
                if len(text)>4500:continue
                intersection=geo_area(make_valid(shp.intersection(poly)))
                purity=intersection/max(geo_area(shp),1e-12)
                coverage=intersection/max(original_area,1e-12)
                if purity>=.999 and coverage>=.995:
                    best={"wkt":text,"poly_area":original_area,
                          "shapefile_component_fidelity":coverage,
                          "polygon_purity":purity,
                          "tolerance_deg":tol}
                    break
            if best:selected.append(best)
        selected_area=sum(x["poly_area"]*x["shapefile_component_fidelity"] for x in selected)
        out[code]={
            "source_multipolygon_parts":len(raw),
            "selected_queryable_parts":len(selected),
            "selected_area_approx_fraction":selected_area/max(total,1e-12),
            "parts":selected,
            "source_has_antimeridian":whole.bounds[2]-whole.bounds[0]>180
        }
    return out

def chunks(geometries):
    result=[];current=[];length=0
    for item in geometries:
        n=len(item["wkt"])+20
        if current and length+n>MAX_WKT_BATCH:
            result.append(current);current=[];length=0
        current.append(item["wkt"]);length+=n
    if current:result.append(current)
    return result

def query_year(taxon,poly_wkts):
    base=[("taxonKey",str(taxon)),("hasCoordinate","true"),
          ("hasGeospatialIssue","false"),("occurrenceStatus","PRESENT"),
          ("year","1700,2026"),("facet","year"),("facetLimit","350"),("limit","0")]
    base.extend(("geometry",g) for g in poly_wkts)
    r=gbif.request("/occurrence/search",base)
    bins=next((f.get("counts",[]) for f in (r.get("facets") or [])
                if str(f.get("field","")).upper()=="YEAR"),None)
    if bins is None:raise RuntimeError("Missing GBIF year facet")
    years=[]
    for z in bins:
        try:y=int(z["name"])
        except (ValueError,TypeError,KeyError):continue
        if 1700<=y<=2026 and int(z.get("count") or 0)>0:years.append(y)
    if int(r.get("count") or 0)>0 and not years:raise RuntimeError("GBIF count positive but dated year facet empty")
    return min(years) if years else None,int(r.get("count") or 0)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-plant-years",type=Path,required=True)
    ap.add_argument("--wgsrpd3-geojson",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    data=readcsv(a.source_plant_years)
    if len(data)!=100:raise RuntimeError("100 frozen case-control host sources changed")
    failed=[r for r in data if r["source_date_status"]=="GBIF_QUERY_FAILED" and "HTTP Error 400" in r["error"]]
    if len(failed)!=11:raise RuntimeError("Original 11 HTTP400 plant region errors changed")
    geometric=polygons(a.wgsrpd3_geojson,{r["island_code"] for r in failed})
    result=[];requests=0
    for row in failed:
        geo=geometric[row["island_code"]]
        selected=geo["parts"]
        valid_polys=[x["wkt"] for x in selected]
        grouped=chunks(selected)
        successes=[];errors=[];records=0
        if geo["source_has_antimeridian"]:
            errors.append("Unsafe antimeridian span")
        elif not grouped:errors.append("No sufficiently accurate POLYGON fragments")
        else:
            for group in grouped:
                if requests>=MAX_NEW_REQUESTS:
                    errors.append("Cap on GBIF requests reached")
                    break
                requests+=1
                try:
                    y,n=query_year(row["gbif_species_key"],group)
                    records+=n
                    if y is not None:successes.append(y)
                except Exception as e:
                    errors.append(str(e)[:180])
                time.sleep(.8)
        fraction=geo["selected_area_approx_fraction"]
        full=not errors and fraction>=.9999
        if successes:
            row["source_date_status"]=STATUS_GOOD
            row["first_gbif_indexed_plant_year"]=str(min(successes))
            row["source_reuse"]="SPLIT_POLYGON_GBIF_DATE"+("_COMPLETE" if full else "_PARTIAL")
            row["error"]=""
        elif full and records==0:
            row["source_date_status"]=STATUS_NO
            row["first_gbif_indexed_plant_year"]=""
            row["source_reuse"]="SPLIT_POLYGON_GBIF_NO_INDEX"
            row["error"]=""
        else:
            row["source_reuse"]="SPLIT_POLYGON_INCOMPLETE_OR_UNRESOLVED"
            row["error"]="; ".join(errors)[:250] or "Source polygon portions not fully covered"
        result.append({
           "plant_id":row["plant_id"],
           "plant":row["wcvp_plant_name"],
           "region":row["island_code"],
           "n_source_polygons":geo["source_multipolygon_parts"],
           "n_queried_polygons":geo["selected_queryable_parts"],
           "area_covered_approximately":fraction,
           "queries":len(grouped),
           "indexed_year_min":min(successes) if successes else None,
           "source_query_record_count":records,
           "all_components_successful":full,
           "status":row["source_date_status"],
           "errors":errors,
        })
        print("GBIF geometry rescue",row["island_code"],row["wcvp_plant_name"],
              row["source_date_status"],"components",geo["source_multipolygon_parts"],
              "covered",round(fraction,5),flush=True)
    rr={ (r["plant_id"],r["island_code"]):r for r in failed}
    full=readcsv(a.source_plant_years)
    for i,r in enumerate(full):
        key=(r["plant_id"],r["island_code"])
        if key in rr:full[i]=rr[key]
    savecsv(a.out_dir/"matched_host_plant_gbif_firstyears_rescued.csv",full)
    audited={
       "schema":"chocho_Q4_GBIF_nonmultipolygon_400_repair_v0.1",
       "status":"GEOGRAPHIC_ERROR_RESCUE_NOT_INTRODUCTION_OR_ESTABLISHMENT",
       "origin":"Original observed-vs-nonindexed matched 100 plant × region GBIF cases",
       "reported_GBIF_400_pairs":len(failed),
       "distinct_official_WGSRPD3_regions":len(geometric),
       "new_GBIF_polygon_batch_queries":requests,
       "new_dated_plant_index_pairs":sum(z["status"]==STATUS_GOOD for z in result),
       "remaining_failed_pairs":sum(z["status"]=="GBIF_QUERY_FAILED" for z in result),
       "source_regions":{k:{"polygon_parts":g["source_multipolygon_parts"],
                           "queryable_parts":g["selected_queryable_parts"],
                           "area_coverage":g["selected_area_approx_fraction"]} for k,g in geometric.items()},
       "outcomes":result,
       "critical_limits":[
         "GBIF occurrence API does not accept MULTIPOLYGON WKT in the documented search predicate; we split only erroring regions into separate CCW POLYGON components.",
         "An indexed year found in any high-purity polygon fragment is real documentary source evidence for part of the WGSRPD3 island area; it is not an exact plant first arrival or guaranteed first record of the entire macroregion.",
         "If summed polygons are partial, source query no-record does not count as botanical absence or a verified negative old host.",
         "Geometries and taxon statuses in previous successful queries stay frozen and are not re-fetched.",
         "Repeated GBIF geometry parameters may be interpreted as an OR in occurrence/search; retrieved records from disjoint components may still be duplicated on shared borders.",
         "The 11 source queries were all HTTP400 before rescue; inability to query some official fragments remains missing."
       ]
    }
    (a.out_dir/"gbif_400_polygon_repair_audit.json").write_text(json.dumps(audited,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GBIF_SPLIT_GEOMETRY_AUDIT",json.dumps({k:v for k,v in audited.items() if k not in ("source_regions","outcomes","critical_limits")},ensure_ascii=False))
if __name__=="__main__":main()
