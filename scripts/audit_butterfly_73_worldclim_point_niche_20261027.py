#!/usr/bin/env python3
"""WorldClim 2.1 source-consistent climate at exact PR72 butterfly points.

Both point climate and the native HOST-PLANT-region climate baseline use
1970–2000 WorldClim 10-arcminute BIO1/BIO12 pixels. No claim about native
butterfly climatic niche, introduction causal history or breeding.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import tempfile
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import rowcol

from audit_butterfly_73_point_source_20261027 import (
    Reader, query_url, source_roster, MAX_QUERY_ROWS, MAX_PANEL_ROWS, PAGE,
)
from analyze_butterfly_host_point_proximity_20261026 import valid_point, BASIS
from analyze_bce_augmented_plant_network_climate_null_20261019 import (
    source_augmented_model,
)
from analyze_butterfly_239_GBIF_occurrence_host_opportunity_20261019 import (
    climate_and_plant_candidates,
)
from audit_butterfly_73_point_tdwg_resource_20261027 import RegionLocator

SCHEMA = "chocho_73_butterfly_worldclim_point_native_host_climate_20261027_v01"
WORLDCLIM_URL = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"
WORLDCLIM_TAG = "WorldClim 2.1, 1970–2000, 10-arc-minute BIO1/BIO12"
SPECIES_REFERENCE_GATE = 50
CATEGORIES = ("inside_5_95", "outside_5_95", "missing_climate",
              "native_host_reference_unavailable")


def weighted_quantile(values, weights, quantiles):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if not len(values) or len(values) != len(weights) or \
       not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)) or \
       np.any(weights <= 0):
        raise ValueError("Invalid weighted climatic distribution")
    ordered = np.argsort(values, kind="stable")
    v, w = values[ordered], weights[ordered]
    frac = (np.cumsum(w) - .5 * w) / sum(w)
    return [float(np.interp(q, frac, v)) for q in quantiles]


def checksummed_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_climate_from_zip(archive, workdir):
    contents = zipfile.ZipFile(archive)
    names = set(contents.namelist())
    wanted = ("wc2.1_10m_bio_1.tif", "wc2.1_10m_bio_12.tif")
    expected = {"wc2.1_10m_bio_%d.tif" % i for i in range(1,20)}
    if not expected <= names or len([n for n in names if n.endswith(".tif")]) != 19:
        raise ValueError("Wrong original WorldClim2.1 10m bundle of 19 BIO rasters")
    for name in names:
        p=Path(name)
        if p.is_absolute() or ".." in p.parts:
            raise ValueError("Unsafe WorldClim source entry")
    climate={}
    source={}
    transform=None
    for key,filename in (("bio1",wanted[0]),("bio12",wanted[1])):
        contents.extract(filename, path=workdir)
        f=workdir/filename
        with rasterio.open(f) as ds:
            if ds.count!=1 or ds.crs.to_epsg()!=4326 or \
               ds.width!=2160 or ds.height!=1080 or \
               not math.isclose(ds.transform.a,1/6,rel_tol=0,abs_tol=1e-7) or \
               not math.isclose(-ds.transform.e,1/6,rel_tol=0,abs_tol=1e-7):
                raise ValueError("Unexpected WorldClim original grid")
            if transform is not None and transform!=ds.transform:
                raise ValueError("BIO1 and BIO12 source grids not identical")
            transform=ds.transform
            data=ds.read(1,masked=True)
        clean=np.asarray(data.filled(np.nan),dtype=np.float32)
        clean[~np.isfinite(clean)] = np.nan
        climate[key]=clean
        source[key+"_sha256"]=checksummed_file(f)
    contents.close()
    valid=(np.isfinite(climate["bio1"]) & np.isfinite(climate["bio12"]) &
           (climate["bio1"] > -100) & (climate["bio1"] < 80) &
           (climate["bio12"] >= 0) & (climate["bio12"] < 25000))
    climate["valid"]=valid
    climate["transform"]=transform
    source.update({
        "zip_sha256":checksummed_file(archive),
        "source_url": WORLDCLIM_URL,
        "model": WORLDCLIM_TAG,
        "bio_grid_dimensions":list(valid.shape),
        "valid_land_climate_grid_cells":int(valid.sum()),
        "climate_units":{"bio1":"°C","bio12":"annual mm"},
    })
    return climate, source


def make_original_region_grid(features, regions, transform, shape_):
    locator=RegionLocator(features, regions)
    items=[(polygon, i+1) for i,polygon in enumerate(locator.polygons)]
    raster=rasterize(items, out_shape=shape_, transform=transform, fill=0,
                     all_touched=False, dtype="uint16")
    if raster.max() != len(regions) or not np.any(raster==1):
        raise ValueError("Original TDWG polygons not rasterized correctly")
    return raster, locator


def native_reference(species, bysp, native_raster, climate):
    native=int(bysp[species]["native"])
    ids=[i+1 for i in range(365) if (native>>i & 1)]
    if len(ids)<10:
        return {"status":"NATIVE_HOST_REGIONAL_SUPPORT_LT10"}
    allmask=np.isin(native_raster,np.asarray(ids,dtype=np.uint16)) & climate["valid"]
    rr,cc=np.nonzero(allmask)
    if len(rr)<SPECIES_REFERENCE_GATE:
        return {"status":"NATIVE_HOST_WORLDCLIM_PIXELS_LT50","pixels":len(rr)}
    t=climate["bio1"][rr,cc].astype(float)
    p=climate["bio12"][rr,cc].astype(float)
    # Equirectangular WGS84 grid cells have different areas by latitude.
    # Each pixel gets an approximate cosine latitude area weight.
    lat=90 + (rr+.5)*climate["transform"].e
    weights=np.cos(np.deg2rad(lat))
    if np.any(weights<=0):
        raise ValueError("WorldClim source land pixel at nonpositive latitude area")
    quantiles=(.05,.10,.90,.95)
    tq=weighted_quantile(t,weights,quantiles)
    pq=weighted_quantile(p,weights,quantiles)
    return {
        "status":"VALID_NATIVE_HOST_RESOURCE_WORLDCLIM_REFERENCE",
        "source_native_host_regions":len(ids),
        "source_native_worldclim_valid_pixels":int(len(rr)),
        "bounds":{
            "area_weighted_5_95":{
                "bio1":[tq[0],tq[3]],"bio12":[pq[0],pq[3]]
            },
            "area_weighted_10_90":{
                "bio1":[tq[1],tq[2]],"bio12":[pq[1],pq[2]]
            },
            "equal_pixel_5_95":{
                "bio1":[float(z) for z in np.quantile(t,[.05,.95])],
                "bio12":[float(z) for z in np.quantile(p,[.05,.95])]
            },
        },
    }


def extract_climate(climate,lat,lon):
    transform=climate["transform"]
    try:
        row,col=rowcol(transform,lon,lat)
    except (ValueError,OverflowError):
        return None
    if not 0<=row<climate["valid"].shape[0] or \
       not 0<=col<climate["valid"].shape[1] or \
       not climate["valid"][row,col]:
        return None
    return (float(climate["bio1"][row,col]),
            float(climate["bio12"][row,col]))


def outside_rectangle(values, rectangle):
    temp,precip=values
    return not (rectangle["bio1"][0] <= temp <= rectangle["bio1"][1] and
                rectangle["bio12"][0] <= precip <= rectangle["bio12"][1])


def score_case(case, reference, raster, locator, bysp, reader, remaining):
    n=0
    receipts=[]
    seen=set()
    year=Counter()
    quality=Counter()
    scores={
        criterion:Counter() for criterion in (
            "area_weighted_5_95", "area_weighted_10_90",
            "equal_pixel_5_95"
        )
    }
    strict=Counter()
    in_added=Counter()
    try:
        first,fdigest=reader.get(query_url(case,0))
        n=int(first["count"])
        if n>MAX_QUERY_ROWS or n>remaining:
            return {
                "status":"NOT_EVALUABLE_QUERY_RECORD_CAP",
                "species":case["species"],"country":case["country"],
                "source_class":case["resource_class"],"live_count":n}
        if n<0:
            raise ValueError("Negative GBIF source count")
        for off in list(range(0,n,PAGE)) or [0]:
            packet,digest=(first,fdigest) if off==0 else reader.get(query_url(case,off))
            if int(packet["count"]) != n or \
               len(packet["results"]) != min(PAGE,max(0,n-off)):
                raise ValueError("GBIF source pagination incomplete or count drift")
            receipts.append({"offset":off,"source_body_sha256":digest})
            for rec in packet["results"]:
                key=rec.get("key")
                if key is None or key in seen:
                    raise ValueError("Missing/duplicate original occurrence key")
                seen.add(key)
                pt,status=valid_point(rec,case["gbif_accepted_species_key"],case["country"])
                if pt is not None and rec.get("basisOfRecord") not in BASIS:
                    pt,status=None,"unaccepted_basis"
                quality[status]+=1
                if pt is None:
                    continue
                b=extract_climate(raster,pt["lat"],pt["lon"])
                if b is None:
                    quality["WorldClim_point_missing"]+=1
                    continue
                quality["WorldClim_point_valid"]+=1
                year["early_2010_2017" if pt["year"]<=2017 else "late_2018_2025"]+=1
                species_data=bysp[case["species"]]
                index,index_status=locator.lookup(pt["lat"],pt["lon"])
                if index_status=="matched":
                    if int(species_data["added"]) >> index & 1:
                        in_added["in_source_introduced_host_region"]+=1
                    else:
                        in_added["not_in_source_introduced_host_region"]+=1
                else:
                    in_added["unmapped_or_ambiguous_WGSRPD_region"]+=1
                for criterion in scores:
                    if reference["status"] != \
                            "VALID_NATIVE_HOST_RESOURCE_WORLDCLIM_REFERENCE":
                        scores[criterion]["native_host_reference_unavailable"]+=1
                        continue
                    o=outside_rectangle(b,reference["bounds"][criterion])
                    scores[criterion]["outside" if o else "inside"]+=1
                    if pt["strict"]:
                        strict[(criterion,"outside" if o else "inside")]+=1
        if len(seen)!=n:
            raise ValueError("Incomplete retrieved GBIF points")
        return {
            "status":"COMPLETE_GBIF_POINT_WORLDCLIM_CLIMATE",
            "species":case["species"],"country":case["country"],
            "source_class":case["resource_class"],
            "original_index_count":case["original_gbif_record_count"],
            "live_count":n,"count_drift":n-case["original_gbif_record_count"],
            "valid_worldclim_indices":quality["WorldClim_point_valid"],
            "source_quality":dict(quality),
            "year":dict(year),
            "point_outside_native_host_rect":{
                k:dict(v) for k,v in scores.items()},
            "point_known_uncertainty_le1000m":{
                k:{z:strict[(k,z)] for z in ("inside","outside")}
                for k in scores
            },
            "point_host_source_region_class":dict(in_added),
            "SHA256_source_pages":receipts,
        }
    except Exception as e:
        return {
            "status":"INCOMPLETE_GBIF_OR_WORLDCLIM_CLIMATE_SOURCE",
            "species":case["species"],"country":case["country"],
            "source_class":case["resource_class"],
            "reason":f"{type(e).__name__}: {str(e)}"[:240],
            "source_page_hashes_prior_failure":receipts,
        }


def group(cases,key):
    complete=[r for r in cases
              if r["status"]=="COMPLETE_GBIF_POINT_WORLDCLIM_CLIMATE"]
    case_ratio=[]
    within=outside=0
    n_valid=0
    for r in complete:
        row=r["point_outside_native_host_rect"][key]
        ins=row.get("inside",0)
        out=row.get("outside",0)
        within+=ins
        outside+=out
        n_valid+=r["valid_worldclim_indices"]
        if ins+out>0:
            case_ratio.append(out/(ins+out))
    return {
        "complete_cases":len(complete),
        "point_valid_total":n_valid,
        "inside":within,"outside":outside,
        "source_point_climate_coverage":within+outside,
        "record_weighted_outside_share":(
            outside/(within+outside) if within+outside>0 else None),
        "equal_species_country_case_mean_outside_share":(
            sum(case_ratio)/len(case_ratio) if case_ratio else None),
        "cases_in_case_mean":len(case_ratio),
    }


def run(args):
    with tempfile.TemporaryDirectory(prefix="chocho-worldclim-") as tmp:
        data=Path(tmp)
        climate,source=load_climate_from_zip(args.worldclim_zip,data)
        model,regions,chelsa,source_meta=source_augmented_model(args.source_root)
        bysp,_,_,_=climate_and_plant_candidates(model,regions,chelsa)
        raw=(args.source_root/"inputs/level3.geojson").read_bytes()
        features=json.loads(raw)["features"]
        reggrid,locator=make_original_region_grid(
            features, regions, climate["transform"],climate["valid"].shape)
        roster=source_roster(args.roster)
        included=sorted({r["species"] for r in roster["roster"]})
        references={sp:native_reference(sp,bysp,reggrid,climate)
                    for sp in included}
        reader=Reader()
        remaining=MAX_PANEL_ROWS
        cases=[]
        for case in sorted(roster["roster"],key=lambda r:(r["species"],r["country"])):
            q=score_case(case,references[case["species"]],climate,
                         locator,bysp,reader,remaining)
            cases.append(q)
            if q["status"]=="COMPLETE_GBIF_POINT_WORLDCLIM_CLIMATE":
                remaining-=q["live_count"]
        n_complete=sum(r["status"]=="COMPLETE_GBIF_POINT_WORLDCLIM_CLIMATE"
                       for r in cases)
        criteria=("area_weighted_5_95","area_weighted_10_90","equal_pixel_5_95")
        report={
            "schema":SCHEMA,
            "status":("COMPLETE_73_SAME_WORLDCLIM_PRODUCT_POINT_AND_HOST_BASELINE"
                      if n_complete==73 and
                      all(r["status"]=="VALID_NATIVE_HOST_RESOURCE_WORLDCLIM_REFERENCE"
                          for r in references.values())
                      else "INCOMPLETE_POINT_OR_NATIVE_WORLDCLIM_SOURCE"),
            "source":source,
            "source_author_only_tdwg_geojson_sha256":hashlib.sha256(raw).hexdigest(),
            "source_73_positive_country_cohort_verified":True,
            "source_butterfly_count":len(included),
            "source_worldclim_native_host_references":references,
            "cases_complete":n_complete,
            "cases_total":73,
            "live_GBIF_network_requests":reader.calls,
            "live_GBIF_retries":reader.retries,
            "cases":cases,
            "by_source_country_resource_class":{
                k:{criterion:group([r for r in cases
                                    if r["source_class"]==k],criterion)
                   for criterion in criteria}
                for k in ("novel","analogue")
            },
            "original_two_largest_source_case_removed":{
                criterion:group([r for r in cases if (r["species"],r["country"])
                   not in {(t["species"],t["country"]) for t in
                           sorted(roster["roster"],
                             key=lambda x:-x["original_gbif_record_count"])[:2]}],
                   criterion) for criterion in criteria
            },
            "no_raw_coordinates_or_full_worldclim_rasters_exported":True,
            "interpretation_gate":(
                "Source-independent raster point climatic novelty relative to "
                "potential native larval-HOST-PLANT regional range, not original "
                "butterfly physiological climatic niche, introduced-host use, "
                "wild establishment or population growth."
            ),
            "computed_utc":datetime.now(timezone.utc).isoformat(),
        }
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
        print(json.dumps({
            "status":report["status"],"complete":n_complete,
            "source_valid_climate_cells":source["valid_land_climate_grid_cells"],
            "native_host_worldclim_references_evaluable":sum(
                r["status"]=="VALID_NATIVE_HOST_RESOURCE_WORLDCLIM_REFERENCE"
                for r in references.values()),
            "by_source_country_class":{
                name:report["by_source_country_resource_class"][name][
                    "area_weighted_5_95"] for name in ("novel","analogue")},
            "network_requests":reader.calls,
        },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worldclim-zip",type=Path,required=True)
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--roster",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
