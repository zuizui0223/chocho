#!/usr/bin/env python3
"""Area-weighted historical-empire exposure: Lenzner GloNAF -> WGSRPD level 3.

Pinned inputs:
  - Lenzner et al. 2022, doi 10.5281/zenodo.6640264,
    full archive MD5 16f5585d2500676951a4b664bbf0ffde
  - TDWG/wgsrpd source commit 52da7828aba9d461dd133c27b3bd7a4407161f54
  - chocho data/frozen/wgsrpd3_level1_lookup_v0.1.json (369 codes).

This is geographic source alignment, NOT a colonial effect on butterflies.
Historical durations are ecological covariates from a published study;
area-weighted exposure does NOT imply a particular voyage, colonizer for every
plant, or the true temporal sequence of introductions.

Never republish source GloNAF geometries or original colonial polygon rows.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import shapefile
from shapely.geometry import shape as shapely_shape
from shapely.ops import transform as geo_transform
from shapely.strtree import STRtree
from shapely.validation import make_valid
from shapely.errors import GEOSException
from pyproj import Transformer

SOURCE_MD5="16f5585d2500676951a4b664bbf0ffde"
TDWG_COMMIT="52da7828aba9d461dd133c27b3bd7a4407161f54"
EMPIRES=("time.emp.GB","time.emp.ESP","time.emp.PRT","time.emp.NL")
ZIP=Path("_temp_historical_sources/lenzner_6640264.zip")
OFFICIAL=Path("data/frozen/wgsrpd3_level1_lookup_v0.1.json")
OUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
TDWG=Path("_external/wgsrpd/geojson/level3.geojson")


def canonical_id(raw):
    if raw is None or str(raw).strip() == "":
        return None
    t=str(raw).strip()
    try:
        return str(int(float(t)))
    except ValueError:
        return t


def load_csv_member(z, name):
    matching=[x for x in z.namelist() if x.endswith("/Datasets/"+name)]
    if len(matching)!=1:
        raise RuntimeError("Expected exactly one source "+name)
    with z.open(matching[0]) as f:
        raw=f.read()
    try:
        decoded=raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        # Some historical region labels contain CP1252 non-breaking spaces.
        # Preserve source characters rather than substituting replacement symbols.
        decoded=raw.decode("cp1252")
    return list(csv.DictReader(io.StringIO(decoded)))


def geom_clean(g):
    if g.is_empty:
        return g
    if not g.is_valid:
        g=make_valid(g)
    return g


def check_source(zip_path):
    h=hashlib.md5()
    with zip_path.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20),b""):
            h.update(chunk)
    if h.hexdigest()!=SOURCE_MD5:
        raise RuntimeError("Colonial archive MD5 drift")


def build(zip_path=ZIP, tdwg_path=TDWG, official_path=OFFICIAL):
    check_source(zip_path)
    original=json.loads(official_path.read_text(encoding="utf8"))
    assert original.get("source_commit")==TDWG_COMMIT
    official_codes=set(original["units"])
    assert len(official_codes)==369
    with zipfile.ZipFile(zip_path) as z:
        regions=load_csv_member(z,"Regions.csv")
        predictors=load_csv_member(z,"Predictors.csv")
        empires=load_csv_member(z,"Empires.csv")
        shp=[k for k in z.namelist() if k.lower().endswith(".shp")]
        shx=[k for k in z.namelist() if k.lower().endswith(".shx")]
        dbf=[k for k in z.namelist() if k.lower().endswith(".dbf")]
        prj=[k for k in z.namelist() if k.lower().endswith(".prj")]
        if not(len(shp)==len(shx)==len(dbf)==len(prj)==1):
            raise RuntimeError("Incomplete source region shapefile")
        proj=z.read(prj[0]).decode("utf8",errors="replace")
        if "WGS_1984" not in proj and "WGS 84" not in proj and "WGS84" not in proj:
            raise RuntimeError("Source GloNAF geometry datum not verified as WGS84")
        source=shapefile.Reader(
            shp=io.BytesIO(z.read(shp[0])),
            shx=io.BytesIO(z.read(shx[0])),
            dbf=io.BytesIO(z.read(dbf[0])))
        assert len(source)==1184

        reg_by_id={}
        for r in regions:
            key=canonical_id(r["OBJIDsic"])
            if key in reg_by_id:
                raise RuntimeError("Duplicate source GloNAF region ID")
            reg_by_id[key]=r
        predictors_by_id={}
        for r in predictors:
            key=canonical_id(r["OBJIDsic"])
            if key in predictors_by_id:
                raise RuntimeError("Duplicate colonial predictor ID")
            v={}
            for empire in EMPIRES:
                raw=(r.get(empire) or "").strip()
                try:
                    val=float(raw)
                except ValueError:
                    raise RuntimeError("Non-numeric duration at source "+key+"/"+empire)
                if not math.isfinite(val) or val<0 or val>600:
                    raise RuntimeError("Impossible colonial duration")
                v[empire]=val
            predictors_by_id[key]=v
        assert len(reg_by_id)==len(predictors_by_id)==1183

        transformer=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
        projected=[]
        identifiers=[]
        geom_issues=Counter()
        for sr in source.iterShapeRecords():
            key=canonical_id(sr.record.as_dict().get("OBJIDsic"))
            if key not in predictors_by_id or key not in reg_by_id:
                geom_issues["shapefile_id_not_in_predictors"]+=1
                continue
            try:
                sg=geom_clean(shapely_shape(sr.shape.__geo_interface__))
                gp=geo_transform(transformer.transform,sg)
                gp=geom_clean(gp)
            except (GEOSException, ValueError, TypeError) as err:
                geom_issues["invalid_source_geometry"]+=1
                continue
            if gp.is_empty or not math.isfinite(gp.area) or gp.area<=0:
                geom_issues["empty_source_polygon"]+=1
                continue
            projected.append(gp)
            identifiers.append(key)
        tree=STRtree(projected)
        if len(projected)<1100:
            raise RuntimeError("Too many rejected source region polygons")

    tdwg_data=json.loads(tdwg_path.read_text(encoding="utf8"))
    features=tdwg_data.get("features",[])
    codes=[str(f["properties"].get("LEVEL3_COD","")).strip() for f in features]
    if len(codes)!=369 or len(set(codes))!=369 or set(codes)!=official_codes:
        raise RuntimeError("Pinned TDWG polygon codes do not match frozen official 369")
    outputs=[]
    sums=Counter()
    area_coverage=[]
    for feature, code in sorted(zip(features,codes),key=lambda x:x[1]):
        geom=geom_clean(shapely_shape(feature["geometry"]))
        g=geo_transform(transformer.transform,geom)
        g=geom_clean(g)
        if g.is_empty or g.area<=0:
            raise RuntimeError("Invalid official polygon: "+code)
        total_area=float(g.area)
        matched_area=0.0
        contributions=Counter()
        contributing_source=set()
        candidate_indices=tree.query(g)
        for candidate in candidate_indices:
            # Shapely 2 returns integer positions, not geometry objects.
            k=int(candidate)
            try:
                overlap=g.intersection(projected[k])
            except GEOSException:
                overlap=make_valid(g).intersection(make_valid(projected[k]))
            if overlap.is_empty or overlap.area<=0:
                continue
            area=float(overlap.area)
            matched_area+=area
            contributing_source.add(identifiers[k])
            for empire in EMPIRES:
                contributions[empire]+=area*predictors_by_id[identifiers[k]][empire]
        cov_ratio=matched_area/total_area
        area_coverage.append(cov_ratio)
        # Sum of overlaps can exceed polygon area if Lenzner regions overlap.
        # Do not mislabel >1 as complete, or silently divide away duplicates.
        if cov_ratio>1.05:
            geom_issues["substantial_region_overlap"]+=1
        row={
            "wgsrpd3_code":code,
            "glonaf_covered_area_fraction_unclipped":cov_ratio,
            "glonaf_overlap_region_count":len(contributing_source),
            "usable_area_fifty_percent":int(cov_ratio>=.5 and cov_ratio<=1.05),
            "usable_area_eighty_percent":int(cov_ratio>=.8 and cov_ratio<=1.05),
        }
        for empire in EMPIRES:
            # Average duration over represented Lenzner footprint, rather
            # than divide by uncovered region area and imply zero colonialism.
            row[empire.replace("time.emp.","years_")]=(
                contributions[empire]/matched_area if matched_area>0 else "")
        outputs.append(row)
    return outputs,{
        "status":"AREA_WEIGHTED_COLONIAL_GEOMETRY_RECONSTRUCTION_NOT_CONSUMER_EFFECT",
        "source_doi":"10.5281/zenodo.6640264",
        "source_md5":SOURCE_MD5,
        "tdwg_source_commit":TDWG_COMMIT,
        "source_glonaf_regions":len(reg_by_id),
        "source_empire_history_rows":len(empires),
        "source_shapefile_records":len(source),
        "source_geometry_accepted":len(projected),
        "source_geometry_issues":dict(geom_issues),
        "official_wgsrpd3_regions":len(outputs),
        "covered_regions_over_50pct":sum(r["usable_area_fifty_percent"] for r in outputs),
        "covered_regions_over_80pct":sum(r["usable_area_eighty_percent"] for r in outputs),
        "regions_overlapping_source_above_105pct":sum(r["glonaf_covered_area_fraction_unclipped"]>1.05 for r in outputs),
        "min_coverage":min(area_coverage),"median_coverage":sorted(area_coverage)[len(area_coverage)//2],
        "max_coverage":max(area_coverage),
        "n_zero_coverage":sum(v==0 for v in area_coverage),
        "interpretation_limits":[
            "Area-weighted source-polygon colonial history is a spatial proxy; it does not show plant transfer dates.",
            "Lenzner GloNAF polygons and TDWG WGSRPD regions have different boundaries and may overlap irregularly.",
            "Colonial durations were not fitted to butterfly or host outcomes.",
            "Only historic GB/ESP/PRT/NL years included; other colonial powers omitted by the original study.",
            "Actual shipping routes, host quality, general-flora control and butterfly observations are separate gates.",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-zip",type=Path,default=ZIP)
    ap.add_argument("--tdwg-geojson",type=Path,default=TDWG)
    ap.add_argument("--official",type=Path,default=OFFICIAL)
    ap.add_argument("--output-dir",type=Path,default=OUT)
    args=ap.parse_args()
    rows,receipt=build(args.source_zip,args.tdwg_geojson,args.official)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/"lenzner_wgsrpd3_colonial_area_crosswalk.csv").open(
            "w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    (args.output_dir/"lenzner_wgsrpd3_colonial_area_audit.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(receipt,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
