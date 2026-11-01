#!/usr/bin/env python3
"""Same-WGSRPD3 historical-exposure comparison: GloNAF vs butterfly resources.

Run only with source-pinned Lenzner 2022 ZIP, exact TDWG Level3 polygons,
Lenzner->TDWG duration crosswalk and a frozen 239-butterfly resource-pair table.

Produces a positive-control check: is shared European empire duration linked
to the co-occurrence of naturalized alien plants AFTER harmonizing geography?
Tests butterflies on the identical pair support, rather than comparing the
original 1,062-region Lenzner association to 284 independent WGSRPD3 sites.

Outputs only derived counts/correlations, not raw GloNAF source species rows,
source polygons, individual site data, or inferred voyages.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
from pyproj import Transformer
import shapefile
from shapely.geometry import shape
from shapely.ops import transform as geo_transform
from shapely.strtree import STRtree
from shapely.validation import make_valid
from shapely.errors import GEOSException

from analyze_colonial_legacy_resource_beta import (
    load as load_colonial_butterfly,
    orthogonal_controls,residual,correlation,shared_years,
)
from build_lenzner_tdwg_colonial_crosswalk import (
    canonical_id,load_csv_member,geom_clean,check_source,
)

BASE=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
ZIP=Path("_temp_historical_sources/lenzner_6640264.zip")
TDWG=Path("_external/wgsrpd/geojson/level3.geojson")
LOCK=BASE/"matched_glonaf_consumer_lock_v1.json"
OUT=BASE/"matched_glonaf_consumer_history_v1.json"
OUTREG=BASE/"matched_glonaf_consumer_region_coverage_v1.csv"
MIN_SOURCE_SHARE=0.5
THRESHOLDS=(0.01,0.05,0.10)


def read_glonaf_taxa(source_zip):
    """Original published alien plants; source site ID, accepted binomial-ish key."""
    with zipfile.ZipFile(source_zip) as z:
        rows,damaged=load_csv_member(z,"GloNAF.csv")
    flora=defaultdict(set)
    failures=Counter()
    for r in rows:
        key=canonical_id(r["OBJIDsic"])
        genus=(r.get("genus") or "").strip()
        epithet=(r.get("Species") or "").strip()
        if not key or not genus or not epithet or "\ufffd" in (genus+epithet):
            failures["missing_or_corrupt_name_or_key"]+=1
            continue
        name=(epithet if epithet.casefold().startswith(genus.casefold()+" ")
              else genus+" "+epithet)
        if len(name)<5 or " sp." in name.casefold():
            failures["not_species_resolvable"]+=1
            continue
        flora[key].add(name.casefold())
    names={name for value in flora.values() for name in value}
    if len(rows)!=470962 or len(names)<15000:
        raise RuntimeError("Original GloNAF flora count/species naming changed")
    index={name:i for i,name in enumerate(sorted(names))}
    bits={key:sum((1<<index[name] for name in values),0)
          for key,values in flora.items()}
    return bits,{"original_rows":len(rows),"distinct_source_site_ids":len(flora),
                 "distinct_species":len(index),
                 "source_text_replacement_characters":damaged,
                 "name_exclusions":dict(failures)}


def read_geometries(source_zip):
    with zipfile.ZipFile(source_zip) as z:
        shp=[name for name in z.namelist() if name.endswith(".shp")]
        shx=[name for name in z.namelist() if name.endswith(".shx")]
        dbf=[name for name in z.namelist() if name.endswith(".dbf")]
        if len(shp)!=len(shx)!=len(dbf):
            raise RuntimeError("Unexpected number of shapefile members")
        if len(shp)!=1 or len(shx)!=1 or len(dbf)!=1:
            raise RuntimeError("Expected one original Lenzner shapefile")
        reader=shapefile.Reader(
            shp=io.BytesIO(z.read(shp[0])),
            shx=io.BytesIO(z.read(shx[0])),
            dbf=io.BytesIO(z.read(dbf[0])),
            encoding="latin1")
        if len(reader)!=1184:
            raise RuntimeError("Original polygon support count mismatch")
        transformer=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
        ids=[]
        geoms=[]
        invalid=Counter()
        for rec in reader.iterShapeRecords():
            key=canonical_id(rec.record.as_dict().get("OBJIDsic"))
            try:
                gg=geom_clean(shape(rec.shape.__geo_interface__))
                g=geom_clean(geo_transform(transformer.transform,gg))
            except (GEOSException,ValueError,TypeError):
                invalid["bad_source_polygon"]+=1
                continue
            if not key or g.is_empty or not math.isfinite(g.area) or g.area<=0:
                invalid["unusable_source_polygon"]+=1
                continue
            ids.append(key)
            geoms.append(g)
    if len(geoms)<1100:
        raise RuntimeError("Too few valid source polygons")
    return ids,geoms,transformer,dict(invalid)


def classify_source_contribution(target_area,source_area,intersection_area,threshold):
    if intersection_area <= 0 or target_area <= 0 or source_area <= 0:
        return False
    # Large source entirely covering a small TDWG unit and small source inside
    # a large TDWG unit are both eligible. Tiny slivers are excluded.
    return (intersection_area/target_area >= threshold
            or intersection_area/source_area >= MIN_SOURCE_SHARE)


def produce_matched_regions(source_zip,tdwg_path,thresholds=THRESHOLDS):
    check_source(source_zip)
    sourcebits,plant_audit=read_glonaf_taxa(source_zip)
    ids,geoms,transformer,invalid=read_geometries(source_zip)
    lookup={code:int(k) for k,code in enumerate(ids)}
    if len(lookup)!=len(ids):
        raise RuntimeError("Source polygon ID duplicated")
    tree=STRtree(geoms)
    sourcearea=[g.area for g in geoms]
    tdwg=json.loads(tdwg_path.read_text(encoding="utf8"))
    polys={}
    for f in tdwg.get("features",[]):
        code=str(f["properties"].get("LEVEL3_COD","")).strip()
        if not code or code in polys:
            raise RuntimeError("Missing or duplicate TDWG Level3 code")
        g=geom_clean(geo_transform(transformer.transform,
                       geom_clean(shape(f["geometry"]))))
        if g.is_empty or g.area<=0:
            raise RuntimeError("Invalid TDWG region polygon: "+code)
        polys[code]=g
    if len(polys)!=369:
        raise RuntimeError("Not official 369 TDWG regions")

    output={t:{} for t in thresholds}
    source_counts={t:{} for t in thresholds}
    threshold_audits=Counter()
    for code in sorted(polys):
        g=polys[code]
        overlap=[]
        for candidate in tree.query(g):
            k=int(candidate)
            try:
                inter=g.intersection(geoms[k])
            except GEOSException:
                inter=make_valid(g).intersection(make_valid(geoms[k]))
            if inter.is_empty or inter.area<=0:
                continue
            overlap.append((k,float(inter.area)))
        for threshold in thresholds:
            mask=0
            matched_source=0
            for k,a in overlap:
                if not classify_source_contribution(g.area,sourcearea[k],a,threshold):
                    continue
                v=sourcebits.get(ids[k],0)
                if v:
                    mask |= v
                    matched_source+=1
            output[threshold][code]=mask
            source_counts[threshold][code]=matched_source
            if matched_source==0:
                threshold_audits[f"no_flora_sources_{threshold}"]+=1
    return output,source_counts,{
        "glonaf":plant_audit,
        "lensed_polygons":len(ids),
        "source_geometry_issues":invalid,
        "tdwg_regions":len(polys),
        "spatial_rule":"source polygon >=50pct within target OR covers >=threshold share of target",
        "thresholds":list(thresholds),
        "matched_source_absence_per_threshold":dict(threshold_audits),
    }


def estimate_on_fixed_pairs(data,masks,threshold,draws=0,seed=31031):
    pair_i,pair_j=data["ia"],data["ib"]
    valid=np.fromiter((bool(masks[int(i)] | masks[int(j)])
                      for i,j in zip(pair_i,pair_j)),bool,count=len(pair_i))
    if int(valid.sum())<20000:
        raise RuntimeError("Too few native/resource region pairs with mapped alien flora")
    i,j=pair_i[valid],pair_j[valid]
    y=np.fromiter((float((masks[int(a)] & masks[int(b)]).bit_count())/
                   (masks[int(a)] | masks[int(b)]).bit_count()
                   for a,b in zip(i,j)),float,count=len(i))
    exp=data["exposure"][valid]
    plant_delta=data["plant_delta"][valid]
    butterfly_delta=data["resource_delta"][valid]
    flora_control=[arr[valid] for arr in data["control_sets"]["colonial_vs_plant_beta"]]
    flora_strict_control=flora_control+[plant_delta]
    butterfly_control=[arr[valid] for arr in
                        data["control_sets"]["colonial_vs_resource_beta_after_flora"]]
    cases={
      "H1_glonaf_alien_plant_legacy":(exp,y,flora_control),
      "H1_strict_after_WCVP_general_flora_delta":(exp,y,flora_strict_control),
      "H2_consumer_resource_delta_after_general_flora":(exp,butterfly_delta,butterfly_control),
      "H3_alien_plant_similarity_vs_consumer_resource_delta":(
          y,butterfly_delta,butterfly_control+[exp]),
    }
    basis={}
    response={}
    observed={}
    exposure={}
    for name,(x,yy,controls) in cases.items():
        z=orthogonal_controls(controls)
        basis[name]=z
        response[name]=residual(yy,z)
        exposure[name]=residual(x,z)
        observed[name]=correlation(exposure[name],response[name])
    if any(v is None for v in observed.values()):
        raise RuntimeError("Matched colonial/flora response unidentifiable")
    out={
       "threshold_share_of_target_region":threshold,
       "eligible_region_nodes":len(masks),
       "evaluated_pair_count":int(valid.sum()),
       "removed_pairs_with_both_alien_plant_lists_empty":int((~valid).sum()),
       "matched_GloNAF_plant_richness_min":min(v.bit_count() for v in masks),
       "matched_GloNAF_plant_richness_median":float(np.median(
           [v.bit_count() for v in masks])),
       "observed_partial_rank_correlations":observed,
       "positive_control_alien_plant_jaccard_mean":float(y.mean()),
       "butterfly_resource_delta_mean_same_pair_set":float(butterfly_delta.mean()),
    }
    if draws==0:
        return out
    sourceyears=data["years"]
    continent=data["continent"]
    groups=defaultdict(list)
    for k,code in enumerate(continent):
        groups[str(code)].append(k)
    fixed=[value for value in groups.values() if len(value)>1]
    rng=np.random.default_rng(seed)
    hist_null={name:[] for name in cases if name!="H3_alien_plant_similarity_vs_consumer_resource_delta"}
    years_min=25
    for draw in range(draws):
        perm=sourceyears.copy()
        for group in fixed:
            perm[group]=sourceyears[rng.permutation(group)]
        x=np.log1p(shared_years(perm[i],perm[j],years_min=years_min))
        for name in hist_null:
            rr=residual(x,basis[name])
            v=correlation(rr,response[name])
            if v is None:
                raise RuntimeError("Null colonial permutation became constant")
            hist_null[name].append(v)
    out["history_permutation_reference"]={
      "draws":draws,"seed":seed,
      "scheme":"whole 4-empire vector shuffle within level1 region",
      "source_continent_groups":{k:len(v) for k,v in sorted(groups.items())},
      "comparisons":{
        key:{
          "observed":observed[key],
          "null_median":float(np.median(v)),
          "null_q025":float(np.quantile(v,.025)),
          "null_q975":float(np.quantile(v,.975)),
          "positive_direction_tail":(1+sum(q>=observed[key] for q in v))/(len(v)+1),
          "abs_two_sided_tail":(1+sum(abs(q)>=abs(observed[key]) for q in v))/(len(v)+1),
        } for key,v in hist_null.items()
      }
    }
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-zip",type=Path,default=ZIP)
    ap.add_argument("--tdwg",type=Path,default=TDWG)
    ap.add_argument("--lock",type=Path,default=LOCK)
    ap.add_argument("--output",type=Path,default=OUT)
    ap.add_argument("--region-output",type=Path,default=OUTREG)
    ap.add_argument("--draws",type=int,default=499)
    args=ap.parse_args()
    lock=json.loads(args.lock.read_text(encoding="utf8"))
    if lock.get("schema")!="chocho_LenznerGloNAF_vs_HOSTS_matchedWGSrpd3_lock_v1":
        raise RuntimeError("Matched comparison lock schema changed")
    if not (99<=args.draws<=999):
        raise ValueError("Require 99-999 permutations")
    data=load_colonial_butterfly()
    assert data["source_info"]["accepted_regions"]==284
    geographies,source_counts,origin=produce_matched_regions(
        args.source_zip,args.tdwg)
    experiments={}
    for t in THRESHOLDS:
        masks=[geographies[t][code] for code in data["ids"]]
        experiments[str(t)]=estimate_on_fixed_pairs(
            data,masks,t,draws=args.draws if t==0.05 else 0)
    result={
      "schema":"chocho_lenzner_glonaf_vs_butterfly_resources_same_WGSRPD3_v1",
      "status":"SAME_GEOGRAPHY_OBSERVATIONAL_POSITIVE_CONTROL_NOT_COLONIAL_CAUSALITY",
      "source_origin":origin,
      "historical_source_zenodo_md5":"16f5585d2500676951a4b664bbf0ffde",
      "butterfly_beta_source":data["source_info"],
      "precomparison_lock_sha256":hashlib.sha256(args.lock.read_bytes()).hexdigest(),
      "primary_threshold":"0.05",
      "comparisons":experiments,
      "interpretation_limits":[
         "Historic exposure is four 2022-study empire-duration variables, not actual colonial-era shipping or import intensity.",
         "GloNAF alien flora present-day assemblage is an outcome LEVEL, while butterfly potential-resource Jaccard DELTA is a change against native-only WCVP counterfactual.",
         "Region maps from source GloNAF polygon union mix parts of source regions; 1pct/10pct thresholds are mandatory sensitivities.",
         "All region pairs share nodes; geographic permutation tail is exploratory and not a causal p-value.",
         "The observed current host catalogue is not date-specific or demonstrated larval use.",
         "GloNAF alien observations and WCVP current range define different plant universes.",
         "This analysis cannot establish host-mediated colonial legacy without independent butterfly occurrence records.",
      ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with args.region_output.open("w",encoding="utf8",newline="") as f:
        w=csv.writer(f)
        w.writerow(["wgsrpd3_code","source_80pct_coverage",
                    "glonaf_alien_species_thr01","glonaf_sources_thr01",
                    "glonaf_alien_species_thr05","glonaf_sources_thr05",
                    "glonaf_alien_species_thr10","glonaf_sources_thr10"])
        for code in data["ids"]:
            w.writerow([code,1,
                        geographies[.01][code].bit_count(),source_counts[.01][code],
                        geographies[.05][code].bit_count(),source_counts[.05][code],
                        geographies[.10][code].bit_count(),source_counts[.10][code]])
    print(json.dumps({k:v for k,v in result.items() if k!="source_origin"},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
