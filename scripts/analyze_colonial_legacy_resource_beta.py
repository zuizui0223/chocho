#!/usr/bin/env python3
"""Does historical shared European empire occupancy match plant/consumer beta shifts?

Exploratory regional association on frozen 239 HOSTS butterflies and 369
WGSRPD3 units. The **source is observational** and the endpoint is potential
larval host opportunity, NOT butterfly arrivals/colonization/abundance.
Pre-outcome controls/eligibility/negative control are in source-locked JSON.
Only pairs with adequate geometry and climate covariates are included.

No ordinary independent-pair p-values. Null resamples whole four-empire
region-history vectors within level-1 region, retaining response and climate.
This is a sensitivity check, NOT randomized exposure or causal mediation.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

BASE=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
COLONIAL=BASE/"lenzner_wgsrpd3_colonial_area_crosswalk.csv"
LOCK=BASE/"colonial_resource_association_lock.json"
LOOKUP=Path("data/frozen/wgsrpd3_level1_lookup_v0.1.json")
BETA=Path("_external/butterfly_beta/provenance/reviewer_defenses/results/"
          "global_colonial_resource_beta_20261031/beta_baseline_v1.pairs.csv.gz")
BETA_RECEIPT=Path("_external/butterfly_beta/provenance/reviewer_defenses/results/"
                  "global_colonial_resource_beta_20261031/beta_baseline_v1.json")
REGIONAL=Path("_external/butterfly_beta/provenance/reviewer_defenses/results/"
              "global_colonial_resource_beta_20261031/beta_baseline_v1.region.csv")
EMPIRES=("GB","ESP","PRT","NL")
R_EARTH=6371.0088


def readcsv(path):
    with Path(path).open("r",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for blk in iter(lambda:f.read(1<<20),b""):
            h.update(blk)
    return h.hexdigest()


def haversine(lat_a,lon_a,lat_b,lon_b):
    a,b=np.radians(lat_a),np.radians(lat_b)
    dl=np.radians(np.asarray(lon_b)-np.asarray(lon_a))
    h=np.sin((b-a)/2)**2+np.cos(a)*np.cos(b)*np.sin(dl/2)**2
    return 2*R_EARTH*np.arcsin(np.sqrt(np.clip(h,0,1)))


def ranked(vals):
    a=np.asarray(vals,dtype=float)
    order=np.argsort(a,kind="stable")
    out=np.empty(len(a),float)
    i=0
    while i<len(a):
        j=i+1
        while j<len(a) and a[order[j]]==a[order[i]]:
            j+=1
        out[order[i:j]]=(i+1+j)/2
        i=j
    return out


def orthogonal_controls(arrays):
    """Projection basis with SVD rank truncation; no fit to outcomes."""
    z=np.column_stack([np.ones(len(arrays[0]))]+[ranked(x) for x in arrays])
    u,s,v=np.linalg.svd(z,full_matrices=False)
    n=np.count_nonzero(s > max(s)*1e-10)
    return u[:,:n]


def residual(value, basis):
    v=ranked(value)
    return v-basis@(basis.T@v)


def correlation(x,y):
    a=x-np.mean(x)
    b=y-np.mean(y)
    den=np.sqrt(np.sum(a*a)*np.sum(b*b))
    return float(np.sum(a*b)/den) if den>1e-12 else None


def shared_years(a,b,years_min=25):
    """No false inference of temporally overlapping colonial rule."""
    aa=np.asarray(a,float)
    bb=np.asarray(b,float)
    mask=(aa>=years_min)&(bb>=years_min)
    return np.sum(np.minimum(aa,bb)*mask,axis=-1)


def load(lock_path=LOCK,colonial_path=COLONIAL,beta_path=BETA,
         look_path=LOOKUP,regional_path=REGIONAL):
    lock=json.loads(lock_path.read_text(encoding="utf8"))
    assert lock["schema"]=="chocho_colonial_global_resource_beta_lock_v1"
    official=json.loads(look_path.read_text(encoding="utf8"))["units"]
    assert len(official)==369
    historical=readcsv(colonial_path)
    assert len(historical)==369 and {r["wgsrpd3_code"] for r in historical}==set(official)
    regional=readcsv(regional_path)
    assert len(regional)==369
    native_richness={r["wgsrpd3_code"]:float(r["plant_native_species"]) for r in regional}
    rows=[]
    climate_fields=("source_temp_mean","source_aridity_mean")
    history={}
    mincov=lock["eligibility"]["min_glonaf_area_coverage"]
    maxcov=lock["eligibility"]["max_glonaf_overlap_coverage"]
    for r in historical:
        code=r["wgsrpd3_code"]
        coverage=float(r["glonaf_covered_area_fraction_unclipped"])
        if not (mincov <= coverage <= maxcov):
            continue
        if not all(r.get(f,"") not in ("","NA","NaN") for f in climate_fields):
            continue
        h=np.asarray([float(r["years_"+e]) for e in EMPIRES],float)
        if not np.all(np.isfinite(h)) or np.any(h<0):
            continue
        v={
            "lat":float(r["representative_lat"]),
            "lon":float(r["representative_lon"]),
            "temp":float(r["source_temp_mean"]),
            "aridity":float(r["source_aridity_mean"]),
            "history":h,
            "level1":official[code],
            "native_richness":native_richness[code],
        }
        if not all(math.isfinite(v[x]) for x in
                   ("lat","lon","temp","aridity","native_richness")):
            continue
        history[code]=v
    if len(history)<180:
        raise RuntimeError("Too few source-matched regions with 80% spatial coverage")
    m=len(history)
    ids=sorted(history)
    index={code:i for i,code in enumerate(ids)}
    region_list=[history[code] for code in ids]
    with gzip.open(beta_path,"rt",encoding="utf8",newline="") as f:
        reader=csv.DictReader(f)
        for r in reader:
            a,b=r["region_a"],r["region_b"]
            if a not in index or b not in index:
                continue
            if not all(r.get(k) not in ("",None) for k in (
                    "resource_delta","plant_delta",
                    "resource_native_jaccard","plant_native_jaccard")):
                continue
            i,j=index[a],index[b]
            if i==j:
                raise RuntimeError("Self pair in butterfly-beta source")
            rows.append((i,j,float(r["resource_delta"]),float(r["plant_delta"]),
                         float(r["resource_native_jaccard"]),
                         float(r["plant_native_jaccard"])))
    if len(rows)<12000:
        raise RuntimeError("Too few usable region pairs for colonial association")
    vals=np.array(rows,float)
    ia=vals[:,0].astype(int)
    ib=vals[:,1].astype(int)
    resource_delta=vals[:,2]
    plant_delta=vals[:,3]
    res_native=vals[:,4]
    flora_native=vals[:,5]
    lat=np.asarray([x["lat"] for x in region_list])
    lon=np.asarray([x["lon"] for x in region_list])
    temps=np.asarray([x["temp"] for x in region_list])
    aridity=np.asarray([x["aridity"] for x in region_list])
    rich=np.asarray([x["native_richness"] for x in region_list])
    continents=np.asarray([x["level1"] for x in region_list])
    years=np.vstack([x["history"] for x in region_list])
    geo=np.log1p(haversine(lat[ia],lon[ia],lat[ib],lon[ib]))
    ct=np.abs(temps[ia]-temps[ib])
    ca=np.abs(aridity[ia]-aridity[ib])
    same=(continents[ia]==continents[ib]).astype(float)
    rich_min=np.log1p(np.minimum(rich[ia],rich[ib]))
    rich_max=np.log1p(np.maximum(rich[ia],rich[ib]))
    if not np.all(np.isfinite([geo,ct,ca,same,rich_min,rich_max])):
        raise RuntimeError("Invalid confounder source")
    years_min=int(lock["exposure"]["minimum_area_weighted_duration_years_per_region"])
    shared=shared_years(years[ia],years[ib],years_min=years_min)
    colonial=np.log1p(shared)
    binary=(shared>0).astype(float)
    assert len(geo)==len(resource_delta)
    baseline=[geo,same,ct,ca,flora_native,res_native,rich_min,rich_max]
    # Three contrasts require matched rank residuals with a fixed pair set:
    plant_controls=[geo,same,ct,ca,flora_native,rich_min,rich_max]
    resource_controls=baseline
    conditional_controls=baseline+[plant_delta]
    # Native-only response tests if "historical" just tags regional biogeography.
    native_controls=[geo,same,ct,ca,flora_native,rich_min,rich_max]
    return {
        "ids":ids,"continent":continents,"years":years,"ia":ia,"ib":ib,
        "shared":shared,"exposure":colonial,"ever":binary,
        "pairs":rows,
        "resource_delta":resource_delta,"plant_delta":plant_delta,
        "resource_native":res_native,"plant_native":flora_native,
        "control_sets":{
            "colonial_vs_plant_beta":plant_controls,
            "colonial_vs_resource_beta":resource_controls,
            "colonial_vs_resource_beta_after_flora":conditional_controls,
            "native_resource_negative_control":native_controls,
        },
        "source_info":{
            "accepted_regions":len(ids),"matched_region_pairs":len(rows),
            "source_butterfly_beta_pairs_sha256":sha256(beta_path),
            "source_colonial_crosswalk_sha256":sha256(colonial_path),
            "source_look_up_sha256":sha256(look_path),
            "source_region_richness_sha256":sha256(regional_path),
            "lock_sha256":sha256(lock_path),
            "shared_any_empire_pairs":int(np.sum(binary>0)),
            "colonial_exposure_varies":float(np.std(colonial)),
            "paired_data_missingness_is_structural":True,
            "eligible_regions_by_continent":dict(Counter(continents)),
        }
    }


def analysis(data,reps=499,seed=31031):
    if reps<99:
        raise ValueError("Minimum 99 exploratory node randomizations")
    x=data["exposure"]
    ia,ib=data["ia"],data["ib"]
    years=data["years"]
    continent=data["continent"]
    base=data["control_sets"]
    responses={
      "colonial_vs_plant_beta":data["plant_delta"],
      "colonial_vs_resource_beta":data["resource_delta"],
      "colonial_vs_resource_beta_after_flora":data["resource_delta"],
      "native_resource_negative_control":data["resource_native"],
    }
    params={}
    for name,y in responses.items():
        basis=orthogonal_controls(base[name])
        yr=residual(y,basis)
        xr=residual(x,basis)
        params[name]={"observed":correlation(xr,yr),
                      "basis":basis,"response_residual":yr}
    grouping=defaultdict(list)
    for i,c in enumerate(continent):
        grouping[str(c)].append(i)
    fixed_groups=[indices for _,indices in sorted(grouping.items())]
    rng=np.random.default_rng(seed)
    null={name:[] for name in responses}
    history_min=25
    for _ in range(reps):
        perm=years.copy()
        for group in fixed_groups:
            if len(group)>1:
                perm[group]=years[rng.permutation(group)]
        xp=np.log1p(shared_years(perm[ia],perm[ib],years_min=history_min))
        for name,obj in params.items():
            rp=residual(xp,obj["basis"])
            v=correlation(rp,obj["response_residual"])
            null[name].append(v)
    results={}
    for name,obj in params.items():
        obs=obj["observed"]
        dist=np.asarray([v for v in null[name] if v is not None],float)
        if obs is None or len(dist)<reps*.99:
            results[name]={"status":"UNIDENTIFIABLE_CONSTANT_RANK_OR_NULL"}
            continue
        results[name]={
            "observed_partial_rank_r":obs,
            "null_median":float(np.median(dist)),
            "null_q025":float(np.quantile(dist,.025)),
            "null_q975":float(np.quantile(dist,.975)),
            "positive_direction_tail":float((1+np.sum(dist>=obs))/(len(dist)+1)),
            "two_sided_abs_tail":float((1+np.sum(np.abs(dist)>=abs(obs)))/(len(dist)+1)),
            "null_draws":int(len(dist)),
        }
    return {
        "status":"COLONIAL_LEGACY_REGION_PAIR_EXPLORATORY_ASSOCIATION_NOT_CAUSAL",
        "source_info":data["source_info"],
        "exposure_definition":"sum min(GB,ESP,PRT,NL area-weighted empire durations when each >=25y), log1p",
        "controls": {
            "plant":"region distance, same continent, temperature and aridity difference, native plant similarity, native plant richness min/max",
            "resource":"plant controls plus native resource similarity",
            "conditional_resource":"resource controls plus general-plant composition delta",
            "negative_control":"geography, climate, prior plant similarity/richness",
        },
        "results":results,
        "randomization":{
            "seed":seed,"draws":reps,
            "method":"shuffle full 4-empire region vectors at region nodes within official WGSRPD level1",
            "continent_group_sizes":dict(Counter(continent)),
            "movable_nodes":sum(len(g) for g in fixed_groups if len(g)>1),
            "not_exact_causal_p_values":True,
        },
        "limits":[
            "Region pairs share nodes, and within-continent exchangeability is not known to be exact.",
            "Source histories may concern distinct colonial dates and do not represent observed bilateral shipping.",
            "Native/current WCVP is a retrospective spatial comparison, not observed colonial-time change.",
            "Lenzner GloNAF, WCVP and HOSTS have distinct observation and taxonomic biases.",
            "Independent GBIF butterfly composition is a separate, not-yet-executed validation.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--lock",type=Path,default=LOCK)
    p.add_argument("--crosswalk",type=Path,default=COLONIAL)
    p.add_argument("--butterfly-pairs",type=Path,default=BETA)
    p.add_argument("--region-richness",type=Path,default=REGIONAL)
    p.add_argument("--official-lookup",type=Path,default=LOOKUP)
    p.add_argument("--draws",type=int,default=499)
    p.add_argument("--output",type=Path,default=BASE/"colonial_global_resource_test_v1.json")
    args=p.parse_args()
    inputs=load(args.lock,args.crosswalk,args.butterfly_pairs,args.official_lookup,args.region_richness)
    result=analysis(inputs,reps=args.draws)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
