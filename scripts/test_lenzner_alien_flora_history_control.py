#!/usr/bin/env python3
"""Test a POSITIVE CONTROL before attributing colonial history to butterflies.

Lenzner et al. (2022) study's own naturalized-alien flora × region lists and
four empire-duration predictors, against within-UN-geographic-group node
permutations. All inputs from MD5-pinned Zenodo 6640264 original ZIP.

This check is deliberately distinct from the new 369-WGSRPD3 native/current
plant and butterfly-potential-resource response. Failure to recover a
colonial-alien flora signal under the chosen simple model is a METHODS/SCOPE
warning, not a biological null for colonial effects on butterflies.

Outputs only analysis aggregates; no original third-party data rows or site
geometries are republished.
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

import numpy as np
import zipfile

ZIP = Path("_temp_historical_sources/lenzner_6640264.zip")
EXPECTED_MD5 = "16f5585d2500676951a4b664bbf0ffde"
OUT = Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/"
           "colonial_alien_flora_positive_control_v1.json")
EMPIRES = ("time.emp.GB", "time.emp.ESP", "time.emp.PRT", "time.emp.NL")
R_EARTH = 6371.0088


def canonical_id(raw):
    if raw is None or str(raw).strip() == "":
        return None
    t = str(raw).strip()
    try:
        return str(int(float(t)))
    except ValueError:
        return t


def read_zip_csv(z, name):
    matched = [v for v in z.namelist() if v.endswith("/Datasets/" + name)]
    if len(matched) != 1:
        raise RuntimeError(f"Cannot uniquely identify published {name}")
    with z.open(matched[0]) as reader:
        # Preserve source numeric ID and ASCII scientific names; malformed
        # diacritics in unused site labels are replaced, not reinterpreted.
        raw = reader.read().decode("utf-8-sig", errors="replace")
    table = list(csv.DictReader(io.StringIO(raw)))
    if not table or "OBJIDsic" not in table[0]:
        raise RuntimeError("Lenzner original column schema changed")
    return table


def verify_archive(path):
    md5 = hashlib.md5()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1<<20), b""):
            md5.update(chunk)
    if md5.hexdigest() != EXPECTED_MD5:
        raise RuntimeError("Lenzner ZIP source checksum mismatch")


def load_nodes(path=ZIP):
    verify_archive(path)
    with zipfile.ZipFile(path) as z:
        predictors = read_zip_csv(z, "Predictors.csv")
        regions = read_zip_csv(z, "Regions.csv")
        glonaf = read_zip_csv(z, "GloNAF.csv")
    regions_by_id = {}
    for r in regions:
        k = canonical_id(r["OBJIDsic"])
        if not k or k in regions_by_id:
            raise RuntimeError("Duplicate/missing Lenzner site identifier")
        regions_by_id[k] = r
    if len(regions_by_id) != 1183:
        raise RuntimeError("Lenzner region source coverage changed")
    preds = {}
    invalid_pred = Counter()
    for r in predictors:
        k = canonical_id(r["OBJIDsic"])
        if not k or k not in regions_by_id:
            invalid_pred["missing_site_key"] += 1
            continue
        if k in preds:
            raise RuntimeError("Duplicate original predictor site")
        try:
            vals = [float(r[p]) for p in EMPIRES]
            lat, lon = float(r["LAT"]), float(r["LON"])
            temp, arid = float(r["temp.mean"]), float(r["aridity.mean"])
            area = float(r["area"])
        except (ValueError, KeyError, TypeError):
            invalid_pred["non_numeric_covariate"] += 1
            continue
        if not (all(math.isfinite(v) and 0<=v<=600 for v in vals)
                and math.isfinite(lat) and math.isfinite(lon)
                and -90<=lat<=90 and -180<=lon<=180
                and math.isfinite(temp) and math.isfinite(arid)
                and math.isfinite(area) and area>0):
            invalid_pred["invalid_range"] += 1
            continue
        geo = (regions_by_id[k].get("UN_geospat") or "").strip()
        preds[k] = {
            "lat":lat, "lon":lon, "temp":temp,"arid":arid,"area":area,
            "empire":tuple(vals),"geounit":geo or "UNRESOLVED",
            "island":str(r.get("island","")).lower() in ("1","true","yes"),
        }
    flora = defaultdict(set)
    keys = set()
    bad = Counter()
    for r in glonaf:
        key = canonical_id(r["OBJIDsic"])
        if key not in preds:
            bad["outside_matched_predictors"] += 1
            continue
        genus = (r.get("genus") or "").strip()
        epithet = (r.get("Species") or "").strip()
        if not genus or not epithet or "\ufffd" in (genus+epithet):
            bad["missing_or_corrupt_taxon"] += 1
            continue
        # Avoid ascribing distinct plant species to the same synonym by case.
        taxon = epithet if epithet.lower().startswith(genus.lower()+" ") else genus+" "+epithet
        if len(taxon)<5 or " sp." in taxon.lower():
            bad["not_species_identifiable"] += 1
            continue
        flora[key].add(taxon.casefold())
        keys.add(taxon.casefold())
    suitable = [k for k in sorted(preds) if len(flora[k])>=10]
    if len(suitable) < 300:
        raise RuntimeError("Too few GloNAF regions with observed alien plants")
    if len(keys) < 1000:
        raise RuntimeError("Too few identifiable alien taxa")
    tax_index = {taxon: i for i,taxon in enumerate(sorted(keys))}
    nodes = []
    for key in suitable:
        mask = 0
        for taxon in flora[key]:
            mask |= 1 << tax_index[taxon]
        node = dict(preds[key])
        node["source_key"] = key
        node["plant_bits"] = mask
        node["alien_species"] = mask.bit_count()
        nodes.append(node)
    return nodes, {
        "source_regions":len(regions),
        "source_predictor_rows":len(predictors),
        "source_glonaf_rows":len(glonaf),
        "source_taxa_after_validated_name":len(keys),
        "eligible_regions":len(nodes),
        "alien_plant_species_records_eligible":sum(x["alien_species"] for x in nodes),
        "unmatched_predictor_or_plant_rows":dict(bad),
        "invalid_predictors":dict(invalid_pred),
        "original_source_zip_md5":EXPECTED_MD5,
    }


def pair_sample(ids, maximum=90000, seed=31031):
    """Hash-defined, response-blind sample of node pairs (both unique IDs).

    The sample and its size cannot change based on alien plant resemblance or
    shared colonial history. Full set may contain hundreds of thousands of
    dependent pairs; node-wise permutations account for repeated endpoints.
    """
    if len(set(ids))!=len(ids):
        raise ValueError("Duplicate region keys")
    total=len(ids)*(len(ids)-1)//2
    n_target=min(total, maximum)
    probability=n_target / total
    sampled=[]
    for i in range(len(ids)):
        for j in range(i+1,len(ids)):
            raw=f"{seed}|{ids[i]}|{ids[j]}".encode("utf8")
            value=int.from_bytes(hashlib.blake2b(raw,digest_size=8).digest(),"big")
            if value / (2**64) < probability:
                sampled.append((i,j))
    if len(sampled)<10000:
        raise RuntimeError("Too few deterministic pairs, check source coverage")
    return np.array(sampled,np.int32),total


def ranks(v):
    ordered=np.argsort(np.asarray(v),kind="stable")
    x=np.asarray(v)
    out=np.empty(len(x),float)
    a=0
    while a<len(x):
        b=a+1
        while b<len(x) and x[ordered[b]]==x[ordered[a]]:
            b+=1
        out[ordered[a:b]]=(a+1+b)/2
        a=b
    return out


def residual(v,basis):
    x=ranks(v)
    return x-basis@(basis.T@x)


def orthogonal(cols):
    z=np.column_stack([np.ones(len(cols[0]))]+[ranks(x) for x in cols])
    u,s,_=np.linalg.svd(z,full_matrices=False)
    return u[:,s>max(s)*1e-10]


def correlation(a,b):
    x=a-a.mean()
    y=b-b.mean()
    den=math.sqrt(float(x@x)*float(y@y))
    return float(x@y/den) if den>0 else None


def years_shared(a,b,min_years=25):
    aa=np.asarray(a,float)
    bb=np.asarray(b,float)
    keep=(aa>=min_years)&(bb>=min_years)
    return np.sum(np.minimum(aa,bb)*keep,axis=-1)


def analysis(nodes, maximum=90000, permutations=199, seed=31031):
    if permutations<99:
        raise ValueError(">=99 node permutations required")
    n=len(nodes)
    pairs,total=pair_sample([v["source_key"] for v in nodes],maximum,seed)
    i,j=pairs[:,0],pairs[:,1]
    y=np.asarray([
        ((nodes[int(a)]["plant_bits"] & nodes[int(b)]["plant_bits"]).bit_count()/
         (nodes[int(a)]["plant_bits"] | nodes[int(b)]["plant_bits"]).bit_count())
        for a,b in pairs
    ])
    years=np.asarray([n["empire"] for n in nodes],float)
    geo=np.asarray([n["geounit"] for n in nodes],object)
    lat=np.radians([n["lat"] for n in nodes])
    lon=np.radians([n["lon"] for n in nodes])
    t=np.asarray([n["temp"] for n in nodes],float)
    dry=np.asarray([n["arid"] for n in nodes],float)
    area=np.asarray([n["area"] for n in nodes],float)
    rich=np.asarray([n["alien_species"] for n in nodes],float)
    dlat=lat[i]-lat[j]
    dlon=lon[i]-lon[j]
    h=np.sin(dlat/2)**2+np.cos(lat[i])*np.cos(lat[j])*np.sin(dlon/2)**2
    distance=np.log1p(2*R_EARTH*np.arcsin(np.sqrt(np.clip(h,0,1))))
    same=(geo[i]==geo[j]).astype(float)
    control=[distance,same,np.abs(t[i]-t[j]),np.abs(dry[i]-dry[j]),
             np.log1p(np.minimum(area[i],area[j])),
             np.log1p(np.maximum(area[i],area[j])),
             np.log1p(np.minimum(rich[i],rich[j])),
             np.log1p(np.maximum(rich[i],rich[j]))]
    basis=orthogonal(control)
    yr=residual(y,basis)
    shared=np.log1p(years_shared(years[i],years[j]))
    xr=residual(shared,basis)
    raw=correlation(ranks(shared),ranks(y))
    partial=correlation(xr,yr)
    groups=defaultdict(list)
    for k,val in enumerate(geo):
        groups[str(val)].append(k)
    fixed=[v for v in groups.values() if len(v)>1]
    movable=sum(len(g) for g in fixed)
    if movable<n*.4:
        raise RuntimeError("Historic source groups too small for geographically controlled node permutations")
    rng=np.random.default_rng(seed)
    null=[]
    for draw in range(permutations):
        perm=np.array(years,copy=True)
        for group in fixed:
            perm[group]=years[rng.permutation(group)]
        xx=np.log1p(years_shared(perm[i],perm[j]))
        value=correlation(residual(xx,basis),yr)
        if value is None:
            raise RuntimeError("Unidentifiable permuted colonial exposure")
        null.append(value)
    nul=np.asarray(null)
    source_sites_shared=int(np.sum(shared>0))
    return {
        "status":"PUBLISHED_STUDY_ALIEN_FLORA_POSITIVE_CONTROL_NOT_NEW_DISCOVERY",
        "sample":{"regions":n,"total_possible_pairs":total,
                  "hashed_response_blind_pairs":len(pairs),
                  "shared_colonial_pairs":source_sites_shared,
                  "shared_colonial_pair_fraction":source_sites_shared/len(pairs),
                  "source_geo_groups":{k:len(v) for k,v in sorted(groups.items())},
                  "movable_source_regions":movable},
        "exposure":"log1p(sum minimum of shared GB ESP PRT NL duration years, both>=25)",
        "response":"Jaccard of published GloNAF recorded naturalized alien plants, original study sites",
        "raw_spearman_rho":raw,
        "distance_climate_area_alien_richness_adjusted_partial_rank_rho":partial,
        "node_permutation":{
            "draws":permutations,"seed":seed,
            "scheme":"shuffle full four-empire vectors within original source UN_geospat group",
            "null_median":float(np.median(nul)),
            "null_q025":float(np.quantile(nul,.025)),
            "null_q975":float(np.quantile(nul,.975)),
            "upper_tail":(1+int(np.sum(nul>=partial)))/(len(nul)+1),
            "two_sided_abs_tail":(1+int(np.sum(abs(nul)>=abs(partial))))/(len(nul)+1),
            "not_a_causal_p_value":True,
        },
        "limitations":[
          "Reanalysis of previously published GloNAF alien-plant data is a POSITIVE CONTROL, never a novel butterfly outcome.",
          "UN_geospat source strata may be finer than biogeographic continents; exchangeability not proven.",
          "Source plant indexing depends on observation effort and original taxonomy.",
          "Empire duration does not establish overlapping calendar years, trade route or plant transport.",
          "Hashed node pairs reuse source regions and cannot be treated as iid samples.",
          "Native WCVP plant/historical butterfly differences are not tested by this source-only positive control.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",type=Path,default=ZIP)
    p.add_argument("--max-pairs",type=int,default=90000)
    p.add_argument("--draws",type=int,default=199)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    nodes,source=load_nodes(args.zip)
    result=analysis(nodes,args.max_pairs,args.draws)
    result["source_data_audit"]=source
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                           encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k!="source_data_audit"},
                     ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
