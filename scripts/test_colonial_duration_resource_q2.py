#!/usr/bin/env python3
"""Q2-0: colonial exposure duration vs island potential host-resource gain.

Published four-empire histories are area-overlaid with WGSRPD3; this is an
observational screening, not a historical introduction mechanism test.
"""
from __future__ import annotations
import argparse, csv, json, math, random
from collections import defaultdict
from pathlib import Path
import numpy as np

EMPIRES=("Great_Britain","Spain","Portugal","Netherlands")

def readcsv(path):
    with path.open(encoding="utf-8",newline="") as h: yield from csv.DictReader(h)

def fit(rows,ys):
    levels=sorted({r["level1"] for r in rows})
    X=np.column_stack([
         np.ones(len(rows)),
         np.array([r["log_legacy"] for r in rows]),
         np.array([math.log(r["area_km2"]) for r in rows]),
         *[np.array([int(r["level1"]==level) for r in rows]) for level in levels[1:]],
    ])
    if np.linalg.matrix_rank(X)!=X.shape[1]:return None
    return float(np.linalg.lstsq(X,np.asarray(ys),rcond=None)[0][1])

def within_permutation(rows,group_by,draws=4999):
    obs=fit(rows,[r["share"] for r in rows])
    inds=defaultdict(list)
    for i,k in enumerate(group_by):inds[k].append(i)
    permutable=[v for v in inds.values() if len({rows[i]["log_legacy"] for i in v})>1]
    effective=len({i for group in permutable for i in group})
    rng=random.Random(20261102)
    null=[]
    for _ in range(draws):
        exposure=[r["log_legacy"] for r in rows]
        for ids in permutable:
            vs=[exposure[i] for i in ids];rng.shuffle(vs)
            for i,value in zip(ids,vs):exposure[i]=value
        perm_rows=[{**r,"log_legacy":exposure[i]} for i,r in enumerate(rows)]
        n=fit(perm_rows,[r["share"] for r in perm_rows])
        if n is not None:null.append(n)
    if obs is None or not null:return {"status":"UNIDENTIFIED","observed":obs,"effective_islands":effective}
    null.sort()
    m=sum(null)/len(null)
    return {"status":"DESCRIPTIVE_PERMUTATION","observed_slope":obs,
            "null_mean":m,
            "observed_minus_null":obs-m,
            "p_high":(1+sum(x>=obs-1e-12 for x in null))/(len(null)+1),
            "p_two_sided":(1+sum(abs(x-m)>=abs(obs-m)-1e-12 for x in null))/(len(null)+1),
            "q025":null[int(.025*(len(null)-1))],"q975":null[int(.975*(len(null)-1))],
            "n_draws":len(null),"effective_islands":effective,
            "effective_blocks":len(permutable)}

def main():
    p=argparse.ArgumentParser()
    for k in ("region-metrics","area-metrics","empire-crosswalk","output"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    reg={r["area_code_l3"]:r for r in readcsv(a.region_metrics) if r["class"]=="island"}
    area={r["code"]:r for r in readcsv(a.area_metrics)}
    hist={r["area_code_l3"]:r for r in readcsv(a.empire_crosswalk)}
    if len(reg)!=125 or len(hist)!=125:raise RuntimeError("Changed fixed island cohort")
    rows=[]
    for code,record in sorted(reg.items()):
        if int(record["contemporary"])<=0 or code not in hist:continue
        h=hist[code];covered=float(h["lenzner_polygon_coverage"])
        if covered<0.8:continue
        years={e:float(h[e+"_weighted_exposure_years"]) for e in EMPIRES}
        maxyrs=max(years.values())
        rows.append({"code":code,"level1":record["level1"],
                    "area_km2":float(area[code]["area_km2"]),
                    "share":int(record["added"])/int(record["contemporary"]),
                    "years":years,"maxyears":maxyrs,
                    "log_legacy":math.log1p(max(0,maxyrs))})
    if len(rows)<30:raise RuntimeError("Not enough published empire geospatial coverage")
    variants={"primary_all_verified":[r for r in rows],
              "excluding_europe":[r for r in rows if r["level1"]!="1"],
              "positive_history_only":[r for r in rows if r["maxyears"]>=50]}
    result={}
    for name, subset in variants.items():
        by_level=[r["level1"] for r in subset]
        median_area=np.median([r["area_km2"] for r in subset])
        by_area=[r["level1"]+"|"+("small" if r["area_km2"]<=median_area else "large") for r in subset]
        result[name]={
          "islands":len(subset),
          "with_empire_50plus_years":sum(r["maxyears"]>=50 for r in subset),
          "mean_introduced_share":float(np.mean([r["share"] for r in subset])),
          "area_median_km2":float(median_area),
          "by_biogeographic_level1":within_permutation(subset,by_level),
          "by_level1_and_area_half":within_permutation(subset,by_area),
        }
    payload={
      "schema":"chocho_colonial_duration_island_resource_screen_v0.1",
      "source":"Lenzner et al 2022 empire history polygons and Q1 exact 239 butterfly host resource reconstruction",
      "status":"OBSERVATIONAL_COLONIAL_EXPOSURE_SCREEN",
      "exposure":"log1p(max_4_empire_area_weighted_years) from published GloNAF, not age since introduction",
      "response":"fraction of contemporary potential butterfly x WGSRPD3 resource opportunity newly supplied by introduced hosts",
      "area_control":"ln(WGSRPD3 area km2) and Level1 biogeographic fixed effects in ordinary least squares",
      "variants":result,
      "caution":[
         "No colonial empire among four does not imply no colonization (France and other empires absent).",
         "Correlations do not demonstrate historical plant imports, ship arrivals, or resident butterflies.",
         "Unmeasured climate and historic human factors, imbalanced islands and area/geospatial uncertainty remain.",
      ],
    }
    a.output.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2,ensure_ascii=False))
if __name__=="__main__":main()
