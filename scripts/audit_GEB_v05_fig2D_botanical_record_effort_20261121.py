#!/usr/bin/env python3
"""Fig 2D isolation x source-registered potential host resources, all-plant indexed effort.
Uses original reproducible 126-island source and exactly 98 successfully indexed
botanical regions. This is sensitivity to a *proxy*, not an observation model.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t as t_dist
import audit_island_port_isolation_239_20261120 as original

EFFORT_BLOB="23f5e5c7f7c875af41ccb5c6094cb4f48f1071fe"
def blob(path):
    b=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+bytes([0])+b).hexdigest()

def check(cond,description):
    if not cond:raise AssertionError(description)

def standardized(values):
    z=np.asarray(values,dtype=float)
    return (z-z.mean())/z.std(ddof=1)

def extract_effect(data,X,response):
    check(np.linalg.matrix_rank(X.to_numpy())==X.shape[1],"Rank-deficient Fig2D design")
    model=sm.OLS(response,X).fit(cov_type="HC3")
    cluster=sm.OLS(response,X).fit(cov_type="cluster",cov_kwds={"groups":data.level1,"use_correction":True})
    degrees=data.level1.nunique()-1
    check(degrees==8,"Unexpected nine TDWG Level1 clusters")
    out={}
    for v in ("distance","ports"):
        b=float(model.params[v]);se=float(cluster.bse[v]);tc=float(t_dist.ppf(.975,degrees))
        out[v]={
          "beta":b,
          "HC3_p":float(model.pvalues[v]),
          "HC3_95":[float(x) for x in model.conf_int().loc[v]],
          "nine_level1_t8_p":float(2*t_dist.sf(abs(b/se),degrees)),
          "nine_level1_t8_95":[b-tc*se,b+tc*se],
        }
    return out

def run_design(frame,effort):
    X=original.design(frame)
    if effort:
        e=standardized(np.log1p(frame[effort].astype(float)))
        check(np.isfinite(e).all(),"Invalid effort covariate")
        X["botanical_index_effort"]=e
    return X

def model_one(frame,y,effort):
    return extract_effect(frame,run_design(frame,effort),y)

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--port-flora",default="results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv")
    a.add_argument("--region-hosts",default="results/island-mainland-resource-q1/region_metrics.csv")
    a.add_argument("--tdwg-geojson",required=True)
    a.add_argument("--effort-table",default="results/three-gate-port-botanical-effort/botanical_effort_flora_island_inputs.csv")
    a.add_argument("--out-dir",default="results/geb_v05_fig2d_effort_sensitivity_20261121")
    z=a.parse_args()
    check(blob(z.effort_table)==EFFORT_BLOB,"Frozen original plant effort table changed")
    f=original.make_frame(z.port_flora,z.region_hosts,z.tdwg_geojson)
    check(len(f)==126 and int(f.native_host_239.sum())==4829,"Not original frozen 126 source")
    effort=pd.read_csv(z.effort_table,dtype={"code":str,"level1":str})
    check(len(effort)==98 and effort.code.is_unique,"Expected source recorded-plant 98 island subset")
    cols=["code","level1","gbif_all_plant_indexed_records","gbif_pre1900_plant_records"]
    s=f.merge(effort[cols],on="code",how="inner",validate="one_to_one",suffixes=("","_effort"))
    check(len(s)==98 and s.code.is_unique,"Sample selection differs from original 98 matched effort cohort")
    check((s.level1.astype(str)==s.level1_effort.astype(str)).all(),"Mismatch region")
    check(s.level1.nunique()==9,"Original coverage of 9 Level1 groups changed")
    check(np.isfinite(s[["gbif_all_plant_indexed_records","gbif_pre1900_plant_records"]]).all().all(),"Invalid effort")
    y_native=standardized(np.log1p(s.native_host_239))
    y_current=standardized(np.log1p(s.current_host_239))
    native_cur=y_current-y_native
    active=s.loc[s.current_host_239>0].copy()
    logodds=np.log((active.introduced_only_host_239+0.5)/(active.native_host_239+0.5))
    targets={"native_host_z":(s,y_native),"contemporary_host_z":(s,y_current),
             "paired_current_minus_native_z":(s,native_cur),
             "introduced_only_share_logodds":(active,logodds)}
    models={}
    for name,(d,y) in targets.items():
        models[name]={
           "same_restricted_cohort_no_effort":model_one(d,y,None),
           "adjusted_all_plant_indexed_records":model_one(d,y,"gbif_all_plant_indexed_records"),
           "adjusted_pre1900_plant_indexed_records":model_one(d,y,"gbif_pre1900_plant_records"),
        }
        models[name]["n_islands"]=len(d)
    output={
      "schema":"chocho_geb_v05_fig2D_98_island_botanical_index_effort_source_20261121_v1",
      "status":"SOURCE_FIXED239_126_ISLANDS;EFFORT_SENSITIVITY_ONLY_98_COVERED_NOT_MAIN_126_CONFIRMATION",
      "source":{"wcvp_port_flora_git_blob":original.FLORA_BLOB,"host_resource_git_blob":original.HOST_REGION_BLOB,
                "GBIF_all_Plant_records_98_git_blob":EFFORT_BLOB,
                "tdwg_original_geometry_sha256":original.sha(z.tdwg_geojson)},
      "cohorts":{"all_source_island_regions":126,"matched_effort_island_regions":len(s),
          "missing_botanical_effort_islands":126-len(s),"active_matching_islands":len(active),
          "level1_regions":int(s.level1.nunique()),
          "excluded_codes":sorted(set(f.code)-set(s.code))},
      "models":models,
      "caveats":[
        "GBIF Plantae indexing over 1700–2026 proxies recording and provider/era intensity; not independent sampling effort and not a WCVP inventory calibration.",
        "The 98 observed islands are a non-random subset of 126: prior source shows excluded botanical islands generally larger and more port-rich; never impute failed GBIF queries as zero effort.",
        "All three fits for each response use the SAME 98-island source panel (or identical 98-based resource-active subset); 126-island coefficient may NOT be compared to a 98-island effort coefficient as if only adjustment changed.",
        "Nine TDWG Level1 geographic clusters give small-cluster t8 sensitivity. Statistical attenuation or persistence does not establish causality.",
        "Potential butterfly known-host presence is NOT observed butterfly presence, host availability quality, historical ship transportation or population rescue.",
      ]
    }
    out=Path(z.out_dir);out.mkdir(parents=True,exist_ok=True)
    (out/"FIG2D_PLANT_INDEX_EFFORT_SENSITIVITY_20261121.json").write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n")
    lines=["# Figure 2D: botanical observation effort sensitivity","",
           "Same original 239-butterfly/HOSTS-WCVP panel, 126 total botanical islands; 98 have GBIF whole-Plantae botanical recording effort.",
           "**98 records are nonrandom. Missing 28 are not zero recording effort.**","",
           "| Outcome | n | No-effort distance β (t8 p) | All-indexed effort distance β (t8 p) | Pre-1900 effort distance β (t8 p) |",
           "|---|---:|---:|---:|---:|"]
    for k,row in models.items():
        rs=[]
        for key in ("same_restricted_cohort_no_effort","adjusted_all_plant_indexed_records","adjusted_pre1900_plant_indexed_records"):
            est=row[key]["distance"]
            rs.append(f"{est['beta']:+.4f} (p={est['nine_level1_t8_p']:.4g})")
        lines.append(f"| {k} | {row['n_islands']} | "+" | ".join(rs)+" |")
    lines.extend(["","These are **observational source sensitivities**. Island-native and contemporary-response models are separately standardized on the same 98 islands. Neither formal mediation nor causal nutrient/plant effects are inferred.",""])
    (out/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print(json.dumps(output,indent=2))
if __name__=="__main__":main()
