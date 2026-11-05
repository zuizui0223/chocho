#!/usr/bin/env python3
"""Offline CIRAD 2022 human N vs cropland N robustness among GBIF-indexed ports.

Source-grain: 3km modelled excretion/supply (NOT effluent delivered to plants);
25/10km buffers at WPI2019 ports; outcome pooled 2001–24 host catalogue.
No network access, files from verified historical GitHub Actions source runs.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
from exploratory_global_wwtp_jones_worldcover import read_source
from exploratory_cirad_n_excretion_ports import models_on_same_coverage

def estimate(df, *, outcome="ratio", with_N=True, with_total=True):
    controls=["log_wwp","log_wwtp","abs_lat","sinlon","coslon"]
    if with_total and outcome=="ratio":controls+=["log_indexed"]
    controls+=["built_fraction_2021","crop_fraction_2021"]
    if with_N:controls+=["log_human_N","log_crop_N"]
    form=outcome+" ~ "+" + ".join([*controls,"C(country)"])
    fit=smf.ols(form,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
    d={"n":len(df),"countries":int(df.country.nunique()),"R2":float(fit.rsquared),
       "formula":form,"status":"DESCRIPTIVE_COUNTRY_CLUSTER_FE_NOT_SEWAGE_CAUSATION"}
    for k in ("log_human_N","log_crop_N","log_wwp","log_wwtp"):
        if k in fit.params:
            d[k]={"beta":float(fit.params[k]),
                  "CI95":list(map(float,fit.conf_int().loc[k])),
                  "p_nominal":float(fit.pvalues[k])}
    return d

def allowed(df):
    a=df.groupby("country").agg(n=("port_id","size"),crop=("log_crop_N","nunique"),
                                human=("log_human_N","nunique"))
    ok=a.index[(a.n>=5)&(a.crop>=2)&(a.human>=2)]
    return df[df.country.isin(ok)].copy()

def run(jones,urban,source):
    src=read_source(jones,urban)
    features=pd.read_csv(source,keep_default_na=False)
    if (len(features),features.port_id.nunique())!=(1596,1596):
        raise ValueError("Source WPI port identity mismatch")
    # A saved CSV stores missing raster cells as empty strings. keep_default_na=False
    # preserves source ISO2 'NA' but requires explicit float parsing for N layers.
    for col in features:
        if col.startswith(("human_","crop_")) or col in ("projected_x","projected_y"):
            features[col]=pd.to_numeric(features[col],errors="coerce")
    r25,p25=models_on_same_coverage(src,features,radius=25,urban=True)
    r10,p10=models_on_same_coverage(src,features,radius=10,urban=True)
    if (len(p25),p25.country.nunique(),len(p10),p10.country.nunique())!=(1041,42,1020,40):
        raise ValueError("Source 25km and 10km fixed cohorts changed")
    original=r25["models"]["both_human_N_and_crop_N"]
    if abs(original["log_crop_N"]["beta"]-.039927909973897)>1e-10:
        raise ValueError("Source regression not reproduced")
    out={
      "schema":"chocho_CIRAD2022_local_N_WPI1596_robustness_20261105_v1",
      "status":"EXPLORATORY_MODELLED_NUTRIENT_POTENTIAL_ASSOCIATION_ONLY",
      "source_CIRAD_run":41923811745,
      "source_port_table_sha256":hashlib.sha256(Path(source).read_bytes()).hexdigest(),
      "matched_1041ports_42countries":True,
      "original_model":original,
      "sensitivity":{},
      "negative_controls":{},
      "same_sample_25_vs_10km":{},
      "leave_one_country_out":{},
      "correlations":{}
    }
    study={
      "primary_all_ports":p25,
      "both_native_and_alien_present":p25[(p25.alien_candidates>0)&(p25.native_candidates>0)],
      "at_least_five_candidate_hosts":p25[(p25.alien_candidates+p25.native_candidates)>=5],
      "min_half_valid_N_source_pixels":p25[
           (p25.human_25km_valid_coverage>=.50)&(p25.crop_25km_valid_coverage>=.50)],
      "crop_N_top_one_percent_excluded":p25[p25.log_crop_N<=p25.log_crop_N.quantile(.99)]
    }
    for name,raw in study.items():
        d=allowed(raw)
        if len(d)<350 or d.country.nunique()<15:
            out["sensitivity"][name]={"status":"INSUFFICIENT","n":len(d),"countries":d.country.nunique()}
            continue
        full=estimate(d);bare=estimate(d,with_N=False);no_total=estimate(d,with_total=False)
        out["sensitivity"][name]={"n":len(d),"countries":int(d.country.nunique()),
             "joint_N":full,"no_N_same_port_cohort":bare,
             "delta_in_sample_R2_N_added":full["R2"]-bare["R2"],
             "without_indexed_total_proxy":no_total}
        print("CIRAD_N_ROBUSTNESS",name,json.dumps({"n":len(d),"crop":full.get("log_crop_N"),
                "human":full.get("log_human_N"),"crop_no_total":no_total.get("log_crop_N")}),flush=True)
    for response in ("alien","native"):
        out["negative_controls"][response]=estimate(p25,outcome=response)
    # The full 10 km sample is smaller, so keep exactly the same port IDs.
    ids=set(p25.port_id)&set(p10.port_id)
    a=allowed(p25[p25.port_id.isin(ids)])
    b=allowed(p10[p10.port_id.isin(ids)])
    common=set(a.port_id)&set(b.port_id)
    for r,d in (("25km",a),("10km",b)):
        same=d[d.port_id.isin(common)].copy()
        out["same_sample_25_vs_10km"][r]=estimate(same)
    out["same_sample_25_vs_10km"]["identical_port_count"] = len(common)
    for country in sorted(p25.country.unique()):
        sample=p25[p25.country!=country]
        f=estimate(sample)
        out["leave_one_country_out"][country]={
            "crop_beta":f["log_crop_N"]["beta"],"human_beta":f["log_human_N"]["beta"]}
    crop=[v["crop_beta"] for v in out["leave_one_country_out"].values()]
    out["leave_one_country_out_crop_beta_range"]=[min(crop),max(crop)]
    out["leave_one_country_out_crop_beta_all_positive"]=all(x>0 for x in crop)
    for pair in (("log_crop_N","crop_fraction_2021"),("log_crop_N","log_wwp"),
                 ("log_human_N","log_wwp"),("log_crop_N","log_human_N")):
        out["correlations"]["_vs_".join(pair)]=float(spearmanr(p25[pair[0]],p25[pair[1]]).statistic)
    out["limitations"]=[
        "Human N raster is national FAOSTAT food-balance protein/excretion × GlobPOP population, not measured sewage discharge into terrestrial port plants.",
        "Cropland N raster is modeled cropland nutrient supply/fertilizer/manure/fixation, not observed soil nitrogen.",
        "Raster values are 10 kg N per ha per year, 3km Equal Earth cells around 2022.",
        "Focal hosts are GBIF candidate host species near ports, not butterflies' realized larval feeding or host plant abundance.",
        "Current WPI is 2019, Jones WWp is 2015, satellite urban cover 2021 and plant catalogue pooled 2001-2024.",
        "The indexed-total control contains response components; removing it is a robustness test, not a fix for all observation bias.",
        "Country fixed effects and country-cluster intervals cannot eliminate residual urban, landuse, observability and cultural/historical confounding.",
        "Post hoc multiple model checks; do not claim independent true effect of agricultural N or establish nutrient causation."
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("jones","urban","grid","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    r=run(a.jones,a.urban,a.grid)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("CIRAD_2022_N_SOURCE_ROBUSTNESS_COMPLETE",json.dumps({
        "negative_controls":r["negative_controls"],
        "sample_sensitivities":{k:{"n":v.get("n"),"crop":v.get("joint_N",{}).get("log_crop_N")} for k,v in r["sensitivity"].items()},
        "crop_LOO_range":r["leave_one_country_out_crop_beta_range"],
        "LOO_crop_all_positive":r["leave_one_country_out_crop_beta_all_positive"],
        "same_port25_vs10":r["same_sample_25_vs_10km"],
        "correlations":r["correlations"],
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
