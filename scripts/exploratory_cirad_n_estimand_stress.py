#!/usr/bin/env python3
"""Estimand/denominator stress test for CIRAD 2022 agricultural-N association.

Reuses the exact source-verified 1,596 WPI port botanical counts, CIRAD 3 km
human/cropland N gridded exposure, Jones 2015 modeled wastewater and ESA 2021
WorldCover. Fits several independent SOURCE-OBSERVATION response definitions.

Response is *GBIF-indexed candidate non-native versus native host species*,
not a count of caterpillars or new introductions. The 2022 modeled cropland N
is not actual sewage N delivered to plant soils. Country-level robust SE do
not identify causation or solve spatial/detection heterogeneity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from exploratory_global_wwtp_jones_worldcover import read_source
from exploratory_cirad_n_excretion_ports import models_on_same_coverage

BASE_TERMS = ["log_wwp", "log_wwtp", "abs_lat", "sinlon", "coslon",
              "built_fraction_2021", "crop_fraction_2021"]
N_TERMS = ["log_human_N", "log_crop_N"]


def country_eligibility(d, *, min_ports=5):
    """Country FE support defined before fitting (no result-sign filtering)."""
    z = d.groupby("country").agg(n=("port_id","size"),
         variation_cropland_N=("log_crop_N","nunique"),
         variation_human_N=("log_human_N","nunique"),
         total_alien=("alien_candidates","sum"),
         total_native=("native_candidates","sum"))
    countries = z.index[(z.n >= min_ports) &
         (z.variation_cropland_N >= 2) &
         (z.variation_human_N >= 2) &
         (z.total_alien > 0) & (z.total_native > 0)]
    return d[d.country.isin(countries)].copy()


def add_responses(d):
    z = d.copy()
    a = z.alien_candidates.astype(float)
    n = z.native_candidates.astype(float)
    if (a<0).any() or (n<0).any():
        raise ValueError("Negative GBIF host candidate counts")
    z["total_alien_native"] = a+n
    z["alien_share_observed"] = np.where(a+n>0,a/(a+n).replace(0,np.nan),np.nan)
    z["ratio_log1p"] = np.log1p(a)-np.log1p(n)
    z["ratio_pc_0p1"] = np.log((a+.1)/(n+.1))
    z["ratio_pc_0p5"] = np.log((a+.5)/(n+.5))
    z["ratio_pc_1"] = np.log((a+1)/(n+1))
    z["ratio_pc_2"] = np.log((a+2)/(n+2))
    return z


def fit_response(z,response,*,include_total_proxy,add_n=True,model="OLS"):
    controls=list(BASE_TERMS)
    if include_total_proxy:
        controls.append("log_indexed")
    if add_n:controls += N_TERMS
    formula=response+" ~ "+" + ".join([*controls,"C(country)"])
    if model=="OLS":
        fitted=smf.ols(formula,data=z).fit(
            cov_type="cluster",cov_kwds={"groups":z.country})
    elif model=="FRAC_LOGIT":
        if z[response].isna().any() or not z[response].between(0,1).all():
            raise ValueError("Fractional outcome must be observed and in [0,1]")
        # One equal-weight unit per *port*, no pseudo-independent trials for
        # the many candidate taxa recorded at one port.
        fitted=smf.glm(formula,data=z,family=sm.families.Binomial()).fit(
            cov_type="cluster",cov_kwds={"groups":z.country})
    else:raise ValueError("Unknown model")
    out={"n_ports":int(fitted.nobs),"n_countries":int(z.country.nunique()),
         "response":response,"model":model,"formula":formula,
         "R2":float(fitted.rsquared) if model=="OLS" else None,
         "loglikelihood":float(fitted.llf),
         "cluster_robust_country":True,
         "port_weighting":"ONE_PORT_ONE_OBSERVATION"}
    for name in N_TERMS:
        if name in fitted.params:
            out[name]={"beta":float(fitted.params[name]),
                       "CI95_country_cluster":list(map(float,fitted.conf_int().loc[name])),
                       "nominal_p":float(fitted.pvalues[name])}
    return out


def study(jones,urban,nfeature):
    source=read_source(jones,urban)
    n=pd.read_csv(nfeature,keep_default_na=False)
    if (len(n),n.port_id.nunique())!=(1596,1596):
        raise ValueError("2022 N source 1596 port identities not unique")
    for col in n.columns:
        if col.startswith(("human_","crop_")) or col in ("projected_x","projected_y"):
            n[col]=pd.to_numeric(n[col],errors="coerce")
    earlier,study_data=models_on_same_coverage(source,n,radius=25,urban=True)
    if (earlier["support"]["n_ports"],earlier["support"]["n_countries"])!=(1041,42):
        raise ValueError("Original WPI N source sample changed")
    original=earlier["models"]["both_human_N_and_crop_N"]
    if abs(original["log_crop_N"]["beta"]-.039927909973897)>1e-8:
        raise ValueError("Reference crop N source coefficient mismatch")

    cohorts={
       "original_1041_all_indexed":study_data,
       "positive_candidate_host_total":study_data.loc[(study_data.alien_candidates+
                                                        study_data.native_candidates)>0],
       "five_plus_candidate_host_species":study_data.loc[(study_data.alien_candidates+
                                                           study_data.native_candidates)>=5],
       "both_native_and_alien_recorded":study_data.loc[
            (study_data.alien_candidates>0)&(study_data.native_candidates>0)],
       "land_N_coverage_ge_50pct":study_data.loc[
            (study_data.human_25km_valid_coverage>=.5)&
            (study_data.crop_25km_valid_coverage>=.5)],
    }
    out={
        "schema":"chocho_WPI1596_CIRAD2022_N_alien_native_denominator_sensitivity_20261105_v1",
        "status":"OBSERVATIONAL_PORT_PLANT_INDEX_ESTIMAND_STRESS_TEST_NOT_NUTRIENT_CAUSATION",
        "source_CIRAD_run":41923811745,
        "n_original_WPI_query_complete":1596,
        "n_source_mapped_ports":1041,
        "original_25km_richness_logratio_OLS":original,
        "cohorts":{},
        "caveats":[
           "GBIF-indexed species richness of candidate host plants is a biased catalogue, not independent biological binomial trials or plant abundance.",
           "Fractional GLM is unweighted at port level, not a binomial likelihood treating candidate species as independent samples.",
           "Source environmental predictors have different dates and pixel radii; 2022 modeled crop-N and human excretion are not direct sewage soil-N.",
           "Conditioning on total host-indexed richness contains outcome components and can induce bias; test without separately.",
           "Restricting to positive GBIF record counts also selects on the outcome and may induce collider bias.",
           "Country cluster intervals are not temporal, spatially blocked, or preregistered confirmation.",
           "N log exposures in source-scaled 10 kg N/ha/year model, not true N deposition measured at host plants.",
           "Multiple analytic choices explored and estimates must be interpreted as diagnostic rather than confirmatory."
        ]
    }
    for cname,raw in cohorts.items():
        d=country_eligibility(add_responses(raw))
        if len(d)<350 or d.country.nunique()<15:
            out["cohorts"][cname]={"status":"INSUFFICIENT_SUPPORT","n":len(d),
                                  "countries":int(d.country.nunique())}
            continue
        models={}
        for resp in ("ratio_pc_0p1","ratio_pc_0p5","ratio_pc_1","ratio_pc_2","ratio_log1p"):
            for idx in (True,False):
                tag=resp+("_with_indexed_total" if idx else "_without_indexed_total")
                models[tag]=fit_response(d,resp,include_total_proxy=idx)
        positive=d[d.total_alien_native>0].copy()
        if len(positive)>=300 and positive.country.nunique()>=15:
            for method in ("OLS","FRAC_LOGIT"):
                for idx in (True,False):
                    tag="alien_share_"+method+("_with_indexed_total" if idx else "_without_indexed_total")
                    models[tag]=fit_response(positive,"alien_share_observed",
                                   include_total_proxy=idx,model=method)
        baselines={}
        for response in ("ratio_pc_0p5","ratio_log1p"):
            base=fit_response(d,response,include_total_proxy=True,add_n=False)
            full=models[response+"_with_indexed_total"]
            baselines[response]={"no_2022_N_R2":base["R2"],
                                 "full_2022_N_R2":full["R2"],
                                 "additional_in_sample_R2":full["R2"]-base["R2"]}
        out["cohorts"][cname]={
          "status":"ESTIMATED_COUNTRY_FE_SOURCE_DESCRIPTIVE",
          "ports":len(d),"countries":int(d.country.nunique()),
          "ports_with_zero_both":int((d.total_alien_native==0).sum()),
          "ports_eligible_fractional":int((d.total_alien_native>0).sum()),
          "models":models,"explanatory_increment":baselines
        }
        compact={x:{"crop_beta":v.get("log_crop_N",{}).get("beta"),
                    "crop_ci":v.get("log_crop_N",{}).get("CI95_country_cluster"),
                    "n":v["n_ports"]} for x,v in models.items()
                 if x in ("ratio_pc_0p5_with_indexed_total","ratio_pc_0p1_with_indexed_total",
                          "ratio_pc_2_with_indexed_total","alien_share_OLS_without_indexed_total",
                          "alien_share_FRAC_LOGIT_without_indexed_total")}
        print("CIRAD_N_OUTCOME_ESTIMAND_CHECK",cname,
              json.dumps({"n":len(d),"country_count":d.country.nunique(),
                          "main":compact},sort_keys=True),flush=True)
    return out


def main():
    ap=argparse.ArgumentParser()
    for key in ("jones","urban","N-features","out"):
        ap.add_argument("--"+key,required=True,type=Path)
    args=ap.parse_args()
    result=study(args.jones,args.urban,args.N_features)
    path=args.out
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("CIRAD_2022_CROPLANDN_ESTIMAND_STRESS_TEST_COMPLETE",json.dumps({
        "n_original":result["n_source_mapped_ports"],
        "cohorts":{k:{"n":v.get("ports"),"countries":v.get("countries"),
            "frac_logit":v.get("models",{}).get("alien_share_FRAC_LOGIT_without_indexed_total",{}).get("log_crop_N"),
            "ratio_pc_05":v.get("models",{}).get("ratio_pc_0p5_with_indexed_total",{}).get("log_crop_N")}
            for k,v in result["cohorts"].items()}},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
