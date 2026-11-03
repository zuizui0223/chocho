#!/usr/bin/env python3
"""Independent HydroWASTE + PANGAEA Jones2015 + WorldCover2021 matched port audit.

Use only frozen successful GitHub Actions artifacts. WWp, WWc, WWt, WWr
are MODELLED million m3/year in 5arc grid neighborhoods, not measured
soil N or in-situ sewage. Outcome is GBIF-indexed relative candidate
butterfly larval-host plant richness, not actual plant or butterfly success.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf


def read_source(jones_path:Path,worldcover_path:Path):
    w=pd.read_csv(worldcover_path,keep_default_na=False)
    j=pd.read_csv(jones_path,keep_default_na=False)
    if (len(w),w.port_id.nunique())!=(1596,1596):
        raise ValueError("Original 1596 WPI + WorldCover source changed")
    if (len(j),j.port_id.nunique())!=(1596,1596):
        raise ValueError("Original 1596 WPI + Jones 2015 source changed")
    for col in ("produced_25km_millionm3yr","collected_25km_millionm3yr",
                "treated_25km_millionm3yr","reused_25km_millionm3yr"):
        if col not in j:raise ValueError("Jones model gridded variable missing "+col)
    # Freeze the host values from both truly independent runs.
    for col in ("country","alien_candidates","native_candidates",
                "all_taxonomically_indexed","WWTP_count_25km"):
        a=j.set_index("port_id")[col]
        b=w.set_index("port_id")[col]
        if col=="country":
            if not a.equals(b):raise ValueError("Country source mismatch")
        elif not np.allclose(a.astype(float).sort_index(),b.astype(float).sort_index()):
            raise ValueError("WWTP/host source mismatch: "+col)
    land=w[["port_id","quality","built_fraction_2021",
            "crop_fraction_2021","tree_fraction_2021",
            "terrestrial_pixels","fraction_nonwater"]].copy()
    matched=j.merge(land,on="port_id",validate="one_to_one")
    if len(matched)!=1596:
        raise ValueError("Original 1596 port grid join incomplete")
    return matched


def model(data,*,outcome="ratio",predictors=(),city_required=False):
    d=data.loc[data.in_within_country_comparison].copy()
    if city_required:d=d.loc[d.quality=="COMPLETE"].copy()
    fields=list(predictors)
    d["abs_lat"]=d.latitude.abs()
    d["sinlon"]=np.sin(np.deg2rad(d.longitude))
    d["coslon"]=np.cos(np.deg2rad(d.longitude))
    d["ln_total"]=np.log1p(d.all_taxonomically_indexed)
    d["ln_wwtp"]=np.log1p(d.WWTP_count_25km)
    d["ln_produced"]=np.log1p(d.produced_25km_millionm3yr)
    d["ln_treated"]=np.log1p(d.treated_25km_millionm3yr)
    d["ln_collected"]=np.log1p(d.collected_25km_millionm3yr)
    d["modeled_not_treated_25km"]=d.produced_25km_millionm3yr-d.treated_25km_millionm3yr
    # Catch incomplete/noncomparable grid source rather than making up negative untreated.
    invalid=(d["modeled_not_treated_25km"]<-.001)
    if invalid.any():
        raise ValueError("Production minus treatment cell-volume inconsistent: negative at source port grid")
    d["ln_not_treated"]=np.log1p(d.modeled_not_treated_25km.clip(lower=0))
    d["modeled_treated_share"]=np.where(
        d.produced_25km_millionm3yr>0,
        d.treated_25km_millionm3yr/d.produced_25km_millionm3yr,
        np.nan)
    d["ratio"]=np.log((d.alien_candidates+.5)/(d.native_candidates+.5))
    d["alien"]=np.log1p(d.alien_candidates)
    d["native"]=np.log1p(d.native_candidates)
    d=d.loc[~d.incomplete_or_missing_grid]
    d=d.dropna(subset=fields)
    # Hold the *same* port set between nested fits for each urban-cohort analysis.
    country=d.groupby("country").agg(n=("port_id","size"),var=("ln_wwtp","nunique"))
    good=set(country.index[(country.n>=5)&(country["var"]>=2)])
    d=d.loc[d.country.isin(good)].copy()
    if len(d)<400 or d.country.nunique()<15:
        return {"status":"INSUFFICIENT_CONDITIONED_COUNTRY_SUPPORT",
                "n":len(d),"countries":int(d.country.nunique())}
    y={"ratio":"ratio","alien":"alien","native":"native"}[outcome]
    controls=["abs_lat","sinlon","coslon","ln_total"]
    if outcome!="ratio":
        # Negative controls native and alien counts reported without
        # conditioning on the sum of both outcome components.
        controls=controls[:-1]
    terms=controls+list(predictors)+["C(country)"]
    f=smf.ols(y+" ~ "+" + ".join(terms),data=d).fit(
        cov_type="cluster",cov_kwds={"groups":d.country})
    result={
      "status":"OBSERVATIONAL_SOURCE_LOCKED_COUNTRY_FE",
      "n_ports":len(d),"n_countries":int(d.country.nunique()),
      "response":outcome,"predictors":list(predictors),
      "controls":controls+["country_FE"],
      "R2":float(f.rsquared)}
    for var in fields:
        low,upper=map(float,f.conf_int().loc[var])
        result[var]={"beta":float(f.params[var]),
                     "cluster_country_95":[low,upper],
                     "nominal_p":float(f.pvalues[var])}
    return result


def run(jones_path:Path,worldcover_path:Path):
    z=read_source(jones_path,worldcover_path)
    if not np.all(z.produced_25km_millionm3yr>=z.treated_25km_millionm3yr-.001):
        raise ValueError("Jones modeled WW production should exceed treatment at the fixed 25km buffer scale")
    if not np.all(z.collected_25km_millionm3yr>=z.treated_25km_millionm3yr-.001):
        raise ValueError("Jones modeled WW collection should exceed treatment at the fixed 25km buffer scale")
    if not np.all(z.treated_25km_millionm3yr>=z.reused_25km_millionm3yr-.001):
        raise ValueError("Jones modeled WW treatment should exceed reuse at the fixed 25km buffer scale")
    summary={
        "schema":"chocho_global_PANGAEA2015_WWTP_worldcover2021_joint_source_port_association_v1",
        "status":"OBSERVATIONAL_PREDICTOR_DECOMPOSITION_NOT_NUTRIENT_CAUSAL_MECHANISM",
        "source_PANGAEA2015_GBIF_ports_run":41563940689,
        "source_ESAWC2021_port_urban_run":41562997326,
        "source_WWTP_port_run":41414513676,
        "source_jones_table_sha256":hashlib.sha256(jones_path.read_bytes()).hexdigest(),
        "source_worldcover_table_sha256":hashlib.sha256(worldcover_path.read_bytes()).hexdigest(),
        "source_years":{"PANGAEA_wastewater":2015,"ESA_WorldCover":2021,
                        "WPI":2019,"indexed_GBIF_hosts":"2001-2024 pooled"},
        "full_WPI_ports":len(z),
        "full_WPI_countries":int(z.country.nunique()),
        "urban_raster_valid_ports":int((z.quality=="COMPLETE").sum()),
        "WWTP_count25_vs_2015_production_Spearman_rho":float(
            spearmanr(z.WWTP_count_25km,z.produced_25km_millionm3yr).statistic),
        "modeled_WW_25km_source_summary":{
            "median_produced_millionm3_per_year":float(z.produced_25km_millionm3yr.median()),
            "median_treated_millionm3_per_year":float(z.treated_25km_millionm3yr.median()),
            "produced_positive_ports":int((z.produced_25km_millionm3yr>0).sum()),
            "treated_positive_ports":int((z.treated_25km_millionm3yr>0).sum()),
            "estimated_model_not_treated_volume_median_millionm3_per_year":
                float((z.produced_25km_millionm3yr-z.treated_25km_millionm3yr).median())
        },
        "same_urban_complete_cohort":{},
        "same_2015_grid_complete_cohort":{},
        "limitations":[
            "These are cross-sectional models relating wastewater infrastructure/production to GBIF-indexed plant candidate lists, not plant biomass, true alien introduction, butterfly resource consumption or colonization.",
            "Jones WWp, WWc, WWt, WWr are MODELLED wastewater volumes in 2015 5arc cells, not N concentration or soil N effects. Jones WWp-WWt is a modelled local not-treated volume proxy, not measured untreated water reaching terrestrial port plants.",
            "Satellite WorldCover 2021 built fraction and 2015 WW volume are retrospective proxies for pooled GBIF 2001-2024 host records; not time-aligned exposures.",
            "WWTP25 count and modelled urban wastewater volume share population drivers. Conditioning on them does NOT identify mediation or causal effect.",
            "Country FE, |latitude|, sin/cos longitude and all-indexed-host count proxy cannot remove all urbanization, plant-recording effort and port activity differences.",
            "Urban complete-case selection differs from full WWTP sample; use within-cohort estimates for interpretation.",
            "Treatment and production grid estimates may inherit common model assumptions; multi-collinearity means joint coefficients are descriptive.",
            "Only source-verified 2019 WPI ports with original complete 5km GBIF host query are included."
        ],
    }
    base=(("ln_wwtp",),
          ("ln_wwtp","ln_produced"),
          ("ln_wwtp","ln_produced","ln_treated"),
          ("ln_wwtp","ln_produced","ln_not_treated"),
          ("ln_produced",),
          ("ln_produced","ln_not_treated"),
          ("ln_wwtp","ln_produced","modeled_treated_share"))
    for name,urban in (("same_2015_grid_complete_cohort",False),("same_urban_complete_cohort",True)):
        output={}
        for pred in base:
            # city adjustment only on urban-complete; all models nested on same
            # cohort except treatment share drops modeled-production=0 cells.
            terms=list(pred)
            if urban:terms+=["built_fraction_2021","crop_fraction_2021"]
            label="+".join(terms)
            output[label]=model(z,outcome="ratio",predictors=terms,city_required=urban)
        for response in ("alien","native"):
            terms=["ln_wwtp","ln_produced"]
            if urban:terms+=["built_fraction_2021","crop_fraction_2021"]
            output[response+"_negative_control"]=model(
                z,outcome=response,predictors=terms,city_required=urban)
        summary[name]=output
        print("GLOBAL_SEWAGE_2015_URBAN_NESTED_MODEL",name,json.dumps(output,sort_keys=True),flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--jones",type=Path,required=True)
    p.add_argument("--worldcover",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    summary=run(a.jones,a.worldcover)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("WWTP_JONES2015_WORLDCOVER2021_GLOBAL_NESTED_SOURCE_PASS",json.dumps({
        "full_WPI":summary["full_WPI_ports"],
        "city_sample":summary["urban_raster_valid_ports"],
        "wwtp_wwp_rho":summary["WWTP_count25_vs_2015_production_Spearman_rho"],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
