#!/usr/bin/env python3
"""Independent national sewage-N inventory as country-level moderator of
port-scale butterfly larval-host GBIF recording: an EXPLORATORY test.

NIES GNINPO-DSM 2025: 2005–2022 COUNTRY × YEAR kgN/year, *modelled*.
National unmanaged pathways (UNM1..UNM5) are not measured site N deposition.
Within-country WPI-port 2015 WW production can be compared against varying
country national N shares via interaction; country fixed effects absorb the
N-share main effect. No invented local kgN/port redistribution.
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
import pycountry
import statsmodels.formula.api as smf

from exploratory_global_wwtp_jones_worldcover import read_source

SOURCE_MD5 = "67db41b6d892c5f1055b9403e194d046"
NIES_DOI = "10.17595/20260227.001"
N_CODES = ["TNE"]+[f"NIS{i}" for i in range(1,6)]+[f"PMN{i}" for i in range(1,4)]+[f"UNM{i}" for i in range(1,6)]


def read_nitrogen(path: Path):
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    if md5 != SOURCE_MD5:
        raise ValueError("NIES official MD5 verification failed")
    d = pd.read_csv(path, dtype={"ISO3":str})
    if len(d) != 3110 or len(d.columns) != 16:
        raise ValueError("Official NIES 3110-country-year source changed")
    if set(d.columns) != {"ISO3","Year",*N_CODES}:
        raise ValueError("Unexpected official NIES nitrogen schema")
    d["ISO3"]=d["ISO3"].str.strip().str.upper()
    if (d.ISO3.str.len()!=3).any():
        raise ValueError("Malformed ISO3 country code")
    for col in ["Year",*N_CODES]:
        d[col]=pd.to_numeric(d[col],errors="coerce")
    if d[["Year",*N_CODES]].isna().any().any() or not np.isfinite(d[N_CODES].to_numpy()).all():
        raise ValueError("Missing or non-finite official NIES N flows")
    if not (d.Year.astype(int)==d.Year).all():
        raise ValueError("Non-integer annual dates")
    d["Year"]=d["Year"].astype(int)
    if (d.Year.min(),d.Year.max())!=(2005,2022):
        raise ValueError("NIES official temporal coverage changed")
    if d.duplicated(["ISO3","Year"]).any():
        raise ValueError("Duplicate ISO3/year N measurements")
    if (d[N_CODES]<0).any().any():
        raise ValueError("Physically invalid negative N flow")
    d["managed_N_kg"]=d[["PMN1","PMN2","PMN3"]].sum(axis=1)
    d["unmanaged_N_kg"]=d[[f"UNM{i}" for i in range(1,6)]].sum(axis=1)
    d["input_NIS_kg"]=d[[f"NIS{i}" for i in range(1,6)]].sum(axis=1)
    # These sums should have shared units and roughly satisfy the source N
    # bookkeeping, but do not force equality or reinterpret "properly managed"
    # nitrogen as removal, or unmanaged as measured river effluent.
    d["input_mass_balance_rel"]=abs(d["input_NIS_kg"]-d.TNE)/d.TNE.replace(0,np.nan)
    d["pathway_mass_balance_rel"]=abs(d["managed_N_kg"]+d["unmanaged_N_kg"]-d.TNE)/d.TNE.replace(0,np.nan)
    d["fraction_unmanaged"]=d["unmanaged_N_kg"]/d.TNE.replace(0,np.nan)
    d["fraction_untreated_sewer"]=d["UNM1"]/d.TNE.replace(0,np.nan)
    d["fraction_plant_managed"]=d["PMN1"]/d.TNE.replace(0,np.nan)
    for col in ["fraction_unmanaged","fraction_untreated_sewer","fraction_plant_managed"]:
        if (d[col].dropna()>1.005).any() or (d[col].dropna()<0).any():
            raise ValueError("Out-of-range national pathway share "+col)
    diag={
        "source_year_count":d.Year.nunique(),
        "source_year_minmax":[int(d.Year.min()),int(d.Year.max())],
        "source_rows":len(d),
        "source_country_count":int(d.ISO3.nunique()),
        "countries_in_2015":int((d.Year==2015).sum()),
        "source_MD5":md5,
        "all_flow_units":"kg N / year (NIES official DataDictionary)",
        "NIS1_to_5_balance_max_fraction":float(d.input_mass_balance_rel.max()),
        "managed_plus_unmanaged_balance_max_fraction":float(d.pathway_mass_balance_rel.max()),
        "unmanaged_fraction_2015_median_countries":float(d.loc[d.Year==2015,"fraction_unmanaged"].median())
    }
    return d,diag


def country_exposure(d,year_window):
    start,stop=year_window
    q=d[(d.Year>=start)&(d.Year<=stop)].copy()
    g=q.groupby("ISO3",as_index=False).agg(
        n_years=("Year","nunique"),kgN_total=("TNE","sum"),
        kgN_unmanaged=("unmanaged_N_kg","sum"),
        kgN_untreated_sewer=("UNM1","sum"),
        kgN_managed_plant=("PMN1","sum"))
    expected=stop-start+1
    g=g[g.n_years==expected].copy()
    if g.duplicated("ISO3").any():raise ValueError("Multiple country values after national aggregation")
    g["unmanaged_fraction"]=g.kgN_unmanaged/g.kgN_total.replace(0,np.nan)
    g["untreated_sewer_fraction"]=g.kgN_untreated_sewer/g.kgN_total.replace(0,np.nan)
    g["plant_managed_fraction"]=g.kgN_managed_plant/g.kgN_total.replace(0,np.nan)
    g["window"]=f"{start}-{stop}"
    return g


def iso3(iso2):
    x=str(iso2).strip().upper()
    if x=="XK":return "XKX"
    obj=pycountry.countries.get(alpha_2=x)
    return obj.alpha_3 if obj else None


def build_port_country_source(jones:Path,worldcover:Path,nitrogen:pd.DataFrame,*,window=(2015,2015),urban=True):
    ports=read_source(jones,worldcover)
    ports["ISO3"]=ports.country.map(iso3)
    exposures=country_exposure(nitrogen,window)
    p=ports.merge(exposures,on="ISO3",how="left",validate="many_to_one",indicator=True)
    count_missing=int((p["_merge"]=="left_only").sum())
    p=p[p.in_within_country_comparison].copy()
    p=p[~p.incomplete_or_missing_grid].copy()
    if urban:p=p[p.quality=="COMPLETE"].copy()
    p=p.dropna(subset=["unmanaged_fraction","untreated_sewer_fraction","plant_managed_fraction",
                      "produced_25km_millionm3yr","built_fraction_2021","crop_fraction_2021"] if urban
               else ["unmanaged_fraction","untreated_sewer_fraction","plant_managed_fraction",
                     "produced_25km_millionm3yr"])
    p["log_wwp"]=np.log1p(p.produced_25km_millionm3yr)
    p["log_wwtp"]=np.log1p(p.WWTP_count_25km)
    p["abs_lat"]=p.latitude.abs()
    p["sin_lon"]=np.sin(np.deg2rad(p.longitude))
    p["cos_lon"]=np.cos(np.deg2rad(p.longitude))
    p["log_indexed"]=np.log1p(p.all_taxonomically_indexed)
    p["y_ratio"]=np.log((p.alien_candidates+.5)/(p.native_candidates+.5))
    p["y_alien"]=np.log1p(p.alien_candidates)
    p["y_native"]=np.log1p(p.native_candidates)
    g=p.groupby("country").agg(n=("port_id","size"),variation=("log_wwp","nunique"),
                                wwtp_variation=("log_wwtp","nunique"))
    ok=g.index[(g.n>=5)&(g.variation>=2)&(g.wwtp_variation>=2)]
    p=p[p.country.isin(ok)].copy()
    if len(p)<450 or p.country.nunique()<20:
        raise ValueError(f"Not enough original ports/countries with NIES coverage: {len(p)}, {p.country.nunique()}")
    p["wwp_country_center"]=p.log_wwp-p.groupby("country").log_wwp.transform("mean")
    # Equal-country, rather than port-number-weighted, standardization of
    # national predictors. Each country has exactly one share per window.
    for raw,label in (("unmanaged_fraction","unmanaged"),
                      ("untreated_sewer_fraction","untreated_sewer"),
                      ("plant_managed_fraction","plant_managed")):
        means=p.groupby("country")[raw].first()
        mu=float(means.mean());sd=float(means.std(ddof=1))
        if sd<.025:raise ValueError("No meaningful cross-country national N-share variation")
        p[label+"_country_z"]=(p[raw]-mu)/sd
        p["interaction_"+label]=p.wwp_country_center*p[label+"_country_z"]
    audit={"n_ports_after_joins":len(p),"countries":int(p.country.nunique()),
           "unmapped_original_ports":count_missing,
           "source_window":list(window),
           "urban_controls":urban,
           "country_N_unmanaged_fraction_range":[float(p.unmanaged_fraction.min()),float(p.unmanaged_fraction.max())],
           "national_fraction_not_used_as_direct_port_exposure":True,
           "country_main_effect_absorbed_by_country_FE":True}
    return p,audit


def clustered_model(p,*,period,urban=True,case="unmanaged",outcome="y_ratio"):
    var="interaction_"+case if case else None
    base=["log_wwp","log_wwtp","abs_lat","sin_lon","cos_lon","C(country)"]
    if outcome=="y_ratio":base.append("log_indexed")
    if urban:base+=["built_fraction_2021","crop_fraction_2021"]
    expression=outcome+" ~ "+" + ".join([*base,var] if case else base)
    fit=smf.ols(expression,data=p).fit(cov_type="cluster",cov_kwds={"groups":p.country})
    res={"status":"ESTIMATED_EXPLORATORY_CROSS_LEVEL_MODERATOR_NOT_SITE_N",
         "national_N_year_window":period,"urban":urban,
         "n":len(p),"countries":int(p.country.nunique()),
         "outcome":outcome,
         "country_cluster_robust":True,"formula":expression,
         "R2":float(fit.rsquared)}
    for col in ["log_wwp","log_wwtp",var]:
        if col in fit.params.index:
            res[col]={"beta":float(fit.params[col]),
                      "95_country_cluster":list(map(float,fit.conf_int().loc[col])),
                      "nominal_p":float(fit.pvalues[col])}
    return res


def run(nies_path:Path,jones_path:Path,urban_path:Path):
    nies,quality=read_nitrogen(nies_path)
    outputs={"schema":"chocho_NIES2025_N_kg_country_year_2015_HydroWASTE_Jones2015_WorldCover2021_interaction_v1",
            "status":"SOURCE_VERIFIED_MODELED_NITROGEN_CROSS_LEVEL_EXPLORATION_NOT_CAUSAL_PORT_N_LOAD",
            "NIES_DOI":NIES_DOI,
            "source_runs":{"nies":41921696817,"jones_port_2015":41563940689,
                           "worldcover_port_2021":41562997326},
            "source_details":quality,"frozen_original_WPI_complete_ports":1596,
            "result_nature":"national N pathway fractions moderate within-country local wastewater-production correlation; no national main effect with country fixed effects",
            "models":{},"population":{},
            "limits":[
                "NIES modeled kgN/year at COUNTRY x YEAR not per-port N concentration, soil N, nitrogen deposition or effluent load to terrestrial vegetation.",
                "PMN1 is nitrogen handled at plants, not nitrogen removed by treatment; UNM group is unmanaged pathways not directly measured aquatic/soil discharge.",
                "Country fixed effects absorb national N fraction main effects; interaction reflects a between-country moderation of within-country port wastewater proxy.",
                "Within-country ports share same national N exposure; cluster errors at country; no pseudoreplicated treatment replicates.",
                "WPI port year 2019, Jones model year 2015, satellite land cover 2021 and butterfly host GBIF pooled 2001–2024, not synchronized.",
                "Global butterfly larval host categories are candidate plants indexed near ports, not locally verified foodplant availability, oviposition or butterfly colonization.",
                "Country sanitation services and wealth co-vary strongly; interactions are not N exposure treatment effects.",
                "National N emissions derived from population and sanitation assumptions, not independent measurements of N delivered to port-adjacent plants.",
                "Country-level N exposure cannot be identified independently with country fixed-effect intercept; no such standalone country N coefficient is reported."
            ]}
    for label,years in (("primary_2015",(2015,2015)),
                        ("prior_2005_2010",(2005,2010)),
                        ("later_2011_2022",(2011,2022))):
        for urban in (True,False):
            p,audit=build_port_country_source(jones_path,urban_path,nies,window=years,urban=urban)
            key=label+("_urban_controlled" if urban else "_full_grid_sample")
            outputs["population"][key]=audit
            models={}
            # Same cohort for each national N moderator, otherwise coefficient
            # changes can be driven by different country sampling.
            models["baseline_no_country_N_interaction"]=clustered_model(p,period=label,urban=urban,case=None)
            for v in ("unmanaged","untreated_sewer","plant_managed"):
                models[v]=clustered_model(p,period=label,urban=urban,case=v)
            if label=="primary_2015" and urban:
                models["native_record_negative_control_unmanaged"]=clustered_model(
                    p,period=label,urban=True,case="unmanaged",outcome="y_native")
                models["alien_record_unmanaged"]=clustered_model(
                    p,period=label,urban=True,case="unmanaged",outcome="y_alien")
            outputs["models"][key]=models
            val=models["unmanaged"]["interaction_unmanaged"]
            print("NIES_N_CROSSLEVEL_JOINT_PORT_HOST",key,
                  json.dumps({"sample":audit,"unmanaged_interaction":val},sort_keys=True),flush=True)
    return outputs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nies",type=Path,required=True)
    p.add_argument("--jones",type=Path,required=True)
    p.add_argument("--worldcover",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.nies,a.jones,a.worldcover)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("NIES_MODEL_COMPLETE",out["population"]["primary_2015_urban_controlled"],
          out["models"]["primary_2015_urban_controlled"]["unmanaged"]["interaction_unmanaged"],flush=True)


if __name__=="__main__":
    main()
