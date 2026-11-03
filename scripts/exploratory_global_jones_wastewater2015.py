#!/usr/bin/env python3
"""Independent PANGAEA 918731 2015 wastewater flow grid × WPI1596 butterfly hosts.

This is observational 2015 WW production/collection/treatment/reuse (million
m3/yr) near WPI port points. NEVER treat any volume as soil nitrogen load.
Country FE and plant-indexing controls are descriptive, not causal.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import statistics
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

SOURCE="https://store.pangaea.de/Publications/JonesE-etal_2020/Global_wastewater.zip"
SHA256="c413f7493fb5c2b5a07d3d974abf3b9af25a32fcadd8eb839cc3fe302838e5eb"
PANGAEA_DOI="10.1594/PANGAEA.918731"
VARIABLES={"WWp":"produced","WWc":"collected","WWt":"treated","WWr":"reused"}
R=6371.0088


def haversine_distance(lat,lon,lats,lons):
    lat0=np.radians(float(lat))
    lat1=np.radians(np.asarray(lats,dtype=float))
    dlat=lat1-lat0
    dlon=np.radians((np.asarray(lons,dtype=float)-float(lon)+180)%360-180)
    x=np.sin(dlat/2)**2+np.cos(lat0)*np.cos(lat1)*np.sin(dlon/2)**2
    return 2*R*np.arcsin(np.sqrt(np.maximum(0,np.minimum(1,x))))


def extract_data(path:Path):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SHA256:
        raise ValueError("PANGAEA 2015 original ZIP SHA256 changed")
    with zipfile.ZipFile(path) as z:
        f=[x for x in z.namelist() if x.endswith(".nc")]
        if len(f)!=1:raise ValueError("Unexpected original Jones et al netCDF archive entries")
        extraction_path=path.with_suffix(".nc")
        extraction_path.write_bytes(z.read(f[0]))
    with Dataset(extraction_path) as ds:
        if (len(ds.dimensions["lat"]),len(ds.dimensions["lon"]))!=(2160,4320):
            raise ValueError("2015 global grid should be 5 arc-min, 2160x4320")
        lats=np.asarray(ds.variables["lat"][:],dtype=float)
        lons=np.asarray(ds.variables["lon"][:],dtype=float)
        arrays={}
        for v in VARIABLES:
            src=ds.variables[v]
            if tuple(src.dimensions)!=("month","lat","lon") or src.shape!=(1,2160,4320):
                raise ValueError("Unexpected Jones global wastewater variable dimensions "+v)
            a=src[0,:,:]
            if np.ma.isMaskedArray(a):
                a=np.ma.filled(a,np.nan)
            a=np.asarray(a,dtype=float)
            a[~np.isfinite(a)]=np.nan
            a[a<-1e10]=np.nan
            arrays[v]=a
        meta={
            "grid_dims":{"lat":len(lats),"lon":len(lons)},
            "latitude_minmax":[float(lats.min()),float(lats.max())],
            "longitude_minmax":[float(lons.min()),float(lons.max())],
            "variables":VARIABLES,
            "unit":"million m3/yr per 5-arcmin cell",
            "reference_year":2015
        }
    return lats,lons,arrays,meta


def source_port_features(path:Path):
    df=pd.read_csv(path,keep_default_na=False)
    if len(df)!=1596 or df.port_id.nunique()!=1596 or df.country.nunique()!=142:
        raise ValueError("WPI HydroWASTE original source population changed")
    if (int(df.in_within_country_comparison.sum()),int(df.loc[df.in_within_country_comparison,"country"].nunique()))!=(1368,52):
        raise ValueError("Source-independent 1368 ports 52 countries FE cohort changed")
    for k in ("latitude","longitude","alien_candidates","native_candidates","all_taxonomically_indexed","WWTP_count_25km"):
        if k not in df.columns:raise ValueError("Original port feature missing "+k)
    return df


def grid_one(lat,lon,lats,lons,arrays,radius_km=25):
    dlat=radius_km/110.5
    dlon=radius_km/(111.15*max(.01,math.cos(math.radians(lat))))
    im=np.flatnonzero(abs(lats-float(lat))<=dlat+.09)
    dx=abs((lons-float(lon)+180)%360-180)
    jm=np.flatnonzero(dx<=dlon+.09)
    if not len(im) or not len(jm):
        return {"cell_centers_in_25km":0,"valid_production_cell_centers":0,
                "incomplete_or_missing_grid":True,
                **{f"{name}_25km_millionm3yr":None for name in VARIABLES.values()}}
    ii,jj=np.meshgrid(im,jm,indexing="ij")
    latgrid=lats[ii]
    longrid=lons[jj]
    distance=haversine_distance(float(lat),float(lon),latgrid,longrid)
    choice=distance<=radius_km
    idxi,idxj=ii[choice],jj[choice]
    production=arrays["WWp"][idxi,idxj]
    valid=np.isfinite(production) & (production>=0)
    output={"cell_centers_in_25km":int(choice.sum()),
            "valid_production_cell_centers":int(valid.sum()),
            "incomplete_or_missing_grid":not bool(valid.any())}
    for v,alias in VARIABLES.items():
        values=arrays[v][idxi,idxj]
        allowed=np.isfinite(values)&(values>=0)&valid
        output[f"{alias}_25km_millionm3yr"]=float(np.sum(values[allowed])) if allowed.any() else None
        output[f"{alias}_valid_grid_cells"]=int(allowed.sum())
    # Do not subtract treatment from production where source describes
    # possible spatial routing and modeled assignment to WWTP locations.
    return output


def raster_exposure(port,lats,lons,arrays):
    data=[]
    for i,row in enumerate(port.itertuples(),1):
        d=grid_one(float(row.latitude),float(row.longitude),lats,lons,arrays)
        data.append({"port_id":int(row.port_id),**d})
        if i%200==0 or i==len(port):
            print("JONES_5ARC_GLOBAL_PORT_SAMPLED",i,len(port),
                "valid",sum(not x["incomplete_or_missing_grid"] for x in data),flush=True)
    x=pd.DataFrame(data)
    assert len(x)==1596 and x.port_id.nunique()==1596
    return port.merge(x,on="port_id",validate="one_to_one")


def clustered(df,terms,name):
    def estimate(controls, label):
        d=df.loc[df.in_within_country_comparison].copy()
        d=d.loc[~d["incomplete_or_missing_grid"]].copy()
        # All WW volumes can be 0; zero denotes modelled grid zero if source valid.
        for v in ("produced","collected","treated","reused"):
            d[f"log_{v}"]=np.log1p(d[f"{v}_25km_millionm3yr"].fillna(0))
        d["log_wwtp"]=np.log1p(d.WWTP_count_25km.astype(float))
        d["abs_lat"]=d.latitude.abs()
        d["sinlon"]=np.sin(np.deg2rad(d.longitude))
        d["coslon"]=np.cos(np.deg2rad(d.longitude))
        d["ln_total"]=np.log1p(d.all_taxonomically_indexed)
        d["ratio"]=np.log((d.alien_candidates+.5)/(d.native_candidates+.5))
        variation=d.groupby("country").log_wwtp.transform("nunique")
        n_country=d.groupby("country").country.transform("size")
        d=d.loc[(variation>=2)&(n_country>=5)].copy()
        if len(d)<500 or d.country.nunique()<18:
            return {"status":"INSUFFICIENT_ELIGIBLE_COUNTRIES","n":len(d),"countries":d.country.nunique()}
        eq="ratio ~ log_wwtp + "+ " + ".join(controls)
        model=smf.ols(eq,data=d).fit(cov_type="cluster",cov_kwds={"groups":d.country})
        b=float(model.params["log_wwtp"])
        lower,upper=[float(x) for x in model.conf_int().loc["log_wwtp"]]
        return {"status":"WWTP_COEFFICIENT_AFTER_MODELED_2015_FLOW_CONTROLS",
                "model_label":label,"n_ports":len(d),"n_countries":d.country.nunique(),
                "WWTP25_log1p_beta":b,"WWTP25_country_cluster_95":[lower,upper],
                "WWTP25_nominal_p":float(model.pvalues["log_wwtp"]),
                "estimated_parameters":len(model.params),
                "model_controls":controls,
                "R2":float(model.rsquared),
                "source_WW_generated_not_N_emissions":True}
    return estimate


def run(port_path,zip_path):
    df=source_port_features(port_path)
    lats,lons,arr,meta=extract_data(zip_path)
    totals={}
    good=np.isfinite(arr["WWp"])
    for name,a in arr.items():
        g=good&np.isfinite(a)&(a>=0)
        totals[name]={"global_valid_cell_count":int(g.sum()),
                     "global_modeled_volume_sum_millionm3yr":float(np.sum(a[g])),
                     "positive_cells":int(np.count_nonzero(a[g]>0))}
    checks={}
    for lhs,rhs in (("WWt","WWp"),("WWc","WWp"),("WWr","WWt")):
        s=good&np.isfinite(arr[lhs])&np.isfinite(arr[rhs])&(arr[lhs]>=0)&(arr[rhs]>=0)
        checks[f"{lhs}_greater_than_{rhs}_grid_cells"]=int(np.count_nonzero(arr[lhs][s]>arr[rhs][s]+.001))
    port=raster_exposure(df,lats,lons,arr)
    goodport=port.loc[~port.incomplete_or_missing_grid]
    x=clustered(port,(),"") if False else None
    # Same cohort for every version; do not compare differing source support.
    estimate=clustered(port,None,"test")
    base=["abs_lat","sinlon","coslon","ln_total","C(country)"]
    models={
        "baseline_wwtp25_country_geo_observation":estimate(base,"baseline"),
        "plus_2015_modeled_wastewater_production25":estimate(base+["log_produced"],"plus_WWp"),
        "plus_2015_production_and_treatment25":estimate(base+["log_produced","log_treated"],"plus_WWp_WWt"),
        "plus_2015_production_collection_and_treatment25":estimate(
            base+["log_produced","log_collected","log_treated"],"plus_WWp_WWc_WWt")
    }
    result={
        "schema":"chocho_PANGAEA_2015_global_wastewater_5arc_HydroWASTE_WPI1596_v1",
        "status":"MODELLED_2015_WASTEWATER_PRODUCTION_TREATMENT_OBSERVATIONAL_NOT_NUTRIENT_CAUSATION",
        "dataset_DOI":PANGAEA_DOI,
        "dataset_sha256_original_zip":SHA256,
        "source_hydrowaste_global_run":41414513676,
        "original_complete_ports":1596,
        "original_country_fixed_effect_ports":1368,
        "grid_metadata":meta,
        "WW_grid_worldwide_totals":totals,
        "grid_same_cell_volume_inequality_diagnostics":checks,
        "ports_with_any_valid_modelled_WW_production_grid_cells":len(goodport),
        "ports_with_positive_modelled_WW_production_25km":int((goodport.produced_25km_millionm3yr>0).sum()),
        "ports_with_positive_modelled_WW_treatment_25km":int((goodport.treated_25km_millionm3yr>0).sum()),
        "median_modelled_production_millionm3yr_near_ports":float(goodport.produced_25km_millionm3yr.median()),
        "median_modelled_treated_millionm3yr_near_ports":float(goodport.treated_25km_millionm3yr.median()),
        "spearman_wwtp25_to_2015_wastewater_production":float(spearmanr(
            goodport.WWTP_count_25km,goodport.produced_25km_millionm3yr).statistic),
        "country_FE_same_source_comparisons":models,
        "limitations":[
          "Jones et al. 2020 PANGAEA data are modeled circa-2015 wastewater flows, not measured sewage composition or nitrogen load.",
          "WPI port positions are 2019, HydroWASTE facility locations circa-2021, GBIF candidate plant record intervals 2001-2024.",
          "5 arcminute approximately 10km global grid with coarse land/ocean coastal uncertainty, 25km radius modeled cell-centre sums.",
          "Treatment, reuse and production are gridded model variables; their same-cell difference cannot be assumed untreated discharge without river routing.",
          "WWTP facility density and wastewater production both covary with urban population/infrastructure; this is a confounding check not causal attribution.",
          "Country FE and plant indexing total do not guarantee correction of port-level urbanization, collection protocols, or true plant growth.",
          "Zero modeled wastewater volume is a model output, not proof of no wastewater; missing cells excluded.",
          "Multiple correlated wastewater predictors in a single model can be unstable; sensitivity only.",
          "No butterfly breeding, larval use, or island stepping stone interaction was measured."
        ]
    }
    return result,port


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ports",type=Path,required=True)
    p.add_argument("--zip",type=Path,required=True)
    p.add_argument("--json",type=Path,required=True)
    p.add_argument("--table",type=Path,required=True)
    a=p.parse_args()
    result,table=run(a.ports,a.zip)
    a.json.parent.mkdir(parents=True,exist_ok=True)
    a.json.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    a.table.parent.mkdir(parents=True,exist_ok=True)
    table.to_csv(a.table,index=False)
    print("JONES2015_GLOBAL_WASTEWATER_WPI_PORTS_MODEL",json.dumps({
        "valid":result["ports_with_any_valid_modelled_WW_production_grid_cells"],
        "WWP_positive":result["ports_with_positive_modelled_WW_production_25km"],
        "global_grid_checks":result["grid_same_cell_volume_inequality_diagnostics"],
        "model":result["country_FE_same_source_comparisons"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
