#!/usr/bin/env python3
"""Source-verified 2022 CIRAD human N excretion vs cropland N supply grids
× 1596 fixed WPI ports with HydroWASTE/Jones wastewater and WorldCover.

These maps give MODELED human-excreted and agricultural nutrient-input
POTENTIAL, NOT measured untreated sewer discharge, plant soil N or butterflies.
Observational source independence is limited: population/land cover are shared
drivers of modeled excretion, 2015 wastewater and 2021 urban covariates.
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
import rasterio
from rasterio.windows import Window
from pyproj import Transformer
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

from exploratory_global_wwtp_jones_worldcover import read_source

FILES={
    "human":("Human_nutrient_excretion2022.gpkg",
             "global_grid_human_N_excretion_2022",
             "b8876ecda3456dbd7cfd5f007683771c"),
    "crop":("Nutrient_supply_to_cropland2022.gpkg",
            "global_grid_N_supply_to_cropland_2022",
            "281be6848535701b1acf65aa4bf87491")
}
RADIUS_MAIN=25000
RADIUS_SENS=10000
DOI="10.18167/DVN1/AWQNSR"


def md5_file(path:Path):
    digest=hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()


def open_source(folder:Path):
    datasets={}
    meta={}
    for key,(name,layer,md5) in FILES.items():
        p=folder/name
        if md5_file(p)!=md5:
            raise ValueError("Original CIRAD/Dataverse GPKG MD5 mismatch "+name)
        ds=rasterio.open(f"GPKG:{p}:{layer}")
        if ds.crs is None or ds.crs.to_string()!="ESRI:54035":
            ds.close()
            raise ValueError("Unexpected original 3km Equal Earth CRS")
        if (round(abs(ds.res[0])),round(abs(ds.res[1])))!=(3000,3000):
            ds.close()
            raise ValueError("Unexpected original 3km geospatial grid resolution")
        if ds.count!=1 or ds.nodata is None:
            ds.close()
            raise ValueError("No original raster mask/no-data gate")
        if ds.width!=11177 or ds.height!=4871:
            ds.close()
            raise ValueError("Unexpected 2022 global source grid dimensions")
        datasets[key]=ds
        meta[key]={"md5":md5,"source_layer":layer,
                   "source_file":name,
                   "crs":ds.crs.to_string(),"width":ds.width,"height":ds.height,
                   "nodata":ds.nodata,"scale":list(ds.scales),
                   "offset":list(ds.offsets),
                   "original_units_from_metadata":"10 kg N/ha/year (as provided); code uses native grid values"}
    if datasets["human"].transform!=datasets["crop"].transform:
        for ds in datasets.values():ds.close()
        raise ValueError("Human N and cropland N are not aligned grids")
    return datasets,meta


def circle_sample(ds,x:float,y:float,radius_m:int):
    tr=ds.transform
    xmin,xmax=x-radius_m,x+radius_m
    ymin,ymax=y-radius_m,y+radius_m
    # Clip window to original source footprint; do not extrapolate ocean or
    # generate made-up terrestrial nitrogen at masked cells.
    ix0=max(0,int(math.floor((xmin-tr.c)/tr.a)))
    ix1=min(ds.width,int(math.ceil((xmax-tr.c)/tr.a)))
    iy0=max(0,int(math.floor((tr.f-ymax)/abs(tr.e))))
    iy1=min(ds.height,int(math.ceil((tr.f-ymin)/abs(tr.e))))
    if ix0>=ix1 or iy0>=iy1:
        return {"valid_cells":0,"circle_cells":0,"mean":None,
                "zero_share_valid":None}
    window=Window(ix0,iy0,ix1-ix0,iy1-iy0)
    arr=ds.read(1,window=window,masked=True)
    rr,cc=np.ogrid[iy0:iy1,ix0:ix1]
    xx=tr.c+(cc+.5)*tr.a
    yy=tr.f+(rr+.5)*tr.e
    circle=((xx-x)**2+(yy-y)**2)<=radius_m**2
    vals=arr[circle]
    unmasked=np.asarray(np.ma.getmaskarray(vals))
    raw=np.asarray(vals.data,dtype=float)
    valid=(~unmasked)&np.isfinite(raw)&(raw>=0)
    keep=raw[valid]
    return {"valid_cells":int(valid.sum()),"circle_cells":int(circle.sum()),
            "mean":float(np.mean(keep)) if len(keep) else None,
            "zero_share_valid":float(np.mean(keep==0)) if len(keep) else None}


def raster_port_features(ports:pd.DataFrame,data_sources,radii=(25000,10000)):
    ds=data_sources["human"]
    trans=Transformer.from_crs("EPSG:4326",ds.crs,always_xy=True)
    output=[]
    for i,row in enumerate(ports.itertuples(index=False),1):
        x,y=trans.transform(float(row.longitude),float(row.latitude))
        d={"port_id":int(row.port_id),
           "projected_x":float(x),"projected_y":float(y)}
        for r in radii:
            for key,layer in data_sources.items():
                value=circle_sample(layer,x,y,r)
                pre=f"{key}_{r//1000}km_"
                for k,v in value.items():
                    d[pre+k]=v
        output.append(d)
        if i%200==0 or i==len(ports):
            print("CIRAD_2022_WPI_PORT_N_GRID_PROGRESS",i,len(ports),flush=True)
    z=pd.DataFrame(output)
    if len(z)!=1596 or z.port_id.nunique()!=1596:
        raise ValueError("CIRAD port/plant source WPI IDs misaligned")
    return z


def models_on_same_coverage(original,z,*,radius=25,urban=True):
    d=original.merge(z,on="port_id",validate="one_to_one")
    if len(d)!=1596:raise ValueError("N grid/GBIF port source join changed")
    d=d[d.in_within_country_comparison & ~d.incomplete_or_missing_grid].copy()
    if urban:d=d[d.quality=="COMPLETE"].copy()
    for prefix in ("human","crop"):
        key=f"{prefix}_{radius}km"
        d[key+"_valid_coverage"]=d[key+"_valid_cells"] / d[key+"_circle_cells"].replace(0,np.nan)
    need=[f"{x}_{radius}km_{y}" for x in ("human","crop") for y in ("mean","valid_coverage")]
    if urban:need+=["built_fraction_2021","crop_fraction_2021"]
    d=d.dropna(subset=need).copy()
    # Complete pixel source is insufficient near coasts; report and exclude
    # near-blank marine masks rather than replacing them with zero nitrogen.
    valid=(d[f"human_{radius}km_valid_coverage"]>=.10)&(d[f"crop_{radius}km_valid_coverage"]>=.10)
    valid &= (d[f"human_{radius}km_valid_cells"]>=8)&(d[f"crop_{radius}km_valid_cells"]>=8)
    d=d[valid].copy()
    d["log_wwtp"]=np.log1p(d.WWTP_count_25km)
    d["log_wwp"]=np.log1p(d.produced_25km_millionm3yr)
    d["log_human_N"]=np.log1p(d[f"human_{radius}km_mean"])
    d["log_crop_N"]=np.log1p(d[f"crop_{radius}km_mean"])
    d["abs_lat"]=d.latitude.abs()
    d["sinlon"]=np.sin(np.deg2rad(d.longitude))
    d["coslon"]=np.cos(np.deg2rad(d.longitude))
    d["log_indexed"]=np.log1p(d.all_taxonomically_indexed)
    d["ratio"]=np.log((d.alien_candidates+.5)/(d.native_candidates+.5))
    d["alien"]=np.log1p(d.alien_candidates)
    d["native"]=np.log1p(d.native_candidates)
    counts=d.groupby("country").agg(n=("port_id","size"),wwtp_variation=("log_wwtp","nunique"),
                                    human_N_variation=("log_human_N","nunique"))
    eligible=counts.index[(counts.n>=5)&(counts.wwtp_variation>=2)&(counts.human_N_variation>=2)]
    d=d[d.country.isin(eligible)].copy()
    if len(d)<400 or d.country.nunique()<15:
        raise ValueError(f"Not enough 3km land-source quality matched ports {len(d)} countries {d.country.nunique()}")
    control=["log_wwp","log_wwtp","abs_lat","sinlon","coslon","log_indexed"]
    if urban:control+=["built_fraction_2021","crop_fraction_2021"]
    common={"n_ports":int(len(d)),"n_countries":int(d.country.nunique()),
            "2022_grid_radius_km":radius,
            "2021_urban_controls":urban,
            "median_human_N_source_grid_raw":float(d[f"human_{radius}km_mean"].median()),
            "median_crop_N_source_grid_raw":float(d[f"crop_{radius}km_mean"].median()),
            "human_N_rank_rho_vs_WWp2015":float(spearmanr(d.log_human_N,d.log_wwp).statistic),
            "human_N_rank_rho_vs_WWTP25km":float(spearmanr(d.log_human_N,d.log_wwtp).statistic),
            "min_human_grid_coverage":float(d[f"human_{radius}km_valid_coverage"].min()),
            "min_crop_grid_coverage":float(d[f"crop_{radius}km_valid_coverage"].min())}
    if urban:
        common["human_N_rank_rho_vs_built2021"]=float(spearmanr(d.log_human_N,d.built_fraction_2021).statistic)
    ret={"support":common,"models":{}}
    cases={
        "baseline_local_WW_and_urban":[],
        "human_excreted_N_addition":["log_human_N"],
        "crop_N_addition":["log_crop_N"],
        "both_human_N_and_crop_N":["log_human_N","log_crop_N"]
    }
    def fit(response,extra,label):
        controls=control if response=="ratio" else [k for k in control if k!="log_indexed"]
        f=smf.ols(response+" ~ "+" + ".join(controls+extra+["C(country)"]),data=d).fit(
            cov_type="cluster",cov_kwds={"groups":d.country})
        out={"status":"STATIC_GRID_ASSOCIATION_NOT_NITROGEN_CAUSAL_IMPACT",
             "source_port_sample_n":len(d),"countries":int(d.country.nunique()),
             "formula":f.model.formula,"response":response,"R2":float(f.rsquared)}
        for term in ["log_human_N","log_crop_N","log_wwp","log_wwtp"]:
            if term in f.params.index:
                out[term]={"beta":float(f.params[term]),
                           "CI95_country_cluster":list(map(float,f.conf_int().loc[term])),
                           "nominal_p":float(f.pvalues[term])}
        return out
    for label,extra in cases.items():
        ret["models"][label]=fit("ratio",extra,label)
    for outcome in ["alien","native"]:
        ret["models"][outcome+"_two_N_sources"]=fit(outcome,["log_human_N","log_crop_N"],outcome)
    return ret,d


def run(jones_path,worldcover_path,source_folder,*,port_table_path):
    source=read_source(jones_path,worldcover_path)
    ds,meta=open_source(source_folder)
    try:
        ports=raster_port_features(source,ds)
    finally:
        for v in ds.values():v.close()
    models={}
    for radius in (25,10):
        for urban in (True,False):
            key=f"{radius}km_"+("urban_2021" if urban else "no_urban")
            result,_=models_on_same_coverage(source,ports,radius=radius,urban=urban)
            models[key]=result
            print("CIRAD_2022_HUMAN_VS_CROP_N_GLOBAL_FE",key,
                  json.dumps({"support":result["support"],
                              "jointN":result["models"]["both_human_N_and_crop_N"]},sort_keys=True),flush=True)
    port_table_path.parent.mkdir(parents=True,exist_ok=True)
    ports.to_csv(port_table_path,index=False)
    return {"schema":"chocho_CIRAD2022_Human_vs_CroplandN_3km_WPI1596_JonesWorldCover_observational_v1",
            "status":"SOURCE_PINNED_LOCAL_N_EXCRETION_POTENTIAL_NOT_WASTEWATER_N_TO_SOIL",
            "citation":"Wassenaar 2026, DOI:10.18167/DVN1/AWQNSR",
            "source_2015_wwp_run":41563940689,"source_worldcover_2021_run":41562997326,
            "raster_meta":meta,"WPI_ports":len(source),
            "N_grid_radii_km":[25,10],"models":models,
            "not_in_repository_original_cirad_data":True,
            "limits":[
                "Human excreted N 2022 from national food balance multiplied by gridded 2022 population; not water-treatment effluent, sewage in soil, or plant-available local nitrogen.",
                "Cropland N from gridded farmland national fertilizer/crop input; not a field soil N exposure and not total regional atmospheric N deposition.",
                "Both modelled N sources share urban/population/landuse drivers with port wastewater production and satellite built cover; full causal separation impossible.",
                "Original N grids sampled in 25km and 10km circles, while GBIF plants sampled 5km WPI; scales and years differ.",
                "No data pixels excluded, never set to zero; valid zeros retained where source provides them.",
                "Marine masks and national projection source errors can bias coastal buffers. Coverage rates retained and minimal support gated.",
                "WPI2019 ports, Jones 2015 produced wastewater, ESA2021 static land cover and pooled 2001-2024 GBIF host lists are not temporally matched.",
                "N densities use as-delivered source numeric scale; no redistribution to port, effluent N flux, butterfly abundance, or larval survival inferred.",
                "OLS nominal p after exploratory multiple model permutations is not confirmatory and is country-clustered."
            ]}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jones",type=Path,required=True)
    ap.add_argument("--worldcover",type=Path,required=True)
    ap.add_argument("--source-dir",type=Path,required=True)
    ap.add_argument("--port-table",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    result=run(args.jones,args.worldcover,args.source_dir,port_table_path=args.port_table)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("CIRAD_3KM_2022_HUMAN_CROP_N_WPI1596_DONE",json.dumps({
       "WPI":result["WPI_ports"],
       "25km_urban_support":result["models"]["25km_urban_2021"]["support"]},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
