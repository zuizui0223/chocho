#!/usr/bin/env python3
"""Join HydroWASTE v1 (58,502 WWTPs) with fixed chocho WPI host-inventory ports.

Global observational co-geography only. HydroWASTE reports treated-water
discharge (m3/day), not N loading; estimated outfalls can be kilometres from
plants. Port/butterfly hosts are GBIF-indexed lists, not true biological absence
or demonstrated ship introductions. No effects of sewage on butterfly fitness.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import time
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

RADIUS_KM=[10,25,50]
EARTH=6371.0088
WWTP_SOURCE="https://raw.githubusercontent.com/rfordatascience/tidytuesday/774fb8cc81e563d6d30efdcb7b48a38565e923cc/data/2022/2022-09-20/HydroWASTE_v10.csv"
WWTP_BLOB="0506f873f6541c5f324478894a9731290d722bce"
WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
HOST_BLOB="fdf73b11ddcb1644f46730123350ae52cc61cae6"
HEADER={"User-Agent":"chocho-global-wwtp-port-host-20261103/1.0"}
KNOWN_EXPOSURES=("wwtp_25km","outfall_25km","wwtp_flow_25km_m3_day")


def git_blob_sha(b:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()


def retry_session_get(url,params=None,timeout=60,attempts=5):
    last=None
    for i in range(attempts):
        try:
            r=requests.get(url,params=params,headers=HEADER,timeout=timeout)
            if r.status_code in (408,425,429,500,502,503,504):
                last=RuntimeError(f"upstream HTTP {r.status_code}")
                time.sleep(min(24,2**i*2))
                continue
            r.raise_for_status()
            return r
        except (requests.Timeout,requests.ConnectionError) as ex:
            last=ex
            time.sleep(min(24,2**i*2))
    raise RuntimeError(f"upstream did not reply successfully: {url[:180]}: {last}")


def download_hydrowaste(cache:Path)->pd.DataFrame:
    cache.parent.mkdir(parents=True,exist_ok=True)
    if not cache.exists():
        r=retry_session_get(WWTP_SOURCE,timeout=100)
        cache.write_bytes(r.content)
    b=cache.read_bytes()
    if git_blob_sha(b)!=WWTP_BLOB:
        raise ValueError(f"HydroWASTE CSV source changed: {git_blob_sha(b)}")
    df=pd.read_csv(io.BytesIO(b),low_memory=False)
    required={"WASTE_ID","LAT_WWTP","LON_WWTP","LAT_OUT","LON_OUT",
              "CNTRY_ISO","STATUS","WASTE_DIS","QUAL_WASTE","QUAL_LOC",
              "COAST_10KM","COAST_50KM","LEVEL"}
    if not required.issubset(df.columns):
        raise ValueError("HydroWASTE required columns unavailable: "+str(required-set(df.columns)))
    if len(df)!=58502 or df.WASTE_ID.nunique()!=58502:
        raise ValueError("Pinned HydroWASTE global source population not 58,502 WWTPs")
    return df


def read_host_inventory(path:Path)->pd.DataFrame:
    b=path.read_bytes()
    if git_blob_sha(b)!=HOST_BLOB:
        raise ValueError("Original complete-port botanical snapshot was modified")
    df=pd.read_csv(io.BytesIO(b),dtype={"country":str})
    if len(df)!=3630 or df.port_id.nunique()!=3630:
        raise ValueError("Original registered WPI 3,630-port population changed")
    if int(df.query_complete.sum())!=1596:
        raise ValueError("Original 1,596 complete botanical catalogue sites changed")
    incomplete=df.loc[~df.query_complete,
        ["contemporary_not_native_in_catalogue","native_in_catalogue"]]
    if incomplete.notna().any().any():
        raise ValueError("Incomplete botanical ports must retain unknown host richness")
    q=df.loc[df.query_complete].copy()
    cols=["contemporary_not_native_in_catalogue","native_in_catalogue",
          "ambiguous_across_regions","no_region_intersection","unknown_catalogue_status"]
    for name in cols:
        if q[name].isna().any() or (q[name]<0).any():
            raise ValueError("Complete port has missing / negative indexed host status "+name)
    if (int(q.contemporary_not_native_in_catalogue.sum()),
        int(q.native_in_catalogue.sum()),
        int((q.contemporary_not_native_in_catalogue>0).sum()))!=(25742,54156,1150):
        raise ValueError("Source original alien/native host pairs changed")
    return q


def wpi_point_lookup(ports:pd.DataFrame):
    ids=ports.port_id.astype(int).tolist()
    parsed=[]
    for start in range(0,len(ids),85):
        subset=ids[start:start+85]
        response=retry_session_get(WPI,{"where":"OBJECTID IN ("+",".join(map(str,subset))+")",
            "outFields":"OBJECTID,INDEX_NO,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
            "returnGeometry":"false","f":"json"},timeout=75).json()
        if "error" in response:
            raise ValueError("WPI identity lookup error: "+str(response["error"])[:300])
        features=response.get("features",[])
        if len(features)!=len(subset):
            raise ValueError(f"WPI mirror ID count {len(features)} != {len(subset)} for batch {start}")
        parsed.extend(v["attributes"] for v in features)
    w=pd.DataFrame(parsed)
    if w.OBJECTID.nunique()!=1596 or set(w.OBJECTID)!=set(ids):
        raise ValueError("World Port Index mirror didn't preserve source port IDs")
    for p in ("LATITUDE","LONGITUDE"):
        w[p]=pd.to_numeric(w[p],errors="coerce")
    if not np.isfinite(w[["LATITUDE","LONGITUDE"]].to_numpy()).all():
        raise ValueError("WPI coordinates must be finite")
    if ((w.LATITUDE.abs()>90)|(w.LONGITUDE.abs()>180)).any():
        raise ValueError("WPI point coordinate out of WGS84 bounds")
    w["country_mirror"]=w.COUNTRY
    w["port_name"]=w.PORT_NAME
    return w[["OBJECTID","INDEX_NO","LATITUDE","LONGITUDE","country_mirror","port_name"]].rename(columns={"OBJECTID":"port_id"})


def xyz(lat,lon):
    lat=np.radians(np.asarray(lat,dtype=float))
    lon=np.radians(np.asarray(lon,dtype=float))
    a=np.cos(lat)
    return np.column_stack((a*np.cos(lon),a*np.sin(lon),np.sin(lat)))


def ball_radius(km):
    return 2*np.sin(km/(2*EARTH))


def exposure(ports:pd.DataFrame,wwtp:pd.DataFrame)->pd.DataFrame:
    ok=wwtp.loc[wwtp.STATUS.isin(["Operational","Not Reported"])].copy()
    ok=ok.loc[np.isfinite(ok[["LAT_WWTP","LON_WWTP","LAT_OUT","LON_OUT"]].to_numpy(dtype=float)).all(axis=1)].copy()
    if len(ok)<30000:
        raise ValueError(f"Too few active/not-reported WWTP points: {len(ok)}")
    pvec=xyz(ports.LATITUDE,ports.LONGITUDE)
    valid_ww=xyz(ok.LAT_WWTP,ok.LON_WWTP)
    valid_out=xyz(ok.LAT_OUT,ok.LON_OUT)
    whole=cKDTree(valid_ww)
    out=cKDTree(valid_out)
    near_ww=whole.query(pvec,k=1)[0]
    near_out=out.query(pvec,k=1)[0]
    result=ports.copy()
    result["wwtp_nearest_km"]=2*EARTH*np.arcsin(np.clip(near_ww/2,0,1))
    result["outfall_nearest_km"]=2*EARTH*np.arcsin(np.clip(near_out/2,0,1))
    flow=ok.WASTE_DIS.astype(float).to_numpy()
    quality=ok.QUAL_WASTE.astype(float).to_numpy()
    quality_loc=ok.QUAL_LOC.astype(float).to_numpy()
    for radius in RADIUS_KM:
        indices=whole.query_ball_point(pvec,r=ball_radius(radius))
        oindices=out.query_ball_point(pvec,r=ball_radius(radius))
        result[f"wwtp_{radius}km"]=[len(q) for q in indices]
        result[f"outfall_{radius}km"]=[len(q) for q in oindices]
        # Positive treated-water volume only: zero and missing remain separate.
        result[f"wwtp_flow_{radius}km_m3_day"]=[
            float(np.nansum(flow[q][(flow[q]>0) & np.isfinite(flow[q])])) for q in indices]
        result[f"outfall_flow_{radius}km_m3_day"]=[
            float(np.nansum(flow[q][(flow[q]>0) & np.isfinite(flow[q])])) for q in oindices]
        result[f"wwtp_reported_flow_{radius}km_m3_day"]=[
            float(np.nansum(flow[q][(quality[q]==1) & (flow[q]>0) & np.isfinite(flow[q])]))
            for q in indices]
        result[f"wwtp_loc_quality12_{radius}km"]=[int(np.isin(quality_loc[q],[1,2]).sum()) for q in indices]
        result[f"wwtp_reported_level_{radius}km"]=[
            int(np.isin(ok.QUAL_LEVEL.to_numpy()[q],[1]).sum()) for q in indices]
    return result


def ranks(values:np.ndarray)->np.ndarray:
    return pd.Series(values).rank().to_numpy(dtype=float)


def descriptive_rank(data:pd.DataFrame):
    result={}
    for x in ("wwtp_10km","wwtp_25km","wwtp_50km",
              "outfall_25km","wwtp_flow_25km_m3_day","wwtp_nearest_km"):
        result[x]={}
        for y in ("contemporary_not_native_in_catalogue","native_in_catalogue","alien_logratio"):
            good=data[[x,y]].dropna()
            val,p=spearmanr(good[x],good[y]) if len(good)>=10 else (np.nan,np.nan)
            result[x][y]={"n":len(good),"rho":float(val) if np.isfinite(val) else None,
                          "nominal_p":float(p) if np.isfinite(p) else None}
    return result


def fit_fixed_effects(data, exposure_column, response,*,quality="main"):
    q=data.copy()
    q=q.loc[(q.native_in_catalogue+q.contemporary_not_native_in_catalogue)>0]
    q["x"]=np.log1p(q[exposure_column])
    q["abs_lat"]=q.LATITUDE.abs()
    q["log_indexed_known_hosts"]=np.log1p(q.native_in_catalogue+q.contemporary_not_native_in_catalogue)
    q["alien_share_logratio"]=np.log((q.contemporary_not_native_in_catalogue+.5)/(q.native_in_catalogue+.5))
    q["alien_logrichness"]=np.log1p(q.contemporary_not_native_in_catalogue)
    q["native_logrichness"]=np.log1p(q.native_in_catalogue)
    country_n=q.groupby("country").port_id.transform("size")
    q=q.loc[country_n>=3].copy()
    # Fixed effect has to exhibit within-country WWTP variation.
    varied=q.groupby("country").x.transform("nunique")
    q=q.loc[varied>=2].copy()
    if len(q)<200 or q.country.nunique()<15:
        return {"status":"INSUFFICIENT_COUNTRY_WITHIN_VARIATION","n":len(q),"countries":q.country.nunique()}
    f={
        "alien_share_logratio":"alien_share_logratio ~ x + log_indexed_known_hosts + abs_lat + C(country)",
        "alien_logrichness":"alien_logrichness ~ x + abs_lat + C(country)",
        "native_logrichness":"native_logrichness ~ x + abs_lat + C(country)",
    }[response]
    model=smf.ols(formula=f,data=q).fit(cov_type="cluster",cov_kwds={"groups":q.country})
    estimate=float(model.params["x"])
    se=float(model.bse["x"])
    lower=estimate-1.96*se
    upper=estimate+1.96*se
    return {"status":"COUNTRY_FE_CLUSTER_ROBUST_DESCRIPTIVE",
            "n":len(q),"countries":int(q.country.nunique()),
            "exposure":exposure_column,"response":response,
            "beta_per_log1p_wwtp_exposure":estimate,
            "cluster_country_se":se,
            "95pct_country_cluster_normal_CI":[lower,upper],
            "pvalue_nominal":float(model.pvalues["x"]),
            "covariates":"country fixed effects + |latitude| + log(total native+alien indexed host species) for alien/native logratio; raw richness models do not include log total",
            "R2":float(model.rsquared),
            "no_causal_identification":True}


def summarize(df,ww):
    df=df.copy()
    df["alien_logratio"]=np.log((df.contemporary_not_native_in_catalogue+.5)/
                                 (df.native_in_catalogue+.5))
    models={}
    for x in ("wwtp_25km","outfall_25km","wwtp_flow_25km_m3_day",
              "wwtp_reported_flow_25km_m3_day","wwtp_loc_quality12_25km"):
        models[x]={}
        for y in ("alien_share_logratio","alien_logrichness","native_logrichness"):
            models[x][y]=fit_fixed_effects(df,x,y)
    status=ww.STATUS.value_counts(dropna=False).to_dict()
    q={
        "schema":"chocho_hydrowaste_v1_global_WPI_host_GBIF_obs_screen_v1",
        "status":"GLOBAL_OBSERVATIONAL_WWTP_PORT_BUTTERFLY_HOST_RECORD_COGEOGRAPHY",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "source_publication_DOI":"10.5194/essd-14-559-2022",
        "source_HydroWASTE_dataset_DOI":"10.6084/m9.figshare.14847786",
        "source_git_HydroWASTE_v10_mirror_blob":WWTP_BLOB,
        "original_port_source_blob":HOST_BLOB,
        "WWTP_population":len(ww),
        "WWTP_operational_or_not_reported":int(ww.STATUS.isin(["Operational","Not Reported"]).sum()),
        "WWTP_status_counts":{str(k):int(v) for k,v in status.items()},
        "WPI_source_registered_ports":3630,
        "WPI_botanically_index_complete_ports":1596,
        "matched_ports":len(df),
        "distinct_port_countries":int(df.country.nunique()),
        "ports_with_25km_operational_or_notreported_WWTP":int((df.wwtp_25km>0).sum()),
        "ports_with_25km_estimated_wwtp_outfall":int((df.outfall_25km>0).sum()),
        "median_nearest_wwtp_distance_km":float(df.wwtp_nearest_km.median()),
        "median_nearest_outfall_distance_km":float(df.outfall_nearest_km.median()),
        "median_wwtp_within25km":float(df.wwtp_25km.median()),
        "registered_positive_host_candidate_ports":int((df.contemporary_not_native_in_catalogue>0).sum()),
        "total_alien_host_pairs":int(df.contemporary_not_native_in_catalogue.sum()),
        "total_native_host_pairs":int(df.native_in_catalogue.sum()),
        "descriptive_Spearman":descriptive_rank(df),
        "country_FE_models":models,
        "country_observations":df.groupby("country").agg(
            n_ports=("port_id","size"),n_ports_with_wwtp25=("wwtp_25km",lambda x:int((x>0).sum())),
            median_wwtp25=("wwtp_25km","median"),
            median_alien_candidates=("contemporary_not_native_in_catalogue","median"),
        ).reset_index().to_dict(orient="records"),
        "limitations":[
            "HydroWASTE reports facility geography and treated wastewater discharge (m3/day), NOT nitrogen concentration, wastewater-to-soil transport, or soil N deposition.",
            "WWTP latitude/longitude represents facilities; LAT_OUT/LON_OUT represents estimated river outfall not direct soil nitrogen to butterfly host plants.",
            "Some facilities are imputed and locations low quality; WWTP country coverage is strongly heterogeneous.",
            "Operational and Not Reported were included, others excluded; Not Reported status assumed operational by HydroWASTE authors.",
            "A missing/zero nearby WWTP in HydroWASTE does NOT imply no sewage facility or wastewater exposure.",
            "WPI geographic point represents retrospective 2019 source, not historical harbor connectivity or shipping introduction.",
            "GBIF candidate introduced hosts are regional contemporary catalogue nonnative status, not local historical introduction or successful butterfly breeding.",
            "Host counts were indexed during 2001-2024, overlapping and unequally recorded 5km port buffers; zeros are nonrecords not verified absences.",
            "Models use country FE and a weak plant-catalogue size proxy, not direct local urbanization/population/sewerage effort matching.",
            "WWTP and ports both correlate with human settlement; positive coefficients cannot identify an independent wastewater effect.",
            "No matched nonport WWTP exposure or temporally prior treatment history in this global WWTP-port screen."
        ]
    }
    return q


def run(port_path,cache,out_json,out_csv):
    host=read_host_inventory(port_path)
    ww=download_hydrowaste(cache)
    wpi=wpi_point_lookup(host)
    j=host.merge(wpi,on="port_id",validate="one_to_one")
    if len(j)!=1596:
        raise ValueError("World port identity join incomplete")
    exposed=exposure(j,ww)
    result=summarize(exposed,ww)
    out_json.parent.mkdir(parents=True,exist_ok=True)
    out_csv.parent.mkdir(parents=True,exist_ok=True)
    out_json.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    exposed.to_csv(out_csv,index=False)
    print("GLOBAL_HYDROWASTE_PORT_1530_SOURCE_GATE",json.dumps({
        "source":result["WWTP_population"],
        "complete_ports":result["matched_ports"],"countries":result["distinct_port_countries"],
        "near_25km_WWTP":result["ports_with_25km_operational_or_notreported_WWTP"],
        "near_25km_outfall":result["ports_with_25km_estimated_wwtp_outfall"],
        "median_nearest_km":result["median_nearest_wwtp_distance_km"],
        "spearman":result["descriptive_Spearman"],
        "FE":result["country_FE_models"]
    },sort_keys=True),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts",type=Path,default=Path("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"))
    p.add_argument("--hydrowaste-cache",type=Path,default=Path("/tmp/HydroWASTE_v10_2022.csv"))
    p.add_argument("--out-json",type=Path,required=True)
    p.add_argument("--out-csv",type=Path,required=True)
    a=p.parse_args()
    run(a.hosts,a.hydrowaste_cache,a.out_json,a.out_csv)


if __name__=="__main__":
    main()
