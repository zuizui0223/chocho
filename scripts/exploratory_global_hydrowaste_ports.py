#!/usr/bin/env python3
"""Global port-host / sewage-treatment infrastructure study (descriptive only).

HydroWASTE v1.0: 58,502 facilities; facility/outfall position, status,
treatment, estimated/reported effluent discharge. WPI 2019: 3,630 ports;
the original chocho plant-query complete 1,596 selected BEFORE model fitting.

The WWTP presence / flow / estimated outfall is NOT measured nitrogen loading,
soil enrichment, actual butterfly host quality, marine-to-land exposure,
nor proof of historical transport or biological dispersal.

Original HydroWASTE copyright CC BY 4.0, Ehalt Macedo et al. (2022).
No raw HydroWASTE dataset is vendored into the chocho Git repository.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import random
import statistics
import time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

import numpy as np
import pandas as pd
import pycountry
from scipy.spatial import cKDTree
from scipy.stats import spearmanr,t as student_t

REPO=Path(__file__).resolve().parents[1]
PORTS=REPO/"provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"
WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
HYDRO_URL="https://raw.githubusercontent.com/rfordatascience/tidytuesday/774fb8cc81e563d6d30efdcb7b48a38565e923cc/data/2022/2022-09-20/HydroWASTE_v10.csv"
PORTS_BLOB="fdf73b11ddcb1644f46730123350ae52cc61cae6"
HYDRO_BLOB="0506f873f6541c5f324478894a9731290d722bce"
R_KM=6371.0088
AGENT="chocho-hydrowaste-wpi-global-observational-pilot/1.0"
TERMINAL={"closed","decommissioned","non-operational","projected","proposed","under construction"}
ACTIVE={"operational","not reported"}


def gitblob_sha(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def iso3(value):
    code=str(value or "").strip().upper()
    if code=="XK":return "XKX"
    try:return pycountry.countries.get(alpha_2=code).alpha_3
    except (KeyError,AttributeError):return None


def source_ports():
    if gitblob_sha(PORTS)!=PORTS_BLOB:
        raise ValueError("Original 2026-10-30 port-flora census SHA changed")
    allrows=pd.read_csv(PORTS,keep_default_na=False)  # ISO2 Namibia 'NA' is not a missing value
    assert len(allrows)==3630 and allrows.port_id.nunique()==3630
    port=allrows[allrows.query_complete.astype(str).str.lower()=="true"].copy()
    if len(port)!=1596 or port.port_id.nunique()!=1596:
        raise ValueError("Source-locked 1596 query-complete ports changed")
    categories=("contemporary_not_native_in_catalogue","native_in_catalogue",
                "ambiguous_across_regions","unknown_catalogue_status",
                "no_region_intersection")
    for field in categories:
        vals=pd.to_numeric(port[field],errors="coerce")
        if vals.isna().any() or (vals<0).any() or np.max(np.abs(vals-vals.round()))>1e-5:
            raise ValueError("Incomplete or invalid count "+field)
        port[field]=vals.astype(int)
    port["alien_candidates"]=port["contemporary_not_native_in_catalogue"]
    port["native_candidates"]=port["native_in_catalogue"]
    port["all_taxonomically_indexed"]=sum(port[f] for f in categories)
    port["known_alien_native"]=port["alien_candidates"]+port["native_candidates"]
    port["country_3"]=port.country.map(iso3)
    if port.country_3.isna().sum()>25:
        raise ValueError("Too many unresolvable source ISO-2 countries")
    return port


def wpi_query(params):
    url=WPI+"?"+urlencode(params)
    for attempt in range(4):
        try:
            req=Request(url,headers={"User-Agent":AGENT,"Accept":"application/json"})
            with urlopen(req,timeout=75) as response:
                data=json.load(response)
            if data.get("error"):
                raise RuntimeError(str(data["error"]))
            return data
        except (HTTPError,URLError,TimeoutError):
            if attempt==3:raise
            time.sleep(1.5*(attempt+1))
    raise AssertionError("unreachable")


def wpi_sites(port):
    ids=sorted(port.port_id.astype(int).tolist())
    features={}
    for a in range(0,len(ids),150):
        selected=ids[a:a+150]
        data=wpi_query({"where":"OBJECTID IN ("+",".join(map(str,selected))+")",
              "outFields":"OBJECTID,INDEX_NO,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
              "returnGeometry":"false","f":"json"})
        for feature in data.get("features",[]):
            row=feature["attributes"]
            pid=int(row["OBJECTID"])
            if pid not in selected or pid in features:
                raise ValueError("Duplicate or unexpected World Port Index OBJECTID")
            if row.get("LATITUDE") is None or row.get("LONGITUDE") is None:
                raise ValueError(f"Missing 2019 WPI coordinate at {pid}")
            lat=float(row["LATITUDE"]);lon=float(row["LONGITUDE"])
            if not(-90<lat<90 and -180<=lon<=180):
                raise ValueError(f"Invalid WPI lat/lon {pid}")
            features[pid]=dict(port_id=pid,latitude=lat,longitude=lon,
                port_name=row.get("PORT_NAME"),wpi_country_label=row.get("COUNTRY"),
                wpi_index_number=row.get("INDEX_NO"))
        print("WPI_GLOBAL_ID_FETCH",a,len(features),flush=True)
    if len(features)!=1596 or set(features)!=set(ids):
        raise ValueError("Not all 1596 source WPI port IDs resolved")
    coords=pd.DataFrame(features.values())
    merged=port.merge(coords,on="port_id",validate="one_to_one")
    if not merged.latitude.notna().all():raise ValueError("Join has missing latitudes")
    # WPI country labels are names rather than ISO2. Check any name
    # that pycountry resolves unambiguously against the source ISO2.
    checked=0;mismatches=[]
    for r in merged.itertuples():
        try:
            recognized=pycountry.countries.lookup(str(r.wpi_country_label)).alpha_2
        except LookupError:
            continue
        checked+=1
        if recognized!=r.country:
            mismatches.append({"port_id":int(r.port_id),
                 "wpi_country":r.wpi_country_label,"original_country":r.country})
    if len(mismatches)>2:
        raise ValueError("WPI 2019 port ID identity/country mismatch: "+str(mismatches[:8]))
    return merged,{"country_names_resolvable":checked,"mismatched_country_samples":mismatches}


def hydrowaste_table(path):
    if gitblob_sha(path)!=HYDRO_BLOB:
        raise ValueError("Global HydroWASTE 2022 TidyTuesday snapshot SHA mismatch")
    rows=pd.read_csv(path,low_memory=False)
    if len(rows)!=58502 or rows.WASTE_ID.nunique()!=58502:
        raise ValueError("HydroWASTE v1.0 row-count or ID mismatch")
    need=["LAT_WWTP","LON_WWTP","LAT_OUT","LON_OUT","WASTE_DIS",
          "QUAL_WASTE","QUAL_LOC","STATUS","CNTRY_ISO","LEVEL"]
    if any(k not in rows.columns for k in need):raise ValueError("Hydro fields missing")
    for col in ["LAT_WWTP","LON_WWTP","LAT_OUT","LON_OUT",
                "WASTE_DIS","QUAL_WASTE","QUAL_LOC","POP_SERVED"]:
        rows[col]=pd.to_numeric(rows[col],errors="coerce")
    if rows[["LAT_WWTP","LON_WWTP"]].isna().any().any():
        raise ValueError("Facilities lack source locations")
    rows["status_norm"]=rows.STATUS.astype(str).str.strip().str.lower()
    rows["active_likely"]=rows.status_norm.isin(ACTIVE)
    rows["qual_loc_1or2"]=rows.QUAL_LOC.isin([1,2])
    rows["flow_m3day_nonnegative"]=rows.WASTE_DIS.where(rows.WASTE_DIS>=0)
    if not rows.active_likely.sum()>50000:
        raise ValueError("Unexpectedly few active/assumed active WWTPs")
    return rows


def xyz(lat,lon):
    lat=np.radians(np.asarray(lat,dtype=float))
    lon=np.radians(np.asarray(lon,dtype=float))
    return np.column_stack((np.cos(lat)*np.cos(lon),
            np.cos(lat)*np.sin(lon),np.sin(lat)))


def chord(km):
    return 2*np.sin((km/R_KM)/2)


def wwtp_exposure(ports,hydro, *,low_quality_excluded=False):
    w=hydro[hydro.active_likely]
    if low_quality_excluded:w=w[w.qual_loc_1or2]
    site=xyz(ports.latitude,ports.longitude)
    p=xyz(w.LAT_WWTP,w.LON_WWTP)
    facility=cKDTree(p)
    d,_=facility.query(site,k=1)
    dist_fac=2*R_KM*np.arcsin(np.clip(d/2,0,1))
    finiteout=np.isfinite(w.LAT_OUT)&np.isfinite(w.LON_OUT)&(w.LAT_OUT.between(-90,90))&w.LON_OUT.between(-180,180)
    wout=w[finiteout]
    outtree=cKDTree(xyz(wout.LAT_OUT,wout.LON_OUT)) if len(wout) else None
    d_out=2*R_KM*np.arcsin(np.clip(outtree.query(site,k=1)[0]/2,0,1)) if outtree else np.full(len(site),np.nan)
    result=pd.DataFrame({"port_id":ports.port_id.astype(int).to_numpy(),
              "nearest_operational_WWTP_km":dist_fac,
              "nearest_estimated_outfall_km":d_out})
    for rad in (5,10,25,50):
        neighbors=facility.query_ball_point(site,r=chord(rad))
        result[f"WWTP_count_{rad}km"]=np.array([len(a) for a in neighbors],dtype=int)
        if rad in (10,25,50):
            flows=w.flow_m3day_nonnegative.to_numpy(dtype=float)
            result[f"WWTP_daily_flow_sum_m3_{rad}km"]=np.array([
                 float(np.nansum(flows[x])) for x in neighbors])
            result[f"WWTP_reported_flow_sites_{rad}km"]=np.array([
                 int(np.isfinite(flows[x]).sum()) for x in neighbors],dtype=int)
        if rad in (10,25,50):
            nearout=outtree.query_ball_point(site,r=chord(rad)) if outtree else [[] for _ in range(len(site))]
            result[f"outfall_count_{rad}km"]=[len(x) for x in nearout]
    return result


def predeclared_eligible(r, h):
    """Select country strata by port sampling availability and HydroWASTE coverage.

    Never select countries by butterfly/plant response, regression direction,
    or nominal p-value. Country-level reporting differences remain.
    """
    countryhydro=Counter(h.loc[h.active_likely,"CNTRY_ISO"].astype(str).str.strip().str.upper())
    r=r.copy()
    r["n_hydrowaste_active_country"]=r.country_3.map(lambda c:countryhydro.get(c,0))
    strata=r.groupby("country").agg(nport=("port_id","size"),
         wwcount_variation=("WWTP_count_25km","nunique"),
         source_wwtp_country=("n_hydrowaste_active_country","max"))
    good=strata[(strata.nport>=5)&(strata.wwcount_variation>=2)&(strata.source_wwtp_country>=5)].index
    r["in_within_country_comparison"]=r.country.isin(good)
    return r, {"countries_in_complete_port_data":int(r.country.nunique()),
       "countries_in_HydroWASTE_covered_population":int((strata.source_wwtp_country>=5).sum()),
       "countries_with_within_country_exposure_variation_and_5plus_ports":len(good),
       "n_within_country_ports":int(r.in_within_country_comparison.sum()),
       "source_countries_eligible":list(sorted(good)),
       "country_census":strata.reset_index().to_dict("records")}


def residual_yx(x,y,df,*,with_sampling_control):
    country=pd.get_dummies(df.country,drop_first=True,dtype=float)
    lat=np.abs(df.latitude.to_numpy(dtype=float))/60.
    lon=np.radians(df.longitude.to_numpy(dtype=float))
    controls=[np.ones(len(x)),lat,np.sin(lon),np.cos(lon)]
    if with_sampling_control:
        controls.append(np.log1p(df.all_taxonomically_indexed.to_numpy(dtype=float)))
    nuisance=np.column_stack(controls+[country.to_numpy()])
    U,s,_=np.linalg.svd(nuisance,full_matrices=False)
    Q=U[:,s>max(s[0],1)*1e-9]
    rx=x-Q@(Q.T@x)
    ry=y-Q@(Q.T@y)
    return rx,ry,Q


def within_country_OLS(df,predictor,response,*,adjust_recorded_total=False,
                       n_perm=299,seed=2217):
    if len(df)<40 or df.country.nunique()<4:
        return {"status":"UNDERPOWERED_COUNTRY_SUPPORT",
                "n_ports":len(df),"n_countries":int(df.country.nunique())}
    x=df[predictor].to_numpy(dtype=float)
    y=df[response].to_numpy(dtype=float)
    rx,ry,Q=residual_yx(x,y,df,with_sampling_control=adjust_recorded_total)
    denom=float(rx@rx)
    if denom<1e-8: return {"status":"NO_WITHIN_COUNTRY_EXPOSURE_VARIATION"}
    beta=float(rx@ry/denom)
    resid=ry-beta*rx
    group=df.country.to_numpy()
    groups=sorted(set(group))
    G=len(groups)
    contributions=np.array([sum(rx[group==c]*resid[group==c]) for c in groups])
    n=len(x);k=Q.shape[1]
    robust_var=(G/(G-1))*((n-1)/(max(1,n-k-1)))*(contributions@contributions)/denom**2
    se=math.sqrt(float(robust_var))
    ci=(beta-student_t.ppf(.975,G-1)*se,beta+student_t.ppf(.975,G-1)*se) if se>0 else (beta,beta)
    indices=[np.flatnonzero(group==country) for country in groups]
    rng=np.random.default_rng(seed)
    null=[]
    for _ in range(n_perm):
        xp=x.copy()
        for ids in indices:
            xp[ids]=rng.permutation(xp[ids])
        pr=xp-Q@(Q.T@xp)
        if float(pr@pr)>1e-10:
            null.append(float(pr@ry/(pr@pr)))
    perm_p=(1+sum(abs(v)>=abs(beta) for v in null))/(len(null)+1) if null else None
    return {"status":"ESTIMATED_EXPLORATORY_COUNTRY_FE",
      "response":response,"predictor":predictor,
      "n_ports":n,"n_countries":G,
      "coefficient_per_log1p_exposure":beta,
      "cluster_country_robust_se":se,
      "cluster_country_robust_approx95":list(ci),
      "within_country_label_permutation_two_sided_p_exploratory":perm_p,
      "within_country_permutations":len(null),
      "percent_change_in_alien_native_ratio_0_to_1_WWTP":(
          (math.exp(beta*math.log(2))-1)*100 if response=="ln_alien_to_native" else None),
      "adjusted_for":["country fixed effects","absolute latitude","sin longitude","cos longitude"]+(
          ["log all port-indexed plant host categories (sampling proxy; includes response)"]
          if adjust_recorded_total else []),
      "spatial_autocorrelation_unadjusted":True}


def run(hydro_path):
    ports=source_ports()
    wpi,identity=wpi_sites(ports)
    ww=hydrowaste_table(hydro_path)
    features=wwtp_exposure(wpi,ww)
    highquality=wwtp_exposure(wpi,ww,low_quality_excluded=True)
    highquality=highquality[["port_id","WWTP_count_25km","nearest_operational_WWTP_km"]].rename(
        columns={"WWTP_count_25km":"WWTP_count_25km_good_location",
                 "nearest_operational_WWTP_km":"nearest_good_location_WWTP_km"})
    df=wpi.merge(features,on="port_id",validate="one_to_one").merge(highquality,on="port_id",validate="one_to_one")
    if len(df)!=1596:raise ValueError("Port census changes after WWTP join")
    df,eligible=predeclared_eligible(df,ww)
    df["ln_alien_to_native"]=np.log((df.alien_candidates+.5)/(df.native_candidates+.5))
    df["ln_alien_count"]=np.log1p(df.alien_candidates)
    df["ln_native_count"]=np.log1p(df.native_candidates)
    df["ln_WWTP_25km"]=np.log1p(df.WWTP_count_25km)
    df["ln_WWTP_10km"]=np.log1p(df.WWTP_count_10km)
    df["ln_WWTP_50km"]=np.log1p(df.WWTP_count_50km)
    df["ln_WWTP_25km_good_location"]=np.log1p(df.WWTP_count_25km_good_location)
    df["ln_WWTP_discharge_25km"]=np.log1p(df.WWTP_daily_flow_sum_m3_25km)
    df["ln_outfall_count_25km"]=np.log1p(df.outfall_count_25km)
    df["WWTP_any_25km"]=df.WWTP_count_25km.gt(0).astype(float)
    subset=df.loc[df.in_within_country_comparison].copy()
    if len(subset)<250 or subset.country.nunique()<10:
        raise ValueError("Insufficient scientifically interpretable country-stratified global port support")
    outcomes=["ln_alien_to_native","ln_alien_count","ln_native_count"]
    predictors=["ln_WWTP_25km","WWTP_any_25km","ln_WWTP_10km",
                "ln_WWTP_50km","ln_WWTP_25km_good_location",
                "ln_WWTP_discharge_25km","ln_outfall_count_25km"]
    models={}
    for pre in predictors:
        models[pre]={}
        for resp in outcomes:
            models[pre][resp]=within_country_OLS(subset,pre,resp,
                        adjust_recorded_total=False,
                        n_perm=399 if pre=="ln_WWTP_25km" else 99,
                        seed=1534)
        if pre in ("ln_WWTP_25km","ln_WWTP_discharge_25km"):
            models[pre]["ln_alien_to_native_plus_observation_proxy"]=within_country_OLS(
                subset,pre,"ln_alien_to_native",
                adjust_recorded_total=True,n_perm=399,seed=551)
        print("HYDROWASTE_PORT_MODEL",pre,json.dumps(models[pre],sort_keys=True),flush=True)
    country_stats=[]
    for country,s in subset.groupby("country"):
        x=s.WWTP_count_25km
        if len(s)>=6 and x.nunique()>=2:
            rho=spearmanr(x,s.ln_alien_to_native).statistic
            country_stats.append({"country":country,"n_ports":len(s),
                                  "rho_wwtp25_vs_alien_native_logratio":float(rho) if np.isfinite(rho) else None,
                                  "n_active_WWTP_catalogued_in_country":int(s.n_hydrowaste_active_country.iloc[0])})
    country_stats.sort(key=lambda x:-x["n_ports"])
    by_near={}
    for k,group in (("near_any_25km",df[df.WWTP_count_25km>0]),
                    ("no_catalogued_25km",df[df.WWTP_count_25km==0])):
        by_near[k]={"n_ports":len(group),
             "median_alien_host_candidates":float(group.alien_candidates.median()) if len(group) else None,
             "median_native_host_candidates":float(group.native_candidates.median()) if len(group) else None,
             "median_known_alien_native_ratio_with_half_smoothing":float(
                ((group.alien_candidates+.5)/(group.native_candidates+.5)).median()) if len(group) else None}
    summary={
       "schema":"chocho_WPI1596_HydroWASTE58502_countryFE_alien_host_observational_v1",
       "status":"GLOBAL_SPATIAL_ASSOCIATION_EXPLORATORY_NO_SEWAGE_NITROGEN_CAUSAL_INFERENCE",
       "source_versions":{
          "ports_source_repo_blob":PORTS_BLOB,
          "HydroWASTE_csv_repo_blob":HYDRO_BLOB,
          "HydroWASTE_csv_pinned_commit":"774fb8cc81e563d6d30efdcb7b48a38565e923cc",
          "original_HydroWASTE_paper":"https://doi.org/10.5194/essd-14-559-2022",
          "WPI_feature_server":WPI,
          "downloaded_at_utc":datetime.now(timezone.utc).isoformat()
        },
       "population":{"WPI2019_ports_full":3630,
            "WPI_ports_query_complete":1596,
            "completed_port_countries":int(df.country.nunique()),
            "HydroWASTE_full_facility_count":58502,
            "HydroWASTE_active_or_not_reported":int(ww.active_likely.sum()),
            "HydroWASTE_operational_reported":int((ww.status_norm=="operational").sum()),
            "HydroWASTE_location_qualities_1_2":int(ww.qual_loc_1or2.sum()),
            "HydroWASTE_estimated_outfall_valid":int((ww.LAT_OUT.notna()&ww.LON_OUT.notna()).sum()),
            "source_original_WPI_identity":identity,
            "predeclared_country_comparison":eligible},
       "unadjusted_port_source_records":by_near,
       "global_port_nearest_WWTP_km_median":float(df.nearest_operational_WWTP_km.median()),
       "global_port_nearest_estimated_outfall_km_median":float(df.nearest_estimated_outfall_km.median()),
       "global_WWTP_any_10km_ports":int((df.WWTP_count_10km>0).sum()),
       "global_WWTP_any_25km_ports":int((df.WWTP_count_25km>0).sum()),
       "global_WWTP_any_50km_ports":int((df.WWTP_count_50km>0).sum()),
       "country_stratified_models":models,
       "within_country_descriptive":country_stats,
       "analytical_scope_and_limits":[
           "WWTP presence is a human wastewater treatment infrastructure proxy and may mainly represent urban size, electricity, sanitation or observation intensity.",
           "HydroWASTE covers known/reportable facilities unevenly across the world. No facility within 25km is not the biological absence of sewage, decentralized treatment or discharges.",
           "WWTP coordinates and modeled river outfalls differ; outfalls do not establish deposition of nutrient pollutants in terrestrial port soils.",
           "WASTE_DIS (m3/day) is treated wastewater discharge, NOT nitrogen load or treatment efficiency/total N concentration.",
           "HydroWASTE and 2019 WPI are modern retrospective sources and cannot explain 2001-2010 botanical arrivals or historical shipping.",
           "The plant response consists of indexed candidate native/non-native *species* counts, not butterfly breeding, host nutritional quality or independent floristic abundance.",
           "Port catalogue country label and facility-country coverage differ; country FE remove only mean between-country heterogeneity.",
           "Latitude, crude longitude and total-indexed-host counts are partial observation/geography controls. Port built-up land cover, port traffic, independent botanical sampling and climate were not fully adjusted.",
           "Within-country shuffled exposures break possible distance correlations and do not preserve spatial autocorrelation; permutation results are exploratory.",
           "Overlapping WWTP search buffers among nearby ports mean that port observations are not independent, even within country.",
           "Intra-country geographic/ecological variation and incomplete plant source availability remain alternative explanations.",
           "WWTP source quality indicators may represent estimated, rather than measured, facility attributes.",
           "No claim of causal wastewater enrichment, nitrogen delivery, ship transport, larval host use, or butterfly colonization can be made."
       ]
    }
    return summary,df


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hydro",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--port-table",type=Path,required=True)
    args=ap.parse_args()
    summary,table=run(args.hydro)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    args.port_table.parent.mkdir(parents=True,exist_ok=True)
    table.to_csv(args.port_table,index=False,float_format="%.8g")
    d=summary["population"]
    print("GLOBAL_HYDROWASTE_PORT_GATE",json.dumps({
        "ports":d["WPI_ports_query_complete"],
        "facilities":d["HydroWASTE_full_facility_count"],
        "active":d["HydroWASTE_active_or_not_reported"],
        "eligible_ports":d["predeclared_country_comparison"]["n_within_country_ports"],
        "eligible_countries":d["predeclared_country_comparison"]["countries_with_within_country_exposure_variation_and_5plus_ports"],
        "near25":summary["global_WWTP_any_25km_ports"],
        "primary_model":summary["country_stratified_models"]["ln_WWTP_25km"]["ln_alien_to_native"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
