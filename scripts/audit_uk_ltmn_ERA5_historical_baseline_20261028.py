#!/usr/bin/env python3
"""ERA5 actual weather vs seasonal baseline on 19 fixed repeatedly visited UKBMS transects.

Experiment: learn explicit UKBMS-visit *detection* during 1976–1995, predict
later 2005–2021 visits with fixed original sites, source weather vs no weather.
This does NOT test local larval host identity or human-maintained plants, which
are absent from these original sources. Visit nondetection != true occupancy.
"""
from __future__ import annotations
import argparse,collections,datetime as dt,hashlib,json,math,time
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from urllib.error import URLError,HTTPError

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from pyproj import Transformer

from audit_uk_ltmn_historical_event_join_20261028 import (
    FIELDS_COUNTS,FIELDS_VISITS,get_source,excel_rows,parse_event)
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk
from audit_uk_ltmn_butterfly_detection_panel_20261028 import (
    CROSSWALK,PERIODS,site_effort,count_at_visits,original_taxonomy_ok)

SCHEMA="chocho_UK_19_LTMN_original_visit_ERA5_1976_2021_outoftime_detection_v01"
ERA5_URL="https://archive-api.open-meteo.com/v1/archive"
ALL_YEARS=(1976,2021)
DAYS_EXPECTED=(dt.date(2021,12,31)-dt.date(1976,1,1)).days+1
SEED=20261028
C=.1
OTHER_BUTTERFLY_GATE=3
MIN_SITES=19

def original_era5_site_weather(lat,lon,request_fn=None):
    """Original ERA5 daily, 46 complete years. No site coords are exported."""
    params={"latitude":f"{lat:.5f}","longitude":f"{lon:.5f}",
            "start_date":"1976-01-01","end_date":"2021-12-31",
            "daily":"temperature_2m_mean,precipitation_sum",
            "models":"era5","timezone":"GMT",
            "temperature_unit":"celsius","precipitation_unit":"mm"}
    url=ERA5_URL+"?"+urlencode(params)
    if request_fn is None:
        def request_fn():
            with urlopen(Request(url,headers={"User-Agent":"chocho-original-LTMN-ERA5-hindcast/1.0"}),
                         timeout=85) as r:return r.read(5_000_001)
    body=request_fn()
    if not body or len(body)>5_000_000:
        raise ValueError("Original site weather response missing/oversized")
    d=json.loads(body)
    if d.get("error"):
        raise ValueError("ERA5 original source rejected query")
    dates=d.get("daily",{}).get("time",[])
    temps=d.get("daily",{}).get("temperature_2m_mean",[])
    rain=d.get("daily",{}).get("precipitation_sum",[])
    if len(dates)!=DAYS_EXPECTED or len(temps)!=DAYS_EXPECTED or len(rain)!=DAYS_EXPECTED:
        raise ValueError(f"ERA5 46-year calendar mismatch {len(dates)}/{DAYS_EXPECTED}")
    if dates[0]!="1976-01-01" or dates[-1]!="2021-12-31":
        raise ValueError("ERA5 source ends changed")
    weather={}
    invalid=0
    for when,t,r in zip(dates,temps,rain):
        if not isinstance(t,(int,float)) or not isinstance(r,(int,float)) or (
              not math.isfinite(t) or not math.isfinite(r) or t< -75 or t>65 or r<0 or r>1200):
            invalid+=1;continue
        weather[dt.date.fromisoformat(when)]=(float(t),float(r))
    if len(weather)<.95*DAYS_EXPECTED:
        raise ValueError("ERA5 annual original weather incomplete")
    return weather,hashlib.sha256(body).hexdigest(),invalid

def prior30(weather,when):
    arr=[weather.get(when-dt.timedelta(days=i)) for i in range(30)]
    if any(x is None for x in arr):return None
    return (sum(x[0] for x in arr)/30,
            math.log1p(sum(x[1] for x in arr)))

def flight_features(taxon,when,lat,lon):
    """No post-training data: fixed source species and original survey calendar."""
    doy=when.timetuple().tm_yday
    z=2*math.pi*doy/365.25
    return [float(taxon==sp) for sp in CROSSWALK]+[
        math.sin(z),math.cos(z),math.sin(2*z),math.cos(2*z),
        lat,lon,(when.year-1990)/10]

def site_weather_covariate_panel(visit_events,strict_sites,observed,cross,weather):
    """Yield only original visited site×species labels, never GBIF pseudoabsences."""
    coord=Transformer.from_crs("EPSG:27700","EPSG:4326",always_xy=True)
    coords={s:coord.transform(cross[s][0],cross[s][1])[::-1] for s in strict_sites}
    rows=[]
    source=collections.Counter()
    sp_order=list(CROSSWALK)
    for site in sorted(strict_sites):
        lat,lon=coords[site]
        for era,(y0,y1) in PERIODS.items():
            events=visit_events[site][era]
            for ev in sorted(events):
                sitekey,yr,month,day,clock=ev
                when=dt.date(yr,month,day)
                temp=prior30(weather[site],when)
                if temp is None:
                    source["no_precise_30_day_original_weather"]+=1
                    continue
                names=observed.get(ev,set())
                for k,sp in enumerate(sp_order):
                    common=CROSSWALK[sp]
                    # The detection source for a focal taxon is qualified by
                    # OTHER butterfly records, not the focal outcome itself.
                    if len(names)-int(common in names)<OTHER_BUTTERFLY_GATE:
                        source["unqualified_visit_by_other_butterflies"]+=1
                        continue
                    flag=int(common in names)
                    source[f"eligible_{era}"]+=1
                    source[f"detected_{era}"]+=flag
                    rows.append((site,sp,era,flag,flight_features(sp,when,lat,lon),temp))
    return rows,dict(source)

def design(raw,idx,with_weather):
    X=np.asarray([list(raw[i][4])+(list(raw[i][5]) if with_weather else []) for i in idx],float)
    y=np.asarray([raw[i][3] for i in idx],dtype=int)
    if not np.isfinite(X).all() or len(set(y))<2:
        raise ValueError("Early training source has incomplete predictor or outcome classes")
    return X,y

def weights_equal_site_species(rows,index):
    c=collections.Counter((rows[i][0],rows[i][1]) for i in index)
    return np.asarray([len(index)/(len(c)*c[rows[i][0],rows[i][1]]) for i in index],float)

def score_site_species(raw,indices,p):
    """Equal per source transect×species, not overcounting common visited sites."""
    errors=collections.defaultdict(list)
    brier=collections.defaultdict(list)
    for i,x in zip(indices,p):
        t=raw[i]
        key=(t[0],t[1])
        prob=float(np.clip(x,1e-7,1-1e-7))
        y=t[3]
        errors[key].append(-y*math.log(prob)-(1-y)*math.log1p(-prob))
        brier[key].append((y-prob)**2)
    if not errors:return {"sites":0,"site_species_cells":0,"logloss":None,"brier":None}
    return {"sites":len({s for s,sp in errors}),
            "site_species_cells":len(errors),
            "n_source_visit_species_labels":len(indices),
            "logloss":sum(sum(v)/len(v) for v in errors.values())/len(errors),
            "brier":sum(sum(v)/len(v) for v in brier.values())/len(brier)}

def train_eval(raw,train,testing,with_weather):
    X,y=design(raw,train,with_weather)
    scaler=StandardScaler()
    X=scaler.fit_transform(X)
    model=LogisticRegression(C=C,max_iter=1000,solver="lbfgs",random_state=SEED)
    model.fit(X,y,sample_weight=weights_equal_site_species(raw,train))
    Xt=np.asarray([list(raw[i][4])+(list(raw[i][5]) if with_weather else []) for i in testing],float)
    if not np.isfinite(Xt).all():raise ValueError("Late test source has missing ERA5 features")
    probs=model.predict_proba(scaler.transform(Xt))[:,1]
    return score_site_species(raw,testing,probs)

def perform(raw,sites):
    tr=[i for i,r in enumerate(raw) if r[2]=="early"]
    te=[i for i,r in enumerate(raw) if r[2]=="late"]
    if len(tr)<1500 or len(te)<1500:
        raise ValueError("Insufficient original UKBMS old/new visited source")
    if {r[0] for r in raw if r[2]=="early"} != set(sites) or {r[0] for r in raw if r[2]=="late"}!=set(sites):
        raise ValueError("Historical two-era site identities not equally supported")
    if len(set(raw[i][3] for i in tr))!=2 or len(set(raw[i][3] for i in te))!=2:
        raise ValueError("Source lacks both explicit count detections and conditional non-detections")
    results={}
    for name,temperature in [("no_weather",False),("ERA5_past30_temperature_precip",True)]:
        results[name]=train_eval(raw,tr,te,temperature)
    results["late_minus_early_outoftime_original_site_score_improvement"]=(
        results["no_weather"]["logloss"]-results["ERA5_past30_temperature_precip"]["logloss"])
    # Five disjoint source-site groups: early training never sees held-out site;
    # late testing also never sees the held-out site in training.
    groups=[[] for _ in range(5)]
    shuffled=sorted(sites,key=lambda s:hashlib.sha256((str(SEED)+s).encode()).hexdigest())
    for i,s in enumerate(shuffled):groups[i%5].append(s)
    folds=[]
    for rank,g in enumerate(groups,1):
        if len(g)<3:raise ValueError("Site fold <3 original UK transects")
        pool=set(sites)-set(g)
        training=[i for i in tr if raw[i][0] in pool]
        testing=[i for i in te if raw[i][0] in g]
        pair={label:train_eval(raw,training,testing,flag) for label,flag
              in (("no_weather",False),("ERA5_past30_temperature_precip",True))}
        folds.append({"fold":rank,"n_heldout_entire_source_sites":len(g),
           "trained_early_source_visit_species_labels":len(training),
           "late_test_original_visit_species_labels":len(testing),
           "no_weather_logloss":pair["no_weather"]["logloss"],
           "ERA5_weather_logloss":pair["ERA5_past30_temperature_precip"]["logloss"],
           "weather_increment_logloss":pair["no_weather"]["logloss"]-pair["ERA5_past30_temperature_precip"]["logloss"],
           "no_weather_brier":pair["no_weather"]["brier"],
           "ERA5_weather_brier":pair["ERA5_past30_temperature_precip"]["brier"]})
    good=sum(f["weather_increment_logloss"]>0 for f in folds)
    results["five_independent_unseen_site_and_era_folds"]=folds
    results["five_fold_positive_weather_gain"]=good
    results["extrapolation_weather_prediction_gate"]=(
        "SOURCE_WEATHER_SPATIOTEMPORAL_HINDCAST_BETTER_THAN_BASELINE"
        if good>=4 and results["late_minus_early_outoftime_original_site_score_improvement"]>0
        else "WEATHER_SPATIOTEMPORAL_HINDCAST_NOT_RELIABLY_BETTER")
    return results

def run(args):
    original=json.loads(args.original_manifest.read_text())
    taxon=original_taxonomy_ok(args.original_manifest)
    cross,crossqc=read_crosswalk()
    original_visits,strict,by_site=site_effort(
        excel_rows(get_source("visits"),FIELDS_VISITS),cross)
    if len(strict)!=19:
        raise ValueError("Original source site panel changed from validated 19")
    counted,qc=count_at_visits(excel_rows(get_source("counts"),FIELDS_COUNTS),original_visits,strict)
    transformer=Transformer.from_crs("EPSG:27700","EPSG:4326",always_xy=True)
    weather={};digests=[];bad_source=0
    for i,site in enumerate(sorted(strict),1):
        east,north,habitat=cross[site]
        lon,lat=transformer.transform(east,north)
        if not(49<lat<60 and -9<lon<3):
            raise ValueError("Geocoded original transect outside UK")
        obj,digest,invalid=original_era5_site_weather(lat,lon)
        weather[site]=obj;digests.append(digest);bad_source+=invalid
        if i%4==0:
            print(json.dumps({"original_ERA5_transect_sites_loaded":i,
                              "original_source_days_each":DAYS_EXPECTED}),flush=True)
        time.sleep(.5)
    raw,sources=site_weather_covariate_panel(by_site,strict,counted,cross,weather)
    results=perform(raw,strict)
    out={"schema":SCHEMA,"status":"SOURCE_VERIFIED_UK_ORIGINAL_TRANSECT_HISTORICAL_CLIMATE_BACKTEST",
      "original_taxon_source":taxon,
      "source":{
        "original_UKBMS_visits_sha256":"c23945feb1568b4b60ac139136c485c234699710eae395a104e0be92a958a84d",
        "original_UKBMS_counts_sha256":"6afcc7d03eba0b657ae0e3a0d9da217b88df3ff39da66d7359df7fbb062ead3f",
        "original_UKBMS_site_metadata_sha256":"e38853b775ecadf39332a8412e2a4a2fa3c83512cf7218e4f1eda172dd83c9cd",
        "ERA5_product":"Open-Meteo historical archive, ERA5 model fixed, 1976–2021 daily tmean and precipitation",
        "n_ERA5_source_transects":len(weather),"ERA5_days_per_site":DAYS_EXPECTED,
        "n_ERA5_dated_invalid_days_all_sites":bad_source,
        "original_ERA5_source_response_sha256_sorted_unlinked_to_locations":sorted(digests),
        "source_original_species":len(CROSSWALK),
        "source_original_transects":len(strict),
        "source_paired_visit_species_labels":len(raw),
        "source_observed_detection_labels":sum(x[3] for x in raw),
        "source_original_time_split_early":[1976,1995],
        "source_original_time_split_late":[2005,2021],
        "early_late_source_labels":sources,
      },
      "model":{
        "name":"L2 logistic C=.1, species/DOY harmonic/latlon/year baseline versus + actual prior30 ERA5 weather",
        "training_only_standardization":True,
        "no_site_identity_predictor":True,
        "original_ERA5_weather_at_visited_historic_date_not_contemporary_climatology":True,
        "same_original_19_site_panel_in_train_and_test_years":True,
        "no_all_taxa_full_checklist_certificate":True,
        "missing_species_name_is_conditional_nondetection_not_true_absence":True,
        "no_larval_host_presence_or_human_management_predictors":True,
        "no_claim_host_mediated_range_expansion":True,
        "no_ERA5_raw_daily_series_or_site_geolocation_exported":True,
      },
      "validation":results,
      "ecological_decision":"CLIMATE_BASELINE_TESTED_ONLY_HOST_EFFECT_NOT_ESTIMATED"}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["status"],"train_source_count":sources.get("eligible_early"),
      "test_source_count":sources.get("eligible_late"),
      "overall_delta_logloss":results["late_minus_early_outoftime_original_site_score_improvement"],
      "five_fold_positive_weather_gain":results["five_fold_positive_weather_gain"],
      "decision":results["extrapolation_weather_prediction_gate"]},sort_keys=True),flush=True)
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--original-manifest",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
