#!/usr/bin/env python3
"""Apply source-authenticated CMIP6 monthly Tmin/Tmax anomalies to original Daymet
weather, then score frozen historical leaf-visit dates with the validated NPN
site-equal temperature/frost logistic model.

Only conditional source-visit leaf probabilities, no predicted real future leaf
calendar, plant biomass, monarch larvae, survival, abundance or generations.
No individual/site coordinates, observer, dates or source records serialized.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import timedelta
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from audit_usanpn_daymet_leaf_prediction_20261028 import (
    C_REG, RANDOM_SEED, WEATHER_DAYS, RETRIEVAL_YEARS, SPECIES,
    source_day_records, choose_sites, weather_for_site, rolling_weather,
    design_rows,
)
from simulate_monarch_cmip6_thermal_window_20261028 import (
    MODELS, SSPS, PERIODS, VARS, archive_name, grid_cell,
    read_sampled_archive, sha256_file, choose_temperature_scale,
    climate_validity_mask,
)

SCHEMA = 'chocho_usanpn_daymet_cmip6_temperature_delta_leaf_visit_v01'
EXPECTED_SITES = 72
EXPECTED_VISITS = 20770
MONTHS = set(range(1,13))
MIN_PHYS_SUPPORT = 0.95
MAX_OUTSIDE_TRAIN_DOMAIN = 0.10
FUTURE_CASES = [(m,s,p) for m in MODELS for s in SSPS for p in PERIODS]
SEASONS = {
    'winter_Nov_Feb': {11,12,1,2},
    'spring_Mar_May': {3,4,5},
    'summer_Jun_Aug': {6,7,8},
    'autumn_Sep_Oct': {9,10},
}


def train_temperature_model(rows):
    """Exact original PR94 baseline + mean30T + frost30 feature and weights."""
    y=np.asarray([r[3] for r in rows],int)
    if len(rows)!=EXPECTED_VISITS or set(y)!={0,1}:
        raise ValueError('Original NPN 20770 historical weather-joined labels not recovered')
    x=np.asarray([list(r[4])+[r[5][-3],r[5][-1]] for r in rows],float)
    sites=[r[0] for r in rows]
    counts=Counter(sites)
    if len(counts)!=EXPECTED_SITES or not np.isfinite(x).all():
        raise ValueError('Original 72-site fixed empirical model cohort changed')
    weights=np.asarray([len(rows)/len(counts)/counts[site] for site in sites],float)
    scaler=StandardScaler()
    scaled=scaler.fit_transform(x)
    fit=LogisticRegression(C=C_REG,max_iter=1000,solver='lbfgs',random_state=RANDOM_SEED)
    fit.fit(scaled,y,sample_weight=weights)
    historical_p=fit.predict_proba(scaled)[:,1]
    return scaler,fit,x,historical_p


def original_dates_for_frozen_rows(records, selected, daily, rows):
    """Replay original source order without new inclusion/host or effort rule."""
    matched=[]
    for (species,state,site,plant,dt),flag in records:
        if site not in selected or site not in daily:
            continue
        if rolling_weather(daily[site],dt) is None:
            continue
        matched.append((site,state,species,int(flag),dt))
    if len(matched)!=len(rows):
        raise ValueError('Source date rejoin differs from original 20770 included outcomes')
    for z,r in zip(matched,rows):
        if z[:4]!=tuple(r[:4]):
            raise ValueError('Source date model order mismatched in station/plant/status')
    return matched


def monthly_anomalies(base, future, valid):
    """Per-month future minus SAME WorldClim baseline, Celsius, no fake daily series."""
    future_delta={}
    for var in VARS:
        future_delta[var]={}
        for month in MONTHS:
            arr=np.asarray(future[var][month],float)-np.asarray(base[var][month],float)
            if len(arr)!=len(valid) or not np.all(np.isfinite(arr[valid])):
                raise ValueError('CMIP6 baseline/future monthly source not aligned')
            if np.any(np.abs(arr[valid])>40.):
                raise ValueError('Physically implausible future monthly warming anomaly')
            future_delta[var][month]=arr
    return future_delta


def shift_history_to_future(daily, delta_tmin_by_month, delta_tmax_by_month):
    """Only shift original daily temperatures; keep same days and precipitation."""
    shifted={}
    for day,(tmin,tmax,precip) in daily.items():
        a=float(tmin)+float(delta_tmin_by_month[day.month])
        b=float(tmax)+float(delta_tmax_by_month[day.month])
        if not (np.isfinite(a) and np.isfinite(b) and a<=b+1e-3):
            raise ValueError('Future anomaly yields inverted Tmin/Tmax day, refuse imputation')
        shifted[day]=(a,b,float(precip))
    return shifted


def climate_features_for_dates(dates,source_weather,delta_by_station):
    """Read past30d data only. Returned meanT + actual shifted daily frost count."""
    out=np.empty((len(dates),2),float)
    for i,(_,_,_,_,when) in enumerate(dates):
        site=dates[i][0]
        original=source_weather[site]
        delt=delta_by_station[site]
        avg=0.;nfrost=0
        for k in range(WEATHER_DAYS):
            day=when-timedelta(days=k)
            tmin,tmax,prcp=original[day]
            lo=float(tmin)+float(delt['tmin'][day.month])
            hi=float(tmax)+float(delt['tmax'][day.month])
            if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi+1e-3):
                raise ValueError('CMIP6 monthly anomaly creates inverted daily tmin/tmax; no silent repair')
            avg+=.5*(lo+hi)
            nfrost+=int(lo<0)
        out[i,0]=avg/WEATHER_DAYS
        out[i,1]=nfrost
    return out


def transport_fraction(future_weather, historical_training_weather):
    """Out-of-training historical 1%–99% range; no CV/fit on future target."""
    past=np.asarray(historical_training_weather,float)
    future=np.asarray(future_weather,float)
    if past.ndim!=2 or past.shape[1]!=2 or future.shape[1]!=2:
        raise ValueError('Training/support weather features invalid')
    quantile=np.quantile(past,[.01,.99],axis=0)
    out=np.any((future<quantile[0]-1e-8)|(future>quantile[1]+1e-8),axis=1)
    return out,quantile


def season_label(month):
    for name,months in SEASONS.items():
        if month in months:return name
    raise ValueError('Month outside year')


def summarize_by_site(dates,p_baseline,p_future,extrapolated):
    """Equal site, not equal repeated individual visit; no raw IDs in output."""
    rows=defaultdict(list)
    for i,(site,state,sp,flag,when) in enumerate(dates):
        for group in ('all_year','season_'+season_label(when.month)):
            rows[(sp,group,site)].append((float(p_baseline[i]),float(p_future[i]),int(extrapolated[i]),int(flag)))
        rows[('ALL_TAXA','all_year',site)].append((float(p_baseline[i]),float(p_future[i]),int(extrapolated[i]),int(flag)))
        rows[('ALL_TAXA','season_'+season_label(when.month),site)].append((float(p_baseline[i]),float(p_future[i]),int(extrapolated[i]),int(flag)))
    groups=defaultdict(list)
    for (species,season,site),v in rows.items():
        groups[(species,season)].append((
            len(v),sum(x[0] for x in v)/len(v),sum(x[1] for x in v)/len(v),
            sum(x[2] for x in v)/len(v),sum(x[3] for x in v)/len(v)))
    answer={}
    for (sp,season),v in sorted(groups.items()):
        nsites=len(v);visits=sum(x[0] for x in v)
        current=sum(x[1] for x in v)/nsites
        future=sum(x[2] for x in v)/nsites
        support=sum(x[3]*x[0] for x in v)/visits
        site_support=sum(x[3] for x in v)/nsites
        qualified=(site_support<=MAX_OUTSIDE_TRAIN_DOMAIN)
        if season=='season_winter_Nov_Feb' and (nsites<10 or visits<100):
            gate='INSUFFICIENT_WINTER_SITE_SAMPLE'
        elif not qualified:
            gate='UNSUPPORTED_EXTRAPOLATION'
        else:
            gate='CONDITIONAL_MODEL_WITHIN_HISTORICAL_FEATURE_SUPPORT'
        answer[f'{sp}::{season}']={
            'distinct_monitored_sites':nsites,'source_plant_visits':visits,
            'observed_leaf_yes_fraction_site_equal':sum(x[4] for x in v)/nsites,
            'predicted_leaf_yes_current_weather_site_equal':current,
            'predicted_leaf_yes_future_delta_weather_site_equal':future,
            'future_minus_current_predicted_probability':future-current,
            'future_weather_outside_source_training_domain_visit_fraction':support,
            'future_weather_outside_source_training_domain_site_equal_fraction':site_support,
            'transfer_gate':gate,
        }
    return answer


def run(args):
    orig=json.loads(args.original_result.read_text())
    driver=json.loads(args.component_result.read_text())
    climate_orig=json.loads(args.cmip6_original.read_text())
    if orig.get('decision')!='HELD_OUT_WEATHER_PREDICTIVE_INCREMENT_REPLICATED' or \
       driver.get('decision')!='PASS_TEMPERATURE_COMPONENT_SPATIAL_AND_TEMPORAL_PREDICTION':
        raise ValueError('Original temperature predictive performance not independently qualified')
    if driver['original_frozen_source']['original_result_sha256'] != hashlib.sha256(args.original_result.read_bytes()).hexdigest():
        raise ValueError('Frozen original full weather source audit changed')
    if climate_orig.get('status')!='SOURCE_LOCKED_CMIP6_THERMAL_SCENARIO_COMPLETE':
        raise ValueError('Original WorldClim CMIP6 result source is not verified')
    accepted, coordinates, effort, npn_qc=source_day_records(args.registry,args.npn_source)
    selected,selection=choose_sites(effort,coordinates)
    if len(selected)!=EXPECTED_SITES:raise ValueError('Source 72-site effort-only cohort changed')
    daily={};digests=[]
    for i,site in enumerate(sorted(selected),1):
        lat,lon=coordinates[site]
        weather,digest=weather_for_site(lat,lon)
        daily[site]=weather
        digests.append(digest)
        if i%12==0:print(json.dumps({'daymet_source_stations_retrieved':i}),flush=True)
        time.sleep(.25)
    sha=hashlib.sha256('\n'.join(sorted(digests)).encode()).hexdigest()
    if sha!=orig['daymet_aggregate_sha256'] or sha!=driver['original_frozen_source']['original_72_weather_site_aggregate_sha256']:
        raise ValueError('Daymet 72 original daily source response hashes changed')
    modelrows,qc=design_rows(accepted,selected,coordinates,daily)
    if len(modelrows)!=EXPECTED_VISITS:
        raise ValueError('Original PR93 frozen 20770 plant visits not recovered')
    dates=original_dates_for_frozen_rows(accepted,selected,daily,modelrows)
    scaler,fit,training_x,p_original=train_temperature_model(modelrows)
    history_weather=training_x[:,-2:]
    # Zero anomaly MUST be numerically exact: seasonal lookback/day count.
    zero_delta={site:{v:{month:0. for month in MONTHS} for v in VARS} for site in selected}
    zero_features=climate_features_for_dates(dates,daily,zero_delta)
    if not np.allclose(zero_features,history_weather,atol=1e-9,rtol=0):
        raise ValueError('Zero-anomaly source weather is not original 30-day feature history')
    rederived=fit.predict_proba(scaler.transform(np.column_stack([training_x[:,:-2],zero_features])))[:,1]
    if not np.allclose(p_original,rederived,atol=1e-12,rtol=0):
        raise ValueError('Zero-anomaly model probability not identical to historical')
    sites=sorted(selected);original_pixel={site:grid_cell(*coordinates[site]) for site in sites}
    if any(x is None for x in original_pixel.values()):raise ValueError('Source NPN source coordinate off WorldClim')
    pairs=sorted(set(original_pixel.values()))
    pixid={pix:i for i,pix in enumerate(pairs)}
    original_archives={k:v['sha256'] for k,v in climate_orig['archives'].items()}
    archive_source=[];panels={}
    for tag,model,ssp,period in [('historical','baseline','baseline','1970-2000')]+[(m+'-'+s+'-'+p,m,s,p) for m,s,p in FUTURE_CASES]:
        vars={}
        for var in VARS:
            name=archive_name(var) if model=='baseline' else archive_name(var,model,ssp,period)
            path=args.climate_dir/name
            if sha256_file(path)!=original_archives[name]:
                raise ValueError('Original source WorldClim SHA256 changed: '+name)
            raw,scale,digest,readpairs=read_sampled_archive(path,pairs)
            if readpairs!=pairs:raise ValueError('Original Source climate sample pixel order drift')
            if var=='tmax':
                factor=choose_temperature_scale(raw[7],scale.values())
                if tag=='historical':stored_scale=factor
                elif not math.isclose(factor,stored_scale,abs_tol=1e-9):
                    raise ValueError('Future CMIP6 temperature storage factor changed')
            vars[var]=(raw,scale)
            archive_source.append(name)
        for var in VARS:
            raw,sc=vars[var]
            vars[var]={month:raw[month]*stored_scale for month in MONTHS}
        mask,bad=climate_validity_mask(vars['tmin'],vars['tmax'])
        panels[tag]={'model':model,'ssp':ssp,'period':period,'data':vars,'physical_mask':mask,'bad':bad}
    common=np.ones(len(pairs),dtype=bool)
    for q in panels.values():common &=q['physical_mask']
    valid_sites={s for s,pix in original_pixel.items() if common[pixid[pix]]}
    if len(valid_sites)<math.ceil(MIN_PHYS_SUPPORT*EXPECTED_SITES):
        raise ValueError('More than 5 percent of original stations lost to physical CMIP6 source')
    kept_dates=[d for d in dates if d[0] in valid_sites]
    kept_rows=[r for r in modelrows if r[0] in valid_sites]
    keep_ids=[i for i,z in enumerate(dates) if z[0] in valid_sites]
    p_qualified=p_original[keep_ids]
    history_qualified=history_weather[keep_ids]
    baseline=panels['historical']['data']
    results=[]
    zero_mask=np.zeros(len(keep_ids),dtype=bool)
    base_groups=summarize_by_site(kept_dates,p_qualified,p_qualified,zero_mask)
    results.append({'gcm':'historical','ssp':'baseline','period':'1970-2000',
                    'groups':base_groups,'climate_common_site_pixels':int(sum(common))})
    for model,ssp,period in FUTURE_CASES:
        tag=model+'-'+ssp+'-'+period
        delta=monthly_anomalies(baseline,panels[tag]['data'],common)
        bysite={site:{var:{month:delta[var][month][pixid[original_pixel[site]]] for month in MONTHS}
                      for var in VARS} for site in valid_sites}
        future_weather=climate_features_for_dates(kept_dates,daily,bysite)
        outside,train_quant=transport_fraction(future_weather,history_weather)
        # fit/scaler ALWAYS the historical original, never refit future.
        predicted=fit.predict_proba(scaler.transform(np.column_stack([
            training_x[keep_ids,:-2],future_weather])))[:,1]
        groups=summarize_by_site(kept_dates,p_qualified,predicted,outside)
        results.append({'gcm':model,'ssp':ssp,'period':period,
                        'groups':groups,'climate_common_site_pixels':int(sum(common))})
    report={
        'schema':SCHEMA,
        'status':'SOURCE_VERIFIED_CMIP6_DELTA_LEAF_MODEL_CONDITIONAL_SCENARIOS',
        'source':{
            'original_daymet_sha256':sha,
            'original_npn_27057_records':len(accepted),
            'original_effort_only_selected_sites':len(selected),
            'original_weather_model_visits':len(modelrows),
            'original_weather_sites':len({r[0] for r in modelrows}),
            'retained_common_climate_sites':len(valid_sites),
            'common_climate_pixels':int(sum(common)),
            'original_source_archives_verified':len(set(archive_source)),
            'source_2018_2025_response_years_verified':npn_qc['source_complete_years'],
            'original_site_coordinates_or_NPN_plant_id_exported':False,
        },
        'model':{'statistical_model':'historical original temp+frost L2 logistic C=.1',
                 'site_equal_train_weights':True,'fitted_on_original_2018to2025_only':True,
                 'temperature_change':'WorldClim CMIP6 future minus historical monthly tmin/tmax each climate grid',
                 'daily_weather':'original 2018–2025 Daymet daily sequences with monthly temp deltas',
                 'frost_count':'recomputed original daily Tmin plus model monthly Tmin delta',
                 'daily_precipitation_source':'original historical Daymet, held fixed',
                 'no_future_realized_plant_leaf_period_inferred':True,
                 'site_equal_historical_source_visit_sample_not_random':True,
                 'zero_anomaly_feature_and_model_prediction_identity_passed':True,
                 'global_training_feature_temp_frost_quantile01_99':train_quant.tolist(),
                 'historical_five_site_fold_temp_gain':driver['validation']['site_fold_aggregate']['temperature_and_frost']['site_equal_heldout_logloss_improvement'],
                 'historical_2024to2025_temp_gain':driver['validation']['time_holdout']['metrics']['temperature_and_frost']['delta_log_loss_vs_reference']},
        'panels':results,
        'interpretation':{'future_leaf_presence_is_a_model_counterfactual':True,
           'observed_historical_leaf_statuses_and_dates_held_fixed':True,
           'future_daily_temperature_sequence_not_directly_observed_or_GCM_supplied':True,
           'causal_warming_effect_or_future_leaf_senescence_identified':False,
           'true_monarch_larval_survival_or_population_forecast':False,
           'no_raw_monitoring_sites_or_daily_climate_exported':True,
           'label':'Exploratory site-equal at-visit predicted NPN leaf-YES probability; NOT future season length'},
    }
    if len(results)!=9:raise ValueError('Incomplete source-locked ensemble')
    output=args.output
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    small={k:results[-1]['groups'].get(k,{}).get('future_minus_current_predicted_probability')
           for k in ('ALL_TAXA::all_year','ALL_TAXA::season_winter_Nov_Feb',
                     'Asclepias syriaca::all_year','Asclepias tuberosa::all_year')}
    print(json.dumps({'status':report['status'],'original_sites':len(selected),
         'qualifying_sites':len(valid_sites),'result_panels':len(results),
         'latest_2090_high_emission_IPSL_conditional_changes':small,
         'qualifying_original_weather_sites':len({d[0] for d in kept_dates})},sort_keys=True),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--npn-source',type=Path,required=True)
    p.add_argument('--original-result',type=Path,required=True)
    p.add_argument('--component-result',type=Path,required=True)
    p.add_argument('--cmip6-original',type=Path,required=True)
    p.add_argument('--climate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    run(args)
