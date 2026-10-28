#!/usr/bin/env python3
"""Posthoc source-locked weather drivers and future-year holdout for NPN leaves.

Requery exactly the PR93 station frame and compare temperature/precipitation
components without original station labels or raw source observations in output.
This is observational prediction, NOT warming causation or future leaf response.
"""
from __future__ import annotations

import argparse,hashlib,json,math,time
from collections import defaultdict
from pathlib import Path

import numpy as np

from audit_usanpn_daymet_leaf_prediction_20261028 import (
 SCHEMA as ORIGINAL_SCHEMA,RETRIEVAL_YEARS,MAX_SITES,SPECIES,
 source_day_records,choose_sites,weather_for_site,design_rows,rolling_weather,
 fold,model_and_predict,score_sites,DAYMET_DOI,
)

SCHEMA='chocho_usanpn_daymet_component_year_holdout_20261028_v01'
YEARS_TRAIN=tuple(range(2018,2024))
YEARS_TEST=(2024,2025)
VARIANTS=('reference','temperature_and_frost','precipitation','all_weather')


def features(row,kind):
    baseline=list(row[4]);weather=list(row[5]);
    if kind=='reference':return baseline
    if kind=='temperature_and_frost':return baseline+[weather[-3],weather[-1]]
    if kind=='precipitation':return baseline+[weather[-2]]
    if kind=='all_weather':return weather
    raise ValueError('Unknown pre-fixed original feature group')


def train_test(train,test):
    y=[r[3] for r in train]
    station=[r[0] for r in train]
    true=[r[3] for r in test]
    held=[r[0] for r in test]
    predictions={}
    for kind in VARIANTS:
        trainx=np.asarray([features(r,kind) for r in train],float)
        testx=np.asarray([features(r,kind) for r in test],float)
        probabilities=model_and_predict(trainx,y,station,testx)
        logloss,brier=score_sites(probabilities,true,held)
        predictions[kind]={'site_equal_log_loss':float(logloss),'site_equal_brier':float(brier)}
    base=predictions['reference']
    for k in VARIANTS[1:]:
        predictions[k]['delta_log_loss_vs_reference']=base['site_equal_log_loss']-predictions[k]['site_equal_log_loss']
        predictions[k]['delta_brier_vs_reference']=base['site_equal_brier']-predictions[k]['site_equal_brier']
    return predictions


def evaluate(rows,years,verify_original=True):
    if len(rows)!=len(years):raise ValueError('Source-year index not aligned with original leaf record')
    byfold=[]
    for fold_id in range(5):
        train=[r for r in rows if fold(r[0])!=fold_id]
        test=[r for r in rows if fold(r[0])==fold_id]
        if len({r[0] for r in test})<5 or len({r[3] for r in test})<2:
            raise ValueError('Original five-fold station QA failed')
        metrics=train_test(train,test)
        byfold.append({'fold':fold_id,'test_sites':len({r[0] for r in test}),
                       'test_plant_visits':len(test),'metrics':metrics})
    nsites=sum(x['test_sites'] for x in byfold)
    if verify_original and nsites!=72:
        raise ValueError('Site holdout frame changed')
    def avg(component):
        return sum(x['test_sites']*x['metrics'][component]['delta_log_loss_vs_reference'] for x in byfold)/nsites
    def positive(component):
        return sum(x['metrics'][component]['delta_log_loss_vs_reference']>0 for x in byfold)
    aggregate={component:{'site_equal_heldout_logloss_improvement':avg(component),
                           'positive_folds':positive(component)}
               for component in VARIANTS[1:]}
    full=aggregate['all_weather']
    # Validate the exact previously reported independent site-fold result.
    if verify_original and not math.isclose(full['site_equal_heldout_logloss_improvement'],
                        .025406028335862998,abs_tol=1e-10):
        raise ValueError('Original source-verified PR93 full weather result not reproduced')
    subset={'train':[], 'test':[]}
    for r,year in zip(rows,years):
        if year in YEARS_TRAIN:subset['train'].append(r)
        elif year in YEARS_TEST:subset['test'].append(r)
        else:raise ValueError('NPN observation year outside frozen source')
    train,test=subset['train'],subset['test']
    no={label:sum(r[3]==label for r in test) for label in (0,1)}
    yr_gate=(len({r[0] for r in test})>=30 and len({r[2] for r in test})>=3
             and len(test)>=1000 and no[0]>=10 and no[1]>=10
             and len({r[0] for r in train})>=30)
    yr={'train_years':list(YEARS_TRAIN),'test_years':list(YEARS_TEST),
        'n_training_plant_visits':len(train),'n_test_plant_visits':len(test),
        'n_test_sites':len({r[0] for r in test}),
        'n_test_species':len({r[2] for r in test}),
        'n_test_explicit_yes':no[1],'n_test_explicit_no':no[0],
        'source_gate':'PASS_TIME_HOLDOUT_INFORMATION' if yr_gate else 'FAIL_TIME_HOLDOUT_INFORMATION',
        'metrics':train_test(train,test) if yr_gate else None,
        'time_holdout_repeats_original_sites_so_not_an_independent_geographic_sample':True}
    temp=aggregate['temperature_and_frost']
    temp_cv=(temp['positive_folds']>=4 and temp['site_equal_heldout_logloss_improvement']>0)
    temp_time=(yr_gate and yr['metrics']['temperature_and_frost']['delta_log_loss_vs_reference']>0)
    permission=bool(temp_cv and temp_time)
    return {
        'site_folds':byfold,'site_fold_aggregate':aggregate,
        'time_holdout':yr,'component_temperature_predictive_cv_gate':temp_cv,
        'temperature_year_holdout_predictive_gate':bool(temp_time),
        'temperature_conditional_scenario_model_gate':(
           'PASS_TEMPERATURE_COMPONENT_SPATIAL_AND_TEMPORAL_PREDICTION'
           if permission else 'FAIL_TEMPERATURE_COMPONENT_SPATIAL_OR_TEMPORAL_PREDICTION'),
        'not_future_leaf_climate_effect_or_monarch_population':True,
    }


def run(registry,receipt,previous,output):
    locked=json.loads(previous.read_text())
    if locked.get('schema')!=ORIGINAL_SCHEMA or locked.get('decision')!='HELD_OUT_WEATHER_PREDICTIVE_INCREMENT_REPLICATED':
        raise ValueError('Prior site-fold original result not source-approved')
    accepted,coord,source_dates,nqc=source_day_records(registry,receipt)
    selected,effort=choose_sites(source_dates,coord)
    if len(selected)!=72:
        raise ValueError('Original 72 selected station frame changed')
    weather={}; digests=[]
    for idx,site in enumerate(sorted(selected),1):
        yy,xx=coord[site]
        met,digest=weather_for_site(yy,xx)
        weather[site]=met;digests.append(digest)
        if idx%12==0:print(json.dumps({'site_weather_reproduced':idx}),flush=True)
        time.sleep(.35)
    totalsha=hashlib.sha256('\n'.join(sorted(digests)).encode()).hexdigest()
    if totalsha!=locked['daymet_aggregate_sha256']:
        raise ValueError('Daymet station source response SHA differs from original PR93 result')
    rows,qc=design_rows(accepted,selected,coord,weather)
    if len(rows)!=20770:
        raise ValueError('Original matched 20770 NPN plant visits changed')
    years=[]
    for (sp,state,site,plant,dt),flag in accepted:
        if site not in weather or rolling_weather(weather[site],dt) is None:continue
        years.append(dt.year)
    assert len(years)==len(rows)
    cv=evaluate(rows,years)
    result={
      'schema':SCHEMA,
      'original_frozen_source':{
        'original_result_sha256':hashlib.sha256(previous.read_bytes()).hexdigest(),
        'original_72_weather_site_aggregate_sha256':totalsha,
        'original_usa_npn_source_status_records':len(accepted),
        'original_usable_plant_visits_with_weather':len(rows),
        'actual_weather_stations':len(weather),
        'actual_host_taxa':len({x[2] for x in rows}),
        'source_years':list(RETRIEVAL_YEARS),
      },
      'models':{'L2_logistic_C':0.1,'baseline_calendar_harmonics':2,
        'weather_window_days':30,'site_equal_weights':True,
        'training_fold_only_standardization':True,
        'Daymet_source_DOI':DAYMET_DOI},
      'validation':cv,
      'inference':{
        'site_folds_are_out_of_station':True,
        'year_holdout_is_out_of_time_same_original_sites':True,
        'weather_prediction_is_not_temperature_causation':True,
        'no_future_leaf_senescence_or_monarch_generation_estimated':True,
        'raw_monitoring_station_coordinates_individuals_and_daily_weather_not_exported':True,
      },
    }
    result['decision']=cv['temperature_conditional_scenario_model_gate']
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'decision':result['decision'],
        'site_components':cv['site_fold_aggregate'],
        'time_holdout_summary':{k:v for k,v in cv['time_holdout'].items() if k not in ('metrics',)},
        'time_metrics':cv['time_holdout'].get('metrics')},sort_keys=True),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--npn-source',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    run(a.registry,a.npn_source,a.previous,a.output)
