"""Synthetic source-independent tests; no actual site or plant outcome used."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / "scripts"))
import numpy as np
import pytest
from datetime import date,timedelta
from simulate_usanpn_daymet_cmip6_leaf_delta_20261028 import (
    MONTHS,monthly_anomalies,climate_features_for_dates,transport_fraction,
    summarize_by_site,season_label,shift_history_to_future,
)


def test_same_product_anomalies_are_month_specific_and_physical():
    baseline={x:{m:np.array([m*1.,m+1.]) for m in MONTHS} for x in ('tmin','tmax')}
    future={x:{m:baseline[x][m]+(2 if x=='tmin' else 3) for m in MONTHS} for x in ('tmin','tmax')}
    delta=monthly_anomalies(baseline,future,np.array([True,True]))
    assert list(delta)==['tmin','tmax']
    assert delta['tmin'][1][1]==pytest.approx(2)
    assert delta['tmax'][12][0]==pytest.approx(3)
    future['tmax'][5][0]+=41
    with pytest.raises(ValueError,match='implausible'):monthly_anomalies(baseline,future,np.array([True,True]))


def test_original_daymet_monthly_delta_recreates_daily_mean_and_frost_no_future_days():
    day=date(2024,1,15)
    weather={day-timedelta(days=k):(-1,9,4) for k in range(30)}
    actual=[(123,'CA','Asclepias tuberosa',1,day)]
    zero={123:{x:{m:0. for m in MONTHS} for x in ('tmin','tmax')}}
    original=climate_features_for_dates(actual,{123:weather},zero)
    assert original.tolist()==[[4.,30.]]
    warming={123:{x:{m:3. for m in MONTHS} for x in ('tmin','tmax')}}
    warm=climate_features_for_dates(actual,{123:weather},warming)
    assert warm.tolist()==[[7.,0.]]
    changed=shift_history_to_future(weather,warming[123]['tmin'],warming[123]['tmax'])
    assert set(changed)==set(weather)
    assert next(iter(changed.values()))==(2.,12.,4.)
    # The **failed original** separate-min/max anomaly would invert this day.
    # The explicitly *post-source revised* mean warming shift must preserve DTR.
    diverging={123:{'tmin':{m:11. for m in MONTHS},'tmax':{m:0. for m in MONTHS}}}
    physical=climate_features_for_dates(actual,{123:weather},diverging)
    assert physical.tolist()==[[9.5,0.]]
    reconstructed=shift_history_to_future(weather,diverging[123]['tmin'],diverging[123]['tmax'])
    assert all((hi-lo)==10. for lo,hi,_ in reconstructed.values())


def test_ood_is_frozen_historical_quantile():
    train=np.array([(x,30-x) for x in np.arange(0,31,dtype=float)])
    future=np.array([[14,16],[32,-1],[-4,34]],float)
    mask,q=transport_fraction(future,train)
    assert mask.tolist()==[False,True,True]
    assert q.shape==(2,2)
    with pytest.raises(ValueError):transport_fraction(np.array([[1,2,3]]),train)


def test_site_equal_not_row_equal_winter_and_no_real_population():
    day=date(2022,12,15)
    dates=[(1,'GA','Asclepias syriaca',1,day)]*100+[(2,'GA','Asclepias syriaca',0,day)]
    base=np.array([.2]*100+[.6]);future=np.array([.4]*100+[.8]);out=np.zeros(101,dtype=bool)
    res=summarize_by_site(dates,base,future,out)
    z=res['Asclepias syriaca::all_year']
    assert z['distinct_monitored_sites']==2
    assert z['source_plant_visits']==101
    assert z['predicted_leaf_yes_current_weather_site_equal']==pytest.approx(.4)
    assert z['predicted_leaf_yes_future_delta_weather_site_equal']==pytest.approx(.6)
    assert z['future_minus_current_predicted_probability']==pytest.approx(.2)
    assert res['Asclepias syriaca::season_winter_Nov_Feb']['transfer_gate']=='INSUFFICIENT_WINTER_SITE_SAMPLE'
    assert season_label(2)=='winter_Nov_Feb'
    assert season_label(8)=='summer_Jun_Aug'


def test_transport_guard_does_not_convert_unsupported_to_population_forecast():
    days=[(i,'GA','Asclepias syriaca',1,date(2023,5,1)) for i in range(1,12)]
    baseline=np.repeat(.4,len(days));future=np.repeat(.7,len(days));out=np.ones(len(days),bool)
    x=summarize_by_site(days,baseline,future,out)['Asclepias syriaca::all_year']
    assert x['transfer_gate']=='UNSUPPORTED_EXTRAPOLATION'
    assert x['future_weather_outside_source_training_domain_site_equal_fraction']==1.
