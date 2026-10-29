"""Synthetic spatial/temporal design and real climate source-fail integrity tests."""
import sys,datetime as dt,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
import pytest
from audit_uk_ltmn_NASA_POWER_historical_baseline_20261028 import (
    prior30,original_nasa_power_site_weather,score_site_species,flight_features,weights_equal_site_species,OTHER_BUTTERFLY_GATE,nasa_site_effort,site_weather_covariate_panel)

def test_30days_only_past_weather():
    day=dt.date(1980,5,31)
    series={day-dt.timedelta(days=k):(10.,2.) for k in range(30)}
    series[day+dt.timedelta(days=1)]=(80.,0.)
    a=prior30(series,day)
    assert a[0]==10.
    assert a[1]==pytest.approx(np.log1p(60))
    del series[day-dt.timedelta(days=29)]
    assert prior30(series,day) is None

def test_source_input_classical_19_original_uk_sites():
    rows=[("S1","Aglais urticae","early",1,[1.0],[1.,2.]),
          ("S1","Aglais urticae","early",0,[1.0],[1.,2.]),
          ("S2","Aglais urticae","early",0,[1.0],[1.,2.])]
    weight=weights_equal_site_species(rows,[0,1,2])
    assert sum(weight)==pytest.approx(3.)
    assert weight[0]==weight[1] and weight[2]>weight[0]

def test_site_species_equal_not_visit_weighted_logloss():
    vals=[("S1","Aglais urticae","late",1,None,None)]*100+[
           ("S2","Aglais urticae","late",0,None,None)]
    pred=np.asarray([.5]*100+[.9])
    s=score_site_species(vals,list(range(101)),pred)
    assert s["site_species_cells"]==2
    assert s["logloss"]==pytest.approx(((-np.log(.5))+(-np.log(1-.9)))/2)

def test_historical_season_features_no_future_outcomes():
    a=flight_features("Aglais urticae",dt.date(1989,5,2),53.,-2.)
    assert len(a)>len([53.,-2.])
    assert a[0]==1.

def test_other_butterfly_visit_gate_constant_and_era_sites():
    assert OTHER_BUTTERFLY_GATE==3
    from datetime import date
    def original_visit(site,year):
        return {"UKBMS_SITENO":site,"YEAR":year,"MONTH":6,"DAY":2,"START_TIME":0.5}
    rows=[original_visit("A",year) for year in (1982,1983,1984,2008,2009,2010)]
    # No two eras survive 10-walk minimum: eligibility is based on visits only.
    source,sites,period=nasa_site_effort(iter(rows),{"A":(200000,300000,("grassland",""))})
    assert sites==set() and len(source)==6
    assert len(period["A"]["early"])==3 and len(period["A"]["late"])==3
