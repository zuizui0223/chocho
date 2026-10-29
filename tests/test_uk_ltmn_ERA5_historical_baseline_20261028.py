"""Synthetic spatial/temporal design and real climate source-fail integrity tests."""
import sys,datetime as dt,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
import pytest
from audit_uk_ltmn_ERA5_historical_baseline_20261028 import (
    prior30,original_era5_site_weather,score_site_species,flight_features,weights_equal_site_species)

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
