"""Source-independent tests: Daymet calendar, NPN effort-only sites, no site CV leakage."""
import sys
import math
from datetime import date, timedelta
from pathlib import Path
from collections import Counter
import numpy as np
import pytest
# The original GEB paper CI does not install optional exploratory ML packages.
# Dedicated Daymet CI installs the fixed version and MUST run all tests.
pytest.importorskip("sklearn")
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_usanpn_daymet_leaf_prediction_20261028 import (
    _safe_daymet_decode,choose_sites,weather_for_site,rolling_weather,fold,
    design_rows,evaluate_prediction,score_sites,model_and_predict,
)


def source_daymet(years=tuple(range(2018,2026))):
    lines=['title: Daymet','year,yday,prcp (mm/day),tmax (deg c),tmin (deg c)']
    for year in years:
        for doy in range(1,366):
            lines.append(f'{year},{doy},2.0,20.0,10.0')
    return ('\n'.join(lines)+'\n').encode()


def test_exact_365_day_calendar_and_weather_features():
    dataset=_safe_daymet_decode(source_daymet())
    assert len(dataset)==8*365
    assert date(2020,12,31) not in dataset
    assert date(2020,2,29) in dataset
    r=rolling_weather(dataset,date(2023,3,15))
    assert r[0]==15 and r[1]==pytest.approx(math.log1p(60)) and r[2]==0
    assert rolling_weather(dataset,date(2018,1,10)) is None
    assert rolling_weather(dataset,date(2021,1,10)) is None  # 2020 Dec31 missing


def test_daymet_fails_closed_on_missing_year():
    with pytest.raises(ValueError,match='complete 8'):
        _safe_daymet_decode(source_daymet(tuple(range(2018,2025))))


def test_site_selection_is_only_effort_not_outcomes():
    points={}
    for site in range(1,61):
        dates={date(2018+(i%3), 1+(i%12), 1+i%25) for i in range(30+site//4)}
        points[site,('Asclepias syriaca' if site%2 else 'Asclepias tuberosa')]=dates
    geo={s:(35.0,-85.0) for s in range(1,61)}
    one,meta=choose_sites(points,geo,max_sites=35)
    two,meta2=choose_sites(points,geo,max_sites=35)
    assert one==two and len(one)==35
    assert meta['source_selected_site_count']==35


def test_scoring_is_site_equal_not_record_equal():
    y=[1,0,1,0]
    sites=[1,1,1,2]
    model=[.9,.1,.9,.9]
    a,b=score_sites(model,y,sites)
    expected=(-math.log(.9)-math.log(.9)-math.log(.9))/3
    assert a>expected  # the single wrong site weighs 1/2, not 1/4
    assert b>0


def test_fold_assignment_is_site_deterministic():
    assert fold(345)==fold(345)
    assert all(0<=fold(x)<5 for x in range(1,100))


def test_small_source_insufficient_not_empty_effect():
    x=[[float(i%2),float(i%5)] for i in range(50)]
    data=[(i%3,'CA','Asclepias syriaca',i%2,x[i],x[i]+[float(i%3)]) for i in range(50)]
    assert evaluate_prediction(data)['status']=='INSUFFICIENT_INDEPENDENT_SITE_TAXON_COVERAGE'


def test_site_model_training_and_prediction_probability():
    x=np.array([[i/10.,i%2] for i in range(30)],dtype=float)
    y=[i%2 for i in range(30)]
    sites=[i%5 for i in range(30)]
    result=model_and_predict(x,y,sites,x)
    assert len(result)==30 and np.isfinite(result).all()
    assert (result>0).all() and (result<1).all()
