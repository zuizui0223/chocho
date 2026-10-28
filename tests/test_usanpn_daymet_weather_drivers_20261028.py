"""Scientific no-network component and temporal comparison controls."""
import math,sys
from pathlib import Path
import numpy as np
import pytest
pytest.importorskip('sklearn')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_usanpn_daymet_weather_drivers_20261028 import features,train_test,evaluate


def make(site,year,label):
    b=[float(site%4),float(year-2021),math.sin(site),math.cos(site),float(site%2),float(site%5),
       float(site%3),float(year-2018),float(site%6),float(site%7),float(site%11),float(site%13),float(site%17)]
    # Baseline 6 species dummies + seven measured/time vars (length13)
    w=b+[float(label+site%3),float(site%5),float(label*2+site%4)]
    return (site,['CA','OR','TX','AZ','MI','FL'][site%6],
            ['Asclepias syriaca','Asclepias tuberosa','Asclepias incarnata'][site%3],label,b,w)


def test_group_features_temperature_and_precipitation_are_distinct():
    row=make(1,2021,1)
    assert features(row,'reference')==row[4]
    assert features(row,'temperature_and_frost')==row[4]+[row[5][-3],row[5][-1]]
    assert features(row,'precipitation')==row[4]+[row[5][-2]]
    assert features(row,'all_weather')==row[5]
    with pytest.raises(ValueError):features(row,'pick_best_after_outcome')


def test_source_limited_temporal_holdout_reports_information_gate():
    rows=[];years=[]
    for site in range(1,73):
        for year in (2018,2020,2023,2024,2025):
            label=(site+year)%2
            rows.append(make(site,year,label));years.append(year)
    result=evaluate(rows,years,verify_original=False)
    assert result['time_holdout']['source_gate']=='FAIL_TIME_HOLDOUT_INFORMATION'
    assert result['temperature_conditional_scenario_model_gate'].startswith('FAIL_')


def test_train_test_does_not_mutate_source_schema():
    train=[make(s,2020,(s+i)%2) for s in range(1,20) for i in range(2)]
    test=[make(s,2024,s%2) for s in range(20,30)]
    score=train_test(train,test)
    assert set(score)=={'reference','temperature_and_frost','precipitation','all_weather'}
    assert all(math.isfinite(score[z]['site_equal_log_loss']) for z in score)
    assert score['reference']['site_equal_log_loss']>=0
