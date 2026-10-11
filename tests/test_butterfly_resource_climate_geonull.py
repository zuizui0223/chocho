import csv
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_butterfly_resource_climate_geonull import load_climate, climate_strata

def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=['wgsrpd3_code','bio1','bio12','climate_sample_points'])
        w.writeheader(); w.writerows(rows)

def test_duplicate_species_rows_do_not_weight_regional_climate(tmp_path):
    p=tmp_path/'climate.csv'
    rows=[{'wgsrpd3_code':'A','bio1':'10','bio12':'100','climate_sample_points':'16'}]*3
    rows += [{'wgsrpd3_code':'B','bio1':'20','bio12':'200','climate_sample_points':'1'}]
    write_csv(p,rows)
    c,receipt=load_climate(p)
    assert c=={'A':(10.0,100.0),'B':(20.0,200.0)}
    assert receipt['unique_valid_regions']==2

def test_conflicting_raw_climate_rejected(tmp_path):
    p=tmp_path/'conflict.csv'
    rows=[{'wgsrpd3_code':'A','bio1':str(v),'bio12':'100','climate_sample_points':'16'} for v in [10,20]]
    write_csv(p,rows)
    with pytest.raises(ValueError,match='conflict'):
        load_climate(p)

def test_missing_climate_frozen_per_region_and_continents_separated():
    c={'A':(10.,100.),'B':(10.,100.),'C':(30.,300.),'D':(30.,300.)}
    continents={k:'1' for k in 'ABCDMN'} | {'X':'2'}
    c['X']=(10.,100.)
    strata,report=climate_strata(continents,c,2)
    assert strata['A']==strata['B']
    assert strata['A']!=strata['C']
    assert strata['A']!=strata['X']
    assert strata['M']!=strata['N']
    assert report['missing_climate_regions']==['M','N']

def test_zero_bins_is_matched_support_benchmark():
    cs={'A':(0.,1.),'B':(50.,10.)}
    labels,_=climate_strata({'A':'1','B':'1','M':'1'},cs,0)
    assert labels['A']==labels['B']
    assert labels['A']!=labels['M']

def test_order_invariance_and_ties_are_not_split():
    co={str(i):'1' for i in range(7)}
    c={str(i):(float(i//2),float(i//3)) for i in range(7)}
    a,_=climate_strata(co,c,3)
    b,_=climate_strata(dict(reversed(list(co.items()))),dict(reversed(list(c.items()))),3)
    assert a==b
    assert a['0']==a['1']
