from datetime import date
from pathlib import Path
import sys,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from simulate_usanpn_milkweed_cmip6_leaf_window_20261028 import masks_from_observed_site_month,seasonal_degree_days,evaluate,WINTER,ONE_GENERATION_GDD


def case(sp='Asclepias syriaca', site=1):
    rec=[]
    for m,year in [(11,2023),(12,2023),(1,2024),(2,2024)]:
        for d in (2,16):
            rec.append(((sp,'CA',site,site,date(year,m,d)),1))
    return rec


def test_leaf_model_order_and_month_mask():
    # Demonstrate physical and logical bounds under synthetic common climate.
    mins={m:np.array([20.0,22.0]) for m in range(1,13)}
    maxs={m:np.array([22.0,24.0]) for m in range(1,13)}
    g=seasonal_degree_days(mins,maxs)
    assert g[11][0]==(21-11.5)*30
    row={'species':'Asclepias syriaca','site':1,'pixel':(2,3),
         'observed_any_yes_month':{m:1 for m in WINTER},
         'observed_all_yes_only_month':{m:0 for m in WINTER}}
    out=evaluate([row],g,{(2,3):0})
    d=out['Asclepias syriaca']
    assert d['site_equal_proportion_threshold_thermal_only']==1
    assert d['site_equal_proportion_threshold_optimistic_monthly_leaf_mask']==1
    assert d['site_equal_proportion_threshold_conservative_monthly_leaf_mask']==0


def test_repeated_winters_from_same_site_are_not_independent():
    mins={m:np.array([20.]) for m in range(1,13)}
    maxs={m:np.array([22.]) for m in range(1,13)}
    g=seasonal_degree_days(mins,maxs)
    all1={m:1 for m in WINTER}; all0={m:0 for m in WINTER}
    rows=[{'species':'Asclepias syriaca','site':1,'pixel':(0,1),
           'observed_any_yes_month':all1,'observed_all_yes_only_month':all1},
          {'species':'Asclepias syriaca','site':1,'pixel':(0,1),
           'observed_any_yes_month':all0,'observed_all_yes_only_month':all0}]
    out=evaluate(rows,g,{(0,1):0})
    assert out['Asclepias syriaca']['n_independent_sites']==1
    assert out['Asclepias syriaca']['n_site_winter_years']==2
    assert out['Asclepias syriaca']['site_equal_proportion_threshold_optimistic_monthly_leaf_mask']==.5
