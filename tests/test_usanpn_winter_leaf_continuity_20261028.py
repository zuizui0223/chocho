from datetime import date
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_usanpn_winter_leaf_continuity_20261028 import complete_winter_gate,site_winter_months

def test_cross_year_complete_winter_coverage_no_interpolation():
    a=[]
    for month,year in [(11,2023),(12,2023),(1,2024),(2,2024)]:
        for d in (4,18):
            a.append((('Asclepias syriaca','CA',1,100,date(year,month,d)),1))
    x=complete_winter_gate(site_winter_months(a))
    r=next(v for v in x['species_fully_monitored_winter_gate'] if v['species']=='Asclepias syriaca')
    assert r['site_winter_years_all_four_months_two_or_more_dates_each']==1
    assert r['strict_winter_siteyears_all_four_months_any_yes']==1
    assert r['strict_winter_siteyears_all_four_months_only_yes_on_sampled_dates']==1
    assert r['information_gate']=='FAIL_FULL_WINTER_SAMPLING_GATE'


def test_missing_january_is_unknown_not_leaf_absence():
    a=[]
    for month,year in [(11,2023),(12,2023),(2,2024)]:
        for d in (4,18):
            a.append((('Asclepias syriaca','AZ',1,1,date(year,month,d)),0))
    x=complete_winter_gate(site_winter_months(a))
    r=next(v for v in x['species_fully_monitored_winter_gate'] if v['species']=='Asclepias syriaca')
    assert r['site_winter_years_all_four_months_observed']==0
    assert r['site_winter_years_all_four_months_two_or_more_dates_each']==0
    assert r['strict_winter_siteyears_at_least_one_month_only_no']==0


def test_same_state_month_requires_replicated_sites_in_both_species():
    a=[]
    for sp,s in [('Asclepias curassavica',1),('Asclepias syriaca',2)]:
        for d in (4,18):
            a.append(((sp,'TX',s,s,date(2024,1,d)),1))
    x=complete_winter_gate(site_winter_months(a))
    assert x['within_state_winter_month_2018_2025']['candidate_paired_state_winter_month_cells']==1
    assert x['within_state_winter_month_2018_2025']['qualifying_cells_min2_repeated_sites_each_species']==0
