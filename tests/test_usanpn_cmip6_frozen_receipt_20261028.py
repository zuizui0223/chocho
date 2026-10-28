"""Offline no-network controls for source-complete 9-scenario NPN-CMIP6 result."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "provenance/archive/monarch/results/USANPN_CMIP6_FIXED_LEAF_MONTH_20261028.json"
CSV = ROOT / "data/exploratory/usanpn_monarch_cmip6_fixed_leaf_site_summary_20261028.csv"
ORIGINAL_SHA256="69737ce9ff14a5c23daa2b78453eaf1795494f2b7dc436814294b6bb8fd7d0ff"
NPN_SOURCE_RECEIPT_SHA256="6684a9082bcf9e90c4c735e0f26b2287cb98971a585515e8036f24c1bcf1bfc4"


def test_original_aggregate_source_untampered_and_years_complete():
    b=SOURCE.read_bytes()
    assert hashlib.sha256(b).hexdigest()==ORIGINAL_SHA256
    d=json.loads(b)
    assert d['schema']=='chocho_usanpn_same_site_cmip6_fixed_monthly_leaf_scenarios_20261028_v01'
    assert d['status']=='EXPLORATORY_REAL_CMIP6_FIXED_OBSERVED_LEAF_MONTH_SITE_OVERLAY'
    assert d['source']['original_receipt_sha256']==NPN_SOURCE_RECEIPT_SHA256
    assert d['source']['source_year_replies_exact_sha_verified']==8
    assert d['source']['unique_source_records']==27057
    assert d['source_site_winter_gate']=={'Asclepias syriaca':30,'Asclepias tuberosa':21}
    assert len(d['CMIP6_original_source_archive_sha256'])==18
    assert len(d['baseline_and_eight_scenarios'])==9
    assert d['physical_common_site_grid_quality']['retained_valid_grid_pixels']==19
    assert d['physical_common_site_grid_quality']['retained_observed_winter_site_years']==51
    assert d['physical_common_site_grid_quality']['excluded_invalid_grid_pixels']==0
    assert d['inference']['leaf_availability_future_projected'] is False
    assert d['inference']['original_monarch_840_hexes_directly_joined'] is False
    assert d['inference']['observed_NPN_month_leaf_status_is_FIXED_counterfactual'] is True
    raw=b.decode('utf-8').lower()
    for restricted in ('"site_id":','"observer":','"individual_id":','"latitude":','"longitude":','"uri":'):
        assert restricted not in raw


def test_every_scenario_is_on_same_sources_and_gated():
    d=json.loads(SOURCE.read_text())
    observed=set()
    for s in d['baseline_and_eight_scenarios']:
        observed.add((s['model'],s['ssp'],s['period']))
        assert set(s['by_milkweed']) == {'Asclepias syriaca','Asclepias tuberosa'}
        for sp,v in s['by_milkweed'].items():
            assert v['n_site_winter_years']==d['source_site_winter_gate'][sp]
            assert v['n_independent_sites']==(14 if sp.endswith('syriaca') else 9)
            assert v['n_distinct_worldclim_10min_pixels']==(12 if sp.endswith('syriaca') else 9)
            a=v['site_equal_proportion_threshold_thermal_only']
            b=v['site_equal_proportion_threshold_optimistic_monthly_leaf_mask']
            c=v['site_equal_proportion_threshold_conservative_monthly_leaf_mask']
            assert 0<=c<=b<=a<=1
    expected={('1970-2000','baseline','historical')}
    expected|={(m,s,t) for m in ('MIROC6','IPSL-CM6A-LR')
                      for s in ('ssp126','ssp585')
                      for t in ('2041-2060','2081-2100')}
    assert observed==expected


def test_original_csv_is_exact_derivation_no_extra_sites():
    d=json.loads(SOURCE.read_text())
    with CSV.open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==18
    keys={}
    for r in rows:
        k=(r['model'],r['ssp'],r['period'],r['host'])
        assert k not in keys
        keys[k]=r
    for s in d['baseline_and_eight_scenarios']:
        for sp,v in s['by_milkweed'].items():
            row=keys[(s['model'],s['ssp'],s['period'],sp)]
            assert int(row['site_winter_years'])==v['n_site_winter_years']
            assert int(row['distinct_sites'])==v['n_independent_sites']
            for label,key in (
              ('thermal_only_fraction','site_equal_proportion_threshold_thermal_only'),
              ('any_leaf_yes_in_month_fraction','site_equal_proportion_threshold_optimistic_monthly_leaf_mask'),
              ('all_sampled_leaf_yes_in_month_fraction','site_equal_proportion_threshold_conservative_monthly_leaf_mask'),
              ('mean_thermal_gdd','site_equal_mean_thermal_gdd')):
                assert abs(float(row[label])-v[key])<1e-12
