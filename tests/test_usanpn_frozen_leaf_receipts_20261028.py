"""Offline original USA-NPN source receipt and policy consistency regression."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'provenance/archive/monarch/results/USANPN_MILKWEED_STATUS_2018_2025.json'
WIN=ROOT/'provenance/archive/monarch/results/USANPN_WINTER_CONTINUITY_2018_2025.json'
SOURCE_SHA='6684a9082bcf9e90c4c735e0f26b2287cb98971a585515e8036f24c1bcf1bfc4'
WIN_SHA='72e65ff5afab84e2c069c5828972b5412ad99f5d2fa859e55d7475bd05b7c87d'

def test_source_identical_and_full_8_year_frame():
    assert hashlib.sha256(SRC.read_bytes()).hexdigest()==SOURCE_SHA
    x=json.loads(SRC.read_text())
    assert x['status']=='COMPLETE_PHENOPHASE_SOURCE_WITH_INFORMATION_GATES'
    assert len(x['requests'])==8
    assert {z['year'] for z in x['requests']}==set(range(2018,2026))
    assert sum(r['source_records'] for r in x['requests'])==30476
    assert x['n_accepted_observed_individual_days']==27057
    assert x['unique_site_months']==5620
    assert x['no_raw_locations_observer_names_or_plant_id_exported']
    assert 'latitude' not in x and 'longitude' not in x and 'site_id' not in x

def test_coupled_gate_matches_winter_coverage_and_no_fabricated_comparison():
    assert hashlib.sha256(WIN.read_bytes()).hexdigest()==WIN_SHA
    x=json.loads(WIN.read_text());y=json.loads(SRC.read_text())
    assert x['retained_source_individual_days']==y['n_accepted_observed_individual_days']
    assert x['source_year_sha256']==[z['response_sha256'] for z in y['requests']]
    w={z['species']:z for z in x['species_fully_monitored_winter_gate']}
    assert len(w)==6
    assert {k for k,v in w.items() if v['information_gate'].startswith('PASS')}=={
        'Asclepias syriaca','Asclepias tuberosa'}
    assert w['Asclepias syriaca']['site_winter_years_all_four_months_two_or_more_dates_each']==30
    assert w['Asclepias tuberosa']['site_winter_years_all_four_months_two_or_more_dates_each']==21
    assert w['Asclepias curassavica']['site_winter_years_all_four_months_two_or_more_dates_each']==1
    cmp=x['within_state_winter_month_2018_2025']
    assert cmp['qualifying_cells_min2_repeated_sites_each_species']==0
    assert cmp['gate']=='INSUFFICIENT_SHARED_STATE_COMPARABILITY'
    assert x['not_directly_merged_with_monarch_840_hexes']
