"""Exact source-verified archived CMIP6 840-hex outcome, no external network."""
import csv
import hashlib
import json
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from freeze_monarch_cmip6_summary_20261028 import verify,rows_from_source

RESULT=ROOT/'provenance/archive/monarch/results/monarch_cmip6_thermal_window_20261028.json'
SUMMARY=ROOT/'data/exploratory/monarch_cmip6_thermal_summary_20261028.csv'
SHA='7942096db3a946e11b17d93a9cdac28fbd9ebf3b5808ddfb207ae359dc770e2d'


def load():
    assert hashlib.sha256(RESULT.read_bytes()).hexdigest()==SHA
    return json.loads(RESULT.read_bytes())


def test_exact_840_hex_source_and_eight_future_scenarios():
    a=load()
    assert verify(a)
    assert a['baseline']['all']['n_100km_hex_localities']==840
    assert a['source']['source_hexes_excluded_due_to_any_physical_climate_qc']==1
    assert a['baseline']['all']['winter_thermal_323_count']==74


def test_source_quality_climate_invalid_one_not_imputed():
    a=load()
    r=a['source']['common_source_climate_quality_audit']
    assert r['excluded_pixel_count']==1
    assert r['excluded_pixel_fraction']<.01
    assert r['invalid_source_pixel_is_never_temperature_imputed']
    assert sum(x['invalid_physical_pixel_count'] for x in r['model_month_source_bad_counts'])==1


def test_frozen_counterfactual_leaf_gate_never_increases_thermal_opportunity():
    a=load()
    for r in [a['baseline']['all']]+[s['groups']['all'] for s in a['scenarios']]:
        assert 0<=r['winter_available_50pct_323_fraction']<=r['thermal_fraction_ge_323']<=1
        assert 0<=r['summer_mean_tmax_over38_screen_fraction']<=1
    high=[s for s in a['scenarios'] if s['ssp']=='ssp585' and s['period']=='2081-2100']
    assert len(high)==2
    for s in high:
        assert s['groups']['all']['thermal_fraction_ge_323']>a['baseline']['all']['thermal_fraction_ge_323']


def test_reported_csv_is_deterministic_derivative_of_source_json():
    source=load()
    with SUMMARY.open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    expected=rows_from_source(source)
    assert len(rows)==len(expected)==9
    assert [r['gcm'] for r in rows]==[r['gcm'] for r in expected]
    for r,e in zip(rows,expected):
        for key,v in e.items():
            assert r[key]==str(v),key
