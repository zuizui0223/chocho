"""Regression holdout to prevent promoting a failed historic-weather result."""
import json,hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'provenance/archive/host_hindcast/results/uk_ltmn_nasa_power_historical_backtest_20261029.json'
DIGEST='18daee3d584983e03822a121601607359da7ed6760d37075bfd33227c895803a'

def test_original_source_checksum_and_negative_predictive_decision():
    blob=SOURCE.read_bytes();d=json.loads(blob)
    assert hashlib.sha256(blob).hexdigest()==DIGEST
    assert d['schema']=='chocho_UK_LTMN_NASA_POWER_1981_2021_outoftime_detection_v01'
    assert d['source']['n_NASA_POWER_source_transects']==19
    assert d['source']['source_original_species']==27
    assert d['source']['source_original_time_split_early']==[1981,1995]
    assert d['source']['source_original_time_split_late']==[2005,2021]
    assert d['source']['NASA_POWER_days_per_site']==14975
    assert d['validation']['five_fold_positive_weather_gain']==1
    assert d['validation']['extrapolation_weather_prediction_gate']=='WEATHER_SPATIOTEMPORAL_HINDCAST_NOT_RELIABLY_BETTER'
    assert d['validation']['late_minus_early_outoftime_original_site_score_improvement']<0
    assert d['ecological_decision']=='CLIMATE_BASELINE_TESTED_ONLY_HOST_EFFECT_NOT_ESTIMATED'

def test_no_historical_host_causal_claim_or_raw_sites():
    d=json.loads(SOURCE.read_text())
    assert d['model']['no_larval_host_presence_or_human_management_predictors']
    assert d['model']['no_claim_host_mediated_range_expansion']
    assert d['model']['missing_species_name_is_conditional_nondetection_not_true_absence']
    s=SOURCE.read_text().lower()
    for bad in ('"ukbms_siteno":','"easting":','"northing":','"latitude":','"longitude":','"observer":'):
        assert bad not in s
