"""Frozen original Daymet weather-driver source result regression tests."""
from pathlib import Path
import json,hashlib,math
R=Path(__file__).resolve().parents[1]/'provenance/archive/monarch/results/USANPN_DAYMET_COMPONENT_YEAR_HOLDOUT_20261028.json'
SHA='4a935e60b4a88d81520acf32d7deaec66fee792b3ebb0a169a3f6cd9add76bc3'

def test_original_response_replicated_and_temperature_skill():
    assert hashlib.sha256(R.read_bytes()).hexdigest()==SHA
    x=json.loads(R.read_text());a=x['validation']
    assert x['original_frozen_source']['actual_weather_stations']==72
    assert x['original_frozen_source']['original_usable_plant_visits_with_weather']==20770
    assert a['temperature_conditional_scenario_model_gate']=='PASS_TEMPERATURE_COMPONENT_SPATIAL_AND_TEMPORAL_PREDICTION'
    assert a['site_fold_aggregate']['temperature_and_frost']['positive_folds']==5
    assert a['site_fold_aggregate']['precipitation']['positive_folds']==1
    assert math.isclose(a['site_fold_aggregate']['all_weather']['site_equal_heldout_logloss_improvement'],.025406028335862998,abs_tol=1e-10)

def test_year_holdout_and_privacy_are_not_assumed():
    x=json.loads(R.read_text());t=x['validation']['time_holdout']
    assert t['train_years']==[2018,2019,2020,2021,2022,2023]
    assert t['test_years']==[2024,2025]
    assert t['n_test_sites']==49 and t['n_test_plant_visits']==5198
    assert t['n_test_explicit_yes']==2521 and t['n_test_explicit_no']==2677
    assert t['metrics']['temperature_and_frost']['delta_log_loss_vs_reference']>0
    lower=R.read_text().lower()
    assert all(z not in lower for z in ('"site_id":','"individual_id":','"latitude":','"longitude":','"observer":','"uri":'))
