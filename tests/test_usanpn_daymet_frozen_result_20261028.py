"""Independent static check of the original NPN×Daymet held-out site source result."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT/'provenance/archive/monarch/results/USANPN_DAYMET_HELDOUT_LEAF_20261028.json'
EXPECTED_SHA256='52a0f857202f79db874c512549af54c05b001beea65d68e67b9f2d81b2e1c014'


def load():
    assert hashlib.sha256(RESULT.read_bytes()).hexdigest() == EXPECTED_SHA256
    return json.loads(RESULT.read_text())


def test_complete_original_source_and_no_reclassified_missing_data():
    d=load()
    assert d['schema']=='chocho_usanpn_daymet_site_holdout_leaf_20261028_v01'
    assert d['study_source']['usanpn']['original_accepted_plant_date_records']==27057
    assert d['daymet_complete_sites']==72 and d['daymet_rejected_sites']==0
    assert d['quality_after_complete_weather']['response_days']==20770
    assert d['quality_after_complete_weather']['sites']==72
    assert d['quality_after_complete_weather']['states']==23
    assert d['quality_after_complete_weather']['host_species']==5
    assert d['quality_after_complete_weather']['yes_statuses']==10742
    assert d['quality_after_complete_weather']['no_statuses']==10028


def test_five_site_folds_and_recomputed_primary_site_equal_score():
    d=load()
    f=d['cv']['source_independent_site_folds']
    assert [a['fold'] for a in f] == [0,1,2,3,4]
    assert sum(a['test_sites'] for a in f)==72
    assert sum(a['test_plant_visits'] for a in f)==20770
    assert sum(a['test_yes_visits'] for a in f)==10742
    assert all(a['heldout_log_loss_lift']>0 for a in f)
    delta=sum(a['heldout_log_loss_lift']*a['test_sites'] for a in f)/72
    assert math.isclose(delta,d['cv']['site_equal_delta_log_loss_baseline_minus_daymet'],abs_tol=1e-12)
    assert math.isclose(delta,.025406028335862998,abs_tol=1e-12)
    assert d['decision']=='HELD_OUT_WEATHER_PREDICTIVE_INCREMENT_REPLICATED'


def test_no_unredacted_locations_individuals_or_source_weather():
    lower=RESULT.read_text().lower()
    assert all(x not in lower for x in ('"site_id":','"individual_id":',
                 '"latitude":','"longitude":','"observer":','"uri":'))
    d=load()
    assert d['daymet_data_not_exported']
    assert d['original_plant_site_individual_identifiers_and_locations_not_exported']
    assert d['cv']['no_statistical_p_value_or_demographic_effect']
