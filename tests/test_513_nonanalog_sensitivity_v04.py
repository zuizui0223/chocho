"""Tests for retained source mask and sample-grain exact climate-extrapolation sensitivity."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_513_2d_host_nonanalog_sensitivity_v04 import (
    count_sign_transitions,source_comparison
)

def test_signed_transitions_need_exactly_identical_species():
    a={'A a':{'delta':-.2},'B b':{'delta':.1},'C c':{'delta':-.1}}
    b={'A a':{'delta':.1},'B b':{'delta':-.05},'C c':{'delta':-.2}}
    assert count_sign_transitions(a,b,'delta')=={
        'decline_to_nondecline':1,'nondecline_to_decline':1,
        'decline_to_decline':1}
    with pytest.raises(ValueError):count_sign_transitions(a,{'A a':{'delta':1}},'delta')

def test_source_already_counted_513_species_and_not_416_plant_artificial_replicates():
    x=(Path(__file__).resolve().parents[1]/'scripts'/
       'audit_513_2d_host_nonanalog_sensitivity_v04.py').read_text()
    assert 'len(groups)!=349' in x
    assert 'len(bysp)!=641' in x
    assert 'source_21ka_predictions_used_only_for_validation' in x

def test_primary_4_and_8_metrics_must_replay_source_with_no_eyeballed_thresholds():
    rows=[{'insect_species':f'Testus species{i}', 'area_delta_km2':-123,
         'component_delta_4':-3,'component_delta_8':-1,
         'largest_fraction_delta_4':-.1,'largest_fraction_delta_8':.05}
         for i in range(513)]
    q=source_comparison(rows,rows,rows)
    assert q['primary_source_replayed']
    assert q['area_direction_transitions']=={'decline_to_decline':513}
    modified=[v.copy() for v in rows]
    modified[11]['area_delta_km2']=0
    with pytest.raises(ValueError):
        source_comparison(modified,rows,rows)
