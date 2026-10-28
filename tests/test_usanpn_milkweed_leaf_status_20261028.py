import json,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_usanpn_milkweed_leaf_status_20261028 import parse_valid,aggregate,get_records

def selected():
    return {'Asclepias syriaca':{'species_id':199,'year_to_leaf_id':{'2024':488}}}

def record(ident=1,day='2024-11-11',state='CA',site=1,ind=1,status=1,spid=199,phid=488):
    return {'observation_id':ident,'observation_date':day,'state':state,'site_id':site,
        'individual_id':ind,'species_id':spid,'phenophase_id':phid,'phenophase_status':status}

def test_real_explicit_no_versus_missing():
    xs=[record(day='2024-11-11',status=0),record(day='2024-11-18',status=1)]
    observed,qc=parse_valid(xs,2024,selected())
    assert len(observed)==2 and sorted([x[1] for x in observed])==[0,1]
    assert qc['accepted_unique_plant_dates']==2
    result=aggregate(observed,qc)
    assert result['winter_data_gate'][0]['two_date_winter_site_years']==1
    assert result['winter_data_gate'][0]['explicit_winter_yes_visits']==1
    assert result['winter_data_gate'][0]['explicit_winter_no_visits']==1
    assert result['winter_data_gate'][0]['status']=='FAIL_NPN_WINTER_SITE_YEAR_INFORMATION_GATE'

def test_absent_month_unobserved_not_imputed_as_zero():
    rows=[record(day='2024-11-11',status=1),record(day='2024-11-20',status=1)]
    accepted,qc=parse_valid(rows,2024,selected())
    result=aggregate(accepted,qc)
    assert len(result['monthly_source'])==1
    assert result['monthly_source'][0]['month']==11
    assert result['monthly_source'][0]['all_yes_site_months']==1


def test_source_nonfocal_species_and_leaf_phase_rejected():
    with pytest.raises(ValueError,match='taxon'):
        parse_valid([record(spid=200)],2024,selected())
    with pytest.raises(ValueError,match='phenophase'):
        parse_valid([record(phid=501)],2024,selected())


def test_unknown_and_conflicting_status_not_negative():
    a=[record(status=1),record(ident=2,status=0),record(ident=3,day='2024-11-18',status='-9999')]
    accepted,qc=parse_valid(a,2024,selected())
    assert accepted==[]
    assert qc['conflicting_same_individual_day']==1
    assert qc['unknown_or_uncertain_status']==1


def test_calendar_cross_year_winter_aggregation():
    entries=[( ('Asclepias syriaca','CA',1,20,__import__('datetime').date(2023,12,12)),1),
             ( ('Asclepias syriaca','CA',1,20,__import__('datetime').date(2024,1,15)),0)]
    x=aggregate(entries,{})
    assert x['winter_data_gate'][0]['two_date_winter_site_years']==1


def test_request_fails_closed_on_unexpected_type():
    with pytest.raises(ValueError,match='JSON list'):
        get_records(2024,selected(),request_fn=lambda:json.dumps({'unexpected':'response'}).encode())
