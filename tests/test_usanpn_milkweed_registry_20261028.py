from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_usanpn_milkweed_registry_20261028 import candidates,registry_match,leaves_id,discover,SPECIES,YEARS

def test_registry_not_substring_or_common_name():
    a=[{'genus':'Asclepias','species':'syriaca','species_id':81},
       {'genus':'Asclepias','species':'incarnata','species_id':82},
       {'scientific_name':'Asclepias syriaca var. x','species_id':83},
       {'species_name':'Asclepias syriaca','species_id':84}]
    assert registry_match(a,'Asclepias syriaca')==81
    assert registry_match(a,'Asclepias incarnata')==82
    assert registry_match(a,'Asclepias curassavica') is None

def test_one_exact_leaves_phase_only():
    row=[{'species_id':81,'phenophases':[
        {'phenophase_id':1,'phenophase_name':'Initial growth','phenophase_category':'Leaves'},
        {'phenophase_id':2,'phenophase_name':'Leaves','phenophase_category':'Leaves'},
        {'phenophase_id':3,'phenophase_name':'Falling leaves','phenophase_category':'Leaves'}]}]
    assert leaves_id(row,81)==2
    assert leaves_id(row,82) is None

def test_ambiguous_phase_fails_closed():
    data=[{'species_id':1,'phenophases':[
        {'phenophase_id':3,'phenophase_name':'Leaves','phenophase_category':'Leaves'},
        {'phenophase_id':8,'phenophase_name':'Leaves','phenophase_category':'Leaves'}]}]
    with pytest.raises(ValueError):leaves_id(data,1)

def test_missing_species_is_unresolved_not_zero_observations():
    registry=[{'genus':'Asclepias','species':'syriaca','species_id':12}]
    def phase(id,year):return [{'species_id':id,'phenophases':[{'phenophase_id':74,'phenophase_name':'Leaves','phenophase_category':'Leaves'}]}]
    out=discover(registry,phase)
    assert out['Asclepias syriaca']['years_leaves_protocol']==len(YEARS)
    assert out['Asclepias incarnata']['species_id'] is None
    assert out['Asclepias incarnata']['years_leaves_protocol']==0
