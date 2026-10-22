"""Neotoma source-access pilot: no fossil age or species occupancy inferred."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from neotoma_fossil_taxon_source_pilot_v01 import (
    preflight_spec,describe_json_shape,species_of_source
)

def openapi():
    keys=['taxonname','ageold','ageyoung','datasettype','limit']
    comps={k+'Query':{'name':k,'in':'query'} for k in keys}
    return {'openapi':'3.0.0','components':{'parameters':comps},
        'paths':{'/v2.0/data/taxa':{'get':{}},
                 '/v2.0/data/occurrences':{'get':{
                     'parameters':[{'$ref':'#/components/parameters/'+k+'Query'} for k in keys]}}}}

def spec():
    return {'Neotoma':{'paths':['/v2.0/data/taxa','/v2.0/data/occurrences'],
       'per_query_limit':5,'ageyoung_calibrated_BP':18000,'ageold_calibrated_BP':24000}}

def test_only_source_schema_accepted_with_exact_age_and_typed_filters():
    q=preflight_spec(openapi(),spec())
    assert set(['taxonname','ageold','ageyoung','datasettype','limit']) <= set(q['occurrence_filter_names'])
    x=openapi();x['paths']['/v2.0/data/occurrences']['get']['parameters']=[]
    with pytest.raises(ValueError):preflight_spec(x,spec())
    z=spec();z['Neotoma']['ageold_calibrated_BP']=26000
    with pytest.raises(ValueError):preflight_spec(openapi(),z)

def test_a_matching_taxon_response_shape_is_not_evidence_of_ancient_occupancy():
    values={'status':'success','data':[{'taxonid':123,'taxonname':'Alnus'}]}
    result=describe_json_shape(values)
    assert result['data_list_length']==1
    assert 'taxonid' in result['data_first_keys']
    assert not any('ancient' in x or 'occupied' in x for x in result.keys())

def test_api_empty_list_cannot_be_interpreted_as_species_absent():
    q=describe_json_shape({'data':[],'success':True})
    assert q['data_list_length']==0
    assert 'absence' not in q and 'occupancy' not in q

def test_source_summary_does_not_contain_genotype_or_fossil_claim():
    source=(Path(__file__).resolve().parents[1]/'scripts'/
        'neotoma_fossil_taxon_source_pilot_v01.py').read_text()
    assert 'fossil_species_occurrence_support_verified_n' in source
    assert 'None' in source
    assert 'no_genetic_response_opened' in source
    assert 'not a species occupancy datum' in source
