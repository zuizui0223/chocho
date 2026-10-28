"""All synthetic: no observational ecology inferred from fixed GBIF facets."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_europe_historical_butterfly_gbif_source_20261028 import (
    ORIGINAL_ROSTER,WINDOWS,COUNTRIES,request_url,validate_facets,summarize,
)


def fake(n=40):
    return {'limit':0,'results':[],'count':n,'facets':[
        {'field':'COUNTRY','counts':[{'name':'GB','count':20},{'name':'DE','count':20}]},
        {'field':'BASIS_OF_RECORD','counts':[
            {'name':'HUMAN_OBSERVATION','count':30},{'name':'PRESERVED_SPECIMEN','count':10}]}
    ]}


def test_source_has_frozen_panel_windows_and_exact_query_filters():
    assert len(ORIGINAL_ROSTER)==len(set(ORIGINAL_ROSTER))==39
    assert len(COUNTRIES)==8 and list(WINDOWS.values())==['1970,1990','1991,2005','2010,2024']
    url=request_url(123,'1970,1990')
    assert 'taxonKey=123' in url and 'limit=0' in url and 'facet=country' in url
    assert 'facet=basisOfRecord' in url and 'occurrenceStatus=PRESENT' in url
    with pytest.raises(ValueError):request_url(-1,'1970,1990')
    with pytest.raises(ValueError):request_url(123,'2000,2020')


def test_valid_source_facets_do_not_imply_biological_absence():
    q=validate_facets(fake())
    assert q['count']==40 and q['countries']=={'GB':20,'DE':20}
    assert q['bases']['PRESERVED_SPECIMEN']==10
    assert validate_facets({'count':0,'limit':0,'results':[],'facets':[]})['countries']=={}


def test_broken_sources_fail_not_zero():
    for changed in [dict(fake(),count=30),dict(fake(),results=[{'key':1}]),
                    dict(fake(),facets=[]),dict(fake(),limit=1)]:
        with pytest.raises(ValueError):validate_facets(changed)
    with pytest.raises(ValueError):validate_facets({'count':1,'limit':0,'results':[],'facets':[
        {'field':'COUNTRY','counts':[{'name':'NA','count':2}]},
        {'field':'BASIS_OF_RECORD','counts':[]}]})


def test_reporting_country_window_eligibility_not_true_colonization():
    panel=[{'species':f'Species{i}','taxon_key':i} for i in range(1,27)]
    r=[]
    for p in panel:
        for w in WINDOWS:
            a=20 if w=='pre_1970_1990' else 50
            c={c:a for c in COUNTRIES[:5]}
            r.append({'species':p['species'],'window':w,'retrieval_status':'COMPLETE',
                      'countries':c,'bases':{'HUMAN_OBSERVATION':a},'count':a*5})
    z=summarize(panel,r,completed=True)
    assert z['historical_source_eligibility_gate']=='PASS_SOURCE_COUNTS_ONLY_NOT_OCCUPANCY'
    assert z['eligible_original_species_count']==26
    assert z['countries_with_eligible_species']==5
    assert summarize(panel,r,completed=False)['historical_source_eligibility_gate']=='FAIL_SOURCE_INFORMATION_OR_INCOMPLETE'
    r[0]['countries']['GB']=0
    assert 'GB' in z['paired_species_by_country']


def test_frozen_39_roster_retains_unresolved_without_guessing(tmp_path,monkeypatch):
    import json,hashlib
    import audit_europe_historical_butterfly_gbif_source_20261028 as module
    originals=[]
    for i,name in enumerate(ORIGINAL_ROSTER,1):
        if name=='Lasiommata megera':
            originals.append({'species':name,'status':'TAXON_UNRESOLVED'})
        else:
            originals.append({'species':name,'status':'PARTIAL_FETCH_NOT_ELIGIBLE',
                              'taxon':{'key':i,'match_type':'EXACT'}})
    originals.extend({'species':f'Other filler {i}','status':'COMPLETE',
                      'taxon':{'key':100000+i}} for i in range(200))
    raw=(json.dumps({'schema':module.MANIFEST_SCHEMA,
         'species_requested':239,'species_receipts':originals})+'\n').encode()
    path=tmp_path/'source.json';path.write_bytes(raw)
    monkeypatch.setattr(module,'MANIFEST_SHA256',hashlib.sha256(raw).hexdigest())
    panel,_=module.read_panel(path)
    assert len(panel)==39
    unresolved=[x for x in panel if x['taxon_key'] is None]
    assert unresolved==[{'species':'Lasiommata megera','taxon_key':None,
                        'taxon_identity':None,'original_retrieval_status':'TAXON_UNRESOLVED'}]
    assert sum(type(x['taxon_key']) is int for x in panel)==38
