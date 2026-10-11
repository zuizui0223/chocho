"""Keep immature occurrence, larval feeding and persistence as different claims."""
import importlib.util
from pathlib import Path
import pytest
from shapely.geometry import box
P = Path(__file__).resolve().parents[1]/'scripts/acquire_butterfly_larval_evidence.py'

def module():
    assert P.exists(), 'Larval evidence acquisition module is not implemented'
    s=importlib.util.spec_from_file_location('larval_evidence',P)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def record(**kw):
    d={'id':1,'quality_grade':'research','captive':False,'obscured':False,
       'geojson':{'coordinates':[0.5,0.5]},'positional_accuracy':100,
       'annotations':[{'controlled_attribute_id':1,'controlled_value_id':6,'vote_score':1}],
       'taxon':{'id':100,'ancestor_ids':[]},'observed_on':'2024-04-10','ofvs':[]}
    d.update(kw);return d

def test_adult_observation_cannot_pass_larval_gate():
    assert module().eligible(record(annotations=[]),100)[0] is False

def test_obscured_or_captive_cannot_pass():
    m=module()
    assert not m.eligible(record(obscured=True),100)[0]
    assert not m.eligible(record(captive=True),100)[0]

def test_larva_in_introduced_host_region_is_not_automatically_feeding():
    m=module(); mapper=m.RegionMapper(['A','B'],[box(0,0,1,1),box(2,0,3,1)])
    row=m.summarize_record(record(), 'Test species',100,mapper,{'B'},{'A','B'})
    assert row['resource_class']=='introduced_only_in_frozen_hosts'
    assert row['evidence_level']=='L1_LARVAL_PRESENCE_NOT_HOST_USE'
    assert row['coordinate_precision']=='known_le_10000m'

def test_boundary_and_no_known_hosts_are_not_assigned_arbitrarily():
    m=module(); mapper=m.RegionMapper(['A','B'],[box(0,0,1,1),box(1,0,2,1)])
    assert mapper.map_point(1,0.5)==(None,'multiple_polygons')
    row=m.summarize_record(record(), 'Test species',100,mapper,set(),set())
    assert row['resource_class']=='no_known_host_mapped'

def test_subspecies_and_exact_taxon_resolution():
    m=module()
    assert m.eligible(record(taxon={'id':101,'ancestor_ids':[100]}),100)[0]
    assert m.resolve('Test species',[{'id':100,'name':'Test species','rank':'species'}])['id']==100
    with pytest.raises(ValueError): m.resolve('Test species',[{'id':200,'name':'Other species','rank':'species'}])

def test_missing_accuracy_is_not_reported_as_high_precision():
    m=module(); mapper=m.RegionMapper(['A'],[box(0,0,1,1)])
    row=m.summarize_record(record(positional_accuracy=None),'Test species',100,mapper,set(),{'A'})
    assert row['coordinate_precision']=='unknown'
