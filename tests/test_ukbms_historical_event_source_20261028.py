"""No ecological claim based on a schema-only source probe."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_ukbms_historical_event_source_20261028 import (
    validate_year_facets,field_presence,DATASET,YEARS
)


def test_fixed_source_and_year_windows():
    assert len(DATASET)==36
    assert list(YEARS.values())==['1970,1990','2010,2024']


def test_real_year_facet_validation():
    data={'limit':0,'count':14,'results':[],
          'facets':[{'field':'YEAR','counts':[{'name':'1980','count':4}, {'name':'1982','count':10}]}]}
    x=validate_year_facets(data)
    assert x['records']==14 and x['years']=={'1980':4,'1982':10}
    assert validate_year_facets({'limit':0,'count':0,'results':[],'facets':[]})['records']==0
    with pytest.raises(ValueError):validate_year_facets(dict(data,limit=5))
    with pytest.raises(ValueError):validate_year_facets(dict(data,count=5))


def test_field_probe_only_exports_counts_not_rows():
    x=field_presence({'results':[
        {'eventID':'a','eventDate':'1980-05-13','locationID':None,'samplingProtocol':'Pollard',
         'occurrenceStatus':'PRESENT','scientificName':'Example species'},
        {'eventID':'b','eventDate':'1980-06-17','locationID':'siteB'}]})
    assert x['returned_sample_rows']==2
    assert x['number_with_field_nonblank']['eventID']==2
    assert x['number_with_field_nonblank']['locationID']==1
    assert '1980-05-13' not in str(x) and 'Pollard' not in str(x)
    assert x['recording_source_samples_are_not_a_site_visit_design'] is True
