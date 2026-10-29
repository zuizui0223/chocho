"""Dryad data catalog and user data safety tests without private fields."""
import sys,json
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_dryad_europe_butterfly_host_tables_20261029 import csv_metadata,discover_file_records,FILES

def test_source_csv_schema_outputs_no_original_biology():
    x=csv_metadata(b'butterfly,host\nAglais urticae,Urtica dioica\n','larval_host_plant.csv')
    assert x['original_rows']==1
    assert x['header']==['butterfly','host']
    assert 'Aglais urticae' not in str(x)

def test_dryad_login_html_not_a_host_dataset():
    with pytest.raises(ValueError):csv_metadata(b'<html>Login required</html>','larval_host_plant.csv')

def test_missing_original_tables_not_taken_as_zero_hosts():
    with pytest.raises(ValueError):discover_file_records({'_embedded':{'stash:files':[]}})
