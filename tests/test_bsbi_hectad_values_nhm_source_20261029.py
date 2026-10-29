"""Coordinate grid checks and archive authenticity, no access to original botanicals."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_bsbi_hectad_values_nhm_source_20261029 import grid10,nested_distributions_scan

def test_british_national_grid_hectads_known():
    assert grid10(529000,180000)=='TQ28'
    assert grid10(451000,205000)=='SP50'
    assert grid10(330000,390000)=='SJ39'
    assert grid10(0,0)=='SV00'
    assert grid10(-1,0) is None
    assert grid10(800000,180000) is None

def test_impostor_bsbi_zip_rejected_by_both_digests():
    with pytest.raises(ValueError,match='original MD5/SHA256'):
        nested_distributions_scan(b'fake')
