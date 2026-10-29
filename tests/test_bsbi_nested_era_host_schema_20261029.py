"""Inner atlas path/CRC and frozen provenance guard without source observations."""
import io,zipfile,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_bsbi_nested_era_host_schema_20261029 import nested_schema

def test_reject_non_source_even_with_expected_nested_filename():
    inner=io.BytesIO()
    with zipfile.ZipFile(inner,'w') as z:z.writestr('hectad.csv','hectad,ddbid,era\nSP10,123,1970-1986\n')
    outer=io.BytesIO()
    with zipfile.ZipFile(outer,'w') as z:z.writestr('hectadDistributions.zip',inner.getvalue())
    with pytest.raises(ValueError,match='source drift'):
        nested_schema(outer.getvalue())
