import io,zipfile,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_island_gshhg_original_source_20261029 import inspect

def test_public_original_sha_missing_fails():
    with pytest.raises(ValueError,match="SHA256"):
        inspect(b"not a real archive")

def test_untouched_source_gates_no_surrogate_polygon():
    with pytest.raises(ValueError,match="SHA256"):
        inspect(bytes(10))
