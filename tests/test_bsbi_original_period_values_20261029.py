"""Only source-period data vocabulary, no host-absence inference."""
from pathlib import Path
import sys,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_original_period_values_20261029 import CLASSES,run

def test_all_periods_kept_distinct():
    assert len(CLASSES)==6
    assert "1970 - 1986" in CLASSES and "1987 - 1999" in CLASSES
    assert CLASSES.index("1970 - 1986") < CLASSES.index("1987 - 1999")

def test_source_hash_missing_stops(tmp_path):
    p=tmp_path/"no.zip"
    with pytest.raises(FileNotFoundError):run(p)
