"""Official OSGB 100km/10km conversion and no-source-zero controls."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_bsbi_original_host_hectad_opportunity_20261029 import to_osgb_hectad,recorded_flag

def test_london_bng_10km_hectad():
    assert to_osgb_hectad(530000,180000)=="TQ38"
    assert to_osgb_hectad(410000,220000)=="SP12"

def test_bsbi_source_blank_is_unrecorded_not_absence():
    assert recorded_flag("") is False
    assert recorded_flag("0") is False
    assert recorded_flag("1") is True
    assert recorded_flag("2") is True
    with pytest.raises(ValueError):recorded_flag("unknown code")
