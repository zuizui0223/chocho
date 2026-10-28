"""No network source calendar and observation-leak test."""
import sys,datetime as dt,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from probe_uk_nasa_power_1981_2021_20261028 import inspect,SOURCE_CALENDAR
import pytest

def test_1981_is_source_start_not_original_1976():
    assert SOURCE_CALENDAR==(dt.date(2021,12,31)-dt.date(1981,1,1)).days+1
    assert SOURCE_CALENDAR>14000

def test_missing_original_variable_fails():
    with pytest.raises(ValueError,match="expected T2M"):
        inspect(json.dumps({"properties":{"parameter":{"T2M":{}}}}).encode())
