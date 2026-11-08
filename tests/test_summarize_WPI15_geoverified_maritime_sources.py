import pytest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from summarize_original_WPI15_geoverified_maritime_sources import (
    SUB,JOIN,MODEL,AIS,WPI,source_15
)

def test_source_only_fifteen_independently_geoverified_original_ports():
    if not (WPI.is_file() and SUB.is_dir() and
            (SUB/"WPI17_independent_PortWatch_source_geography_LOCODE.csv").is_file()):
        pytest.skip("WPI source details omitted from anonymous review snapshot")
    m,original=source_15()
    assert len(original)==194
    assert len(m)==15
    assert len({x["original_island_id"] for x in m.values()})==14
    assert all(x["verified_PortWatch_LOCODE"] for x in m.values())
    assert "758" not in m and "764" not in m
    assert {"28","428","817","1416"}.issubset(set(m))
    assert JOIN.is_file() and MODEL.is_file() and AIS.is_file()
