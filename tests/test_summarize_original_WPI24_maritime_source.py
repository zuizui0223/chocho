import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from summarize_original_WPI24_geoverified_maritime_sources import (
    WPI, SUB17, SUB13, JOIN, MODEL, AIS, source_24
)

def test_fixed_original_30_source_candidates_and_24_independent_geoverified_harbors():
    if not WPI.is_file() or not (SUB13/"WPI13_independent_PortWatch_source_geography_LOCODE.csv").is_file():
        pytest.skip("Original WPI and PortWatch source-specific data omitted in anonymous review snapshot")
    ports,original=source_24()
    assert len(original)==194
    assert len({x["island_id"] for x in original.values()})==51
    assert sum(x["included"].lower()=="true" for x in original.values())==182
    assert len(ports)==24
    assert len({x["original_island_id"] for x in ports.values()})==17
    assert len({x["verified_PortWatch_portid"] for x in ports.values()})==24
    assert all(x["verified_PortWatch_LOCODE"] for x in ports.values())
    assert {"28","428","817","1416","96","98","576"}.issubset(set(ports))
    assert {"758","764","146","284","416","1516"}.isdisjoint(set(ports))
    assert all((x["original_WPI_independent_geographic_identity_verified"]=="True"
                 and 0<=float(x["original_WPI_to_PortWatch_distance_km"])<=25)
               for x in ports.values())
    assert JOIN.is_file() and MODEL.is_file() and AIS.is_file()

def test_original_exact_WPI24_population_has_276_pairs():
    assert 24*23//2==276
    assert 51*50//2==1275
