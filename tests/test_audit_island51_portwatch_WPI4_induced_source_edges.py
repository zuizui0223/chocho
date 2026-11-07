import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_portwatch_WPI4_induced_source_edges import (
  fixed_samples,SOURCE_LAYER,PROVIDER_IDS,WPI_IDS,MAX_SAMPLE_ROWS,PROVIDER_SAMPLE
)

def test_four_original_WPI_ports_are_source_verified_before_any_edges():
    assert WPI_IDS==("28","428","817","1416")
    assert PROVIDER_IDS==("port141","port1153","port647","port254")
    assert MAX_SAMPLE_ROWS==16
    assert SOURCE_LAYER.endswith("Spillover_Simulator_Maritime_Connections/FeatureServer/0")
    if not PROVIDER_SAMPLE.exists():
        pytest.skip("Original PortWatch source pilot omitted from anonymous review snapshot")
    rows=fixed_samples()
    assert len(rows)==4
    assert len({x["original_island_group"] for x in rows})==4
    assert [x["portwatch_id"] for x in rows]==list(PROVIDER_IDS)
    assert all(x["source_registry_geodesic_km"]<=25 for x in rows)
