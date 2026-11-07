import pytest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from audit_island51_portwatch_five_WPI_source_pilot import (
   source,normalize,distance,IDS,EXPECTED,ISO3,PORTS,SOURCE_BLOB,
   REGISTRY,MAX_ROWS_PER_BOX,HALF_BOX_DEGREES,MAX_CENTROID_KM
)

def test_source_only_original_enriched_pilot_unchanged():
    assert IDS==("28","428","817","1416","1525")
    assert len(EXPECTED)==5 and set(EXPECTED)==set(IDS)
    assert ISO3=={"GB":"GBR","ES":"ESP","CY":"CYP","LK":"LKA","US":"USA"}
    assert MAX_ROWS_PER_BOX==15 and HALF_BOX_DEGREES==.35
    assert MAX_CENTROID_KM==25
    if not PORTS.is_file():
        pytest.skip("Frozen original WPI file is intentionally absent from the anonymous-review bundle")
    data=source()
    assert len(data)==194 and len({x["island_id"] for x in data.values()})==51
    for pid,expected in EXPECTED.items():
        x=data[pid]
        assert (x["port_name"],x["country_code"],x["island_id"])==expected

def test_exact_frozen_canonical_name_prevents_fuzzy_inflation():
    assert normalize("Santa Cruz de Tenerife")=="SANTA CRUZ DE TENERIFE"
    assert normalize("SANTA CRUZ DE TENERIFE")=="SANTA CRUZ DE TENERIFE"
    assert normalize("Belfast")=="BELFAST"
    assert normalize("Hull")!="KINGSTON UPON HULL"
    assert normalize("Harwich")!="FELIXSTOWE"

def test_bounded_terrain_geography_and_portwatch_registry_only():
    assert REGISTRY.endswith("PortWatch_ports_database/FeatureServer/0")
    assert abs(distance(54.6,-5.9,54.6,-5.9))<1e-10
    assert 10<distance(26,-80,26.1,-80)<15
    assert abs(distance(26,-80,26.1,-80)-distance(26.1,-80,26,-80))<1e-10
