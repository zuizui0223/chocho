import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from plan_native_host_sdm_v01 import spatial_thin,stable_priority

def test_source_only_thin_keeps_points_farther_than_ten_km():
    rows=[{"gbif_id":str(i),"latitude":0.0,"longitude":float(i)*.2,
           "native_wgsrpd3_code":"AAA"} for i in range(12)]
    sel=spatial_thin(rows,"123",kilometers=10,cap=200)
    assert len(sel)==len(rows)

def test_repeated_coordinates_and_duplicates_do_not_bloat_native_presences():
    rows=[{"gbif_id":str(i),"latitude":45.0,"longitude":10.0,
           "native_wgsrpd3_code":"ABC"} for i in range(50)]
    assert len(spatial_thin(rows,"123",kilometers=10,cap=200))==1

def test_same_taxon_thinning_is_order_invariant():
    rows=[{"gbif_id":str(i),"latitude":float(i)/10,
           "longitude":float(i)/5,"native_wgsrpd3_code":"ABC"} for i in range(250)]
    q=spatial_thin(rows,"123",cap=200)
    r=spatial_thin(list(reversed(rows)),"123",cap=200)
    assert [x["gbif_id"] for x in q]==[x["gbif_id"] for x in r]
    assert len(q)==200

def test_different_taxa_use_independent_reproducible_sampling():
    row={"gbif_id":"abc123","latitude":1.0,"longitude":2.0}
    assert stable_priority(row,"123")==stable_priority(row,"123")
    assert stable_priority(row,"123")!=stable_priority(row,"124")

def test_invalid_coordinate_fails_closed():
    with pytest.raises(ValueError):
        spatial_thin([{"gbif_id":"abc","latitude":92.0,"longitude":0}],"123")
