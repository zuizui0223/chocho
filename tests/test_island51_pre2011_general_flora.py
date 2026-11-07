import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from acquire_island51_pre2011_public_flora import original_source_rows,query_params,QUERY_CONST
from run_island51_pre2011_general_flora_prediction import (
    ROOT,ROSTER,PAIRS,MODELS,CONTRASTS,SCOPES,general_flora_features
)
from exploratory_island51_circuit_network_holdout import source

def test_all_194_original_source_port_coordinates_and_early_archive_locked():
    port,older=original_source_rows()
    assert len(port)==194 and len(older)==194
    assert len({r["island_id"] for r in port.values()})==51
    assert sum(r["marine_eligible"] for r in port.values())==182
    assert sum(older.values())==155002
    assert sum(v>0 for v in older.values())==141
    assert sum(v==0 for v in older.values())==53
    for pid,source_port in port.items():
        params=query_params(source_port)
        assert params["geoDistance"].endswith(",5km")
        assert params["year"]=="2001,2010"
        assert params["taxonKey"]=="6"
        assert params["hasGeospatialIssue"]=="false"
        assert params["occurrenceStatus"]=="PRESENT"
        assert params["facet"]=="speciesKey"
        assert params["facetLimit"]=="1000"
        assert params["limit"]=="0"

def test_exact_51_original_islands_and_1275_source_pairs():
    q,ids,_=source(PAIRS,ROSTER)
    assert len(ids)==51 and len(q)==1275
    assert len(MODELS)==6 and len(CONTRASTS)==6 and len(SCOPES)==4

def test_new_early_general_flora_jaccard_missingness_controls():
    q,ids,_=source(PAIRS,ROSTER)
    plants={k:set() for k in ids}
    plants[ids[0]]={3,5}
    plants[ids[1]]={5,10}
    controls,jaccard,masks=general_flora_features(q,ids,plants)
    assert all(np.isfinite(x).all() and x.shape==(1275,) for x in [jaccard,*controls.values()])
    assert sum(len(s)>0 for s in plants.values())==2
    assert int(masks["both_early_species_indexed"].sum())==1
    loc=q.index[((q.island_a==ids[0])&(q.island_b==ids[1])) |
                ((q.island_a==ids[1])&(q.island_b==ids[0]))]
    assert len(loc)==1
    assert np.isclose(jaccard[int(loc[0])],1/3)
    assert controls["both_endpoints_no_2001_2010_plant_specieskey"].sum()>0
    assert controls["only_one_endpoint_no_2001_2010_plant_specieskey"].sum()>0
