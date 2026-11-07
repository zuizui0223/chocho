import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from run_island51_AIS_node_null import (
    ROOT,PAIRS,ROSTER,GEOGRAPHY,SEED,NULL_REPS,freeze_latitude_bands,
    generate_source_node_assignments,ais_source_matrix,COMPARE
)
from exploratory_island51_circuit_network_holdout import source
from run_island51_host_circuit_joint import (
    INPUT_PATHS,load_original_port_counts,load_resistance_matrix
)
from run_island51_latitude_controls import make_features

def test_all_frozen_source_latitudes_and_hemisphere_bands():
    q,ids,_=source(PAIRS,ROSTER)
    groups,table,latitude=freeze_latitude_bands(ids)
    assert len(ids)==51 and len(q)==1275 and len(groups)==10
    assert [len(x) for x in groups]==[5,5,6,5,5,5,5,5,5,5]
    assert (table.hemisphere.eq("S")).sum()==10
    assert (table.hemisphere.eq("N")).sum()==41
    assert len(set(table.island_id))==51
    assert set(latitude)==set(ids)

def test_all_199_surrogates_are_unique_symmetric_node_labelings():
    q,ids,_=source(PAIRS,ROSTER)
    bands,table,_=freeze_latitude_bands(ids)
    maps,attempts=generate_source_node_assignments(ids,bands)
    assert SEED==20261107 and NULL_REPS==199
    assert len(maps)==200 and np.array_equal(maps[0],np.arange(51))
    assert len({tuple(p) for p in maps})==200
    assert all(sorted(x.tolist())==list(range(51)) for x in maps)
    assert all(set(p[group])==set(group) for p in maps for group in bands)
    assert 199<=attempts<=100000

def test_original_source_AIS_circuit_reconstruction_symmetric():
    q,ids,_=source(PAIRS,ROSTER)
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    ports=load_original_port_counts(ROSTER,ids,prep)
    matrices={key:load_resistance_matrix(ROOT/INPUT_PATHS[key],ids)
              for key in ("uniform","ais","wind")}
    features,names,_=make_features(q,ids,ports,matrices)
    initial=features["new_butterfly_jaccard"]["HRE_AIS"][
        :,names["new_butterfly_jaccard"]["HRE_AIS"].index("log1p_normalized_ais_circuit")]
    mat=ais_source_matrix(q,ids,initial)
    i=q.i.to_numpy(dtype=int);j=q.j.to_numpy(dtype=int)
    assert mat.shape==(51,51)
    assert np.allclose(mat,mat.T,rtol=0,atol=0)
    assert np.array_equal(mat[i,j],initial)
    assert len(COMPARE)==2
