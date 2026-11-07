import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from run_island51_AIS_local750_node_null import (
    ROOT,PAIRS,ROSTER,GEOGRAPHY,SEED,NULL_REPS,MAX_SOURCE_SWAP_DISTANCE_KM,
    freeze_local_swap_pairs,generate_local_vertex_assignments,ais_source_matrix,
    source_haversine_km,COMPARE
)
from exploratory_island51_circuit_network_holdout import source
from run_island51_host_circuit_joint import (
    INPUT_PATHS,load_original_port_counts,load_resistance_matrix
)
from run_island51_latitude_controls import make_features

def test_original_51_island_750km_distance_capped_pairs():
    q,ids,_=source(PAIRS,ROSTER)
    groups,table,unchanged=freeze_local_swap_pairs(ids)
    assert len(ids)==51 and len(q)==1275
    assert MAX_SOURCE_SWAP_DISTANCE_KM==750.0
    assert len(groups)==18 and len(table)==18 and len(unchanged)==15
    assert set(table.hemisphere)=={"N","S"}
    assert table.centroid_greatcircle_distance_km.max()<=750.
    assert table.centroid_greatcircle_distance_km.max()>700.
    assert len({idx for pair in groups for idx in pair})==36
    for row in table.itertuples(index=False):
        assert source_haversine_km(row.latitude_a,row.longitude_a,row.latitude_b,row.longitude_b)<=750.
        assert row.latitude_a*row.latitude_b>0
    assert set(unchanged)==set(ids)-set(table.island_a)-set(table.island_b)

def test_exact_199_source_local_node_swap_maps():
    q,ids,_=source(PAIRS,ROSTER)
    group,table,unchanged=freeze_local_swap_pairs(ids)
    maps,attempts=generate_local_vertex_assignments(ids,group)
    assert SEED==20261107 and NULL_REPS==199
    assert len(maps)==200 and attempts>=199
    assert len({tuple(x.tolist()) for x in maps})==200
    assert np.array_equal(maps[0],np.arange(51))
    assert all(sorted(x.tolist())==list(range(51)) for x in maps)
    assert all(set(p[pair].tolist())==set(pair.tolist()) for p in maps for pair in group)
    positions={name:idx for idx,name in enumerate(ids)}
    assert all(int(x[positions[s]])==positions[s] for x in maps for s in unchanged)

def test_original_source_native_AIS_graph_is_unchanged():
    q,ids,_=source(PAIRS,ROSTER)
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    ports=load_original_port_counts(ROSTER,ids,prep)
    mats={key:load_resistance_matrix(ROOT/INPUT_PATHS[key],ids)
          for key in ("uniform","ais","wind")}
    feature,names,_=make_features(q,ids,ports,mats)
    k="log1p_normalized_ais_circuit"
    v=feature["new_butterfly_jaccard"]["HRE_AIS"][:,names["new_butterfly_jaccard"]["HRE_AIS"].index(k)]
    m=ais_source_matrix(q,ids,v)
    assert m.shape==(51,51)
    assert np.array_equal(m,m.T)
    assert np.array_equal(m[q.i.to_numpy(dtype=int),q.j.to_numpy(dtype=int)],v)
    assert len(COMPARE)==2
