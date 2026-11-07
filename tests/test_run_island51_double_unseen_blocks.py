import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_island51_double_unseen_blocks import (
  partition_geographic_blocks,holdout_plan,GEOGRAPHY,BLOCK_SIZES,CATEGORIES,
  ROOT,PAIRS,ROSTER)
from exploratory_island51_circuit_network_holdout import source

def test_group_rule_fixed_and_all_islands():
    q,ids,_=source(PAIRS,ROSTER)
    blocks,geo=partition_geographic_blocks(GEOGRAPHY,ids)
    assert len(blocks)==51 and len(geo)==51
    assert [int((geo.longitude_block==k).sum()) for k in range(5)]==list(BLOCK_SIZES)
    for k in range(5):
        longitudes=geo.loc[geo.longitude_block.eq(k),"representative_longitude"].to_numpy()
        assert np.all(np.diff(longitudes)>=0)

def test_all_pairs_double_unseen_no_train_endpoint():
    q,ids,_=source(PAIRS,ROSTER)
    block,_=partition_geographic_blocks(GEOGRAPHY,ids)
    folds,desc=holdout_plan(q,ids,block)
    assert len(folds)==15 and int(desc.test_pairs.sum())==1275
    assert int(desc.loc[desc.pair_class.eq("within_longitude_block"),"test_pairs"].sum())==235
    assert int(desc.loc[desc.pair_class.eq("between_longitude_blocks"),"test_pairs"].sum())==1040
    tested=np.zeros(len(q),dtype=int)
    for fid,train,test in folds:
        held=set()
        for pair in q[test].itertuples(index=False):
            held.add(pair.island_a);held.add(pair.island_b)
        assert all((row.island_a not in held and row.island_b not in held) for row in q[train].itertuples(index=False))
        tested+=test.astype(int)
    assert np.array_equal(tested,np.ones(len(q),dtype=int))
