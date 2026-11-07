import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_island51_top_port_training_influence import (
    ROSTER,PAIRS,ROOT,SCENARIO_DROPS,MODEL_NAMES,COMPARE,
    verify_port_roster,scenario_plan
)
from run_island51_latitude_controls import partition_geographic_blocks,holdout_plan,GEOGRAPHY
from exploratory_island51_circuit_network_holdout import source
from run_island51_host_circuit_joint import load_original_port_counts,INPUT_PATHS
import json

def original_sources():
    q,ids,_=source(PAIRS,ROSTER)
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    counts=load_original_port_counts(ROSTER,ids,prep)
    return q,ids,counts

def test_original_top_two_source_port_groups_and_four_unselected_scenarios():
    q,ids,counts=original_sources()
    ranks=verify_port_roster(ids,counts)
    assert ranks.iloc[0].island_id=="15" and int(ranks.iloc[0].original_port_count)==39
    assert ranks.iloc[1].island_id=="9" and int(ranks.iloc[1].original_port_count)==21
    assert sum(counts.values())==194
    assert SCENARIO_DROPS=={"FULL":(),"DROP_15":("15",),
                            "DROP_9":("9",),"DROP_15_9":("15","9")}

def test_original_test_pairs_in_all_scenarios_and_extra_training_exclusions():
    q,ids,counts=original_sources()
    block,_=partition_geographic_blocks(GEOGRAPHY,ids)
    folds,spec=holdout_plan(q,ids,block)
    masks,ledger=scenario_plan(q,ids,folds,spec)
    assert len(masks)==60 and len(ledger)==60
    for scenario in SCENARIO_DROPS:
        subset=ledger[ledger.scenario.eq(scenario)]
        assert len(subset)==15 and int(subset.test_pairs.sum())==1275
        assert int(subset[subset.pair_class.eq("within_longitude_block")].test_pairs.sum())==235
        assert int(subset[subset.pair_class.eq("between_longitude_blocks")].test_pairs.sum())==1040
    ids_of_pairs=q[["island_a","island_b"]]
    for scenario,fid,train,test in masks:
        assert int((train&test).sum())==0
        blocked=set(SCENARIO_DROPS[scenario])
        held_test=set(ids_of_pairs.loc[test,"island_a"])|set(ids_of_pairs.loc[test,"island_b"])
        assert all((a not in blocked and b not in blocked and
                    a not in held_test and b not in held_test)
                   for a,b in ids_of_pairs.loc[train].itertuples(index=False,name=None))

def test_exact_latitude_adjusted_parent_models_and_five_source_contrasts():
    assert MODEL_NAMES==("HRE_LAT","HRE_AIS_LAT","HRE_WIND_LAT","HCRE_LAT")
    assert len(COMPARE)==5 and all(name.endswith("_LAT") for name,*_ in COMPARE)
