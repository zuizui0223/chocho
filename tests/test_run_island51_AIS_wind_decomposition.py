import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_island51_AIS_wind_decomposition import (
    ROOT,MODELS,COMPARISONS,make_features,partition_geographic_blocks,holdout_plan,
    GEOGRAPHY,ROSTER,PAIRS)
from exploratory_island51_circuit_network_holdout import source
from run_island51_host_circuit_joint import INPUT_PATHS,load_original_port_counts,load_resistance_matrix

def test_fixed_marine_2_by_2():
    assert MODELS==("HRE","HRE_AIS","HRE_WIND","HCRE")
    assert [x[0] for x in COMPARISONS]==[
        "HRE_to_HRE_AIS","HRE_to_HRE_WIND","HRE_AIS_to_HCRE",
        "HRE_WIND_to_HCRE","HRE_to_HCRE"]

def test_composition_and_historical_effort_held_fixed_in_all_four_models():
    q,ids,_=source(PAIRS,ROSTER)
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    ports=load_original_port_counts(ROSTER,ids,prep)
    mats={key:load_resistance_matrix(ROOT/INPUT_PATHS[key],ids)
          for key in ("uniform","ais","wind")}
    designs,names,history=make_features(q,ids,ports,mats)
    for target in designs:
        assert set(names[target]["HRE"])==set(names[target]["HRE_AIS"])-{"log1p_normalized_ais_circuit"}
        assert set(names[target]["HRE"])==set(names[target]["HRE_WIND"])-{"log1p_normalized_wind_circuit"}
        assert set(names[target]["HRE"])==set(names[target]["HCRE"])-{
            "log1p_normalized_ais_circuit","log1p_normalized_wind_circuit"}
        assert all(len(x)==1275 and np.isfinite(x).all() for x in designs[target].values())
    block,_=partition_geographic_blocks(GEOGRAPHY,ids)
    folds,manifest=holdout_plan(q,ids,block)
    assert len(folds)==15 and int(manifest.test_pairs.sum())==1275
