import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from exploratory_island51_circuit_network_holdout import source
from run_island51_cross_community_marine import (
    ROOT,PAIRS,ROSTER,SOURCE,MODEL_NAMES,PANELS,COMPARISONS,SCOPES,
    validate_source,build_common,model_designs
)
from run_island51_host_circuit_joint import INPUT_PATHS,load_original_port_counts,load_resistance_matrix

def original():
    q,ids,_=source(PAIRS,ROSTER)
    prep=json.loads((ROOT/INPUT_PATHS["preparation"]).read_text())
    counts=load_original_port_counts(ROSTER,ids,prep)
    matrices={name:load_resistance_matrix(ROOT/INPUT_PATHS[name],ids)
        for name in ("uniform","ais","wind")}
    return q,ids,counts,matrices

def test_all_historical_source_records_and_target_independent_cohort():
    q,ids,ports,matrices=original()
    plants,coverage,island,receipt,cov,jac,good,high=validate_source(q,ids)
    assert len(q)==1275 and len(ids)==51 and sum(ports.values())==194
    assert sum(len(plants[k])>0 for k in ids)==45
    assert int(good.sum())==990 and int(high.sum())==741
    assert int((~good).sum())==285
    assert len(coverage)==194 and len(island)==51 and receipt["total_port_speciesKey_rows"]==27598
    assert np.all(jac[good]>=0) and np.all(jac[good]<=1)
    assert len(MODEL_NAMES)==4 and len(PANELS)==4 and len(COMPARISONS)==5 and len(SCOPES)==4

def test_plant_response_not_predicted_with_own_community_similarity_or_future_hosts():
    q,ids,ports,matrices=original()
    plants,coverage,island,receipt,cov,jac,good,high=validate_source(q,ids)
    designs,names,targets,meta=model_designs(q,ids,ports,matrices,island,cov,jac)
    assert set(designs)=={"early_general_flora","earlier_butterflies",
           "later_butterflies","later_butterflies_retro_hosts"}
    assert np.array_equal(targets["early_general_flora"],jac)
    assert np.array_equal(targets["earlier_butterflies"],q.earlier_butterfly_jaccard)
    assert np.array_equal(targets["later_butterflies"],q.new_butterfly_jaccard)
    assert np.array_equal(targets["later_butterflies_retro_hosts"],q.new_butterfly_jaccard)
    for panel in designs:
        for key in MODEL_NAMES:
            assert designs[panel][key].shape[0]==1275
            assert np.isfinite(designs[panel][key]).all()
        assert all("early_general_flora_jaccard" not in str(x).lower()
                   for x in names[panel]["BASE"])
    assert "previous_2001_2010_butterfly_jaccard" not in names["early_general_flora"]["BASE"]
    assert "previous_2001_2010_butterfly_jaccard" not in names["earlier_butterflies"]["BASE"]
    assert "previous_2001_2010_butterfly_jaccard" in names["later_butterflies"]["BASE"]
    assert "retrospective_host_jaccard_pooled_2001_2024" not in names["later_butterflies"]["BASE"]
    assert "retrospective_host_jaccard_pooled_2001_2024" in names["later_butterflies_retro_hosts"]["BASE"]
    assert "native_AIS_circuit_resistance" not in names["early_general_flora"]["BASE"]
    assert "native_AIS_circuit_resistance" in names["early_general_flora"]["AIS"]
    assert "mean_windspeed_circuit_resistance" not in names["early_general_flora"]["AIS"]
    assert "mean_windspeed_circuit_resistance" in names["early_general_flora"]["BOTH"]
