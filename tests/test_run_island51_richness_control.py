import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_island51_richness_control import richness_features, host_counts, compare_folds, MODELS, COMPARISONS

def test_symmetric_pairwise_host_richness_terms():
    q = pd.DataFrame({"i":[0,0,1],"j":[1,2,2]})
    r = richness_features(q, ["1","2","3"], {"1":0,"2":4,"3":8})
    v = np.log1p([0,4,8])
    assert np.allclose(r["sum_log1p_candidate_host_taxa"],
                       [v[0]+v[1], v[0]+v[2], v[1]+v[2]])
    assert np.allclose(r["abs_difference_log1p_candidate_host_taxa"],
                       [v[1]-v[0],v[2]-v[0],v[2]-v[1]])

def test_richness_missing_island_rejected():
    with pytest.raises(ValueError, match="identity"):
        richness_features(pd.DataFrame({"i":[0],"j":[1]}),["1","2"],{"1":4,"3":5})

def test_richness_does_not_depend_on_butterfly_outcome():
    q = pd.DataFrame({"i":[0], "j":[1], "new_butterfly_species":[3]})
    a = richness_features(q, ["1","2"], {"1":3,"2":4})
    q["new_butterfly_species"] = 10000
    b = richness_features(q, ["1","2"], {"1":3,"2":4})
    assert all(np.array_equal(a[k],b[k]) for k in a)

def test_fixed_comparisons_and_model_support():
    assert len(MODELS) == 6
    assert len(COMPARISONS) == 7
    assert [c[0] for c in COMPARISONS if c[3].startswith("primary")] == ["CR_to_HCR"]
    with pytest.raises(ValueError, match="Incomplete"):
        compare_folds(pd.DataFrame({"target":[],"model":[],"heldout_island_id":[]}))
