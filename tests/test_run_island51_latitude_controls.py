import pytest
np=pytest.importorskip("numpy")
pd=pytest.importorskip("pandas")
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_island51_latitude_controls import (
    MODELS,BASE_MODELS,COMPARISONS,latitude_pair_controls,PAIRS,ROSTER
)
from exploratory_island51_circuit_network_holdout import source

def test_latitude_control_names_and_original_51_island_pair_support():
    q,ids,_=source(PAIRS,ROSTER)
    f=latitude_pair_controls(q,ids)
    assert tuple(f)==(
        "sum_signed_latitude_degrees",
        "sum_absolute_latitude_degrees",
        "abs_latitude_separation_degrees",
        "opposite_hemisphere_indicator",
    )
    assert len(q)==1275 and len(ids)==51
    assert all(np.isfinite(v).all() and v.shape==(1275,) for v in f.values())
    assert set(np.unique(f["opposite_hemisphere_indicator"]))<={0.,1.}

def test_pair_order_does_not_change_geographic_controls():
    q,ids,_=source(PAIRS,ROSTER)
    fixed=latitude_pair_controls(q,ids)
    other=q.copy()
    other["i"]=q["j"].copy()
    other["j"]=q["i"].copy()
    swapped=latitude_pair_controls(other,ids)
    assert all(np.array_equal(fixed[k],swapped[k]) for k in fixed)

def test_original_and_latitude_adjusted_models_and_comparisons_frozen():
    assert BASE_MODELS==("HRE","HRE_AIS","HRE_WIND","HCRE")
    assert len(MODELS)==8 and len(COMPARISONS)==12
    assert [c[0] for c in COMPARISONS][:5]==[
       "HRE_to_HRE_AIS","HRE_to_HRE_WIND","HRE_AIS_to_HCRE",
       "HRE_WIND_to_HCRE","HRE_to_HCRE"]
    assert [c[0] for c in COMPARISONS][5:10]==[
       c[0]+"_LAT" for c in COMPARISONS[:5]]
    assert COMPARISONS[-2][0]=="HRE_to_HRE_LAT"
    assert COMPARISONS[-1][0]=="HCRE_to_HCRE_LAT"
