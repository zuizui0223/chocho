"""Source-free checks for independent GBIF-recording outcome comparison."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_gbif_31_observed_host_beta.py"


def load():
    for requirement in ("pyproj","shapefile"):
        pytest.importorskip(requirement)
    sys.path.insert(0,str(ROOT/"scripts"))
    spec=importlib.util.spec_from_file_location("gbif31_independent_beta",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_pairwise_real_record_assemblage_Jaccard_and_defined_zero_policy():
    m=load()
    r={"A":{"Vanessa cardui","Danaus plexippus"},"B":{"Danaus plexippus"},
       "C":{"Aglais urticae"}}
    bits,names=m.jaccard_from_sets(r,["A","B","C"])
    assert len(names)==3
    values,zeros=m.ji(bits,np.array([0,0,1]),np.array([1,2,2]))
    assert zeros==0
    assert values.tolist()==pytest.approx([.5,0,0])
    empty, names=m.jaccard_from_sets({"A":set(),"B":set()},["A","B"])
    z,undefined=m.ji(empty,np.array([0]),np.array([1]))
    assert z[0]==0 and undefined==1


def test_heldout_region_scoring_preserves_all_pair_nodes():
    m=load()
    n=6
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    ia=np.array([i for i,j in pairs])
    ib=np.array([j for i,j in pairs])
    host=np.array([((i+j)%3)/2 for i,j in pairs])
    observed=.15+.6*host
    control=np.ones(len(pairs))
    models={
      "M0_general_flora_and_effort":[control],
      "M1_known_host":[control,host],
      "M2_other_aliens":[control,host*0],
      "M3_both":[control,host,host*0],
    }
    scoring,comparisons=m.evaluate_loio(
        observed,models,ia,ib,node_count=n,lam=0.01)
    assert set(scoring)==set(models)
    assert scoring["M1_known_host"]["island_or_region_equal_mse"] < scoring["M0_general_flora_and_effort"]["island_or_region_equal_mse"]
    assert comparisons["M1_known_host"]["regions_improved"]==n
    assert comparisons["M2_other_aliens"]["mean_region_equal_mse_difference"]==pytest.approx(0,abs=1e-12)
    assert comparisons["M1_known_host"]["mse_reduction_percent"]>0


def test_whole_world_division_holdout_excludes_all_held_region_pairs():
    m=load()
    n=9
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    ia=np.array([i for i,j in pairs])
    ib=np.array([j for i,j in pairs])
    y=np.array([.2+(((i+j)%3)/10) for i,j in pairs])
    host=np.array([float((i*j)%7)/7 for i,j in pairs])
    worlds=np.array(["1"]*3+["3"]*3+["7"]*3)
    features={
      "M0_general_flora_and_effort":[np.ones(len(ia))],
      "M1_known_host":[np.ones(len(ia)),host],
      "M2_other_aliens":[np.ones(len(ia)),1-host],
      "M3_both":[np.ones(len(ia)),host,1-host],
    }
    result=m.evaluate_division_heldout(y,features,ia,ib,worlds,lam=10.)
    assert set(result["models"])==set(features)
    folds=result["models"]["M0_general_flora_and_effort"]["folds"]
    assert len(folds)==3
    assert all(x["heldout_region_count"]==3 for x in folds)
    assert all(x["training_region_pair_count"]==15 for x in folds)
    assert all(x["test_region_pair_count"]==21 for x in folds)
    assert result["comparisons"]["M1_known_host"]["divisions_improved"]<=3
    assert "no outcome" in result["method"]
    with pytest.raises(ValueError,match="labels"):
        m.evaluate_division_heldout(y,features,ia,ib,worlds[:-1],10.)

