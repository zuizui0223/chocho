"""Offline test of region/year source-to-independent GBIF beta methods."""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analyze_gbif_global224_observed_beta.py"


def load():
    pytest.importorskip("pyproj")
    pytest.importorskip("shapefile")
    pytest.importorskip("pandas")
    sys.path.insert(0,str(ROOT/"scripts"))
    spec=importlib.util.spec_from_file_location("global224_observed",SCRIPT)
    v=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v)
    return v


def test_observed_species_bitsets_do_not_infer_absence():
    m=load()
    seen={"r1":{"A species","B species"},"r2":{"B species"},
          "r3":{"C species"}}
    masks,taxa=m.species_masks(seen,["r1","r2","r3"])
    assert taxa==["A species","B species","C species"]
    v,undef=m.jaccard(masks,np.asarray([0,0,1]),np.asarray([1,2,2]))
    assert undef==0
    assert v.tolist()==pytest.approx([.5,0.,0.])
    assert m.PERIODS["prior"].start==2001
    assert m.PERIODS["later"].stop==2025


def test_entire_region_holdout_gain_not_pairwise_random_split():
    m=load()
    n=36
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    ia=np.asarray([a for a,b in pairs])
    ib=np.asarray([b for a,b in pairs])
    signal=np.asarray([((a+b)%7)/7 for a,b in pairs])
    response=.1+.6*signal
    models={
      "M0_all_flora_effort_prior":[np.ones(len(ia)),np.asarray([(a+b)%5 for a,b in pairs],float)],
      "M1_known_host":[np.ones(len(ia)),np.asarray([(a+b)%5 for a,b in pairs],float),signal],
      "M2_other_aliens":[np.ones(len(ia)),np.asarray([(a+b)%5 for a,b in pairs],float),np.zeros(len(ia))],
      "M3_both":[np.ones(len(ia)),np.asarray([(a+b)%5 for a,b in pairs],float),signal,np.zeros(len(ia))],
    }
    ids=[f"X{i:02d}" for i in range(n)]
    labels=["1"]*12+["2"]*12+["3"]*12
    folds=m.score_held_out(response,models,ia,ib,ids,labels,lam=.1)
    summary,contrasts=m.consolidate_scores(folds)
    assert all(len(v)==n for v in folds.values())
    assert all(r["training_pairs"]==595 and r["test_pairs"]==35 for r in folds["M1_known_host"])
    assert summary["M1_known_host"]["region_equal_MSE"]<summary["M0_all_flora_effort_prior"]["region_equal_MSE"]
    assert contrasts["M1_known_host"]["regions_improved"]==n
    world=m.world_holdout(response,models,ia,ib,ids,labels,lam=.1)
    assert world["world_divisions_evaluable"]==3
    assert all(v["training_pairs"]==276 for v in world["by_world_division"].values())
    assert all(v["test_pairs"]==354 for v in world["by_world_division"].values())


def test_source_missing_not_confused_with_indexed_absence(tmp_path):
    m=load()
    with pytest.raises(FileNotFoundError):
        m.load_annual_data(audit=tmp_path/"never_there.json")


def test_region_effort_gates_are_predeclared_and_primary_unchanged():
    m=load()
    assert m.REGION_GATES=={
        "primary":(5,10,20),
        "relaxed":(5,5,10),
        "stringent":(10,20,50),
    }
    with pytest.raises(ValueError,match="Unknown predeclared"):
        m.analyse(eligibility="best-p-value-after-scan")
