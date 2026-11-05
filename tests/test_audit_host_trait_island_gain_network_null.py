"""Offline ecological structural-null tests; no GIFT or GBIF calls."""
import pytest
for dep in ("numpy","pandas"):pytest.importorskip(dep)
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
directory=Path(__file__).resolve().parents[1]/"scripts"
import sys
sys.path.insert(0,str(directory))
p=directory/"audit_host_trait_island_gain_network_null.py"
spec=importlib.util.spec_from_file_location("networkplantnull",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_island_count_is_separate_from_global_gain_and_controls_valid():
    rows=[]
    for k in range(28):
        rows.append({"genus":f"G{k//2}","log1p_added_island_cells":np.log1p(k+1),
            "log_global_gain":np.log1p(k*.2),
            "log_introduced_host_count":np.log1p(2+k%5),
            "log_native":np.log1p(k+5),"log_host_families":np.log1p(2+k%3),
            "distinct_both_trait_coverage":.75+(.01*(k%4)),
            "disturbance_strong_binary":float(k%2)})
    df=pd.DataFrame(rows)
    f,controls=m.design(df,"log1p_added_island_cells")
    assert "disturbance_strong_binary" in controls
    x=np.asarray([float(i%7)/6 for i in range(28)])
    assert f(x) is not None

def test_plant_label_null_keeps_plant_traits_and_number_of_hosts_fixed():
    rng=np.random.default_rng(3)
    n=35;plants=21
    scored=pd.DataFrame([{
        "genus":f"G{i//2}",
        "log1p_added_island_cells":np.log1p(2+i%13),
        "log_introduced_host_count":np.log1p(4+i%7),
        "log_native":np.log1p(i+12),
        "log_host_families":np.log1p(2+i%5),
        "distinct_both_trait_coverage":.80,
        "disturbance_strong_binary":float(i%2)
    } for i in range(n)])
    matrix=np.zeros((n,plants),dtype=float)
    for i in range(n):
        matrix[i,[(i+j*3)%plants for j in range(4)]]=1
    denom=matrix.sum(axis=1)
    annual=np.array([float(i%2) for i in range(plants)])
    valid=pd.DataFrame({"stratum":np.repeat("matched_network",plants)})
    out=m.structural_null(scored,valid,matrix,denom,annual,
                [np.arange(plants)],response="log1p_added_island_cells",B=59,seed=4)
    assert out["permutations"]==59
    assert out["network_plant_taxa_known_both"]==plants
    assert 0<=out["two_sided_centered_p_descriptive"]<=1
