"""Offline arithmetic checks for global report-prevalence control."""
from __future__ import annotations
import importlib.util
import math
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_nutrient_global_prevalence.py"
spec=importlib.util.spec_from_file_location("nutrient_global_prevalence_test",SCRIPT)
assert spec and spec.loader
import sys
sys.path.insert(0,str(SCRIPT.parent))
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_residualization_removes_linear_global_reporting_trend():
    y=[4+2*v for v in range(15)]
    x=[v for v in range(15)]
    assert max(abs(a) for a in mod.residualize(y,[x]))<1e-10


def test_partial_spearman_no_control_matches_raw_spearman():
    data=[{"ports":i+2,"N":i%7+1,
           "number_countries_with_indexed_observations":i//3+1,
           "global_indexed_records_2001_2024":i*i+12}
          for i in range(1,21)]
    assert abs(mod.partial_spearman(data,[])-
               mod.associations(data)["raw"])<1e-12
    r=mod.associations(data)
    assert all(-1<=v<=1 for v in r.values())
    assert math.isfinite(r["partial_both"])


def test_genus_bootstrap_is_reproducible():
    sample=[
        {"genus":str(i//4),"ports":1+(i*7)%23,"N":1+i%8,
         "number_countries_with_indexed_observations":1+(i*3)%9,
         "global_indexed_records_2001_2024":2+i*i}
        for i in range(28)
    ]
    a=mod.bootstrap_by_genus(sample,n_draws=49,seed=54)
    b=mod.bootstrap_by_genus(sample,n_draws=49,seed=54)
    assert a==b
    assert all(v["draws"]==49 for v in a.values())
