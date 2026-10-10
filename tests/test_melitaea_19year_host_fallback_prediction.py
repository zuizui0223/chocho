"""Test leak-free chronological features and frozen model mechanics."""
import importlib.util
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
p=ROOT/"scripts/analyze_melitaea_19year_host_fallback_prediction.py"
spec=importlib.util.spec_from_file_location("predictor",p)
assert spec and spec.loader
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_predictor_uses_past_current_not_future_hosts():
    import audit_melitaea_19year_fallback_support as src
    raw=("year\tpatch\tpopulation\tplantago\tveronica\tprevious_population\n"
         "2005\ta\t3\t2\t3\t3\n"
         "2006\ta\t2\t2\t1\t3\n"
         "2007\ta\t0\t0\t0\t2\n").encode()
    rows,_=src.ingest(raw)
    net={"a":{"area":1000.0,"x":120000.0,"y":6700000.0}}
    panel=mod.build_panel(rows,net)
    assert len(panel)==1
    r=panel[0]
    assert r["y"]==1 and r["decline"]==1 and r["backup_high"]==1
    assert r["interaction"]==1
    assert r["base"][1:5]==[2.0,1.0,2.0,3.0]
    assert len(r["base"])==len(mod.PRIMARY_FIELDS)
    # The 2007 HOST columns are not allowed into features.
    rows[("a",2007)]["plantago"]=3.0
    rows[("a",2007)]["veronica"]=3.0
    assert mod.build_panel(rows,net)[0]["base"]==r["base"]

def test_static_area_coordinates_join_schema():
    tsv=b"patch\tx\ty\tarea\na\t1000\t2000\t500\n"
    d=mod.read_static(tsv)
    assert d["a"]=={"x":1000.0,"y":2000.0,"area":500.0}

def test_group_bootstrap_deterministic():
    dif=np.array([-0.1,0.2,0.3,-0.05])
    groups=["x","x","y","z"]
    a=mod.cluster_interval(dif,groups,20261010,draws_not_defined) if False else None
    x=mod.cluster_interval(dif,groups,seed=20261010,iterations=50)
    y=mod.cluster_interval(dif,groups,seed=20261010,iterations=50)
    assert x==y and x[0]<=x[1]
