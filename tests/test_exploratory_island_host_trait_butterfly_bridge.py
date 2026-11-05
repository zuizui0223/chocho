"""Offline checks of exact introduced-host species and genus-level inference."""
import pytest
for d in ("numpy","pandas"): pytest.importorskip(d)
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
p=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_island_host_trait_butterfly_bridge.py"
s=importlib.util.spec_from_file_location("hostbridge",p)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def test_distinct_plant_taxa_vs_reused_island_edges():
    x=pd.DataFrame([
        {"species":"A a","island_code":"ALD","plant_id":"1","known_lifecycle":True,"known_both":True,"annual":True,"annual_herb":True},
        {"species":"A a","island_code":"CAN","plant_id":"1","known_lifecycle":True,"known_both":True,"annual":True,"annual_herb":True},
        {"species":"A a","island_code":"CAN","plant_id":"2","known_lifecycle":False,"known_both":False,"annual":False,"annual_herb":False},
        {"species":"A a","island_code":"CAN","plant_id":"3","known_lifecycle":True,"known_both":True,"annual":False,"annual_herb":False}
    ])
    panel=pd.DataFrame([{"species":"A a","native_resource_units":10,"host_family_count":2,
                         "log_resource_expansion":.3,"host_taxonomy_lower_bound_adequate":True}])
    desc=pd.DataFrame([{"species":"A a","disturbance_affinity":"Disturbance-associated (strong)"}])
    r=m.plant_level_species_summary(x,panel,desc).iloc[0]
    assert r.island_added_cells==2
    assert (r.distinct_total,r.distinct_known_both)==(3,2)
    assert abs(r.distinct_annual_herb_fraction-.5)<1e-10
    assert r.edge_total==4 and abs(r.edge_annual_herb_fraction-2/3)<1e-10
    assert abs(r.distinct_both_trait_coverage-2/3)<1e-10

def test_no_source_data_cannot_produce_positive_claim():
    rows=[{"genus":f"G{i}","response":float(i),"predict":float(i%2),"control":1.0} for i in range(12)]
    x=pd.DataFrame(rows)
    assert m.adjusted_beta(x,"response","predict",["control"]) is None
    assert m.genus_bootstrap(x,"response","predict",["control"],B=9)["status"]=="UNIDENTIFIABLE_FIXED_SOURCE"

def test_genus_cluster_bootstrap_returns_valid_interval():
    rows=[]
    for i in range(24):
        for rep in range(2):
            z=float((i+rep)%2)
            rows.append({"genus":f"G{i}","response":.23*z+.017*i+.002*rep,
                         "predict":z,"control":(i+rep)/10.0})
    d=pd.DataFrame(rows)
    x=m.genus_bootstrap(d,"response","predict",["control"],B=99,seed=3)
    assert x["status"]=="EXPLORATORY_NONCAUSAL_GENUS_BOOTSTRAP"
    assert x["successes"]==99 and len(x["beta_boot95"])==2 and x["beta"]>0
