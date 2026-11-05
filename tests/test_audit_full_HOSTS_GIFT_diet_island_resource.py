"""Original WCVP-HOSTS diet-trait estimator tests, including zero island gain."""
import pytest
for pkg in ("numpy","pandas"):pytest.importorskip(pkg)
import sys,importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
base=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(base))
sp=importlib.util.spec_from_file_location("dietstudy",base/"audit_full_HOSTS_GIFT_diet_island_resource.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_WCVP_id_only_lossless_integer_normalization():
    assert m.id_key("12345.0")=="12345"
    assert m.id_key("12345")=="12345"
    assert m.id_key("KEW-12345")=="KEW-12345"
    assert m.id_key("123.500")=="123.500"

def test_zero_island_gain_butterfly_kept_in_independent_diet_model(tmp_path):
    p=pd.DataFrame([
        {"species":f"G{i} s","native_resource_units":20+i,
         "host_family_count":2,"log_resource_expansion":.1*i,
         "host_taxonomy_lower_bound_adequate":True} for i in range(239)
    ])
    desc=pd.DataFrame([{"species":f"G{i} s","disturbance_affinity":
         ("Disturbance-associated (strong)" if i%2 else "Disturbance-avoidant (strong)")}
         for i in range(239)])
    hosts=[]
    for i in range(239):
        for j in range(2):
            hosts.append({"insect_species":f"G{i} s","accepted_plant_name_id":str(10000+i*2+j),
                          "accepted_name":f"Xy{i} planta{j}","family":"Fabaceae"})
    c=pd.DataFrame([{"plant_name_id":str(10000+i*2+j),
                     "lifecycle":"annual" if (i+j)%2 else "perennial",
                     "growth_form":"herb","taxon_key":f"x{i}{j} p"}
         for i in range(239) for j in range(2)])
    edges=pd.DataFrame([{"species":f"G{i} s","island_code":"ALD","plant_id":"10000"}
         for i in range(100)])
    # Source input must reject when the Q2 frozen source is different.
    for name,z in [("h",pd.DataFrame(hosts)),("c",c),("p",p),("d",desc),("e",edges)]:
        z.to_csv(tmp_path/f"{name}.csv",index=False)
    with pytest.raises(ValueError,match="Original corrected Q2 island add source"):
        m.load(*(tmp_path/f"{x}.csv" for x in ("h","c","p","d","e")))

def test_noncausal_species_regression_on_all_host_diet():
    rows=[]
    for i in range(36):
        ann=(i%7)/6
        rows.append({"genus":f"G{i//2}","log_island_added":.3*ann+.07*i,
                     "annual_herb_fraction":ann,
                     "disturbance_binary":i%2,"log_diet_hosts":np.log1p(3+i%3),
                     "log_native":np.log1p(10+i),"log_host_families":np.log1p(2+i%4),
                     "coverage":.7+.01*(i%3)})
    out=m.genus_bootstrap(pd.DataFrame(rows),"log_island_added","annual_herb_fraction",
      ["disturbance_binary","log_diet_hosts","log_native","log_host_families","coverage"],B=59,seed=5)
    assert out["status"]=="EXPLORATORY_NONCAUSAL_GENUS_BOOTSTRAP"
    assert len(out["beta_boot95"])==2 and out["successes"]==59
