"""Country-cluster estimator and source support tests; no internet."""
import importlib.util, sys, numpy as np, pandas as pd
from pathlib import Path
S=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(S))
sp=importlib.util.spec_from_file_location("n_audit",S/"exploratory_cirad_n_robustness.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_ratio_formula_and_country_cluster():
    rng=np.random.default_rng(22)
    rows=[]
    for i in range(25):
        for j in range(10):
            crop=(j*.1)+(i*.03)
            human=j*.08+i*.02+((j*j)%5)*.02
            wwp=j*.09+i*.03+((j*11)%7)*.02
            rows.append({"country":f"C{i:02d}","port_id":i*10+j,
                "log_crop_N":crop,"log_human_N":human,
                "log_wwp":wwp,"log_wwtp":.4*wwp+.03*j*j,
                "abs_lat":i/10,"sinlon":j/10,"coslon":i*j/100,
                "log_indexed":j*.1+i*.02,
                "built_fraction_2021":((i+j)%5)/10,
                "crop_fraction_2021":((i+3*j)%7)/10,
                "ratio":1.1*crop+.08*rng.normal(),
                "alien":.2*crop+.05*rng.normal(),
                "native":.1*crop+.05*rng.normal()})
    d=pd.DataFrame(rows)
    z=m.estimate(d)
    assert z["n"]==250 and z["countries"]==25
    assert z["log_crop_N"]["beta"]>0
    assert len(z["log_crop_N"]["CI95"])==2
    assert m.estimate(d,with_N=False).get("log_crop_N") is None
    assert m.estimate(d,with_total=False).get("log_crop_N") is not None

def test_identifiability_support_screen():
    rows=[]
    for i in range(10):
        for j in range(3 if i<4 else 6):
            rows.append({"country":str(i),"port_id":i*10+j,
                         "log_crop_N":float(j),"log_human_N":float(j)})
    good=m.allowed(pd.DataFrame(rows))
    assert len(good)==36 and good.country.nunique()==6
