"""Offline tests for candidate-host response definition and country cluster models."""
from __future__ import annotations
import pytest
for dep in ("numpy","pandas","statsmodels","scipy","rasterio","pyproj"):
    pytest.importorskip(dep)
import importlib.util,math,sys
import numpy as np
import pandas as pd
from pathlib import Path

S=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(S))
spec=importlib.util.spec_from_file_location("estimand",S/"exploratory_cirad_n_estimand_stress.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_no_zero_as_absence_and_pseudo_count_sensitivity():
    d=pd.DataFrame({"country":["AA","BB","CC"],"port_id":[1,2,3],
                    "alien_candidates":[0,2,5],"native_candidates":[0,0,10]})
    q=m.add_responses(d)
    assert np.isnan(q.loc[0,"alien_share_observed"])
    assert np.isclose(q.loc[1,"alien_share_observed"],1)
    assert np.isclose(q.loc[2,"alien_share_observed"],1/3)
    assert np.isclose(q.loc[0,"ratio_pc_0p5"],0)
    assert q.loc[1,"ratio_pc_0p1"]>q.loc[1,"ratio_pc_2"]
    assert not d.columns.isin(["alien_share_observed"]).any()


def synthetic(ncountries=23,n_per=10):
    rng=np.random.default_rng(1299)
    rows=[]
    for k in range(ncountries):
        for i in range(n_per):
            x=(i-4.5)*.16+ rng.normal(scale=.13)
            h=rng.normal(size=1)[0]*.5 + .2*i
            ww=rng.normal(size=1)[0]*.4 + i*.1
            built=rng.uniform()
            crop=rng.uniform()
            frac=1/(1+np.exp(-(-.4+.45*x+.1*h+.04*ww+.1*built)))
            tot=rng.integers(4,35)
            a=rng.binomial(tot,frac)
            rows.append({"port_id":k*100+i,"country":f"C{k:02d}",
                "alien_candidates":a,"native_candidates":tot-a,
                "log_crop_N":x+.07*k,"log_human_N":h+.11*k,
                "log_wwp":ww,"log_wwtp":ww*.3+rng.normal(scale=.2),
                "abs_lat":k*.3,"sinlon":np.sin(k),"coslon":np.cos(k),
                "log_indexed":np.log1p(tot+9),
                "built_fraction_2021":built,"crop_fraction_2021":crop})
    return m.add_responses(pd.DataFrame(rows))


def test_country_FE_and_port_level_fractional_likelihood():
    z=m.country_eligibility(synthetic())
    assert (len(z),z.country.nunique())==(230,23)
    ols=m.fit_response(z,"ratio_pc_0p5",include_total_proxy=True)
    share=m.fit_response(z,"alien_share_observed",include_total_proxy=False,
                         model="FRAC_LOGIT")
    assert ols["n_ports"]==230 and share["n_ports"]==230
    assert share["port_weighting"]=="ONE_PORT_ONE_OBSERVATION"
    assert share["log_crop_N"]["beta"]>0
    assert len(share["log_crop_N"]["CI95_country_cluster"])==2
    assert "log_indexed" not in share["formula"]
    assert "log_indexed" in ols["formula"]


def test_low_country_counts_are_excluded_without_sign_selection():
    z=synthetic(ncountries=2,n_per=8)
    z=z.iloc[:].copy()
    z.loc[z.country=="C01","native_candidates"]=0
    q=m.country_eligibility(z)
    assert q.country.nunique()==1
    assert set(q.country)=={"C00"}
