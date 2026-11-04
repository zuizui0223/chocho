"""Source-constrained national N country-year parsing and FE model tests."""
import pytest
for dep in ("numpy","pandas","statsmodels","scipy","rasterio","pyproj","pycountry"):
    pytest.importorskip(dep)
import importlib.util
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
spec=importlib.util.spec_from_file_location("nies",ROOT/"scripts"/"exploratory_nies_country_port_interaction.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

import numpy as np
import pandas as pd
import pytest


def test_kgN_fraction_country_aggregation_without_synthetic_port_downscale():
    f=pd.DataFrame([
      {"ISO3":"AAA","Year":2005,"TNE":100,"UNM1":20,"UNM2":20,"UNM3":0,"UNM4":0,"UNM5":0,"PMN1":60},
      {"ISO3":"AAA","Year":2006,"TNE":200,"UNM1":20,"UNM2":20,"UNM3":0,"UNM4":0,"UNM5":0,"PMN1":160}
    ])
    f["unmanaged_N_kg"]=f[[f"UNM{i}" for i in range(1,6)]].sum(axis=1)
    a=m.country_exposure(f,(2005,2006))
    assert len(a)==1
    assert math.isclose(a.unmanaged_fraction.iloc[0],80/300)
    assert math.isclose(a.untreated_sewer_fraction.iloc[0],40/300)
    assert math.isclose(a.plant_managed_fraction.iloc[0],220/300)
    assert a.n_years.iloc[0]==2


def test_2015_country_N_interaction_is_not_country_main_effect():
    rng=np.random.default_rng(3)
    rows=[]
    for i in range(25):
        c=f"C{i:02d}"
        for j in range(8):
            x=(j-3.5)/2
            country_N=i/24
            y=.13*x +.27*x*(country_N-.5) + .03*rng.normal()
            rows.append({"country":c,"log_wwp":x,"wwp_country_center":x,
                         "interaction_unmanaged":x*(country_N-.5),
                         "log_wwtp":x*.3,
                         "abs_lat":float(i),"sin_lon":0.2,"cos_lon":.4,
                         "log_indexed":.7,"built_fraction_2021":.2,"crop_fraction_2021":.1,
                         "y_ratio":y})
    d=pd.DataFrame(rows)
    z=m.clustered_model(d,period="2015",urban=True,case="unmanaged")
    assert z["countries"]==25
    assert z["n"]==200
    assert z["interaction_unmanaged"]["beta"]>0
    assert "unmanaged_fraction" not in z["formula"]
    assert "C(country)" in z["formula"]
    assert z["status"].startswith("ESTIMATED")


def test_national_iso3_conversion_preserves_namibia_NA_code():
    assert m.iso3("NA")=="NAM"
    assert m.iso3("JP")=="JPN"
    assert m.iso3("XK")=="XKX"
    assert m.iso3("NO_SUCH") is None
