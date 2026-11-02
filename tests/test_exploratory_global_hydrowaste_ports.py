"""Offline source unit and model checks for global sewage treatment proximity."""
import importlib.util
import math
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_global_hydrowaste_ports.py"
spec=importlib.util.spec_from_file_location("hydro_global",SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_github_frozen_port_inventory():
    ports=m.source_ports()
    assert len(ports)==1596
    assert ports.port_id.nunique()==1596
    assert ports.country.nunique()==142
    assert ports.alien_candidates.ge(0).all()
    assert ports.native_candidates.ge(0).all()
    assert ports.all_taxonomically_indexed.ge(ports.known_alien_native).all()


def test_country_resolution_and_spherical_distance():
    assert m.iso3("US")=="USA"
    assert m.iso3("GB")=="GBR"
    assert m.iso3("DE")=="DEU"
    assert m.iso3("JP")=="JPN"
    assert np.isclose(m.chord(25),2*np.sin(25/2/m.R_KM))
    a=m.xyz([0,0],[0,1])
    d=np.linalg.norm(a[1]-a[0])
    dist=2*m.R_KM*np.arcsin(d/2)
    assert 110<dist<112


def test_wwtp_proximity_counts_facility_and_outfall_not_equivalent():
    p=pd.DataFrame({"port_id":[1,2],
                    "latitude":[0,0],"longitude":[0,2]})
    h=pd.DataFrame({"active_likely":[True,True,False],
        "qual_loc_1or2":[True,False,True],
        "LAT_WWTP":[0,0,0],"LON_WWTP":[0.03,0.1,2],
        "LAT_OUT":[0,0,0],"LON_OUT":[0.3,0.5,2],
        "flow_m3day_nonnegative":[120.,45.,1000.]})
    f=m.wwtp_exposure(p,h)
    assert f.loc[0,"WWTP_count_5km"]==1
    assert f.loc[0,"WWTP_count_25km"]==2
    assert f.loc[0,"WWTP_daily_flow_sum_m3_25km"]==165
    assert f.loc[0,"outfall_count_25km"]==0
    assert f.loc[1,"WWTP_count_25km"]==0
    limited=m.wwtp_exposure(p,h,low_quality_excluded=True)
    assert limited.loc[0,"WWTP_count_25km"]==1


def test_country_fixed_effects_regression_sign_and_permutation():
    rng=np.random.default_rng(1234)
    n=240
    country=np.repeat(["US","CA","GB","FR","DK","JP"],40)
    x=np.tile(np.linspace(0,2,40),6)
    lat=np.tile(np.linspace(26,53,40),6)
    lon=np.tile(np.linspace(-25,25,40),6)
    labels=np.repeat(np.arange(6)*.6,40)
    y=1.5*x+labels+rng.normal(scale=.07,size=n)
    df=pd.DataFrame({"country":country,"latitude":lat,"longitude":lon,
                    "all_taxonomically_indexed":100,
                    "ln_WWTP_25km":x,"ln_alien_to_native":y})
    val=m.within_country_OLS(df,"ln_WWTP_25km","ln_alien_to_native",n_perm=99)
    assert val["status"]=="ESTIMATED_EXPLORATORY_COUNTRY_FE"
    assert 1.1<val["coefficient_per_log1p_exposure"]<1.9
    assert val["n_countries"]==6 and val["n_ports"]==n
    assert val["within_country_label_permutation_two_sided_p_exploratory"]<.05


def test_nonoperational_facilities_are_never_main_exposure():
    assert m.ACTIVE=={"operational","not reported"}
    assert "projected" not in m.ACTIVE
    assert "decommissioned" not in m.ACTIVE
