"""Offline checks: urban source and butterfly host background N alignment."""
import importlib.util,sys,math
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_port_urban_n_support.py"
s=importlib.util.spec_from_file_location("urban_n_support",f)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def test_urban_rural_strata_defined_without_observed_plants():
    p={"built_fraction_among_land_pixels":.6}
    c={"built_fraction_among_land_pixels":.001}
    assert m.classify_land_use(p,c)=="PORT_BUILT20_REMOTE_UNDER5"
    c["built_fraction_among_land_pixels"]=.25
    assert m.classify_land_use(p,c)=="BOTH_BUILT20"
    p["built_fraction_among_land_pixels"]=.02
    assert m.classify_land_use(p,c)=="MIXED_OR_MODERATE"

def test_rank_with_ties():
    assert m.ranks([4,2,4,1])==[3.5,2,3.5,1]
    assert m.rho([1,2,3,4],[4,3,2,1])==-1
    assert m.rho([1,2,3],[1,2,3]) is None

def test_n_contrast_excludes_undersampled_sources():
    source={("1",2,"port","later_2011_2024"):{"plant_records":40,
       "N_scored_candidate_host_species":2,"N_scored_other_plant_species":5,
       "mean_N_candidate_hosts":7,"mean_N_other_scored_plants":6}}
    assert m.n_contrast(source,("1",2,"port","later_2011_2024"))==1
    source[("1",2,"port","later_2011_2024")]["plant_records"]=20
    assert m.n_contrast(source,("1",2,"port","later_2011_2024")) is None
