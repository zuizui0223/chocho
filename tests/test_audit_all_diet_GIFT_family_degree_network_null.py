"""Unit checks: conservative full larval-diet family/degree null."""
import pytest
for p in ("numpy","pandas","pyreadr"):pytest.importorskip(p)
import importlib.util,sys
from pathlib import Path
import numpy as np
import pandas as pd
src=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(src))
sp=importlib.util.spec_from_file_location("all_host_trait_null",src/"audit_all_diet_GIFT_family_degree_network_null.py")
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def toy():
    n=40;plants=28
    d=pd.DataFrame([{
        "species":f"G{i} s","genus":f"G{i//2}",
        "disturbance_binary":float(i%2),
        "log_island_added":np.log1p(i%11),
        "log_global_gain":np.log1p(i%15)/2,
        "log_diet_hosts":np.log1p(4+i%6),
        "log_native":np.log1p(10+i),
        "log_host_families":np.log1p(2+i%4),
        "coverage":.65+.01*(i%5)
        } for i in range(n)])
    pl=pd.DataFrame([{
        "plant_id":str(i),
        "annual_herb":float(i%2),
        "family":"Fabaceae" if i<15 else "Poaceae",
        "degree_bin":"degree1" if i%3 else "degree2to3",
        "native_bin":"native0to3" if i%5 else "native4to20",
        } for i in range(plants)])
    M=np.zeros((n,plants),dtype=float)
    for i in range(n):
        M[i,[(i+2*j)%plants for j in range(5)]]=1
    return d,pl,M,M.sum(axis=1)

def test_full_host_plant_shuffle_is_not_synthetic_butterfly_repetition():
    d,p,M,den=toy()
    a=m.pnull(d,p,M,den,outcome="log_island_added",stratify="family_degree",B=69)
    b=m.pnull(d,p,M,den,outcome="log_island_added",stratify="family_degree_native",B=69)
    assert a["permutations"]==69
    assert 0<=a["p_two_centered_descriptive"]<=1
    assert 0<=b["p_two_centered_descriptive"]<=1
    assert a["network_plant_taxa_known_both"]==28 if "network_plant_taxa_known_both" in a else a["known_plant_species_count"]==28

def test_permutations_reject_wrong_model():
    d,p,M,den=toy()
    with pytest.raises(ValueError,match="Unexpected"):
        m.pnull(d,p,M,den,outcome="log_island_added",stratify="unrelated",B=3)
