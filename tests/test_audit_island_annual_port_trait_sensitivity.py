"""Offline six-primary-outcome multiplicity and incomplete-trait checks."""
import pytest
for dep in ("numpy","pandas","pyreadr"):
    pytest.importorskip(dep)
import importlib.util,sys
from pathlib import Path
f=Path(__file__).resolve().parents[1]/"scripts"/"audit_island_annual_port_trait_sensitivity.py"
sys.path.insert(0,str(f.parent))
spec=importlib.util.spec_from_file_location("annual_audit",f)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def test_BH_q_not_interpreted_as_new_p():
    v=m.b_h([.0085,.019,.0285,.0645,.4085,.4595])
    assert len(v)==6
    assert v[0]>=.0085 and v[1]>=.019
    assert v[2]>=.0285 and all(0<=x<=1 for x in v)
    assert sorted(v)==v
def test_small_cohort_does_not_become_supported_result():
    import pandas as pd
    i=[]
    for x in range(5):
        i.append({"code":f"S{x}","level1":"1","introduced_lifeform_both_known":10,
            "native_lifeform_both_known":10,"introduced_both_trait_coverage":.8,
            "native_both_trait_coverage":.8,
            "introduced_annual_herb_fraction_both_known":.5,
            "native_annual_herb_fraction_both_known":.4,
            "GBIF_allplant_indexed":200,
            "log2ports":x,"log_area":x+1,"introduced_N":20,"native_N":30})
    q=m.eval_subset(pd.DataFrame(i),min_cov=.2,min_species=20)
    assert q["cohort_n"]==0
    assert q["introduced_annual_herb"]["status"]=="INSUFFICIENT_EXCHANGEABLE_GLOBAL_TRAIT_SUPPORT"
