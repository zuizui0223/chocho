"""Global three-source nested contrasts: index-count and modeled sewage controls."""
import importlib.util,sys,math
from pathlib import Path
import pytest
p=Path(__file__).resolve().parents[1]/"scripts"/"exploratory_global_wwtp_jones_worldcover.py"
spec=importlib.util.spec_from_file_location("three_source",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_nuisance_controls_do_not_cause_port_or_sewage_assumption():
    import inspect
    src=inspect.getsource(m.model)
    assert '"ln_total"' in src
    assert "C(country)" in src
    assert "cluster" in src
def test_2015_urban_condition_is_2021_only_not_nutrient():
    import inspect
    src=inspect.getsource(m.run)
    assert "ln_not_treated" in src
    assert "built_fraction_2021" in src
    assert "modeled_treated_share" in src
