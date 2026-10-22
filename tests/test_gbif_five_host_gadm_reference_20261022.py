"""No-network exact hypergeometric and conservative candidate-plant checks."""
from __future__ import annotations
import importlib.util
from fractions import Fraction
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=ROOT/"scripts/analyze_gbif_five_host_gadm_reference_20261022.py"
spec=importlib.util.spec_from_file_location("local_gadm",s)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def row(N,H,j,O):
    return {"all_butterfly_gadm1_records":N,
            "all_butterfly_records_in_plant_gadm1":H,
            "immature_source_labels":j,
            "immature_labels_in_plant_detected_gadm1":O}

def test_exact_conditional_hypergeom_matches_toy():
    r=[row(4,2,1,1)]
    d=m.exact_reference(r)
    assert d==[Fraction(1,2),Fraction(1,2)]
    v=m.summarize_reference(r)
    assert v["null_expected_labels"]==.5
    assert v["observed_minus_reference_pp"]==50
    assert v["conditional_finite_label_reassignment_upper_tail_NOT_ecological_p"]==.5

def test_two_stage_finite_margins_are_independently_conditional():
    r=[row(4,2,1,1),row(4,2,1,0)]
    d=m.exact_reference(r)
    assert d==[Fraction(1,4), Fraction(1,2),Fraction(1,4)]
    v=m.summarize_reference(r)
    assert v["observed_minus_reference_pp"]==0
    assert sum(d)==1

def test_reject_ineligible_sources_and_invalid_hypergeom():
    with pytest.raises(ValueError):
        m.analyse({"schema":"changed"},{},{})
    with pytest.raises(ValueError):
        m.exact_reference([row(3,4,1,1)])
    with pytest.raises(ValueError):
        m.exact_reference([row(3,2,4,1)])
